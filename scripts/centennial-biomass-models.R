
library(ncdf4)
library(ggplot2); library(gridExtra); library(scales)
library(mgcv)
library(plyr); library(parallel)
library(maps)
library(corrplot)
library(RCurl)
library(sp)
library(udunits2)
library(fields)

#### Christy's Functions with knots changed to 100
calc.derivs <- function(model.gam, newdata, vars, n=100, eps=1e-7, alpha=0.05, lwr=NULL, upr=NULL, return.sims=F){
  # Note: this function can be used to generate a 95% CI on the full model.gam OR terms
  
  # -----------
  # Getting the derivative of a GAM and the confidence around that derivative to statistically detect periods of change
  # This is following Gavin Simpson's post here: 
  # http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/
  # His handy-dandy functions can be found here: https://github.com/gavinsimpson/random_code/
  #      This script is based off of derivFun.R
  # -----------
  library(MASS)
  set.seed(321)
  
  # If the model.gam is a mixed model.gam (gamm) rather than a normal gam, extract just the gam portion
  if(class(model.gam)[[1]]=="gamm") model.gam <- model.gam$gam
  
  # If we're interested in an even split to our upper & lower bounds based on our alpha (upper/lower not specified), 
  # calculate them here
  if(is.null(lwr) & is.null(upr)){ lwr=alpha/2; upr= 1-alpha/2 }
  
  # Determining what the terms in our model are
  m.terms <- attr(terms(model.gam), "term.labels")
  
  # finding which columns are numeric
  df.model <- model.frame(model.gam)
  cols.num <- vector()
  for(j in 1:ncol(df.model)){
    if(is.numeric(df.model[,j])) cols.num <- c(cols.num, names(df.model)[j])
  }
  
  # Generate a random distribution of betas using the covariance matrix
  coef.gam <- coef(model.gam)
  Rbeta <- mvrnorm(n=n, model.gam$coefficients, model.gam$Vp)
  
  # From derivFun.R, "Deriv"
  # Create the prediction matrices
  X0 <- predict(model.gam, newdata=newdata, type="lpmatrix")
  
  # Create a new data frame where the numbers are shifted just a bit 
  #  so we can look at the difference (i.e. slope)
  newD <- newdata
  newD[,m.terms[m.terms %in% cols.num]] <- newdata[,m.terms[m.terms %in% cols.num]]+eps
  
  X1 <- predict(model.gam, newdata=newD, type="lpmatrix")
  
  # Getting the "finite difference approximation of the first derivative"
  Xp <- (X1 - X0) / eps # Change in Y per unit X
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  
  for(v in vars) {
    Xi <- Xp * 0 # zeroing out our data frame 
    want <- which(substr(names(coef.gam),1,(nchar(v)+3))==paste0("s(",v,")")) # Finding which columns belong to this factor
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(model.gam)
    
    # Generating a distribution of simulated derivatives
    sim.tmp <- data.frame(Xp[,want] %*% t(Rbeta[,want]) )
    sim.mean <- apply(sim.tmp, 1, mean)
    sim.lwr <- apply(sim.tmp, 1, quantile, lwr)
    sim.upr <- apply(sim.tmp, 1, quantile, upr)
    sig <- as.factor(ifelse(sim.lwr*sim.upr>0, "*", NA))
    
    df.tmp <- data.frame(newdata, 
                         mean=sim.mean,
                         lwr=sim.lwr,
                         upr=sim.upr,
                         sig=sig,
                         var=as.factor(v))
    
    sim.tmp$var <- as.factor(v)
    
    if(v == vars[1]){ 
      df.out <- df.tmp 
      df.sim <- sim.tmp
    } else {
      df.out <- rbind(df.out, df.tmp)
      df.sim <- rbind(df.sim, sim.tmp)
    }
    
  }
  if(return.sims==T){
    out <- list()
    out[["ci"]]   <- df.out
    out[["sims"]] <- df.sim
  } else {
    out <- df.out
  }
  
  return(out)
}
calc.stability <- function(x, width){
  if(length(which(is.na(x)))==length(x)){
    return(NA)
  }else{
    dat.tmp <- data.frame(Y=x, Year=1:length(x))
    k.use = round(length(x)/width, 0)
    mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
    mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
    return(mod.deriv)
  }

}

#### ALL BIOMASS FILES IN kgC/m^2
ed.agb <- readRDS("Data/ED2.AGB.rds")
lpj.agb <- readRDS("Data/LPJ-GUESS.AGB.rds")
lpj.wsl.agb <- readRDS("Data/LPJ-WSL.v1.TotLivBiom.rds")
triffid.agb <- readRDS("Data/TRIFFID.TotLivBio_PFT.rds")
linkages.agb <- readRDS("Data/PalEON_regional_LINKAGES.AGB.rds")

load("Data/PalEON_siteInfo_all.RData")
paleon <- read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/EcosystemStability/master/data/paleon_site_info.csv'))


ed.df <- data.frame(time = seq(1,13932,by = 1161),ed = ed.agb[seq(1,13932,by = 1161),])
lpj.df <- data.frame(time = seq(1,1161,100),lpj = lpj.agb[seq(1,1161,100),,13])

lpj.wsl.df <- data.frame(time = seq(1,1161,100),lpj = lpj.wsl.agb[seq(1,1161,100),])

triffid.agb[which(triffid.agb<0)] <- NA
trif<-apply(triffid.agb,2,rowSums,na.rm=TRUE)
trif.sum <- apply(triffid.agb[seq(1,13932,by = 1161),,1:5],2,rowSums)
trif.sum[which(trif.sum<0)] <- NA
triffid.df <- data.frame(time = seq(1,13932,by = 1161),triffid = trif.sum)

linkages.agb[which(linkages.agb<0)] <- NA
linkages.df <- data.frame(time = seq(1,1161,100),linkages = linkages.agb[seq(1,1161,100),])


linkages.s <- ed.s <- lpj.s <- lpj.wsl.s <- triffid.s <- list()

### Takes awhile to run
for(i in 1:ncol(linkages.agb)){
    linkages.s[[i]] <- calc.stability(linkages.agb[,i])
    ed.s[[i]] <- calc.stability(ed.agb[,i])
    lpj.s[[i]] <- calc.stability(lpj.agb[,i,13])
    lpj.wsl.s[[i]] <- calc.stability(lpj.wsl.agb[,i])
    triffid.s[[i]] <- calc.stability(trif[,i])
  print(i)
}

stability.models = list(ed.s = abs(ed.s), lpj.s = abs(lpj.s), lpj.wsl.s = abs(lpj.wsl.s), triffid.s = abs(triffid.s), linkages.s = abs(linkages.s))
save(stability.models,
     file = 'stability.models.Rdata')

keep.derivs <- keep.signif <- matrix(NA,254,5)

colnames(keep.derivs) <- c('ed.mean.deriv','lpj.mean.deriv','lpj.wsl.mean.deriv','triffid.mean.deriv','linkages.mean.deriv')
colnames(keep.signif) <- c('ed.nsig','lpj.nsig','lpj.wsl.nsig','triffid.nsig','linkages.nsig')


stability.df <- cbind(paleon$lat, paleon$lon, keep.derivs, keep.signif)

write.csv(stability.df, file='models.biomass.stability.csv')




