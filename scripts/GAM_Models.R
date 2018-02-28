# -------------------------------------------
# Assessing the scales and magnitude of climate change variability over the past millennium in the MIP drivers
# Author: Christy Rollinson, crollinson@gmail.com

# 1. Stastical detection of significant change
#    - Use loess or TP regression splines
# 2. Comparison with paleoclimate reconstructions
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra); library(scales)
library(mgcv)
library(plyr); library(parallel)

# Setting the path to this repository
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_Regional_Extract/"
path.data2 <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"

# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
path.gamm.func <- "~/Desktop/Research/R_Functions/"
# mip.utils <- "~/Dropbox/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

path.co2 <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/co2/paleon_annual_co2.nc"
# -------------------------------------------

# -------------------------------------------
# Define some useful variables
# -------------------------------------------
# Set up some time variables just to help with indexing
yrs <- 850:2010
mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)

us <- map_data("state")
# -------------------------------------------

# -------------------------------------------
# Covert raw monthly drivers to annual
# Base dataset = phase2_met_regional_v2_monthly (available on iPlant)
# -------------------------------------------
load(file.path(path.data, "PalEON_siteInfo_all.RData"))
paleon.models <- paleon
paleon.models <- paleon.models[,1:4]
paleon.models$latlon <- as.factor(paleon.models$latlon)
summary(paleon.models)

# Loading in the model output
# Vars (the cascade): GPP --> NPP --> NEE --> LAI --> AGB
ed.gpp1  <- readRDS(file.path(path.data, "ED2/ED2.v1.2017-06-23.GPP.rds"))
ed.npp1  <- readRDS(file.path(path.data, "ED2/ED2.v1.2017-06-23.NPP.rds"))
ed.nee1  <- readRDS(file.path(path.data, "ED2/ED2.v1.2017-06-23.NEE.rds"))
ed.lai1  <- readRDS(file.path(path.data, "ED2/ED2.v1.2017-06-23.LAI.rds"))
ed.bm1 <- readRDS(file.path(path.data, "ED2/ED2.v1.2017-06-23.AGB.rds"))

lpjg.gpp1 <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.GPP.rds"))
lpjg.npp1 <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.NPP.rds"))
lpjg.nee1 <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.NEE.rds"))
lpjg.lai1 <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.LAI.rds"))
lpjg.bm <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.AGB.rds"))

lpjw.gpp1 <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.GPP.rds"))
lpjw.npp1 <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.NPP.rds"))
lpjw.nee1 <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.NEE.rds"))
# lpjw.lai1 <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.LAI.rds")) # Doesn't Exist
lpjw.bm <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.TotLivBiom.rds"))

# link.gpp <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.GPP.rds")) # LINKAGES doens't do GPP
link.npp <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.NPP.rds"))
link.nee <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.NEE.rds"))
link.lai <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.LAI.rds"))
link.bm <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.AGB.rds"))

# triff.gpp1 <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.GPP.rds")) # Not Extracted
triff.npp1 <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.NPP.rds"))
# triff.nee1 <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.NEE.rds"))  # Not Extracted
triff.lai1 <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.LAI.rds"))
triff.bm1 <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.TotLivBio_PFT.rds"))

# Calculating GPP & NEE
triff.autoresp <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.AutoResp.rds"))
triff.hetrresp <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.HeteroResp.rds"))

triff.gpp1 <- triff.npp1 + triff.autoresp
triff.nee1 <- triff.hetrresp - triff.npp1

rm(triff.autoresp, triff.hetrresp)

# Unit Check
mean(ed.bm1[ed.bm1>0], na.rm=T)
mean(lpjg.bm[lpjg.bm>0], na.rm=T)
mean(lpjw.bm[lpjw.bm>0], na.rm=T)
mean(link.bm[link.bm>0], na.rm=T)
mean(triff.bm1[triff.bm1>0], na.rm=T)

mean(ed.npp1[ed.npp1>0], na.rm=T)
mean(lpjg.npp1[lpjg.npp1>0], na.rm=T)
mean(lpjw.npp1[lpjw.npp1>0], na.rm=T)
mean(link.npp[link.npp>0], na.rm=T)
mean(triff.npp1[triff.npp1>0], na.rm=T)

# LINKAGES BM & NPP look like they need to be divided by 10, but Ann says no
# link.npp <- link.npp * 0.1
# link.bm <- link.bm * 0.1

# Doing a bit of formatting
lpjg.bm <- lpjg.bm[,,dim(lpjg.bm)[3]] # Pull total Biomass

# triff.npp1[triff.npp1==-9999] <- NA
# triff.bm1[triff.bm1==-9999] <- NA
triff.lai1 <- apply(triff.lai1, c(1,2), sum)

triff.gpp1 <- triff.gpp1[1:nrow(lpjw.npp1),]
triff.npp1 <- triff.npp1[1:nrow(lpjw.npp1),]
triff.nee1 <- triff.nee1[1:nrow(lpjw.npp1),]
triff.lai1 <- triff.lai1[1:nrow(lpjw.npp1),]
triff.bm1 <- apply(triff.bm1[1:nrow(lpjw.npp1),,], c(1,2), sum, na.rm=T)

# load in the paleon domain info;
# This got generatlpjw using domain_environment_extraction.R
paleon <- read.csv(file.path(path.repo, "data/paleon_models_environment_master.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
summary(paleon)

# Load the driver info
tair    <- readRDS(file.path(path.data2, "Met/tair_all.rds"))
precipf <- readRDS(file.path(path.data2, "Met/precipf_all.rds"))
pdsi    <- readRDS(file.path(path.data2, "Met/pdsi_calib_1931-1990_all.rds"))



# aggregate to annual resolution
pdsi.ann <- tair.ann <- tair.jja <- precip.ann <- precip.jja <- matrix(ncol=ncol(tair), nrow=length(yrs))
ed.gpp    <- ed.npp    <- ed.nee    <- ed.lai    <- ed.bm <- matrix(ncol=ncol(lpjw.npp1), nrow=length(yrs))
lpjg.gpp  <- lpjg.npp  <- lpjg.nee  <- lpjg.lai  <- matrix(ncol=ncol(lpjw.npp1), nrow=length(yrs))
lpjw.gpp  <- lpjw.npp  <- lpjw.nee  <-              matrix(ncol=ncol(lpjw.npp1), nrow=length(yrs))
triff.gpp <- triff.npp <- triff.nee <- triff.lai <- triff.bm <- matrix(ncol=ncol(lpjw.npp1), nrow=length(yrs))

dimnames(tair.ann)[[1]] <- yrs
for(i in 1:length(yrs)){
  # Generating indices for the cells we want to aggregate across
  rows.yrs <- which(time.mos$year==yrs[i])
  rows.jja <- which(time.mos$year==yrs[i] & time.mos$month %in% c(6:8))
  
  # doing the aggregation
  pdsi.ann  [i,] <- colMeans(pdsi   [rows.jja,])
  tair.ann  [i,] <- colMeans(tair   [rows.yrs,])
  precip.ann[i,] <- colMeans(precipf[rows.yrs,])
  
  ed.gpp   [i,] <- colMeans(ed.gpp1   [rows.yrs,])
  ed.npp   [i,] <- colMeans(ed.npp1   [rows.yrs,])
  ed.nee   [i,] <- colMeans(ed.nee1   [rows.yrs,])
  ed.lai   [i,] <- colMeans(ed.lai1   [rows.yrs,])
  ed.bm    [i,] <- colMeans(ed.bm1    [rows.yrs,])
  
  lpjg.gpp[i,] <- colMeans(lpjg.gpp1[rows.yrs,])
  lpjg.npp[i,] <- colMeans(lpjg.npp1[rows.yrs,])
  lpjg.nee[i,] <- colMeans(lpjg.nee1[rows.yrs,])
  lpjg.lai[i,] <- colMeans(lpjg.lai1[rows.yrs,])

  lpjw.gpp[i,] <- colMeans(lpjw.gpp1[rows.yrs,])
  lpjw.npp[i,] <- colMeans(lpjw.npp1[rows.yrs,])
  lpjw.nee[i,] <- colMeans(lpjw.nee1[rows.yrs,])
  # lpjw.lai[i,] <- colMeans(lpjw.lai1[rows.yrs,])
  
  triff.gpp  [i,] <- colMeans(triff.gpp1  [rows.yrs,])
  triff.npp  [i,] <- colMeans(triff.npp1  [rows.yrs,])
  triff.nee  [i,] <- colMeans(triff.nee1  [rows.yrs,])
  triff.lai  [i,] <- colMeans(triff.lai1  [rows.yrs,])
  triff.bm   [i,] <- colMeans(triff.bm1   [rows.yrs,])
}

# -------------------------------------------

# -------------------------------------------
# 1. Calculating and Mapping met stability
# -------------------------------------------
# source("R/0_TimeAnalysis.R")
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))
source(file.path(path.gamm.func, "Calculate_GAMM_Posteriors.R"))

calc.stability <- function(x, width){
  dat.tmp <- data.frame(Y=x, Year=1:length(x))
  k.use=round(length(x)/width, 0)
  mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
  
  yrs.cent <- rev(seq(length(x), 1, by=-width)) # Go backwards so everything lines up in 1850
  
  mod.out <- list()
  
  mod.out$gam.post <- post.distns(mod.gam, newdata=dat.tmp[yrs.cent, ], vars="Year", return.sims=T)$sims # Note col1=X (index); col2=year
  
  mod.out$mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
  
  return(mod.out)
}


# -------------------
# 1A. Met 
# -------------------
pdsi.list1 <- list()
tair.list1 <- list()
precip.list1 <- list()
for(i in 1:ncol(pdsi.ann)){
  # for(i in 1:10){
  pdsi.list1  [[i]] <- pdsi.ann  [which(yrs<=1850), i]
  tair.list1  [[i]] <- tair.ann  [which(yrs<=1850), i]
  precip.list1[[i]] <- precip.ann[which(yrs<=1850), i]
}

pdsi.out1   <- mclapply(pdsi.list1  , calc.stability, mc.cores=8, width=100)
tair.out1   <- mclapply(tair.list1  , calc.stability, mc.cores=8, width=100)
precip.out1 <- mclapply(precip.list1, calc.stability, mc.cores=8, width=100)

# Plugging in the mean absolute value of the derivative
for(i in 1:length(pdsi.out1)){
  pdsi.diff   <- apply(pdsi.out1  [[i]]$gam.post[,3:ncol(pdsi.out1  [[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  tair.diff   <- apply(tair.out1  [[i]]$gam.post[,3:ncol(tair.out1  [[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  precip.diff <- apply(precip.out1[[i]]$gam.post[,3:ncol(precip.out1[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  
  paleon[i,"pdsi.diff"  ] <- mean(abs(pdsi.diff))
  paleon[i,"tair.diff"  ] <- mean(abs(tair.diff))
  paleon[i,"precip.diff"] <- mean(abs(precip.diff))
  
  paleon[i,"pdsi.deriv"  ] <- mean(abs(pdsi.out1  [[i]]$mod.deriv$mean))
  paleon[i,"tair.deriv"  ] <- mean(abs(tair.out1  [[i]]$mod.deriv$mean))
  paleon[i,"precip.deriv"] <- mean(abs(precip.out1[[i]]$mod.deriv$mean))
  paleon[i,"pdsi.nyr"    ] <- length(pdsi.out1[[i]]$mod.deriv[!is.na(pdsi.out1[[i]]$mod.deriv$sig),"sig"])
  paleon[i,"tair.nyr"    ] <- length(tair.out1[[i]]$mod.deriv[!is.na(tair.out1[[i]]$mod.deriv$sig),"sig"])
  paleon[i,"precip.nyr"  ] <- length(precip.out1[[i]]$mod.deriv[!is.na(precip.out1[[i]]$mod.deriv$sig),"sig"])
}

method.comp <- stack(paleon[,c("pdsi.diff", "tair.diff", "precip.diff")])
names(method.comp) <- c("diff", "var")
method.comp$var <- as.factor(unlist(lapply(stringr::str_split(method.comp$var, "[.]"), function(x) {x[1]})))
method.comp$deriv <- stack(paleon[,c("pdsi.deriv", "tair.deriv", "precip.deriv")])[,1]
summary(method.comp)

png(file.path(path.google, "Current Figures/Stability_GAMs", "Diff_v_Deriv_ClimateDrivers.png"), height=6, width=7, unit="in", res=320)
ggplot(data=method.comp) +
  facet_wrap(~var, scales="free", ncol=2) +
  geom_point(aes(x=diff, y=deriv, color=var)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", col="black") + theme_bw()
dev.off()

ggplot(data=paleon) +
  geom_tile(aes(x=lon, y=lat, fill=log(pdsi.deriv))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(paleon$lon), ylim=range(paleon$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(paleon$pdsi.deriv))) + 
  theme_bw()

ggplot(data=paleon) +
  geom_tile(aes(x=lon, y=lat, fill=tair.nyr)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(paleon$lon), ylim=range(paleon$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(paleon$pdsi.deriv))) + 
  theme_bw()

ggplot(data=paleon) +
  geom_tile(aes(x=lon, y=lat, fill=precip.nyr)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(paleon$lon), ylim=range(paleon$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(paleon$pdsi.deriv))) + 
  theme_bw()


ggplot(data=paleon) +
  geom_tile(aes(x=lon, y=lat, fill=log(tair.deriv))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(paleon$lon), ylim=range(paleon$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(paleon$tair.deriv))) + 
  theme_bw()

ggplot(data=paleon) +
  geom_tile(aes(x=lon, y=lat, fill=log(precip.deriv))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(paleon$lon), ylim=range(paleon$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(paleon$precip.deriv))) + 
  theme_bw()
# -------------------

# -------------------
# 1B. Models
# -------------------
# Tricking missing points into having something just to not have things crash
ed.gpp[is.na(ed.gpp)] <- -9999
ed.npp[is.na(ed.npp)] <- -9999
ed.nee[is.na(ed.nee)] <- -9999
ed.lai[is.na(ed.lai)] <- -9999
ed.bm[is.na(ed.bm)]  <- -9999

triff.gpp[is.na(triff.gpp)] <- -9999
triff.npp[is.na(triff.npp)] <- -9999
triff.nee[is.na(triff.nee)] <- -9999
triff.lai[is.na(triff.lai)] <- -9999
triff.bm[is.na(triff.bm)]  <- -9999

ed.gpp.list <- list()
ed.npp.list <- list()
ed.nee.list <- list()
ed.lai.list <- list()
ed.bm.list <- list()

lpjg.gpp.list <- list()
lpjg.npp.list <- list()
lpjg.nee.list <- list()
lpjg.lai.list <- list()
lpjg.bm.list <- list()

lpjw.gpp.list <- list()
lpjw.npp.list <- list()
lpjw.nee.list <- list()
# lpjw.nlai.list <- list()
lpjw.bm.list <- list()

# link.gpp.list <- list()
link.npp.list <- list()
link.nee.list <- list()
link.lai.list <- list()
link.bm.list <- list()

triff.gpp.list <- list()
triff.npp.list <- list()
triff.nee.list <- list()
triff.lai.list <- list()
triff.bm.list <- list()

for(i in 1:ncol(ed.npp)){
  ed.gpp.list  [[i]] <- ed.gpp  [which(yrs<1850), i]
  ed.npp.list  [[i]] <- ed.npp  [which(yrs<1850), i]
  ed.nee.list  [[i]] <- ed.nee  [which(yrs<1850), i]
  ed.lai.list  [[i]] <- ed.lai  [which(yrs<1850), i]
  ed.bm.list   [[i]] <- ed.bm   [which(yrs<1850), i]

  lpjg.gpp.list  [[i]] <- lpjg.gpp  [which(yrs<1850), i]
  lpjg.npp.list  [[i]] <- lpjg.npp  [which(yrs<1850), i]
  lpjg.nee.list  [[i]] <- lpjg.nee  [which(yrs<1850), i]
  lpjg.lai.list  [[i]] <- lpjg.lai  [which(yrs<1850), i]
  lpjg.bm.list   [[i]] <- lpjg.bm   [which(yrs<1850), i]
  
  lpjw.gpp.list  [[i]] <- lpjw.gpp  [which(yrs<1850), i]
  lpjw.npp.list  [[i]] <- lpjw.npp  [which(yrs<1850), i]
  lpjw.nee.list  [[i]] <- lpjw.nee  [which(yrs<1850), i]
  # lpjw.lai.list  [[i]] <- lpjw.lai  [which(yrs<1850), i]
  lpjw.bm.list   [[i]] <- lpjw.bm   [which(yrs<1850), i]
  
  # link.gpp.list  [[i]] <- link.gpp  [which(yrs<1850), i]
  link.npp.list  [[i]] <- link.npp  [which(yrs<1850), i]
  link.nee.list  [[i]] <- link.nee  [which(yrs<1850), i]
  link.lai.list  [[i]] <- link.lai  [which(yrs<1850), i]
  link.bm.list   [[i]] <- link.bm   [which(yrs<1850), i]
  
  triff.gpp.list  [[i]] <- triff.gpp  [which(yrs<1850), i]
  triff.npp.list  [[i]] <- triff.npp  [which(yrs<1850), i]
  triff.nee.list  [[i]] <- triff.nee  [which(yrs<1850), i]
  triff.lai.list  [[i]] <- triff.lai  [which(yrs<1850), i]
  triff.bm.list   [[i]] <- triff.bm   [which(yrs<1850), i]
}

ed.gpp.out   <- mclapply(ed.gpp.list, calc.stability, mc.cores=8, width=100)
ed.npp.out   <- mclapply(ed.npp.list, calc.stability, mc.cores=8, width=100)
ed.nee.out   <- mclapply(ed.nee.list, calc.stability, mc.cores=8, width=100)
ed.lai.out   <- mclapply(ed.lai.list, calc.stability, mc.cores=8, width=100)
ed.bm.out    <- mclapply(ed.bm.list , calc.stability, mc.cores=8, width=100)

lpjg.gpp.out   <- mclapply(lpjg.gpp.list, calc.stability, mc.cores=8, width=100)
lpjg.npp.out   <- mclapply(lpjg.npp.list, calc.stability, mc.cores=8, width=100)
lpjg.nee.out   <- mclapply(lpjg.nee.list, calc.stability, mc.cores=8, width=100)
lpjg.lai.out   <- mclapply(lpjg.lai.list, calc.stability, mc.cores=8, width=100)
lpjg.bm.out    <- mclapply(lpjg.bm.list , calc.stability, mc.cores=8, width=100)

lpjw.gpp.out   <- mclapply(lpjw.gpp.list, calc.stability, mc.cores=8, width=100)
lpjw.npp.out   <- mclapply(lpjw.npp.list, calc.stability, mc.cores=8, width=100)
lpjw.nee.out   <- mclapply(lpjw.nee.list, calc.stability, mc.cores=8, width=100)
# lpjw.lai.out   <- mclapply(lpjw.lai.list, calc.stability, mc.cores=8, width=100)
lpjw.bm.out    <- mclapply(lpjw.bm.list , calc.stability, mc.cores=8, width=100)

# link.gpp.out   <- mclapply(link.gpp.list, calc.stability, mc.cores=8, width=100)
link.npp.out   <- mclapply(link.npp.list, calc.stability, mc.cores=8, width=100)
link.nee.out   <- mclapply(link.nee.list, calc.stability, mc.cores=8, width=100)
link.lai.out   <- mclapply(link.lai.list, calc.stability, mc.cores=8, width=100)
link.bm.out    <- mclapply(link.bm.list , calc.stability, mc.cores=8, width=100)

triff.gpp.out   <- mclapply(triff.gpp.list, calc.stability, mc.cores=8, width=100)
triff.npp.out   <- mclapply(triff.npp.list, calc.stability, mc.cores=8, width=100)
triff.nee.out   <- mclapply(triff.nee.list, calc.stability, mc.cores=8, width=100)
triff.lai.out   <- mclapply(triff.lai.list, calc.stability, mc.cores=8, width=100)
triff.bm.out    <- mclapply(triff.bm.list , calc.stability, mc.cores=8, width=100)

# Plugging in the mean absolute value of the derivative
# triff.npp2 <- triff.npp
# triff.npp2[triff.npp2==-9999] <- NA
for(i in 1:length(ed.npp.out)){
  
  if(min(ed.npp[,i])>-9999){
    ed.gpp.diff <- apply(ed.gpp.out[[i]]$gam.post[,3:ncol(ed.gpp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    ed.npp.diff <- apply(ed.npp.out[[i]]$gam.post[,3:ncol(ed.npp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    ed.nee.diff <- apply(ed.nee.out[[i]]$gam.post[,3:ncol(ed.nee.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    ed.lai.diff <- apply(ed.lai.out[[i]]$gam.post[,3:ncol(ed.lai.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    ed.bm.diff  <- apply(ed.bm.out [[i]]$gam.post[,3:ncol(ed.bm.out [[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    
    paleon.models[i,"ed.gpp.mean"   ] <- mean(ed.gpp.list[[i]])
    paleon.models[i,"ed.npp.mean"   ] <- mean(ed.npp.list[[i]])
    paleon.models[i,"ed.nee.mean"   ] <- mean(ed.nee.list[[i]])
    paleon.models[i,"ed.lai.mean"   ] <- mean(ed.lai.list[[i]])
    paleon.models[i,"ed.bm.mean"    ] <- mean(ed.bm.list[[i]])

    paleon.models[i,"ed.gpp.diff"   ] <- mean(abs(ed.gpp.diff))
    paleon.models[i,"ed.npp.diff"   ] <- mean(abs(ed.npp.diff))
    paleon.models[i,"ed.nee.diff"   ] <- mean(abs(ed.nee.diff))
    paleon.models[i,"ed.lai.diff"   ] <- mean(abs(ed.lai.diff))
    paleon.models[i,"ed.bm.diff"    ] <- mean(abs(ed.bm.diff))

    paleon.models[i,"ed.gpp.deriv"  ] <- mean(abs(ed.gpp.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"ed.npp.deriv"  ] <- mean(abs(ed.npp.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"ed.nee.deriv"  ] <- mean(abs(ed.nee.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"ed.lai.deriv"  ] <- mean(abs(ed.lai.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"ed.bm.deriv"   ] <- mean(abs(ed.bm.out   [[i]]$mod.deriv$mean))

    paleon.models[i,"ed.gpp.nyr"    ] <- length(ed.gpp.out[[i]]$mod.deriv[!is.na(ed.gpp.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"ed.npp.nyr"    ] <- length(ed.npp.out[[i]]$mod.deriv[!is.na(ed.npp.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"ed.nee.nyr"    ] <- length(ed.nee.out[[i]]$mod.deriv[!is.na(ed.nee.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"ed.lai.nyr"    ] <- length(ed.lai.out[[i]]$mod.deriv[!is.na(ed.lai.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"ed.bm.nyr"     ] <- length(ed.bm.out [[i]]$mod.deriv[!is.na(ed.bm.out [[i]]$mod.deriv$sig),"sig"])
  }
  
  # -----------
  lpjg.gpp.diff <- apply(lpjg.gpp.out[[i]]$gam.post[,3:ncol(lpjg.gpp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  lpjg.npp.diff <- apply(lpjg.npp.out[[i]]$gam.post[,3:ncol(lpjg.npp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  lpjg.nee.diff <- apply(lpjg.nee.out[[i]]$gam.post[,3:ncol(lpjg.nee.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  lpjg.lai.diff <- apply(lpjg.lai.out[[i]]$gam.post[,3:ncol(lpjg.lai.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  lpjg.bm.diff  <- apply(lpjg.bm.out [[i]]$gam.post[,3:ncol(lpjg.bm.out [[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  
  paleon.models[i,"lpjg.gpp.mean"   ] <- mean(lpjg.gpp.list[[i]])
  paleon.models[i,"lpjg.npp.mean"   ] <- mean(lpjg.npp.list[[i]])
  paleon.models[i,"lpjg.nee.mean"   ] <- mean(lpjg.nee.list[[i]])
  paleon.models[i,"lpjg.lai.mean"   ] <- mean(lpjg.lai.list[[i]])
  paleon.models[i,"lpjg.bm.mean"    ] <- mean(lpjg.bm.list[[i]])
  
  paleon.models[i,"lpjg.gpp.diff"   ] <- mean(abs(lpjg.gpp.diff))
  paleon.models[i,"lpjg.npp.diff"   ] <- mean(abs(lpjg.npp.diff))
  paleon.models[i,"lpjg.nee.diff"   ] <- mean(abs(lpjg.nee.diff))
  paleon.models[i,"lpjg.lai.diff"   ] <- mean(abs(lpjg.lai.diff))
  paleon.models[i,"lpjg.bm.diff"    ] <- mean(abs(lpjg.bm.diff))
  
  paleon.models[i,"lpjg.gpp.deriv"  ] <- mean(abs(lpjg.gpp.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjg.npp.deriv"  ] <- mean(abs(lpjg.npp.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjg.nee.deriv"  ] <- mean(abs(lpjg.nee.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjg.lai.deriv"  ] <- mean(abs(lpjg.lai.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjg.bm.deriv"   ] <- mean(abs(lpjg.bm.out   [[i]]$mod.deriv$mean))
  
  paleon.models[i,"lpjg.gpp.nyr"    ] <- length(lpjg.gpp.out[[i]]$mod.deriv[!is.na(lpjg.gpp.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjg.npp.nyr"    ] <- length(lpjg.npp.out[[i]]$mod.deriv[!is.na(lpjg.npp.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjg.nee.nyr"    ] <- length(lpjg.nee.out[[i]]$mod.deriv[!is.na(lpjg.nee.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjg.lai.nyr"    ] <- length(lpjg.lai.out[[i]]$mod.deriv[!is.na(lpjg.lai.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjg.bm.nyr"     ] <- length(lpjg.bm.out [[i]]$mod.deriv[!is.na(lpjg.bm.out [[i]]$mod.deriv$sig),"sig"])
  # -----------
  
  # -----------
  lpjw.gpp.diff <- apply(lpjw.gpp.out[[i]]$gam.post[,3:ncol(lpjw.gpp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  lpjw.npp.diff <- apply(lpjw.npp.out[[i]]$gam.post[,3:ncol(lpjw.npp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  lpjw.nee.diff <- apply(lpjw.nee.out[[i]]$gam.post[,3:ncol(lpjw.nee.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  # lpjw.lai.diff <- apply(lpjw.lai.out[[i]]$gam.post[,3:ncol(lpjw.lai.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  lpjw.bm.diff  <- apply(lpjw.bm.out [[i]]$gam.post[,3:ncol(lpjw.bm.out [[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
  
  paleon.models[i,"lpjw.gpp.mean"   ] <- mean(lpjw.gpp.list[[i]])
  paleon.models[i,"lpjw.npp.mean"   ] <- mean(lpjw.npp.list[[i]])
  paleon.models[i,"lpjw.nee.mean"   ] <- mean(lpjw.nee.list[[i]])
  paleon.models[i,"lpjw.lai.mean"   ] <- NA # mean(lpjw.lai.list[[i]])
  paleon.models[i,"lpjw.bm.mean"    ] <- mean(lpjw.bm.list[[i]])
  
  paleon.models[i,"lpjw.gpp.diff"   ] <- mean(abs(lpjw.gpp.diff))
  paleon.models[i,"lpjw.npp.diff"   ] <- mean(abs(lpjw.npp.diff))
  paleon.models[i,"lpjw.nee.diff"   ] <- mean(abs(lpjw.nee.diff))
  paleon.models[i,"lpjw.lai.diff"   ] <- NA # mean(abs(lpjw.lai.diff))
  paleon.models[i,"lpjw.bm.diff"    ] <- mean(abs(lpjw.bm.diff))
  
  paleon.models[i,"lpjw.gpp.deriv"  ] <- mean(abs(lpjw.gpp.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjw.npp.deriv"  ] <- mean(abs(lpjw.npp.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjw.nee.deriv"  ] <- mean(abs(lpjw.nee.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjw.lai.deriv"  ] <- NA # mean(abs(lpjw.lai.out  [[i]]$mod.deriv$mean))
  paleon.models[i,"lpjw.bm.deriv"   ] <- mean(abs(lpjw.bm.out   [[i]]$mod.deriv$mean))
  
  paleon.models[i,"lpjw.gpp.nyr"    ] <- length(lpjw.gpp.out[[i]]$mod.deriv[!is.na(lpjw.gpp.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjw.npp.nyr"    ] <- length(lpjw.npp.out[[i]]$mod.deriv[!is.na(lpjw.npp.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjw.nee.nyr"    ] <- length(lpjw.nee.out[[i]]$mod.deriv[!is.na(lpjw.nee.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjw.lai.nyr"    ] <- NA # length(lpjw.lai.out[[i]]$mod.deriv[!is.na(lpjw.lai.out[[i]]$mod.deriv$sig),"sig"])
  paleon.models[i,"lpjw.bm.nyr"     ] <- length(lpjw.bm.out [[i]]$mod.deriv[!is.na(lpjw.bm.out [[i]]$mod.deriv$sig),"sig"])
  # -----------
  
  if(min(link.npp[,i])>-9999){
    # link.gpp.diff <- apply(link.gpp.out[[i]]$gam.post[,3:ncol(link.gpp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    link.npp.diff <- apply(link.npp.out[[i]]$gam.post[,3:ncol(link.npp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    link.nee.diff <- apply(link.nee.out[[i]]$gam.post[,3:ncol(link.nee.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    link.lai.diff <- apply(link.lai.out[[i]]$gam.post[,3:ncol(link.lai.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    link.bm.diff  <- apply(link.bm.out [[i]]$gam.post[,3:ncol(link.bm.out [[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    
    paleon.models[i,"link.gpp.mean"   ] <- NA # mean(link.gpp.list[[i]])
    paleon.models[i,"link.npp.mean"   ] <- mean(link.npp.list[[i]])
    paleon.models[i,"link.nee.mean"   ] <- mean(link.nee.list[[i]])
    paleon.models[i,"link.lai.mean"   ] <- mean(link.lai.list[[i]])
    paleon.models[i,"link.bm.mean"    ] <- mean(link.bm.list[[i]])
    
    paleon.models[i,"link.gpp.diff"   ] <- NA # mean(abs(link.gpp.diff))
    paleon.models[i,"link.npp.diff"   ] <- mean(abs(link.npp.diff))
    paleon.models[i,"link.nee.diff"   ] <- mean(abs(link.nee.diff))
    paleon.models[i,"link.lai.diff"   ] <- mean(abs(link.lai.diff))
    paleon.models[i,"link.bm.diff"    ] <- mean(abs(link.bm.diff))
    
    paleon.models[i,"link.gpp.deriv"  ] <- NA # mean(abs(link.gpp.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"link.npp.deriv"  ] <- mean(abs(link.npp.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"link.nee.deriv"  ] <- mean(abs(link.nee.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"link.lai.deriv"  ] <- mean(abs(link.lai.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"link.bm.deriv"   ] <- mean(abs(link.bm.out   [[i]]$mod.deriv$mean))
    
    paleon.models[i,"link.gpp.nyr"    ] <- NA # length(link.gpp.out[[i]]$mod.deriv[!is.na(link.gpp.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"link.npp.nyr"    ] <- length(link.npp.out[[i]]$mod.deriv[!is.na(link.npp.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"link.nee.nyr"    ] <- length(link.nee.out[[i]]$mod.deriv[!is.na(link.nee.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"link.lai.nyr"    ] <- length(link.lai.out[[i]]$mod.deriv[!is.na(link.lai.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"link.bm.nyr"     ] <- length(link.bm.out [[i]]$mod.deriv[!is.na(link.bm.out [[i]]$mod.deriv$sig),"sig"])
  }
  
  if(max(triff.npp[,i])>-9999){
    triff.gpp.diff <- apply(triff.gpp.out[[i]]$gam.post[,3:ncol(triff.gpp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    triff.npp.diff <- apply(triff.npp.out[[i]]$gam.post[,3:ncol(triff.npp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    triff.nee.diff <- apply(triff.nee.out[[i]]$gam.post[,3:ncol(triff.nee.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    triff.lai.diff <- apply(triff.lai.out[[i]]$gam.post[,3:ncol(triff.lai.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    triff.bm.diff  <- apply(triff.bm.out [[i]]$gam.post[,3:ncol(triff.bm.out [[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
    
    paleon.models[i,"triff.gpp.mean"   ] <- mean(triff.gpp.list[[i]])
    paleon.models[i,"triff.npp.mean"   ] <- mean(triff.npp.list[[i]])
    paleon.models[i,"triff.nee.mean"   ] <- mean(triff.nee.list[[i]])
    paleon.models[i,"triff.lai.mean"   ] <- mean(triff.lai.list[[i]])
    paleon.models[i,"triff.bm.mean"    ] <- mean(triff.bm.list[[i]])
    
    paleon.models[i,"triff.gpp.diff"   ] <- mean(abs(triff.gpp.diff))
    paleon.models[i,"triff.npp.diff"   ] <- mean(abs(triff.npp.diff))
    paleon.models[i,"triff.nee.diff"   ] <- mean(abs(triff.nee.diff))
    paleon.models[i,"triff.lai.diff"   ] <- mean(abs(triff.lai.diff))
    paleon.models[i,"triff.bm.diff"    ] <- mean(abs(triff.bm.diff))
    
    paleon.models[i,"triff.gpp.deriv"  ] <- mean(abs(triff.gpp.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"triff.npp.deriv"  ] <- mean(abs(triff.npp.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"triff.nee.deriv"  ] <- mean(abs(triff.nee.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"triff.lai.deriv"  ] <- mean(abs(triff.lai.out  [[i]]$mod.deriv$mean))
    paleon.models[i,"triff.bm.deriv"   ] <- mean(abs(triff.bm.out   [[i]]$mod.deriv$mean))
    
    paleon.models[i,"triff.gpp.nyr"    ] <- length(triff.gpp.out[[i]]$mod.deriv[!is.na(triff.gpp.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"triff.npp.nyr"    ] <- length(triff.npp.out[[i]]$mod.deriv[!is.na(triff.npp.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"triff.nee.nyr"    ] <- length(triff.nee.out[[i]]$mod.deriv[!is.na(triff.nee.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"triff.lai.nyr"    ] <- length(triff.lai.out[[i]]$mod.deriv[!is.na(triff.lai.out[[i]]$mod.deriv$sig),"sig"])
    paleon.models[i,"triff.bm.nyr"     ] <- length(triff.bm.out [[i]]$mod.deriv[!is.na(triff.bm.out [[i]]$mod.deriv$sig),"sig"])
  }
  
}

mods <- c("ed", "lpjg", "lpjw", "link", "triff")
model.stability <- data.frame(paleon.models[,c("lon", "lat", "latlon")],
                              Model = rep(c("ED2", "LPJ-GUESS", "LPJ-WSL", "LINKAGES", "TRIFFID"), each=nrow(paleon.models)),
                              mean.gpp  = stack(paleon.models[,paste0(mods, ".gpp.mean" )])[,1],
                              mean.npp  = stack(paleon.models[,paste0(mods, ".npp.mean" )])[,1],
                              mean.nee  = stack(paleon.models[,paste0(mods, ".nee.mean" )])[,1],
                              mean.lai  = stack(paleon.models[,paste0(mods, ".lai.mean" )])[,1],
                              mean.bm   = stack(paleon.models[,paste0(mods, ".bm.mean"  )])[,1],
                              diff.gpp  = stack(paleon.models[,paste0(mods, ".gpp.diff" )])[,1],
                              diff.npp  = stack(paleon.models[,paste0(mods, ".npp.diff" )])[,1],
                              diff.nee  = stack(paleon.models[,paste0(mods, ".nee.diff" )])[,1],
                              diff.lai  = stack(paleon.models[,paste0(mods, ".lai.diff" )])[,1],
                              diff.bm   = stack(paleon.models[,paste0(mods, ".bm.diff"  )])[,1],
                              deriv.gpp = stack(paleon.models[,paste0(mods, ".gpp.deriv")])[,1],
                              deriv.npp = stack(paleon.models[,paste0(mods, ".npp.deriv")])[,1],
                              deriv.nee = stack(paleon.models[,paste0(mods, ".nee.deriv")])[,1],
                              deriv.lai = stack(paleon.models[,paste0(mods, ".lai.deriv")])[,1],
                              deriv.bm  = stack(paleon.models[,paste0(mods, ".bm.deriv" )])[,1],
                              nyr.gpp   = stack(paleon.models[,paste0(mods, ".gpp.nyr"  )])[,1],
                              nyr.npp   = stack(paleon.models[,paste0(mods, ".npp.nyr"  )])[,1],
                              nyr.nee   = stack(paleon.models[,paste0(mods, ".nee.nyr"  )])[,1],
                              nyr.lai   = stack(paleon.models[,paste0(mods, ".lai.nyr"  )])[,1],
                              nyr.bm    = stack(paleon.models[,paste0(mods, ".bm.nyr"   )])[,1]
)

# We have some lingering NA issues
model.stability[!is.na(model.stability$mean.nee) & model.stability$mean.nee==-9999,c("mean.nee", "diff.nee", "deriv.nee", "nyr.nee")] <- NA
summary(model.stability)
# -------------------

summary(paleon)
summary(model.stability)

write.csv(paleon, file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"), row.names=F)
write.csv(model.stability, file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"), row.names=F)
# -------------------------------------------
