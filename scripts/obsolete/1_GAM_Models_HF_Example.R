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
hf.lat <- 42.54
hf.lon <- -72.18
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
  
  yrs.cent <- rev(seq(length(x), 1, by=-100)) # Go backwards so everything lines up in 1850
  
  mod.out <- list()
  
  mod.out$gam.post <- post.distns(mod.gam, newdata=dat.tmp[, ], vars="Year", return.sims=T) # Note col1=X (index); col2=year
  # mod.out$gam.post <- post.distns(mod.gam, newdata=dat.tmp[yrs.cent, ], vars="Year", return.sims=T) # Note col1=X (index); col2=year
  mod.out$mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
  
  return(mod.out)
}

# -------------------
# 1A. Met 
# -------------------
# Finding which grid cell corresponds to HF
hf.dist.reg <- sqrt((paleon$lat-hf.lat)^2 + (paleon$lon-hf.lon)^2)
hf.ind.reg <- which(hf.dist.reg==min(hf.dist.reg))

# Creating a met list
met.out <- list()
met.out$pdsi   <- pdsi.ann  [which(yrs<=1850), hf.ind.reg]
met.out$precip <- precip.ann[which(yrs<=1850), hf.ind.reg]
met.out$tair   <- tair.ann  [which(yrs<=1850), hf.ind.reg]

yrs.cent <- rev(seq(length(met.out$pdsi), 1, by=-100))

met.out1   <- mclapply(met.out  , calc.stability, mc.cores=min(8, length(met.out)), width=100)

plot(met.out$tair, type="l")
lines(met.out1$tair$gam.post$ci$mean, type="l", lwd=5, col="red")
lines(met.out1$tair$gam.post$ci$lwr, type="l", lwd=2, col="red", lty="dashed")
lines(met.out1$tair$gam.post$ci$upr, type="l", lwd=2, col="red", lty="dashed")
points(met.out1$tair$gam.post$ci$mean[yrs.cent]~yrs.cent, pch=19, col="blue")

# lines(met.out1$tair$mod.deriv$mean + mean(met.out1$tair$mod.deriv$Y), type="l", lwd=2, col="red")
# lines(met.out1$tair$mod.deriv$lwr + mean(met.out1$tair$mod.deriv$Y), type="l", lwd=2, col="red", lty="dashed")
# lines(met.out1$tair$mod.deriv$upr + mean(met.out1$tair$mod.deriv$Y), type="l", lwd=2, col="red", lty="dashed")
# -------------------

# -------------------
# Models
# -------------------
# Finding which grid cell corresponds to HF
hf.dist.mod <- sqrt((paleon.models$lat-hf.lat)^2 + (paleon.models$lon-hf.lon)^2)
hf.ind.mod <- which(hf.dist.mod==min(hf.dist.mod))

# Zeroing out anythign that didn't get run
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


ed.list <- list()
ed.list$gpp <- ed.gpp  [which(yrs<=1850), hf.ind.mod]
ed.list$npp <- ed.npp  [which(yrs<=1850), hf.ind.mod]
ed.list$nee <- ed.nee  [which(yrs<=1850), hf.ind.mod]
ed.list$lai <- ed.lai  [which(yrs<=1850), hf.ind.mod]
ed.list$bm  <- ed.bm   [which(yrs<=1850), hf.ind.mod]

lpjg.list <- list()
lpjg.list$gpp <- lpjg.gpp  [which(yrs<=1850), hf.ind.mod]
lpjg.list$npp <- lpjg.npp  [which(yrs<=1850), hf.ind.mod]
lpjg.list$nee <- lpjg.nee  [which(yrs<=1850), hf.ind.mod]
lpjg.list$lai <- lpjg.lai  [which(yrs<=1850), hf.ind.mod]
lpjg.list$bm  <- lpjg.bm   [which(yrs<=1850), hf.ind.mod]

lpjw.list <- list()
lpjw.list$gpp <- lpjw.gpp  [which(yrs<=1850), hf.ind.mod]
lpjw.list$npp <- lpjw.npp  [which(yrs<=1850), hf.ind.mod]
lpjw.list$nee <- lpjw.nee  [which(yrs<=1850), hf.ind.mod]
lpjw.list$lai <- lpjw.lai  [which(yrs<=1850), hf.ind.mod]
lpjw.list$bm  <- lpjw.bm   [which(yrs<=1850), hf.ind.mod]

link.list <- list()
link.list$gpp <- link.gpp  [which(yrs<=1850), hf.ind.mod]
link.list$npp <- link.npp  [which(yrs<=1850), hf.ind.mod]
link.list$nee <- link.nee  [which(yrs<=1850), hf.ind.mod]
link.list$lai <- link.lai  [which(yrs<=1850), hf.ind.mod]
link.list$bm  <- link.bm   [which(yrs<=1850), hf.ind.mod]

triff.list <- list()
triff.list$gpp <- triff.gpp  [which(yrs<=1850), hf.ind.mod]
triff.list$npp <- triff.npp  [which(yrs<=1850), hf.ind.mod]
triff.list$nee <- triff.nee  [which(yrs<=1850), hf.ind.mod]
triff.list$lai <- triff.lai  [which(yrs<=1850), hf.ind.mod]
triff.list$bm  <- triff.bm   [which(yrs<=1850), hf.ind.mod]

ed.out   <- mclapply(ed.list, calc.stability, mc.cores=min(8, length(ed.list)), width=100)
lpjg.out   <- mclapply(lpjg.list, calc.stability, mc.cores=min(8, length(lpjg.list)), width=100)
lpjw.out   <- mclapply(lpjw.list, calc.stability, mc.cores=min(8, length(lpjw.list)), width=100)
link.out   <- mclapply(link.list, calc.stability, mc.cores=min(8, length(link.list)), width=100)
triff.out   <- mclapply(triff.list, calc.stability, mc.cores=min(8, length(triff.list)), width=100)

plot(ed.list$bm, type="l")
lines(ed.out$bm$gam.post$ci$mean, type="l", lwd=5, col="red")
lines(ed.out$bm$gam.post$ci$lwr, type="l", lwd=2, col="red", lty="dashed")
lines(ed.out$bm$gam.post$ci$upr, type="l", lwd=2, col="red", lty="dashed")
points(ed.out$bm$gam.post$ci$mean[yrs.cent]~yrs.cent, pch=19, col="blue")

plot(ed.out$bm$gam.post$ci$mean[yrs.cent]~yrs.cent, pch=19, col="blue")
lines(ed.list$bm, type="l")
lines(ed.out$bm$gam.post$ci$mean, type="l", lwd=5, col="red")
lines(ed.out$bm$gam.post$ci$lwr, type="l", lwd=2, col="red", lty="dashed")
lines(ed.out$bm$gam.post$ci$upr, type="l", lwd=2, col="red", lty="dashed")
points(ed.out$bm$gam.post$ci$mean[yrs.cent]~yrs.cent, pch=19, col="blue")


plot(lpjg.list$bm, type="l")
lines(lpjg.out$bm$gam.post$ci$mean, type="l", lwd=5, col="red")
lines(lpjg.out$bm$gam.post$ci$lwr, type="l", lwd=2, col="red", lty="dashed")
lines(lpjg.out$bm$gam.post$ci$upr, type="l", lwd=2, col="red", lty="dashed")
points(lpjg.out$bm$gam.post$ci$mean[yrs.cent]~yrs.cent, pch=19, col="blue")

save(ed.list , ed.out   , file=file.path(path.google, "Current Data/Stability_GAMs", "Stability_Example_ED2_HarvardForest.RData"))
save(lpjg.list , lpjg.out , file=file.path(path.google, "Current Data/Stability_GAMs", "Stability_Example_LPJ-GUESS_HarvardForest.RData"))
save(lpjw.list , lpjw.out , file=file.path(path.google, "Current Data/Stability_GAMs", "Stability_Example_LPJ-WSL_HarvardForest.RData"))
save(link.list , link.out , file=file.path(path.google, "Current Data/Stability_GAMs", "Stability_Example_LINKAGES_HarvardForest.RData"))
save(triff.list, triff.out, file=file.path(path.google, "Current Data/Stability_GAMs", "Stability_Example_TRIFFID_HarvardForest.RData"))
# -------------------


# -------------------
# Making a figure that approximates pollen data
# -------------------
df.ed <- data.frame(Year=850:1850, Biomass = ed.list$bm, 
                    gam.mean = ed.out$bm$gam.post$ci$mean, 
                    gam.lwr = ed.out$bm$gam.post$ci$lwr,
                    gam.upr = ed.out$bm$gam.post$ci$upr)
summary(df.ed)

png(file.path(path.google, "Current Figures/GAM_Methods", "HF_TimeSeries1_Points.png"), height=9, width=10, units="in", res=220)
ggplot(data=df.ed) +
  geom_point(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, y=gam.mean), size=10, color="blue") +
  geom_linerange(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, ymin=gam.lwr, ymax=gam.upr), size=5, color="blue") +
  scale_y_continuous(limits=range(df.ed$Biomass), name="Biomass") +
  theme_bw() +
  theme(axis.text=element_text(size=rel(2)),
        axis.title=element_text(size=rel(2.5)))
dev.off()

png(file.path(path.google, "Current Figures/GAM_Methods", "HF_TimeSeries2_Actual.png"), height=9, width=10, units="in", res=220)
ggplot(data=df.ed) +
  geom_line(data=df.ed[,], aes(x=Year, y=Biomass), size=0.5, color="black") +
  geom_point(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, y=gam.mean), size=10, color="blue") +
  geom_linerange(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, ymin=gam.lwr, ymax=gam.upr), size=5, color="blue") +
  scale_y_continuous(limits=range(df.ed$Biomass), name="Biomass") +
  theme_bw() +
  theme(axis.text=element_text(size=rel(2)),
        axis.title=element_text(size=rel(2.5)))
dev.off()

png(file.path(path.google, "Current Figures/GAM_Methods", "HF_TimeSeries3_GAM.png"), height=9, width=10, units="in", res=220)
ggplot(data=df.ed) +
  geom_line(aes(x=Year, y=Biomass), size=0.5, color="black") +
  geom_ribbon(aes(x=Year, ymin=gam.lwr, ymax=gam.upr), alpha=0.5, fill="blue") +
  geom_line(aes(x=Year, y=gam.mean), size=5, color="blue") +
  geom_point(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, y=gam.mean), size=10, color="blue4") +
  geom_linerange(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, ymin=gam.lwr, ymax=gam.upr), size=5, color="blue4") +
  scale_y_continuous(limits=range(df.ed$Biomass), name="Biomass") +
  theme_bw() +
  theme(axis.text=element_text(size=rel(2)),
        axis.title=element_text(size=rel(2.5)))
dev.off()

# -------------------
