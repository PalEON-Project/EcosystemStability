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
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"

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
# Vars: GPP, Biomass
ed.npp1  <- readRDS(file.path(path.data, "ED2/ED2.NPP.rds"))
ed.bm1 <- readRDS(file.path(path.data, "ED2/ED2.AGB.rds"))
lpjg.npp1 <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.NPP.rds"))
lpjg.bm <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.AGB.rds"))
lpjw.npp1 <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.NPP.rds"))
lpjw.bm <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.TotLivBiom.rds"))
link.npp <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.NPP.rds"))
link.bm <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.AGB.rds"))

lpjg.bm <- lpjg.bm[,,dim(lpjg.bm)[3]]

# load in the paleon domain info;
# This got generated using domain_environment_extraction.R
paleon <- read.csv(file.path(path.repo, "data/paleon_models_environment_master.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
summary(paleon)

# Load the driver info
tair    <- readRDS(file.path(path.data, "Met/tair_all.rds"))
precipf <- readRDS(file.path(path.data, "Met/precipf_all.rds"))
pdsi    <- readRDS(file.path(path.data, "Met/pdsi_calib_1931-1990_all.rds"))



# aggregate to annual resolution
pdsi.ann <- tair.ann <- tair.jja <- precip.ann <- precip.jja <- matrix(ncol=ncol(tair), nrow=length(yrs))
ed.npp <- ed.bm <- lpjg.npp <- lpjw.npp <- matrix(ncol=ncol(ed.npp1), nrow=length(yrs))

dimnames(tair.ann)[[1]] <- yrs
for(i in 1:length(yrs)){
  # Generating indices for the cells we want to aggregate across
  rows.yrs <- which(time.mos$year==yrs[i])
  rows.jja <- which(time.mos$year==yrs[i] & time.mos$month %in% c(6:8))
  
  # doing the aggregation
  pdsi.ann  [i,] <- colMeans(pdsi   [rows.yrs,])
  tair.ann  [i,] <- colMeans(tair   [rows.yrs,])
  precip.ann[i,] <- colMeans(precipf[rows.yrs,])
  
  ed.npp  [i,] <- colMeans(ed.npp1  [rows.yrs,])
  ed.bm   [i,] <- colMeans(ed.bm1   [rows.yrs,])
  lpjg.npp[i,] <- colMeans(lpjg.npp1[rows.yrs,])
  lpjw.npp[i,] <- colMeans(lpjw.npp1[rows.yrs,])
}

# -------------------------------------------

# -------------------------------------------
# 1. Calculating and Mapping met stability
# -------------------------------------------
# source("R/0_TimeAnalysis.R")
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))

calc.stability <- function(x){
  dat.tmp <- data.frame(Y=x, Year=1:length(x))
  k.use=round(length(x)/25, 0)
  mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
  mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
  return(mod.deriv)
}


# -------------------
# 1A. Met 
# -------------------
pdsi.list1 <- list()
tair.list1 <- list()
precip.list1 <- list()
for(i in 1:ncol(pdsi.ann)){
  pdsi.list1  [[i]] <- pdsi.ann  [which(yrs<1850), i]
  tair.list1  [[i]] <- tair.ann  [which(yrs<1850), i]
  precip.list1[[i]] <- precip.ann[which(yrs<1850), i]
}

pdsi.out1   <- mclapply(pdsi.list1  , calc.stability, mc.cores=8)
tair.out1   <- mclapply(tair.list1  , calc.stability, mc.cores=8)
precip.out1 <- mclapply(precip.list1, calc.stability, mc.cores=8)

# Plugging in the mean absolute value of the derivative
for(i in 1:length(pdsi.out1)){
  paleon[i,"pdsi.deriv"  ] <- mean(abs(pdsi.out1  [[i]]$mean))
  paleon[i,"tair.deriv"  ] <- mean(abs(tair.out1  [[i]]$mean))
  paleon[i,"precip.deriv"] <- mean(abs(precip.out1[[i]]$mean))
  paleon[i,"pdsi.nyr"    ] <- length(pdsi.out1[[i]][!is.na(pdsi.out1[[i]]$sig),"sig"])
  paleon[i,"tair.nyr"    ] <- length(tair.out1[[i]][!is.na(tair.out1[[i]]$sig),"sig"])
  paleon[i,"precip.nyr"  ] <- length(precip.out1[[i]][!is.na(precip.out1[[i]]$sig),"sig"])
}

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
# Tricking missing points into having somethign just to not have things crash
ed.npp[is.na(ed.npp)] <- -9999
ed.bm[is.na(ed.bm)]  <- -9999

ed.npp.list <- list()
ed.bm.list <- list()
lpjg.npp.list <- list()
lpjg.bm.list <- list()
lpjw.npp.list <- list()
lpjw.bm.list <- list()
link.npp.list <- list()
link.bm.list <- list()
for(i in 1:ncol(ed.npp)){
  ed.npp.list  [[i]] <- ed.npp  [which(yrs<1850), i]
  ed.bm.list   [[i]] <- ed.bm   [which(yrs<1850), i]
  lpjg.npp.list[[i]] <- lpjg.npp[which(yrs<1850), i]
  lpjg.bm.list [[i]] <- lpjg.bm [which(yrs<1850), i]
  lpjw.npp.list[[i]] <- lpjw.npp[which(yrs<1850), i]
  lpjw.bm.list [[i]] <- lpjw.bm [which(yrs<1850), i]
  link.npp.list[[i]] <- link.npp[which(yrs<1850), i]
  link.bm.list [[i]] <- link.bm [which(yrs<1850), i]
}

ed.npp.out   <- mclapply(ed.npp.list, calc.stability, mc.cores=8)
ed.bm.out    <- mclapply(ed.bm.list , calc.stability, mc.cores=8)
lpjg.npp.out <- mclapply(lpjg.npp.list, calc.stability, mc.cores=8)
lpjg.bm.out  <- mclapply(lpjg.bm.list , calc.stability, mc.cores=8)
lpjw.npp.out <- mclapply(lpjw.npp.list, calc.stability, mc.cores=8)
lpjw.bm.out  <- mclapply(lpjw.bm.list , calc.stability, mc.cores=8)
link.npp.out <- mclapply(link.npp.list, calc.stability, mc.cores=8)
link.bm.out  <- mclapply(link.bm.list , calc.stability, mc.cores=8)

# Plugging in the mean absolute value of the derivative
for(i in 1:length(ed.npp.out)){
  
  if(min(ed.npp[,i])>-9999){
    paleon.models[i,"ed.npp.deriv"  ] <- mean(abs(ed.npp.out  [[i]]$mean))
    paleon.models[i,"ed.bm.deriv"   ] <- mean(abs(ed.bm.out   [[i]]$mean))
    paleon.models[i,"ed.bm.nyr"     ] <- length(ed.bm.out [[i]][!is.na(ed.bm.out [[i]]$sig),"sig"])
    paleon.models[i,"ed.npp.nyr"    ] <- length(ed.npp.out[[i]][!is.na(ed.npp.out[[i]]$sig),"sig"])
  }
  paleon.models[i,"lpjg.npp.deriv"] <- mean(abs(lpjg.npp.out[[i]]$mean))
  paleon.models[i,"lpjg.bm.deriv" ] <- mean(abs(lpjg.bm.out [[i]]$mean))
  paleon.models[i,"lpjg.bm.nyr"     ] <- length(lpjg.bm.out [[i]][!is.na(lpjg.bm.out [[i]]$sig),"sig"])
  paleon.models[i,"lpjg.npp.nyr"    ] <- length(lpjg.npp.out[[i]][!is.na(lpjg.npp.out[[i]]$sig),"sig"])
  
  paleon.models[i,"lpjw.npp.deriv"] <- mean(abs(lpjw.npp.out[[i]]$mean))
  paleon.models[i,"lpjw.bm.deriv" ] <- mean(abs(lpjw.bm.out [[i]]$mean))
  paleon.models[i,"lpjw.bm.nyr"     ] <- length(lpjw.bm.out [[i]][!is.na(lpjw.bm.out [[i]]$sig),"sig"])
  paleon.models[i,"lpjw.npp.nyr"    ] <- length(lpjw.npp.out[[i]][!is.na(lpjw.npp.out[[i]]$sig),"sig"])
  
  if(min(link.npp[,i])>-9999){
    paleon.models[i,"link.npp.deriv"] <- mean(abs(link.npp.out[[i]]$mean))
    paleon.models[i,"link.bm.deriv" ] <- mean(abs(link.bm.out [[i]]$mean))
    paleon.models[i,"link.bm.nyr"     ] <- length(link.bm.out [[i]][!is.na(link.bm.out [[i]]$sig),"sig"])
    paleon.models[i,"link.npp.nyr"    ] <- length(link.npp.out[[i]][!is.na(link.npp.out[[i]]$sig),"sig"])
    
  }
}

mods <- c("ed", "lpjg", "lpjw", "link")
model.stability <- data.frame(paleon.models[,c("lon", "lat", "latlon")],
                              Model = rep(c("ED2", "LPJ-GUESS", "LPJ-WSL", "LINKAGES"), each=nrow(paleon.models)),
                              deriv.bm  = stack(paleon.models[,paste0(mods, ".bm.deriv" )])[,1],
                              deriv.npp = stack(paleon.models[,paste0(mods, ".npp.deriv")])[,1],
                              bm.nyr    = stack(paleon.models[,paste0(mods, ".bm.nyr" )])[,1],
                              npp.nyr   = stack(paleon.models[,paste0(mods, ".npp.nyr" )])[,1]
                              )
summary(model.stability)
# -------------------

summary(paleon)
summary(model.stability)

write.csv(paleon, file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers.csv"), row.names=F)
write.csv(model.stability, file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models.csv"), row.names=F)
# -------------------------------------------
