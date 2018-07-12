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
# 1. Load and Format Data
# -------------------------------------------
load(file.path(path.data, "PalEON_siteInfo_all.RData"))
paleon.models <- paleon
paleon.models <- paleon.models[,1:4]
paleon.models$latlon <- as.factor(paleon.models$latlon)
summary(paleon.models)

pft.table <- read.csv(file.path(path.google, "Current Data/PFT_Conversions", "ModelPFT_Table.csv"), na.strings="")
pft.table

# Loading in the model output
# Vars (the cascade): GPP --> NPP --> NEE --> LAI --> AGB
ed.fcomp    <- readRDS(file.path(path.data, "ED2/ED2.v1.2017-06-23.Fcomp.rds"))
lpjg.fcomp  <- readRDS(file.path(path.data, "LPJ-GUESS/LPJ-GUESS.Fcomp.rds"))
lpjw.fcomp  <- readRDS(file.path(path.data, "LPJ-WSL/LPJ-WSL.v1.Fcomp.rds"))
link.fcomp  <- readRDS(file.path(path.data, "LINKAGES/PalEON_regional_LINKAGES.Fcomp.rds"))
triff.fcomp <- readRDS(file.path(path.data, "TRIFFID/TRIFFID.Fcomp.rds"))
dim(ed.fcomp)
dim(lpjg.fcomp)
dim(lpjw.fcomp)
dim(link.fcomp)
dim(triff.fcomp)

# ED & TRIFFID need months converted to Years
ed.fcomp.mo <- ed.fcomp
triff.fcomp.mo <- triff.fcomp[1:dim(ed.fcomp.mo)[1],,]

# Aggregate to annual resolution if necessary
ed.fcomp    <- array(dim=c(length(yrs), dim(ed.fcomp.mo)[2:3]))
triff.fcomp <- array(dim=c(length(yrs), dim(triff.fcomp.mo)[2:3]))
for(i in 1:length(yrs)){
  # Generating indices for the cells we want to aggregate across
  rows.yrs <- which(time.mos$year==yrs[i])
  ed.fcomp   [i,,] <- apply(ed.fcomp.mo[rows.yrs,,], c(2,3), mean)
  triff.fcomp[i,,] <- apply(triff.fcomp.mo[rows.yrs,,], c(2,3), mean)
}
dimnames(ed.fcomp)[[1]] <- dimnames(lpjg.fcomp)[[1]] <- dimnames(lpjw.fcomp)[[1]] <- dimnames(link.fcomp)[[1]] <- dimnames(triff.fcomp)[[1]] <- yrs

dimnames(ed.fcomp   )[[3]] <- pft.table$ED2[!is.na(pft.table$ED2)]
dimnames(lpjg.fcomp )[[3]] <- pft.table$LPJ.GUESS[!is.na(pft.table$LPJ.GUESS)]
dimnames(lpjw.fcomp )[[3]] <- pft.table$LPJ.WSL[!is.na(pft.table$LPJ.WSL)]
dimnames(link.fcomp )[[3]] <- pft.table$LINKAGES[!is.na(pft.table$LINKAGES)]
dimnames(triff.fcomp)[[3]] <- pft.table$JULES[!is.na(pft.table$JULES)]

# Subset to the pre-1850 window for simplicity
ed.fcomp    <- ed.fcomp   [which(yrs<1850),,]
lpjg.fcomp  <- lpjg.fcomp [which(yrs<1850),,]
lpjw.fcomp  <- lpjw.fcomp [which(yrs<1850),,]
link.fcomp  <- link.fcomp [which(yrs<1850),,]
triff.fcomp <- triff.fcomp[which(yrs<1850),,]

# Tricking missing points into having something just to not have things crash
ed.fcomp   [is.na(ed.fcomp   )] <- -9999
lpjg.fcomp [is.na(lpjg.fcomp )] <- -9999
lpjw.fcomp [is.na(lpjw.fcomp )] <- -9999
link.fcomp [is.na(link.fcomp )] <- -9999
triff.fcomp[is.na(triff.fcomp)] <- -9999
# -------------------------------------------


# -------------------------------------------
# 2. Calculating Stability
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
# Set up a blank dataframe where everything will go
stab.fcomp <- data.frame()
mods <- c("ed", "lpjg", "lpjw", "link", "triff")

# Need to calculate stability for each PFT in each grid cell

# ------
# ED2 
# ------
pb <- txtProgressBar(min=0, max=dim(ed.fcomp)[[3]], style=3)
for(PFT in 1:dim(ed.fcomp)[3]){
  setTxtProgressBar(pb, PFT)
  
  if(max(ed.fcomp[,,PFT])<=0) next # If this PFT isn't present, skip it!
  
  # Set up a blank data frame 
  pft.out <- data.frame(paleon.models, Model="ED2", PFT=dimnames(ed.fcomp)[[3]][PFT])
  
  # Store things in a list for parallelization
  dat.pft <- list()
  for(i in 1:dim(ed.fcomp)[2]){
    dat.pft[[i]] <- ed.fcomp[,i,PFT]
  }
  
  # Calculate stability for each grid cell in parallel
  fcomp.out <- mclapply(dat.pft, calc.stability, mc.cores=8, width=100)
  
  # Save the output
  for(i in 1:length(fcomp.out)){
    if(min(ed.fcomp[,i,PFT])>-9999){
      fcomp.diff <- apply(fcomp.out[[i]]$gam.post[,3:ncol(fcomp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
      
      pft.out[i,"pft.mean"     ] <- mean(dat.pft[[i]])
      pft.out[i,"pft.diff.abs" ] <- mean(abs(fcomp.diff))
      pft.out[i,"pft.deriv.abs"] <- mean(abs(fcomp.out[[i]]$mod.deriv$mean))
    }  
  }
  
  # Add what we just calculated to everything we've already done
  stab.fcomp <- rbind(stab.fcomp, pft.out)
}
# ------

# ------
# LPJ-GUESS 
# ------
pb <- txtProgressBar(min=0, max=dim(lpjg.fcomp)[[3]], style=3)
for(PFT in 1:dim(lpjg.fcomp)[3]){
  setTxtProgressBar(pb, PFT)
  
  if(max(lpjg.fcomp[,,PFT])<=0) next # If this PFT isn't present, skip it!
  
  # Set up a blank data frame 
  pft.out <- data.frame(paleon.models, Model="LPJ-GUESS", PFT=dimnames(lpjg.fcomp)[[3]][PFT])
  
  # Store things in a list for parallelization
  dat.pft <- list()
  for(i in 1:dim(lpjg.fcomp)[2]){
    dat.pft[[i]] <- lpjg.fcomp[,i,PFT]
  }
  
  # Calculate stability for each grid cell in parallel
  fcomp.out <- mclapply(dat.pft, calc.stability, mc.cores=8, width=100)
  
  # Save the output
  for(i in 1:length(fcomp.out)){
    if(min(lpjg.fcomp[,i,PFT])>-9999){
      fcomp.diff <- apply(fcomp.out[[i]]$gam.post[,3:ncol(fcomp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
      
      pft.out[i,"pft.mean"     ] <- mean(dat.pft[[i]])
      pft.out[i,"pft.diff.abs" ] <- mean(abs(fcomp.diff))
      pft.out[i,"pft.deriv.abs"] <- mean(abs(fcomp.out[[i]]$mod.deriv$mean))
    }  
  }
  
  # Add what we just calculated to everything we've already done
  stab.fcomp <- rbind(stab.fcomp, pft.out)
}
# ------

# ------
# LPJ-WSL 
# ------
pb <- txtProgressBar(min=0, max=dim(lpjw.fcomp)[[3]], style=3)
for(PFT in 1:dim(lpjw.fcomp)[3]){
  setTxtProgressBar(pb, PFT)
  
  if(max(lpjw.fcomp[,,PFT])<=0) next # If this PFT isn't present, skip it!
  
  # Set up a blank data frame 
  pft.out <- data.frame(paleon.models, Model="LPJ-WSL", PFT=dimnames(lpjw.fcomp)[[3]][PFT])
  
  # Store things in a list for parallelization
  dat.pft <- list()
  for(i in 1:dim(lpjw.fcomp)[2]){
    dat.pft[[i]] <- lpjw.fcomp[,i,PFT]
  }
  
  # Calculate stability for each grid cell in parallel
  fcomp.out <- mclapply(dat.pft, calc.stability, mc.cores=8, width=100)
  
  # Save the output
  for(i in 1:length(fcomp.out)){
    if(min(lpjw.fcomp[,i,PFT])>-9999){
      fcomp.diff <- apply(fcomp.out[[i]]$gam.post[,3:ncol(fcomp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
      
      pft.out[i,"pft.mean"     ] <- mean(dat.pft[[i]])
      pft.out[i,"pft.diff.abs" ] <- mean(abs(fcomp.diff))
      pft.out[i,"pft.deriv.abs"] <- mean(abs(fcomp.out[[i]]$mod.deriv$mean))
    }  
  }
  
  # Add what we just calculated to everything we've already done
  stab.fcomp <- rbind(stab.fcomp, pft.out)
}
# ------

# ------
# LINKAGES 
# ------
pb <- txtProgressBar(min=0, max=dim(link.fcomp)[[3]], style=3)
for(PFT in 1:dim(link.fcomp)[3]){
  setTxtProgressBar(pb, PFT)
  
  if(max(link.fcomp[,,PFT])<=0 | toupper(dimnames(link.fcomp)[[3]][PFT])=="TOTAL") next # If this PFT isn't present, skip it!
  
  # Set up a blank data frame 
  pft.out <- data.frame(paleon.models, Model="LINKAGES", PFT=dimnames(link.fcomp)[[3]][PFT])
  
  # Store things in a list for parallelization
  dat.pft <- list()
  for(i in 1:dim(link.fcomp)[2]){
    dat.pft[[i]] <- link.fcomp[,i,PFT]
  }
  
  # Calculate stability for each grid cell in parallel
  fcomp.out <- mclapply(dat.pft, calc.stability, mc.cores=8, width=100)
  
  # Save the output
  for(i in 1:length(fcomp.out)){
    if(min(link.fcomp[,i,PFT])>-9999){
      fcomp.diff <- apply(fcomp.out[[i]]$gam.post[,3:ncol(fcomp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
      
      pft.out[i,"pft.mean"     ] <- mean(dat.pft[[i]])
      pft.out[i,"pft.diff.abs" ] <- mean(abs(fcomp.diff))
      pft.out[i,"pft.deriv.abs"] <- mean(abs(fcomp.out[[i]]$mod.deriv$mean))
    }  
  }
  
  # Add what we just calculated to everything we've already done
  stab.fcomp <- rbind(stab.fcomp, pft.out)
}
# ------

# ------
# TRIFFID 
# ------
pb <- txtProgressBar(min=0, max=dim(triff.fcomp)[[3]], style=3)
for(PFT in 1:dim(triff.fcomp)[3]){
  setTxtProgressBar(pb, PFT)
  
  if(max(triff.fcomp[,,PFT])<=0) next # If this PFT isn't present, skip it!
  
  # Set up a blank data frame 
  pft.out <- data.frame(paleon.models, Model="TRIFFID", PFT=dimnames(triff.fcomp)[[3]][PFT])
  
  # Store things in a list for parallelization
  dat.pft <- list()
  for(i in 1:dim(triff.fcomp)[2]){
    dat.pft[[i]] <- triff.fcomp[,i,PFT]
  }
  
  # Calculate stability for each grid cell in parallel
  fcomp.out <- mclapply(dat.pft, calc.stability, mc.cores=8, width=100)
  
  # Save the output
  for(i in 1:length(fcomp.out)){
    if(min(triff.fcomp[,i,PFT])>-9999){
      fcomp.diff <- apply(fcomp.out[[i]]$gam.post[,3:ncol(fcomp.out[[i]]$gam.post)], 2, function(x) diff(x, na.rm=TRUE)/100)
      
      pft.out[i,"pft.mean"     ] <- mean(dat.pft[[i]])
      pft.out[i,"pft.diff.abs" ] <- mean(abs(fcomp.diff))
      pft.out[i,"pft.deriv.abs"] <- mean(abs(fcomp.out[[i]]$mod.deriv$mean))
    }  
  }
  
  # Add what we just calculated to everything we've already done
  stab.fcomp <- rbind(stab.fcomp, pft.out)
}
# ------


# -------------------
stab.fcomp <- stab.fcomp[!toupper(stab.fcomp$PFT)=="TOTAL",]
summary(stab.fcomp)

write.csv(stab.fcomp, file.path(path.google, "Current Data/Stability", "Stability_Models_FCOMP_100.csv"), row.names=F)
# -------------------------------------------

