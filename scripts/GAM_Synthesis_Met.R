# -------------------------------------------
# Comparing stability of annually resolved climate in models and Data
# Author: Christy Rollinson, crollinson@gmail.com
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra); library(scales); library(grid)
library(mgcv)
library(plyr); library(parallel)

# Setting the path to this repository
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

# Path to where the raw output is
path.models <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"

# Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
path.gamm.func <- "~/Desktop/Research/R_Functions/"
# Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" 

# -------------------------------------------

# -------------------------------------------
# -------------------------------------------
stab.met  <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers.csv"))
stab.lbda <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_LBDA.csv"))
summary(stab.lbda)
summary(stab.met)

# Standardizing the derivatives relative to the mean

# -------------------------------------------
