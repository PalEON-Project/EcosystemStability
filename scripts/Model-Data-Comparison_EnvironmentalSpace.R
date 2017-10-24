# -------------------------------------------
# Model-Data comparisons: Comparisons of Model & Data Stability in Climate & Environmental Space

# Here we will:
# 2. Look at patterns of stability in environmental space
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# Load Libaries & Set file paths
# -------------------------------------------
library(ggplot2); library(gridExtra); library(scales)
library(sp); library(raster)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"
# -------------------------------------------

# -------------------------------------------
# Read in & align predictor layers
# -------------------------------------------
# -----------
# Empirical Drought Reconstruction (LBDA, annually-resolved)
# -----------
lbda <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_LBDA_100.csv"))
coordinates(lbda) <- lbda[,c("lon", "lat")]

# Set up a blank raster
lbda.rast <- raster(ncol=length(unique(lbda$lon)), nrow=length(unique(lbda$lat)))
extent(lbda.rast) <- extent(lbda)
projection(lbda.rast) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Rasterize the layer we want
lbda.deriv <- rasterize(lbda, lbda.rast, lbda$deriv.abs, fun=mean)
# projection(lbda.deriv) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
names(lbda.deriv) <- "deriv.pdsi"
lbda.deriv

# -----------

# -----------
# Model Drivers
# -----------
drivers.all <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
coordinates(drivers.all) <- drivers.all[,c("lon", "lat")]
drivers.all$whc.tot[drivers.all$whc.tot>1e3] <- NA
summary(drivers.all)

# Set up a blank raster
domain <- raster(ncol=length(unique(drivers.all$lon)), nrow=length(unique(drivers.all$lat)))
extent(domain) <- extent(drivers.all)
projection(domain) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Rasterize the layer we want
drivers <- rasterize(drivers.all, domain, drivers.all$depth, fun=mean)
# projection(drivers) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
names(drivers) <- "soil.depth"
drivers <- stack(drivers)
drivers

drivers$whc <- rasterize(drivers.all, domain, drivers.all$whc.tot, fun=mean)
drivers$tair.sett <- rasterize(drivers.all, domain, drivers.all$tair.yr.set, fun=mean)
drivers$precip.sett <- rasterize(drivers.all, domain, drivers.all$precip.yr.set, fun=mean)
drivers$deriv.pdsi <- rasterize(drivers.all, domain, drivers.all$pdsi.deriv, fun=mean)
drivers$deriv.tair <- rasterize(drivers.all, domain, drivers.all$tair.deriv, fun=mean)
drivers$deriv.precip <- rasterize(drivers.all, domain, drivers.all$precip.deriv, fun=mean)
drivers
plot(drivers)
# -----------
# -------------------------------------------

# -------------------------------------------
# Read in, format, & extract the point data
# -------------------------------------------
# -----------
# STEPPS (empirical composition, centennially-resolved)
# -----------
stepps <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_STEPPS.csv"))
coordinates(stepps) <- stepps[,c("lon", "lat")]
projection(stepps) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(stepps)

plot(drivers$whc)
plot(stepps, add=T, col="black", cex=0.2)
# -----------

# -----------
# ReFAB (empirical biomass, centennially-resolved)
# # NOTE: Derivative needs to be the mean ABSOLUTE VALUE!  NEEDS TO BE RE-RUN!!
# -----------
refab <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_ReFAB.csv"))
refab$deriv.abs <- abs(refab$refab.mean.slope) # NOTE: SHOULDN"T HAVE TO DO THIS! FIX IN FUTURE VERSIONS
names(refab)[which(names(refab)=="long")] <-"lon"
coordinates(refab) <- refab[,c("lon", "lat")]
projection(refab) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(refab)
# -----------

# -----------
# Model Output
# -----------
models1 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))
coordinates(models1) <- models1[,c("lon", "lat")]
projection(models1) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(models1)
# -----------

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Spatial_Extent.png"), height=5, width=10, units="in", res=320)
plot(drivers$tair.sett-273.15)
plot(models1, add=T, col="gray50", pch=20, cex=0.5)
plot(stepps, add=T, col="black", cex=0.2)
plot(refab, add=T, col="dodgerblue", pch=20, cex=1)
dev.off()
# -------------------------------------------
