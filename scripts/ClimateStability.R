# --------------------------------------------
# Calculating stability indices for climate
# Stabilty Index = sum of squares of second derivative
# Climate = Temperature, Precipitation
# Script Author: Christy Rollinson (crollinson@mortonarb.org)

# Workflow
# 0. Define file paths etc
# 1. Read in & format climate data (temp, precip)
#    1.a. Annual Means
#    1.b. Summer Mean (June, July, August)
#    1.c. 100-year averages of annual & summer temp/precip
# 2. Calculate Stability Index for all grid cells
#    -  Call calc.second.deriv (calc_second_deriv.R)
# 3. Save output, generate & save a couple figures
# --------------------------------------------

# --------------------------------------------
# 0. Define file paths etc
# --------------------------------------------
# Path to github repository/working directory
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Set up some time variables just to help with indexing
yrs <- 850:2010
mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)
# --------------------------------------------

# --------------------------------------------
# 1. Read in & format climate data (temp, precip)
#    1.a. Annual Means
#    1.b. Summer Mean (June, July, August)
# 
# Notes: Since we're only interested in climate at places we have model 
#        runs, I'm just going to pull T & P from LPJ-GUESS (assuming that 
#        those are fine, which should be a safe 
#        screw ED up and it's running with the right met... I think this is
#        a safe assumption)
# --------------------------------------------
# load in the paleon domain info;
# This got generated using domain_environment_extraction.R
paleon <- read.csv(file.path(path.repo, "data/paleon_models_environment_master.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
summary(paleon)

# Load the temp & precip data
tair    <- readRDS(file.path(path.data, "Met/tair_all.rds"))
precipf <- readRDS(file.path(path.data, "Met/precipf_all.rds"))
pdsi <- readRDS(file.path(path.data, "Met/pdsi_calib_1931-1990_all.rds"))

# aggregate to annual resolution
pdsi.ann <- tair.ann <- tair.jja <- precip.ann <- precip.jja <- matrix(ncol=ncol(tair), nrow=length(yrs))

for(i in 1:length(yrs)){
  # Generating indices for the cells we want to aggregate across
  rows.yrs <- which(time.mos$year==yrs[i])
  rows.jja <- which(time.mos$year==yrs[i] & time.mos$month %in% c(6:8))
  
  # doing the aggregation
  pdsi.ann  [i,] <- colMeans(pdsi   [rows.yrs,])
  tair.ann  [i,] <- colMeans(tair   [rows.yrs,])
  precip.ann[i,] <- colMeans(precipf[rows.yrs,])
  tair.jja  [i,] <- colMeans(tair   [rows.jja,])
  precip.jja[i,] <- colMeans(precipf[rows.jja,])
}

# Getting the 100-year averages that are comparable to STEPPS & REFAB
yrs.cent <- seq(900, 1800, by=100)
pdsi.ann2 <- tair.ann2 <- tair.jja2 <- precip.ann2 <- precip.jja2 <- matrix(ncol=ncol(tair), nrow=length(yrs.cent))
for(i in 1:length(yrs.cent)){
  # Generating indices for the cells we want to aggregate across
  rows.yrs <- which(yrs>=yrs.cent[i]-50 & yrs<=yrs.cent[i]+50)

  # doing the aggregation
  pdsi.ann2  [i,] <- colMeans(pdsi.ann  [rows.yrs,])
  tair.ann2  [i,] <- colMeans(tair.ann  [rows.yrs,])
  precip.ann2[i,] <- colMeans(precip.ann[rows.yrs,])
  tair.jja2  [i,] <- colMeans(tair.jja  [rows.yrs,])
  precip.jja2[i,] <- colMeans(precip.jja[rows.yrs,])
}

summary(tair.ann[,1:10])
summary(tair.ann2[,1:10])
# --------------------------------------------


# --------------------------------------------
# 2. Calculate Stability Index for all grid cells
#    -  Call calc.second.deriv (calc_second_deriv.R)
# --------------------------------------------
source(file.path(path.repo, "scripts/calc_second_deriv.R"))

# Doing the actual calculation on each cell
paleon$stab.pdsi.ann      <- apply(pdsi.ann[which(yrs<=1850),] , 2, calc.second.deriv, h=1, H=1)
paleon$stab.pdsi.ann.cent <- apply(pdsi.ann2, 2, calc.second.deriv, h=1, H=100)

paleon$stab.tair.ann      <- apply(tair.ann[which(yrs<=1850),] , 2, calc.second.deriv, h=1, H=1)
paleon$stab.tair.ann.cent <- apply(tair.ann2, 2, calc.second.deriv, h=1, H=100)
paleon$stab.tair.jja      <- apply(tair.jja[which(yrs<=1850),] , 2, calc.second.deriv, h=1, H=1)
paleon$stab.tair.jja.cent <- apply(tair.jja2, 2, calc.second.deriv, h=1, H=100)

paleon$stab.precip.ann      <- apply(precip.ann[which(yrs<=1850),] , 2, calc.second.deriv, h=1, H=1)
paleon$stab.precip.ann.cent <- apply(precip.ann2, 2, calc.second.deriv, h=1, H=100)
paleon$stab.precip.jja      <- apply(precip.jja[which(yrs<=1850),] , 2, calc.second.deriv, h=1, H=1)
paleon$stab.precip.jja.cent <- apply(precip.jja2, 2, calc.second.deriv, h=1, H=100)
# --------------------------------------------


# --------------------------------------------
# 3. Save output, generate & save a couple figures
# --------------------------------------------
summary(paleon)

# Save the file
write.csv(paleon, file.path(path.repo, "data/PalEON_ClimateStability.csv"), row.names=F, eol="\r\n")


# Do some graphing
library(ggplot2)

summary(paleon)
# Stacking things together
paleon2 <- stack(paleon[,26:ncol(paleon)])
paleon2[,c("lon", "lat", "latlon", "umw", "domain.paleon", "x", "y")] <- paleon[,c("lon", "lat", "latlon", "umw", "domain.paleon", "x", "y")]
summary(paleon2)

for(i in 1:nrow(paleon2)){
  paleon2[i, "var"       ] <- strsplit(paste(paleon2$ind[i]), "[.]")[[1]][2] 
  paleon2[i, "season"    ] <- strsplit(paste(paleon2$ind[i]), "[.]")[[1]][3] 
  paleon2[i, "resolution"] <- ifelse(is.na(strsplit(paste(paleon2$ind[i]), "[.]")[[1]][4]), "annual", "century")
}
paleon2$var        <- as.factor(paleon2$var)
paleon2$season     <- as.factor(paleon2$season)
paleon2$resolution <- as.factor(paleon2$resolution)
summary(paleon2)

us <- map_data("state")
for(v in unique(paleon2$var)){
  pdf(file.path(path.repo, "figures", paste0("Stability_Met_Region_", v, ".pdf")))
  print(
    ggplot(data=paleon2[paleon2$var==v,]) +
      facet_grid(season~resolution, scales="free_x") +
      geom_histogram(aes(values)) + 
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="annual",]) +
      ggtitle("Annual Resolution") +
      facet_grid(season~.) +
      geom_tile(aes(x=lon, y=lat, fill=values)) + 
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      coord_equal(xlim=range(paleon2$lon), ylim=range(paleon2$lat)) +
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="century",]) +
      facet_grid(season~.) +
      geom_tile(aes(x=lon, y=lat, fill=values)) + 
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      coord_equal(xlim=range(paleon2$lon), ylim=range(paleon2$lat)) +
      ggtitle("Centennial Resolution") +
      theme_bw()
  )
  dev.off()
}

for(v in unique(paleon2$var)){
  pdf(file.path(path.repo, "figures", paste0("Stability_Met_Models_", v, ".pdf")))
  print(
    ggplot(data=paleon2[paleon2$var==v & !is.na(paleon2$umw),]) +
      facet_grid(season~resolution, scales="free_x") +
      geom_histogram(aes(values)) + 
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="annual",]) +
      ggtitle("Annual Resolution") +
      facet_grid(season~.) +
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      geom_point(aes(x=lon, y=lat, color=values)) + 
      coord_equal(xlim=range(paleon2$lon), ylim=range(paleon2$lat)) +
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="century",]) +
      facet_grid(season~.) +
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      geom_point(aes(x=lon, y=lat, color=values)) + 
      coord_equal(xlim=range(paleon2$lon), ylim=range(paleon2$lat)) +
      ggtitle("Centennial Resolution") +
      theme_bw()
  )
  dev.off()
}

for(v in unique(paleon2$var)){
  pdf(file.path(path.repo, "figures", paste0("Stability_Met_UMW_", v, ".pdf")))
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$umw=="y" & !is.na(paleon2$umw),]) +
      facet_grid(season~resolution, scales="free_x") +
      geom_histogram(aes(values)) + 
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="annual",]) +
      facet_grid(season~.) +
      geom_tile(aes(x=lon, y=lat, fill=values)) + 
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      coord_equal(xlim=range(paleon2[paleon2$umw=="y", "lon"], na.rm=T), ylim=range(paleon2[paleon2$umw=="y", "lat"], na.rm=T)) +
      ggtitle("Annual Resolution") +
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="century",]) +
      facet_grid(season~.) +
      geom_tile(aes(x=lon, y=lat, fill=values), size=2) + 
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      coord_equal(xlim=range(paleon2[paleon2$umw=="y", "lon"], na.rm=T), ylim=range(paleon2[paleon2$umw=="y", "lat"], na.rm=T)) +
      ggtitle("Centennial Resolution") +
      theme_bw()
  )
  dev.off()
}

for(v in unique(paleon2$var)){
  pdf(file.path(path.repo, "figures", paste0("Stability_Met_Models_UMW_", v, ".pdf")))
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$umw=="y" & !is.na(paleon2$umw),]) +
      facet_grid(season~resolution, scales="free_x") +
      geom_histogram(aes(values)) + 
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="annual" & paleon2$umw=="y" & !is.na(paleon2$umw),]) +
      facet_grid(season~.) +
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      geom_point(aes(x=lon, y=lat, color=values), size=2) + 
      coord_equal(xlim=range(paleon2[paleon2$umw=="y", "lon"], na.rm=T), ylim=range(paleon2[paleon2$umw=="y", "lat"], na.rm=T)) +
      ggtitle("Annual Resolution") +
      theme_bw()
  )
  print(
    ggplot(data=paleon2[paleon2$var==v & paleon2$resolution=="century" & paleon2$umw=="y" & !is.na(paleon2$umw),]) +
      facet_grid(season~.) +
      geom_path(data=us,aes(x=long, y=lat, group=group)) + 
      geom_point(aes(x=lon, y=lat, color=values), size=2) + 
      coord_equal(xlim=range(paleon2[paleon2$umw=="y", "lon"], na.rm=T), ylim=range(paleon2[paleon2$umw=="y", "lat"], na.rm=T)) +
      ggtitle("Centennial Resolution") +
      theme_bw()
  )
  dev.off()
}

# --------------------------------------------
