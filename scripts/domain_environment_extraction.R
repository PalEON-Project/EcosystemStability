# --------------------------------------------
# Script Description
# --------------------------------------------
# Purpose: Extract all of the abiotic info for the paleon modeling grid cells and
#          add info as to what subdomains each grid cell is a part of
#
# Analysis Details:
# - Environment (Driver) Variables: 
#     - Water Holding Capactity (calculated from depth, % sand, % clay)
#        - could break down and run these 3 factors independently rather than combining
#     - Temperature: Mean annual, summer (June-July-August) 
#     - Precipitation: Mean annual, summer (JJA)
# - Spatial Extent:
#    - All sites in ecosystem model runs.
#    * maybe also run with full spatial coverage (includes sites many models didn't do)
# - Temporal Extent (for climate):
#    - Settlement Era: 1800-1850 mean
#    - Empirical Climate: 1901-1930 mean
#
#
# Workflow
# 0. Define file paths etc
# 1. Read in & format environemnt data (temp, precip, WHC)
#    - Annual/Summer Mean, 1800-1850/1901-1930
# 2. Perform Oridination (PCA)
#    - Extract loadings & PCs
# 3. Save output, 
#    - generate & save a couple figures
# --------------------------------------------

# --------------------------------------------
# 0. Define file paths etc
# --------------------------------------------
# Path to github repository/working directory
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

# Path to where data are
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"
path.met  <- "~/Desktop/Research/PalEON_MIP_Region/phase2_met_regional_v2_monthly/"
path.soil <- "~/Dropbox/PalEON_CR/env_regional/phase2_env_drivers_v2/soil/"

# This is the path to my pdsi code is, which has the formula for calculating WHC
# This is currently part of my met ensemble code and code 
# and can be found here: https://github.com/PalEON-Project/modeling_met_ensemble.git
path.pdsi <- "~/Dropbox/PalEON_CR/met_ensemble/scripts/"

# Lets just save processed to the Google Drive folder
path.out <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data"

# Set up some time variables just to help with indexing
yrs <- 850:2010
mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)
# --------------------------------------------


# --------------------------------------------
# Reading in & extracting data
# --------------------------------------------
# loading the key to dimensions in the .RDS files; this is Christy's sheet generated
# this loads an object called "paleon"
paleon <- read.csv(file.path(path.repo, "data/paleon_site_info.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
paleon <- paleon[,c("lon", "lat", "latlon", "umw", "x", "y")]
summary(paleon)

# Extracting soil data for the appropriate sites 
# Note: We need this for the "upper" and "lower" soil layers
library(ncdf4)
tair.nc <- nc_open(file.path(path.met, "tair.nc"))
precipf.nc <- nc_open(file.path(path.met, "precipf.nc"))
sand.t <- nc_open(file.path(path.soil, "paleon_soil_t_sand.nc"))
sand.s <- nc_open(file.path(path.soil, "paleon_soil_s_sand.nc"))
clay.t <- nc_open(file.path(path.soil, "paleon_soil_t_clay.nc"))
clay.s <- nc_open(file.path(path.soil, "paleon_soil_s_clay.nc"))
depth  <- nc_open(file.path(path.soil, "paleon_soil_soil_depth.nc"))

# Getting the lat/lon index for each point to extract soil & met
# Note: for some reason the met & the soil are not in the same order
lon <- ncvar_get(sand.t, "longitude")
lat <- ncvar_get(sand.t, "latitude")
lon2 <- ncvar_get(tair.nc, "lon")
lat2 <- ncvar_get(tair.nc, "lat")


# Creating a dataframe to plug all sites into
paleon.all <- data.frame(lon=rep(lon, length=length(lon)*length(lat)), lat=rep(lat, each=length(lon)))
paleon.all$latlon <- as.factor(paste0("lat", paleon.all$lat, "lon", paleon.all$lon))
summary(paleon.all)


# Merging the subset and the full region together; we can always separate this back out later
paleon.all <- merge(paleon.all, paleon, all=T)
summary(paleon.all)

# Determining where we should be pulling to get the settlement-era & empirical climates
set.all <- which(time.mos$year>=1800 & time.mos$year<=1850)
set.jja <- which(set.all %in% which(time.mos$year>=1800 & time.mos$year<=1850 & time.mos$month>=6 & time.mos$month<=8))
cru.all <- which(time.mos$year>=1901 & time.mos$year<=1930)
cru.jja  <- which(cru.all %in% which(time.mos$year>=1901 & time.mos$year<=1930 & time.mos$month>=6 & time.mos$month<=8))

# Lets convert precip (kg/m2/s = mm/s) to total mm per yr or summer
library(lubridate)
sec_2_yr <- 60*60*24*365.25
sec_2_smr <- 60*60*24*sum(days_in_month(6:8))

# This is the ugly way of doing it; could do it prettier since we're doing the whole region, but whatever
for(i in 1:nrow(paleon.all)){
  # Extracting soil info
  x.ind  <- which(lon == paleon.all[i,"lon"])
  y.ind  <- which(lat == paleon.all[i,"lat"])
  paleon.all[i,"sand.t"] <- ncvar_get(sand.t, "t_sand", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"sand.s"] <- ncvar_get(sand.s, "s_sand", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"clay.t"] <- ncvar_get(clay.t, "t_clay", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"clay.s"] <- ncvar_get(clay.s, "s_clay", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"depth"]  <- ncvar_get(depth , "soil_depth", c(x.ind, y.ind), c(1,1)) # cm
  
  # Extracting met data
  # Convert to 
  x.ind2 <- which(lon2 == paleon.all[i,"lon"])
  y.ind2 <- which(lat2 == paleon.all[i,"lat"])
  
  tair.set   <- ncvar_get(tair.nc   , "tair"   , c(x.ind2, y.ind2, min(set.all)), c(1,1,length(set.all)))
  precip.set <- ncvar_get(precipf.nc, "precipf", c(x.ind2, y.ind2, min(set.all)), c(1,1,length(set.all)))
  paleon.all[i,"tair.yr.set"   ] <- mean(tair.set)
  paleon.all[i,"tair.jja.set"  ] <- mean(tair.set[set.jja])
  paleon.all[i,"precip.yr.set" ] <- mean(precip.set)
  paleon.all[i,"precip.jja.set"] <- mean(precip.set[set.jja])
  
  tair.cru   <- ncvar_get(tair.nc   , "tair"   , c(x.ind2, y.ind2, min(cru.all)), c(1,1,length(cru.all)))
  precip.cru <- ncvar_get(precipf.nc, "precipf", c(x.ind2, y.ind2, min(cru.all)), c(1,1,length(cru.all)))
  paleon.all[i,"tair.yr.cru"   ] <- mean(tair.cru)
  paleon.all[i,"tair.jja.cru"  ] <- mean(tair.cru[cru.jja])
  paleon.all[i,"precip.yr.cru" ] <- mean(precip.cru)
  paleon.all[i,"precip.jja.cru"] <- mean(precip.cru[cru.jja])
}

# Get rid of the cells that are water
paleon.all <- paleon.all[!is.na(paleon.all$sand.t),]
summary(paleon.all)

# To convert precip (kg/m2/s = mm/s) to total mm per yr or summer
library(lubridate)
sec_2_yr <- 60*60*24*365.25
sec_2_smr <- 60*60*24*sum(days_in_month(6:8))



# Calculating Water Holding Capacity for each site from texture
source(file.path(path.pdsi, "calc.awc.R"))
paleon.all$awc.t <- calc.awc(paleon.all$sand.t, paleon.all$clay.t)
paleon.all$awc.s <- calc.awc(paleon.all$sand.s, paleon.all$clay.s)

# Calculating water hodling capacity **in cm**
# If our soil depth is greater than 30 cm, we assume a 30 cm topsoil; if our soil depth is less than 
# 30 cm, we assume all but 1 cm of that is topsoil (yeah, that's not great, but it's what I came up with)
paleon.all$whc.t   <- ifelse(paleon.all$depth > 30, paleon.all$awc.t * 30, paleon.all$awc.t * paleon.all$depth-1) # 30 cm depth
paleon.all$whc.s   <- ifelse(paleon.all$depth > 30, paleon.all$awc.s * (paleon.all$depth-30), paleon.all$awc.s*1) # depth - 30 cm depth 
paleon.all$whc.tot <- paleon.all$whc.t + paleon.all$whc.s
summary(paleon.all)


# Mapping our environment 
{
  library(ggplot2); library(grid)
  plot.clay.t   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=clay.t)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "coral2"), name="Clay 1")
  plot.clay.s   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=clay.s)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "coral2"), name="Clay 2")
  plot.sand.t   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=sand.t)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "goldenrod2"), name="Sand 1")
  plot.sand.s   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=sand.s)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "goldenrod2"), name="Sand 2")
  plot.wcap.t   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=whc.t)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "cyan3"), name="WHC 1")
  plot.wcap.s   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=whc.s)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "cyan3"), name="WHC 2", limits=c(0,1000))
  plot.depth   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=depth)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "tan3"), name="Depth")
  plot.whc     <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=whc.tot)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue"), name="Total WHC", limits=c(0,1000))
  
  tair.yr.set   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.yr.set-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red3"), name="Tair.yr")
  tair.gs.set   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.jja.set-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red2"), name="Tair.jja")
  precip.yr.set <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.yr.set*sec_2_yr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.yr")
  precip.gs.set <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.jja.set*sec_2_smr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.jja")
  
  tair.yr.cru   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.yr.cru-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red3"), name="Tair.yr")
  tair.gs.cru   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.jja.cru-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red2"), name="Tair.jja")
  precip.yr.cru <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.yr.cru*sec_2_yr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.yr")
  precip.gs.cru <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.jja.cru*sec_2_smr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.jja")
  
  
  # ---------------
  # Some exploratory graphs to see if we'd expect many differences; 
  # in general, summer shows stronger spatial patterns
  # ---------------
  pdf(file.path(path.out, "Current Figures/Data_Raw", "Environmental_Drivers.pdf"))
  # The non-climate environment
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(4, 2)))
  print(plot.clay.t, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot.clay.s, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(plot.sand.t, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot.sand.s, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(plot.wcap.t, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(plot.wcap.s, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
  print(plot.depth , vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
  print(plot.whc   , vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
  
  
  # Compare year vs. summer: settlement era
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(tair.yr.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(tair.gs.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(precip.yr.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(precip.gs.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  # Compare year vs. summer: CRUNCEP
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(tair.yr.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(tair.gs.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(precip.yr.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(precip.gs.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  # Compare Settlement vs CRU: Year
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(tair.yr.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(tair.yr.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(precip.yr.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(precip.yr.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  
  # Compare Settlement vs CRU: Summer
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))
  print(tair.gs.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(tair.gs.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(precip.gs.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(precip.gs.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  dev.off()
  # ---------------
}
# --------------------------------------------

# --------------------------------------------
# Adding some more definitions of regions
# --------------------------------------------
library(rgdal); library(raster)
# Loading the Paleon 1 SetVeg Domain:
domain.paleon <- raster("data/Paleon_Domain_Masks/paleon_full_ll_v0.1.tif")
plot(domain.paleon)

models.sp <- SpatialPointsDataFrame(coords=c(paleon.all[,c("lon", "lat")]), data = paleon.all)

paleon.all$domain.paleon <- as.factor(extract(domain.paleon, models.sp))

# Recoding the numeric domains to the factors listed here: 
# https://paleon.geography.wisc.edu/doku.php/data_and_products;rasters
library(car)
paleon.all$domain.paleon <- recode(paleon.all$domain.paleon, "'1'='New England'; '2'='IL'; '3'='IN'; '4'='ME'; '5'='MI - lower'; '6'='MN'; '7'='NJ'; '8'='NY'; '9'='OH'; '10'='PA'; '11'='WI'; '12'='MI - upper'")
summary(paleon.all)

write.csv(paleon.all, file.path(path.out, "Current Data", "paleon_models_environment_master.csv"), row.names=F)
# --------------------------------------------
