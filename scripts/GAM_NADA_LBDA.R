# -------------------------------------------
# Assessing the scales and magnitude of climate change variability over the past millennium in the NADA
# Author: Christy Rollinson, crollinson@gmail.com

# 1. Stastical detection of significant change
#    - Use loess or TP regression splines
# 2. Comparison with paleoclimate reconstructions
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# 0. Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra); library(scales); library(grid)
library(mgcv)
library(plyr); library(parallel)

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"

# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
path.gamm.func <- "~/Desktop/Research/R_Functions/"

# Path to github repository/working directory
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

# Path to where the raw output is
path.data <- "~/Desktop/Research/PalEON_MIP_Region/NADA/"

# Lets just save processed to the Google Drive folder
path.out <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/Current Data/Stability_Index/"
path.fig <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/Current Figures/Stability_Index/"

# Set up some time variables just to help with indexing
yrs <- 850:2010
mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)

us <- map_data("state")
# -------------------------------------------


# --------------------------------------------
# 1. Load & Extract PDSI data from NADA & LBDA
# 
# All files downloaded from here on 9 June 2017
# https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring/north-american-drought-variability
# --------------------------------------------
library(ncdf4)

# Load the Paleon site data to get a bounding box
paleon <- read.csv(file.path(path.repo, "data/paleon_models_environment_master.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
summary(paleon)

# This is the v2 NADA at 2.5 degree resolution
nada.nc <- nc_open(file.path(path.data, "NADAv2-2008.nc"))

# This is the living blended drought atlas PDSI
lbda.nc <- nc_open(file.path(path.data, "nada_hd2_cl.nc"))


nada.lat <- ncvar_get(nada.nc, "lat")
nada.lon <- ncvar_get(nada.nc, "lon")
nada.time <- ncvar_get(nada.nc, "time")
nada.nc$dim$time$units # Time goes present to past; starting in 2006-0

lbda.lat <- ncvar_get(lbda.nc, "lat")
lbda.lon <- ncvar_get(lbda.nc, "lon")
lbda.time <- ncvar_get(lbda.nc, "time")
lbda.nc$dim$time$units # Goes from past to present: 0-2005

nada.res <- mean(diff(nada.lon))
nada.lon.ind <- which(nada.lon+nada.res/2>=min(paleon$lon) & nada.lon-nada.res/2<=max(paleon$lon))
nada.lat.ind <- which(nada.lat+nada.res/2>=min(paleon$lat) & nada.lat-nada.res/2<=max(paleon$lat))
nada.time.ind <- which(nada.time>=850 & nada.time<=1850)
nada.lon  <- nada.lon [nada.lon.ind]
nada.lat  <- nada.lat [nada.lat.ind]
nada.time <- nada.time[nada.time.ind]


lbda.res <- mean(diff(lbda.lon))
lbda.lon.ind <- which(lbda.lon+lbda.res/2>=min(paleon$lon) & lbda.lon-lbda.res/2<=max(paleon$lon))
lbda.lat.ind <- which(lbda.lat+lbda.res/2>=min(paleon$lat) & lbda.lat-lbda.res/2<=max(paleon$lat))
lbda.time.ind <- which(lbda.time>=850 & lbda.time<=1850)
lbda.lon  <- lbda.lon [lbda.lon.ind]
lbda.lat  <- lbda.lat [lbda.lat.ind]
lbda.time <- lbda.time[lbda.time.ind]

nada.raw <- ncvar_get(nada.nc, "PDSI")[nada.lon.ind, nada.lat.ind, nada.time.ind]
lbda.raw <- ncvar_get(lbda.nc, "pdsi")[lbda.time.ind, lbda.lat.ind, lbda.lon.ind]
# lbda.raw <- ncvar_get(lbda.nc, "pdsi")
dim(nada.raw)
dim(lbda.raw)

# Don't forget to close the netcdf files!
nc_close(nada.nc); nc_close(lbda.nc)

# Re-arrange LBDS to be lon, lat, time
lbda.raw <- aperm(lbda.raw, c(3,2,1))

# Arrange everything to go from least to most in all dimensions
nada.raw <- nada.raw[order(nada.lon), order(nada.lat), order(nada.time)]
lbda.raw <- lbda.raw[order(lbda.lon), order(lbda.lat), order(lbda.time)]

# Getting our dimension indices right now
nada.lon <- nada.lon[order(nada.lon)]
nada.lat <- nada.lat[order(nada.lat)]
nada.time <- nada.time[order(nada.time)]
lbda.lon <- lbda.lon[order(lbda.lon)]
lbda.lat <- lbda.lat[order(lbda.lat)]
lbda.time <- lbda.time[order(lbda.time)]



# There are some weird values that aren't real, so lets get rid of them
lbda.raw[lbda.raw<=-99] <- NA

# Getting the 100-year averages that are comparable to STEPPS & REFAB
yrs.cent <- seq(900, 1800, by=100)
# --------------------------------------------


# --------------------------------------------
# 2. Run a GAM to look at stability for all grid cells
#    -  Call my GAM helper functions
#    2.0. Set up the GAM & source the helper functions
#    2.1. NADA
#    2.2. LBDA
# --------------------------------------------

# ------------------
# 2.0. Set up the GAM & source the helper functions
# ------------------
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))

calc.stability <- function(x, width=25){
  dat.tmp <- data.frame(Y=x, Year=1:length(x))
  k.use=round(length(x[!is.na(x)])/width, 0)
  mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
  mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
  return(mod.deriv)
}
# ------------------

# ------------------
# 2.1. Calculate NADA stability
# ------------------
nada.df <- data.frame(lon=rep(nada.lon), lat=rep(nada.lat, each=length(nada.lon)))
nada.deriv <- array(dim=dim(nada.raw))
nada.sig   <- array(dim=dim(nada.raw))

pb <- txtProgressBar(min=0, max=dim(nada.raw)[1] * dim(nada.raw)[2], style=3)
pb.ind=1
for(i in 1:dim(nada.raw)[1]){
  for(j in 1:dim(nada.raw)[2]){
    row.ind <- which(nada.df$lon==nada.lon[i] & nada.df$lat==nada.lat[j])
    
    setTxtProgressBar(pb, pb.ind)
    pb.ind=pb.ind+1

    if(all(is.na(nada.raw[i,j,]))) next
    
    stab.tmp <- calc.stability(nada.raw[i,j,], width=25)
    stab.tmp[is.na(stab.tmp$Y), c("mean", "lwr", "upr", "sig")] <- NA

    nada.deriv[i,j,] <- stab.tmp$mean
    nada.sig  [i,j,] <- stab.tmp$sig
    
    # Store some summaryinfo
    nada.df[row.ind, "deriv.mean"] <- mean(stab.tmp$mean, na.rm=T)
    nada.df[row.ind, "deriv.abs"] <- mean(abs(stab.tmp$mean), na.rm=T)
    nada.df[row.ind, "n.yrs"] <- length(which(!is.na(stab.tmp$Y)))
    nada.df[row.ind, "n.yrs.sig"] <- length(which(!is.na(stab.tmp$sig)))
    
  }
}
nada.df$fract.sig <- nada.df$n.yrs.sig/nada.df$n.yrs 
write.csv(nada.df, file.path(path.google, "Current Data/Stability_GAMs", "Stability_NADA.csv"), row.names=F)
summary(nada.df)

nada.deriv <- ggplot(data=nada.df) +
  geom_tile(aes(x=lon, y=lat, fill=deriv.abs)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(nada.df$lon), ylim=range(nada.df$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(nada.df$deriv.abs, na.rm=T)) +
  theme_bw() +
  ggtitle("Mean Absolute Rate of Change")

nada.sig <- ggplot(data=nada.df) +
  geom_tile(aes(x=lon, y=lat, fill=n.yrs.sig)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(nada.df$lon), ylim=range(nada.df$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(nada.df$n.yrs.sig, na.rm=T)) +
  theme_bw() +
  ggtitle("Number of Years Showing Change")

png(file.path(path.google, "Current Figures/Stability_GAMs", "Stability_NADA.png"), height=4, width=6, units="in", res=320)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(nada.deriv, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(nada.sig, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

# ------------------

# ------------------
# 2.2. Calculate LBDA stability
# ------------------
lbda.df <- data.frame(lon=rep(lbda.lon), lat=rep(lbda.lat, each=length(lbda.lon)))
lbda.deriv <- array(dim=dim(lbda.raw))
lbda.sig   <- array(dim=dim(lbda.raw))

pb <- txtProgressBar(min=0, max=dim(lbda.raw)[1] * dim(lbda.raw)[2], style=3)
pb.ind=1
for(i in 1:dim(lbda.raw)[1]){
  for(j in 1:dim(lbda.raw)[2]){
    row.ind <- which(lbda.df$lon==lbda.lon[i] & lbda.df$lat==lbda.lat[j])
    
    setTxtProgressBar(pb, pb.ind)
    pb.ind=pb.ind+1
    
    if(all(is.na(lbda.raw[i,j,]))) next
    
    stab.tmp <- calc.stability(lbda.raw[i,j,], width=25)
    stab.tmp[is.na(stab.tmp$Y), c("mean", "lwr", "upr", "sig")] <- NA
    
    lbda.deriv[i,j,] <- stab.tmp$mean
    lbda.sig  [i,j,] <- stab.tmp$sig
    
    # Store some summaryinfo
    lbda.df[row.ind, "deriv.mean"] <- mean(stab.tmp$mean, na.rm=T)
    lbda.df[row.ind, "deriv.abs"] <- mean(abs(stab.tmp$mean), na.rm=T)
    lbda.df[row.ind, "n.yrs"] <- length(which(!is.na(stab.tmp$Y)))
    lbda.df[row.ind, "n.yrs.sig"] <- length(which(!is.na(stab.tmp$sig)))
    
  }
}
lbda.df$fract.sig <- lbda.df$n.yrs.sig/lbda.df$n.yrs 
summary(lbda.df)
write.csv(lbda.df, file.path(path.google, "Current Data/Stability_GAMs", "Stability_LBDA.csv"), row.names=F)

lbda.deriv <- ggplot(data=lbda.df) +
  geom_tile(aes(x=lon, y=lat, fill=deriv.abs)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(lbda.df$lon), ylim=range(lbda.df$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(lbda.df$deriv.abs, na.rm=T)) +
  theme_bw() +
  ggtitle("Mean Absolute Rate of Change")

lbda.sig <- ggplot(data=lbda.df) +
  geom_tile(aes(x=lon, y=lat, fill=n.yrs.sig)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(lbda.df$lon), ylim=range(lbda.df$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(lbda.df$n.yrs.sig, na.rm=T)) +
  theme_bw() +
  ggtitle("Number of Years Showing Change")


# lbda.frac <- ggplot(data=lbda.df) +
#   geom_tile(aes(x=lon, y=lat, fill=fract.sig)) +
#   geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
#   coord_equal(xlim=range(lbda.df$lon), ylim=range(lbda.df$lat), expand=0) +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(lbda.df$fract.sig, na.rm=T)) +
#   theme_bw() +
#   ggtitle("Fraction of Years Showing Change")


png(file.path(path.google, "Current Figures/Stability_GAMs", "Stability_LBDA.png"), height=4, width=6, units="in", res=320)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(lbda.deriv, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(lbda.sig, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
# print(lbda.frac, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
dev.off()

# ------------------

# --------------------------------------------
