# --------------------------------------------
# Extracting NADA PDSI for our region and calculating the stability index
# Stabilty Index = sum of squares of second derivative
# Climate = Temperature, Precipitation
# Script Author: Christy Rollinson (crollinson@mortonarb.org)

# Workflow
# 0. Define file paths etc
# 1. Load & Extract PDSI data from NADA
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
path.data <- "~/Desktop/Research/PalEON_MIP_Region/NADA/"

# Set up some time variables just to help with indexing
yrs <- 850:2010
# mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)
# --------------------------------------------


# --------------------------------------------
# 1. Load & Extract PDSI data from NADA
# 
# All files downloaded from here on 9 June 2017
# https://www.ncdc.noaa.gov/data-access/paleoclimatology-data/datasets/tree-ring/north-american-drought-variability
# --------------------------------------------
library(ncdf4)

# Load the Paleon site data to get a bounding box
paleon <- read.csv(file.path(path.repo, "data/paleon_site_info.csv"))
summary(paleon)

# This is the v2 NADA at 2.5 degree resolution
nada.nc <- nc_open(file.path(path.data, "NADAv2-2008.nc"))

# This is the living blended drought atlas PDSI
lbda.nc <- nc_open(file.path(path.data, "nada_hd2_cl.nc"))


nada.lat <- ncvar_get(nada.nc, "lat")
nada.lon <- ncvar_get(nada.nc, "lon")
nada.time <- ncvar_get(nada.nc, "time")
nada.nc$dim$time$units # Time goes present to past; starting in 2007

lbda.lat <- ncvar_get(lbda.nc, "lat")
lbda.lon <- ncvar_get(lbda.nc, "lon")
lbda.time <- ncvar_get(nada.nc, "time")
lbda.nc$dim$time$units # Time goes present to past; starting in 2007

nada.lon.ind <- which(nada.lon>=min(paleon$lon) & nada.lon<=max(paleon$lon))
nada.lat.ind <- which(nada.lat>=min(paleon$lat) & nada.lat<=max(paleon$lat))

lbda.lon.ind <- which(lbda.lon>=min(paleon$lon) & lbda.lon<=max(paleon$lon))
lbda.lat.ind <- which(lbda.lat>=min(paleon$lat) & lbda.lat<=max(paleon$lat))

nada.raw <- ncvar_get(nada.nc, "PDSI")[nada.lon.ind,nada.lat.ind,]

dim(nada.raw)

# --------------------------------------------


# --------------------------------------------
# 2. Calculate Stability Index for all grid cells
#    -  Call calc.second.deriv (calc_second_deriv.R)
# --------------------------------------------
source(file.path(path.repo, "scripts/calc_second_deriv.R"))

nada.stability <- data.frame(lat=rep(nada.lat[nada.lat.ind], each=length(nada.lon.ind)),
                             lon=rep(nada.lon[nada.lon.ind], length.out=length(nada.lon.ind)*length(nada.lat.ind)))

for(i in 1:length(nada.lon.ind)){
  for(j in 1:length(nada.lat.ind)){
    # Finding the site index for our data frame
    site.ind <- which(nada.stability$lat==nada.lat[nada.lat.ind[j]] & nada.stability$lon==nada.lon[nada.lon.ind[i]])

    # Pulling the raw data
    ts.raw <- nada.raw[i,j,]
    
    if(length(ts.raw[!is.na(ts.raw)])==0) next
    
    # Finding out what years actually have data
    yr1 <- min(which(!is.na(ts.raw)))
    yr2 <- max(which(!is.na(ts.raw)))
    
   
    # Saving the time frame in our stability info
    nada.stability[site.ind,"yr.end"]   <- nada.time[yr1] 
    nada.stability[site.ind,"yr.start"] <- nada.time[yr2] 
    nada.stability[site.ind,"stability"] <- calc.second.deriv(ts.raw[yr2:yr1], h=1, H=1)
    
  }
}

nada.stability <- nada.stability[complete.cases(nada.stability),]
nada.stability$n.yrs <- nada.stability$yr.end - nada.stability$yr.start +1
summary(nada.stability)
write.csv(nada.stability, file.path(path.repo, "data/NADA_Stability.csv"), row.names=F, eol="\r\n")

library(ggplot2)

us <- map_data("state")

pdf(file.path(path.repo, "figures/NADA_Stability.pdf"))
print(
ggplot(data=nada.stability) +
  geom_tile(aes(x=lon, y=lat, fill=stability)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(nada.stability$lon), ylim=range(nada.stability$lat)) +
  theme_bw()
)
# print(
# ggplot(data=nada.stability) +
#   geom_tile(aes(x=lon, y=lat, fill=n.yrs)) +
#   geom_path(data=us,aes(x=long, y=lat, group=group), color="gray30") + 
#   coord_equal(xlim=range(nada.stability$lon), ylim=range(nada.stability$lat)) +
#   coord_equal()
# )
print(
ggplot(data=nada.stability) +
  geom_tile(aes(x=lon, y=lat, fill=stability/n.yrs)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(nada.stability$lon), ylim=range(nada.stability$lat)) +
  theme_bw()
)
print(
ggplot(data=nada.stability) +
  geom_histogram(aes(stability/n.yrs))
)
dev.off()

# --------------------------------------------
