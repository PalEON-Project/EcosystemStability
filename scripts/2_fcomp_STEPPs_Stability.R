library(sp); library(rgdal)
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"

# Load data from Google Drive
stepps <- readRDS(file.path(path.google, 'Current Data/raw_data/STEPPS/STEPPS_predictions_median_sd.RDS'))
stepps$taxon <- as.factor(stepps$taxon)
stepps$year <- 1950-stepps$time
stepps <- stepps[stepps$year>=800,] # Subset to just the time slice we're interested for this
summary(stepps)

stepps.sp <- stepps[stepps$taxon=="OAK" & stepps$year==1750,c("cell", "x", "y")]
# stepps.sp[,c("x", "y")] <- stepps.sp[,c("x", "y")]*1e6 # Andria had convereted to megameters
stepps.sp <- SpatialPointsDataFrame(coords=stepps.sp[,c("x", "y")]*1e6, stepps.sp, proj4string = CRS("+init=epsg:3175")) # x,y * 1e6 to convert to meters
stepps.sp <- spTransform(stepps.sp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
stepps.sp <- data.frame(stepps.sp)[,1:5]
names(stepps.sp)[which(names(stepps.sp) %in% c("x.1", "y.1"))] <- c("lon", "lat")
summary(stepps.sp)

# Add lat/lon in
stepps <- merge(stepps, stepps.sp)
summary(stepps)

# Adding in the stability just for the period for which we have climate data
lbda.df <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_LBDA_100.csv"))
summary(lbda.df)

# Setting up a blank data frame to store stability metrics in
stab.stepps <- data.frame(cell=rep(unique(stepps$cell), each=length(unique(stepps$taxon))),
                          taxon=unique(stepps$taxon))
head(stab.stepps)

# Looping through each grid cell
for(CELL in unique(stab.stepps$cell)){
  # Extracting the stepps lat/lon
  stepps.lon <- unique(stepps[stepps$cell==CELL,"lon"])
  stepps.lat <- unique(stepps[stepps$cell==CELL,"lat"])
  
  stab.stepps[stab.stepps$cell==CELL, "lon"] <- stepps.lon
  stab.stepps[stab.stepps$cell==CELL, "lat"] <- stepps.lat
  
  # finidng the closest LBDA Grid cell
  lon.dist <- stepps.lon - lbda.df$lon
  lat.dist <- stepps.lat - lbda.df$lat
  ind.lbda <- which(sqrt(lon.dist^2 + lat.dist^2)==min(sqrt(lon.dist^2 + lat.dist^2)))
  
  # Get the number of years to look at by rounding to the nearest 100
  yr.min <- 1850-mean(lbda.df$n.yrs[ind.lbda])
  
  stab.stepps[stab.stepps$cell==CELL,"lbda.min"]   <- yr.min
  
  # n.yrs <- round(mean(lbda.df$n.yrs[ind.lbda]), -2) # Adding mean() in case a site is equidistant between lbda points
  # stab.stepps[stab.stepps$cell==CELL,"n.yrs.lbda"] <- n.yrs
  # time.bin <- 99:(99-n.yrs/100) # Updating our time bin based on number of years
  
  
  # Subset our cell data & order it through time  to speed things up
  dat.cell <- stepps[stepps$cell==CELL,]
  dat.cell <- dat.cell[order(dat.cell$year),]

  # Looping through each PFT 
  for(TAX in unique(stepps[stepps$cell==CELL, "taxon"])){
    df.tmp <- dat.cell[dat.cell$cell==CELL & dat.cell$taxon==TAX,]
    
    # Looking at stability over the full 1k period (Model Sim period)
    diff.1k <- diff(df.tmp$median)
    stab.stepps[stab.stepps$cell==CELL & stab.stepps$taxon==TAX, "stepps.mean.1k"    ] <- mean(df.tmp$median)
    stab.stepps[stab.stepps$cell==CELL & stab.stepps$taxon==TAX, "stepps.diff.1k"    ] <-  mean(diff.1k, na.rm = F)
    stab.stepps[stab.stepps$cell==CELL & stab.stepps$taxon==TAX, "stepps.diff.abs.1k"] <-  mean(abs(diff.1k), na.rm = F)
    
    if(is.na(yr.min)) next # Skip anything that's not in an LBDA grid cell for some reason
    
    # Looking at stability just where we have empirical drought data
    diff.lbda <- diff(df.tmp[df.tmp$year>=yr.min,"median"])
    stab.stepps[stab.stepps$cell==CELL & stab.stepps$taxon==TAX, "stepps.mean.lbda"    ] <- mean(df.tmp[df.tmp$year>=yr.min,"median"])
    stab.stepps[stab.stepps$cell==CELL & stab.stepps$taxon==TAX, "stepps.diff.lbda"    ] <-  mean(diff.lbda, na.rm = F)
    stab.stepps[stab.stepps$cell==CELL & stab.stepps$taxon==TAX, "stepps.diff.abs.lbda"] <-  mean(abs(diff.lbda), na.rm = F)
  } # End Taxon Loop
} # End Grid Cell Loop
summary(stab.stepps)

write.csv(stab.stepps, file = file.path(path.google, "Current Data/Stability", 'Stability_STEPPS2.csv'))
