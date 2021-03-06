library(sp)
library(rgdal)
library(foreach)
library(doParallel)
library(reshape2)

path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"

load(file.path(path.google, "Current Data/raw_data/STEPPS_NEUS_2k", 'stepps_neus_2k_input_data.rdata'))
taxa <- colnames(y)

stepps <- readRDS(file.path(path.google, "Current Data/raw_data/STEPPS_NEUS_2k", 'STEPPS_iterations_NEUS.RDS'))
dim(stepps)

iter.keep <- sample(251:dim(stepps)[4], 100) # For speed, parsing down to a random subset of 100 iterations (we don't need 1746!)

stepps = melt(stepps[,,,iter.keep])
colnames(stepps) = c("cell", "taxon", "time", "iter", "value")
dim(stepps)

stepps$taxon = taxa[stepps$taxon]
stepps$taxon <- as.factor(stepps$taxon)
stepps$year <- 1850-stepps$time*100 # for NEUS most recent reconstruction is 1750 CE
stepps <- stepps[stepps$year>=800,] # Subset to just the time slice we're interested for this
stepps$x = coord.agg.final[stepps$cell, 'east']
stepps$y = coord.agg.final[stepps$cell, 'north']


stepps.sp <- stepps[stepps$taxon=="Oak" & stepps$year==1750,c("cell", "x", "y")]
# stepps.sp[,c("x", "y")] <- stepps.sp[,c("x", "y")]*1e6 # Andria had convereted to megameters
stepps.sp <- SpatialPointsDataFrame(coords=stepps.sp[,c("x", "y")], stepps.sp, proj4string = CRS("+init=epsg:3175")) # x,y * 1e6 to convert to meters
stepps.sp <- spTransform(stepps.sp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
stepps.sp <- data.frame(stepps.sp)[,1:5]
names(stepps.sp)[which(names(stepps.sp) %in% c("x.1", "y.1"))] <- c("lon", "lat")
#summary(stepps.sp)

# Add lat/lon in
# stepps <- merge(stepps, stepps.sp, by=c("x", "y"), all.x=TRUE)
stepps$lon = stepps.sp$lon[match(stepps$cell, stepps.sp$cell)]
stepps$lat = stepps.sp$lat[match(stepps$cell, stepps.sp$cell)]

# summary(stepps)

# Adding in the stability just for the period for which we have climate data
lbda.df <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_LBDA_100.csv"))
summary(lbda.df)

ggplot2::ggplot(data=lbda.df) +
  ggplot2::coord_equal() +
  ggplot2::geom_tile(ggplot2::aes(x=lon, y=lat, fill=n.yrs))

# Setting up a blank data frame to store stability metrics in
stab.stepps <- data.frame(cell=rep(unique(stepps$cell), each=length(unique(stepps$taxon))), 
                          taxon=unique(stepps$taxon))

N_taxa = length(unique(stepps$taxon))
N_iter = length(unique(stepps$iter))
N_cell = length(unique(stepps$cell))

taxa = as.vector(unique(stepps[, "taxon"]))


stab.stepps <- data.frame(cell = rep(rep(seq(1, N_cell), each=N_taxa), ntimes=N_iter),
                          lat = rep(NA, N_iter*N_cell*N_taxa),
                          lon = rep(NA, N_iter*N_cell*N_taxa),
                          taxon = taxa[rep(seq(1, N_taxa), ntimes=N_iter*N_cell)],
                          iter = rep(seq(1, N_iter), each=N_cell*N_taxa),
                          stepps.mean.1k = rep(NA, N_iter*N_cell*N_taxa),
                          stepps.diff.1k = rep(NA, N_iter*N_cell*N_taxa),
                          stepps.diff.abs.1k = rep(NA, N_iter*N_cell*N_taxa),
                          stepps.mean.lbda = rep(NA, N_iter*N_cell*N_taxa),
                          stepps.diff.lbda = rep(NA, N_iter*N_cell*N_taxa),
                          stepps.diff.abs.lbda = rep(NA, N_iter*N_cell*N_taxa),
                          ldba.min = rep(NA, N_iter*N_cell*N_taxa))

# head(stab.stepps)

# Looping through each grid cell
#setup parallel backend to use many processors
# cores=detectCores()
# cores <- 8
# cl <- makeCluster(cores) #not to overload your computer
# registerDoParallel(cl)
pb <- txtProgressBar(min=0, max=N_taxa*N_iter*N_cell, style=3)
pb.ind = 0
for(CELL in 1:N_cell){#N_cell){
# foreach(CELL=1:N_cell)%dopar%{
# out <- foreach(CELL=1:N_cell, .combine=rbind) %dopar% {
  
  # # print(paste0("Cell ", CELL, " of ", N_cell))
  # stab.stepps <- data.frame(cell = rep(rep(CELL, each=N_taxa), ntimes=N_iter),
  #                           lat = rep(NA, N_iter*N_taxa),
  #                           lon = rep(NA, N_iter*N_taxa),
  #                           taxon = taxa[rep(seq(1, N_taxa), ntimes=N_iter)],
  #                           iter = rep(seq(1, N_iter), each=N_taxa),
  #                           stepps.mean.1k = rep(NA, N_iter*N_taxa),  
  #                           stepps.diff.1k = rep(NA, N_iter*N_taxa),
  #                           stepps.diff.abs.1k = rep(NA, N_iter*N_taxa),
  #                           stepps.mean.lbda = rep(NA, N_iter*N_taxa),
  #                           stepps.diff.lbda = rep(NA, N_iter*N_taxa),
  #                           stepps.diff.abs.lbda = rep(NA, N_iter*N_taxa),
  #                           lbda.min = rep(NA, N_iter*N_taxa))
  
  idx.cell = stepps$cell==CELL
  
  # Extracting the stepps lat/lon
  stepps.lon <- unique(stepps[idx.cell,"lon"])
  stepps.lat <- unique(stepps[idx.cell,"lat"])
  
  idx.cell.stab = stab.stepps$cell==CELL
  stab.stepps[idx.cell.stab, "lon"] <- stepps.lon
  stab.stepps[idx.cell.stab, "lat"] <- stepps.lat
  
  # finidng the closest LBDA Grid cell
  lon.dist <- stepps.lon - lbda.df$lon
  lat.dist <- stepps.lat - lbda.df$lat
  ind.lbda <- which(sqrt(lon.dist^2 + lat.dist^2)==min(sqrt(lon.dist^2 + lat.dist^2)))
  
  # Get the number of years to look at by rounding to the nearest 100
  yr.min <- 1850-mean(lbda.df$n.yrs[ind.lbda])
  
  stab.stepps[idx.cell.stab,"lbda.min"]   <- yr.min
  
  # n.yrs <- round(mean(lbda.df$n.yrs[ind.lbda]), -2) # Adding mean() in case a site is equidistant between lbda points
  # stab.stepps[stab.stepps$cell==CELL,"n.yrs.lbda"] <- n.yrs
  # time.bin <- 99:(99-n.yrs/100) # Updating our time bin based on number of years
  
  for (iter in 1:N_iter){
    # if ((iter %% 5)==0) {print(paste0("Iteration --> ", iter))} 
    
    # Subset our cell data & order it through time  to speed things up
    dat.cell <- stepps[which(idx.cell&(stepps$iter==iter)),]
    dat.cell <- dat.cell[order(dat.cell$year),]
    
    # Looping through each PFT 
    for(TAX in taxa){
      setTxtProgressBar(pb, pb.ind); pb.ind=pb.ind+1
      # print(TAX)
      df.tmp <- dat.cell[dat.cell$cell==CELL & dat.cell$taxon==TAX,]
      
      idx.stab = which(idx.cell.stab & stab.stepps$taxon==TAX & stab.stepps$iter==iter)
      
      # Looking at stability over the full 1k period (Model Sim period)
      diff.1k <- diff(df.tmp$value)
      stab.stepps[idx.stab, "stepps.mean.1k"    ] <- mean(df.tmp$value)
      stab.stepps[idx.stab, "stepps.diff.1k"    ] <-  mean(diff.1k, na.rm = F)
      stab.stepps[idx.stab, "stepps.diff.abs.1k"] <-  mean(abs(diff.1k), na.rm = F)
      
      if(is.na(yr.min)) next # Skip anything that's not in an LBDA grid cell for some reason
      
      # Looking at stability just where we have empirical drought data
      df.tmp.sub = df.tmp[df.tmp$year>=yr.min,"value"]
      diff.lbda <- diff(df.tmp.sub)
      stab.stepps[idx.stab, "stepps.mean.lbda"    ] <- mean(df.tmp.sub)
      stab.stepps[idx.stab, "stepps.diff.lbda"    ] <-  mean(diff.lbda, na.rm = F)
      stab.stepps[idx.stab, "stepps.diff.abs.lbda"] <-  mean(abs(diff.lbda), na.rm = F)
      
    } # End Taxon Loop
  } # End Iter Loop
  # stab.stepps[idx.cell.stab,]
}
# stopCluster(cl)

# stab.stepps[which((stab.stepps$cell == 1) & (stab.stepps$iter==1)),]
# 
# summary(stab.stepps)

# write.csv(out, file = "data/Stability_STEPPS2.csv")
write.csv(stab.stepps, file = file.path(path.google, "Current Data/Stability", 'Stability_STEPPS2_NEUS_iter.csv'))

# stab.stepps.agg <- aggregate(stab.stepps[,"stepps.mean.1l", "stepps.diff.abs.1k"], 
#                              by=stab.stepps[,c("cell", "lat", "lon", "taxon")],
#                              FUN=mean, na.rm=T)
library(ggplot2)
ggplot(data=stab.stepps[stab.stepps$iter==100 & stab.stepps$taxon=="Hemlock",]) +
  coord_equal() +
  geom_point(aes(x=lon, y=lat, color=stepps.diff.abs.1k))
