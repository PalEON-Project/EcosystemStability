path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"

# Load data from Google Drive
load(file.path(path.google, 'Current Data/raw_data/ReFAB/ReFAB.all.samps.list_v2.Rdata'))
load(file.path(path.google, 'Current Data/raw_data/ReFAB/refab.sites.lat.lon_v2.Rdata'))

### Andria's function for calculating significance
prob_sig <- function(x, prob){
  n_side = max(sum(x>0), sum(x<0))
  
  if (n_side/length(x) > prob){
    return(mean(x))
  } else {
    return(NA)
  }
}

### loop for getting mean differences and significance for all refab sites
mean.mat.list <- diff.mat.list <- list()
refab.mean <- refab.diff.mean <- sig_vals <- matrix(NA, length(all.samps.list), 99)
for(i in 1:length(all.samps.list)){
  if(length(all.samps.list[[i]]) == 0){
    mean.mat.list[[i]] <- NULL #My dataset is currently missing billy's lake which is in position 2. I can remove this when Billy's Lake is added.
    diff.mat.list[[i]] <- NULL #My dataset is currently missing billy's lake which is in position 2. I can remove this when Billy's Lake is added.
  } else {
    mean.mat.list[[i]] <- apply(all.samps.list[[i]], 1, mean, na.rm=FALSE) 
    diff.mat.list[[i]] <- apply(all.samps.list[[i]], 1, function(x) diff(rev(x)/100, na.rm=FALSE)) 
    # dividing by 100 because estimates are cenntennial. 
    # rev() because they are not sequential in time b[1] is present b[100] is 10,000 years BP
    # From here out values go from past to present
    
    refab.mean[i,] <- mean(mean.mat.list[[i]])
    refab.diff.mean[i,] <- rowMeans(diff.mat.list[[i]])
    
    sig_vals[i,] <- apply(t(diff.mat.list[[i]]), c(2), prob_sig, prob=.85)
    
  }
}

### checking calculations blue and red should match
boxplot(t(diff.mat.list[[1]]))
points(diff(rev(colMeans(all.samps.list[[1]]))/100), col='blue', pch=19)
points(refab.diff.mean[1,] ,col='red')

### summing accross time for significance
time.bin <- 90:99 ## 850 - 1850 AD : 1000 - 100 years before present
mean.all = rowMeans(refab.mean[,time.bin],na.rm = F)
diff.mean.all = rowMeans(refab.diff.mean[,time.bin],na.rm = F)
diff.mean.abs.all = rowMeans(abs(refab.diff.mean[,time.bin]),na.rm = F)
n.signif.all = apply(sig_vals[,time.bin], 1, function(x) sum(!is.na(x)))

lat.lon.df.missing <- matrix(NA,length(all.samps.list),2)
# lat.lon.df.missing[1,] <- as.numeric(lat.lon.df[1,])
# lat.lon.df.missing[3:62,] <- as.matrix(lat.lon.df[2:61,])

refab.mean.slope = data.frame(lat = lat.lon.df[,1], lon = lat.lon.df[,2],
                              refab.mean.1k = mean.all,
                              refab.mean.slope.1k =diff.mean.all,
                              refab.mean.slope.abs.1k = diff.mean.abs.all,
                              n.signif.1k = n.signif.all)

lat.lon.df[which(is.na(refab.mean.slope$refab.mean.1k)),]

# Adding in the stability just for the period for which we have climate data
lbda.df <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_LBDA_100.csv"))
summary(lbda.df)

# Loop through each point and get the number of points we should look at
for(i in 1:nrow(refab.mean.slope)){
  if(is.na(refab.mean.slope$lat[i])) next 
  
  # Find the closest LDBA grid cell
  lon.refab <- refab.mean.slope$lon[i]
  lat.refab <- refab.mean.slope$lat[i]
  
  lon.dist <- lon.refab - lbda.df$lon
  lat.dist <- lat.refab - lbda.df$lat
  ind.lbda <- which(sqrt(lon.dist^2 + lat.dist^2)==min(sqrt(lon.dist^2 + lat.dist^2)))
  
  # if(length(ind.lbda)>1) ind.lbda <- ind.libda[1] # If a point is equidistant
  
  # Get the number of years to look at by rounding to the nearest 100
  n.yrs <- round(mean(lbda.df$n.yrs[ind.lbda]), -2) # Adding mean() in case a site is equidistant between lbda points
  
  if(is.na(n.yrs)) next # Skip anythign that's not in an LBDA grid cell for some reason
  
  time.bin <- 99:(99-n.yrs/100) # Updating our time bin based on number of years
  
  # Extracting/storing the info for the LBDA period
  refab.mean.slope[i,"n.yrs.lbda"] <- n.yrs
  refab.mean.slope[i, "refab.mean.lbda"] <-  mean(refab.mean[i,time.bin],na.rm = F)
  refab.mean.slope[i, "refab.mean.slope.lbda"    ] <-  mean(refab.diff.mean[i,time.bin],na.rm = F)
  refab.mean.slope[i, "refab.mean.slope.abs.lbda"] <-  mean(abs(refab.diff.mean[i,time.bin]),na.rm = F)
}
summary(refab.mean.slope)

write.csv(refab.mean.slope, file = file.path(path.google, "Current Data/Stability", 'Stability_ReFAB.csv'))
