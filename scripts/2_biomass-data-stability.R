
load('/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/Current Data/ReFAB.all.samps.list.Rdata')
load('/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/Current Data/refab.sites.lat.lon.Rdata')

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
diff.mat.list <- list()
refab.diff.mean <- sig_vals <- matrix(NA, 62, 99)
for(i in 1:62){
  if(i == 2){
    diff.mat.list[[i]] <- NULL #My dataset is currently missing billy's lake which is in position 2. I can remove this when Billy's Lake is added.
  } else {
    diff.mat.list[[i]] <- apply(all.samps.list[[i]], 1, function(x) diff(rev(x)/100, na.rm=TRUE)) 
    # dividing by 100 because estimates are cenntennial. 
    # rev() because they are not sequential in time b[1] is present b[100] is 10,000 years BP
    
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
diff.mean = rowMeans(refab.diff.mean[,time.bin],na.rm = TRUE)
diff.mean.abs = rowMeans(abs(refab.diff.mean[,time.bin]),na.rm = TRUE)
n.signif = apply(sig_vals[,time.bin], 1, function(x) sum(!is.na(x)))

lat.lon.df.missing <- matrix(NA,62,2)
lat.lon.df.missing[1,] <- as.numeric(lat.lon.df[1,])
lat.lon.df.missing[3:62,] <- as.matrix(lat.lon.df[2:61,])

refab.mean.slope = data.frame(lat = lat.lon.df.missing[,1], lon = lat.lon.df.missing[,2],
                            refab.mean.slope =diff.mean,
                            refab.mean.slope.abs = diff.mean.abs,
                            n.signif = n.signif)

write.csv(refab.mean.slope, file = '~/Google Drive/PalEON_ecosystem-change_models-vs-data/refab.mean.slope.csv')
