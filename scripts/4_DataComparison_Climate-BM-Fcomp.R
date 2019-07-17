# -------------------------------------------
# Data-Climate Comparisons: 
# 
# H0: Climate change causes ecosystem change, as seen through shifts in 
#     biomass and composition (consistent & proportional responses; data)
#        Prediction: Locations that have experienced the most changes in 
#                    climate show the most changes in biomass and composition
#        Prediction: Places with the greatest change in composition also 
#                    have the greatest change in biomass (Zhang 2018 Nature Climate Change)
# H1: Changes in composition facilitate biomass stability or vice versa
#        Prediction: There is a negative correlation between stability of
#                    composition and biomass (stable biomass because comp shift)
#        Prediction: Sites with higher diversity have higher biomass stability
# -------------------------------------------
# Steps:
# -------------------------------------------
# 0. Load libraries, set file paths etc.
# 1. Read in Datasets; Standardize Stability Metrics
# 2. Compare Biomass & Composition to Climate (climate window; do they show the same pattern?)
# 3. Compare Biomass to Composition (Full time?; are they correlated at all?)
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# Load Libaries & Set file paths
# -------------------------------------------
library(ggplot2); library(gridExtra); library(scales)
library(grid)

library(sp)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"

# Storing the paths for US state outlines
us <- map_data("state")

# Set up a table with colors for models & data
# Using a CB palette from here: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
dat.colors <- data.frame(model=c("LBDA", "STEPPS", "ReFAB", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"),
                         type=c(rep("emipirical", 3), rep("model", 7)),
                         color=c(rep("#000000", 3), rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
dat.colors
# -------------------------------------------

# -------------------------------------------
# 1. Read in Datasets; Standardize Stability Metrics
# -------------------------------------------
# -----------
# Empirical Drought Reconstruction (LBDA, annually-resolved)
# -----------
lbda <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_LBDA_100.csv"))
lbda$model <- as.factor("LBDA")
lbda$class <- as.factor("climate")
lbda$var <- as.factor("pdsi")
lbda$type <- as.factor("empirical")
lbda$resolution <- as.factor("annual")
names(lbda)[which(names(lbda)=="n.yrs.sig")] <- "n.sig"
# lbda$stability.lbda <- -log(lbda$diff.abs/abs(mean(lbda$lbda.mean, na.rm=T)))  # Note: Positive numbers mean MORE stable
lbda$variability.lbda <- lbda$diff.abs/abs(mean(lbda$lbda.mean, na.rm=T))  # Note: Positive numbers mean MORE stable
summary(lbda)


# png(file.path(path.google, "Current Figures/Stability_Synthesis", "PDSI_LBDA_TimeSeriesDepth.png"), height=5, width=8, units="in", res=320)
# ggplot(data=lbda) +
#   geom_tile(aes(x=lon, y=lat, fill=n.yrs)) +
#   geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
#   ggtitle("Time Series Length: LBDA") +
#   coord_equal(xlim=range(lbda$lon), ylim=range(lbda$lat), expand=0) +
#   theme_bw() 
# dev.off()
# -----------

# -------------------------------------------
# STEPPS (empirical composition, centennially-resolved)
# -------------------------------------------
# stepps <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_STEPPS2_iters.csv"))
stepps.umw <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_STEPPS2_iters.csv"))
summary(stepps.umw)

stepps.neus <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_STEPPS2_NEUS_iter.csv"))
summary(stepps.neus)

vars.stepps <- c("stepps.mean.1k", "stepps.diff.1k", "stepps.diff.abs.1k", "stepps.mean.lbda", "stepps.diff.lbda", "stepps.diff.abs.lbda", "lbda.min")
stepps.umw <- aggregate(stepps.umw[,vars.stepps],
                        by=stepps.umw[,c("lat", "lon", "taxon")],
                        FUN=mean)

# NEUS has iterations at this point
stepps.neus <- aggregate(stepps.neus[,vars.stepps],
                         by=stepps.neus[,c("lat", "lon", "taxon")],
                         FUN=mean)
stepps.umw$model <- as.factor("STEPPS-UMW")
stepps.neus$model <- as.factor("STEPPS-NEUS")

# stepps.umw$variability.1k.region <- stepps.umw$stepps.diff.abs.1k/abs(mean(stepps.umw$stepps.mean.1k, na.rm=T))
# stepps.umw$variability.lbda.region <- stepps.umw$stepps.diff.abs.lbda/abs(mean(stepps.umw$stepps.mean.lbda, na.rm=T))
# stepps.neus$variability.1k.region <- stepps.neus$stepps.diff.abs.1k/abs(mean(stepps.neus$stepps.mean.1k, na.rm=T))
# stepps.neus$variability.lbda.region <- stepps.neus$stepps.diff.abs.lbda/abs(mean(stepps.neus$stepps.mean.lbda, na.rm=T))


# Bind UMW & NEUS together
stepps <- rbind(stepps.umw[,names(stepps.neus)], stepps.neus)
summary(stepps)

# Change names to match up with drivers
stepps$taxon <- toupper(stepps$taxon)
# stepps$model <- as.factor("STEPPS")
stepps$class <- as.factor("composition")
stepps$var <- stepps$taxon
stepps$type <- as.factor("empirical")
stepps$resolution <- as.factor("centennial")
# stepps$stability.1k <- -log(stepps$stepps.diff.abs.1k/abs(mean(stepps$stepps.mean.1k, na.rm=T)))
# stepps$stability.lbda <- -log(stepps$stepps.diff.abs.lbda/abs(mean(stepps$stepps.mean.lbda, na.rm=T)))
stepps$variability.1k <- stepps$stepps.diff.abs.1k/abs(mean(stepps$stepps.mean.1k, na.rm=T))
stepps$variability.lbda <- stepps$stepps.diff.abs.lbda/abs(mean(stepps$stepps.mean.lbda, na.rm=T))
summary(stepps)
unique(stepps$taxon)

write.csv(stepps, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_STEPPS_aggregated.csv"), row.names=F)


ggplot(stepps[stepps$taxon=="OAK",]) +
  coord_equal() +
  geom_point(aes(x=lon, y=lat, color=stepps.mean.1k))

# Parsing down STEPPS
dat.sites.stepps <- stepps[stepps$taxon=="OAK", c("lon", "lat")]
# summary(dat.sites.stepps)
for(i in 1:nrow(dat.sites.stepps)){
  # ------------
  # STEPPS Composition
  # ------------
  # Find the closest stepps site
  # stepps.dist <- sqrt((stepps$lon - dat.sites.stepps$lon[i])^2 + (stepps$lat - dat.sites.stepps$lat[i])^2)
  stepps.ind <- which(stepps$lat==dat.sites.stepps$lat[i] & stepps$lon==dat.sites.stepps$lon[i] & !stepps$taxon %in% c("Deciduous", "Evergreen")) # Don't include generic Decid/Evergreen
  
  # Subset to just that site
  fcomp.tmp.all <- stepps[stepps.ind,]
  fcomp.tmp <- fcomp.tmp.all[fcomp.tmp.all$stepps.mean.1k>1e-3,]
  
  # Find the dominant PFT 
  ind.fcomp.1k <- which(fcomp.tmp$stepps.mean.1k==max(fcomp.tmp$stepps.mean.1k))
  
  dat.sites.stepps[i, "richness.1k"   ] <- nrow(fcomp.tmp)
  dat.sites.stepps[i, "H.prime.1k"    ] <- - sum(fcomp.tmp$stepps.mean.1k * log(fcomp.tmp$stepps.mean.1k)) # calculate shannon-weiner index
  dat.sites.stepps[i, "dom.pft.1k"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.1k])
  dat.sites.stepps[i, "dom.mean.1k"   ] <- fcomp.tmp$stepps.mean.1k[ind.fcomp.1k]
  dat.sites.stepps[i, "diff.abs.1k"   ] <- fcomp.tmp$stepps.diff.abs.1k[ind.fcomp.1k]
  # dat.sites.stepps[i, "var.stepps.1k"] <- fcomp.tmp$stability.1k[ind.fcomp.1k]
  dat.sites.stepps[i, "var.1k.all"] <- fcomp.tmp$variability.1k[ind.fcomp.1k]
  
  fcomp.tmp <- fcomp.tmp.all[fcomp.tmp.all$stepps.mean.lbda>1e-3,]
  ind.fcomp.lbda <- which(fcomp.tmp$stepps.mean.lbda==max(fcomp.tmp$stepps.mean.lbda))
  if(length(ind.fcomp.lbda)>0){
    dat.sites.stepps[i, "richness.lbda"   ] <- nrow(fcomp.tmp)
    dat.sites.stepps[i, "H.prime.lbda"    ] <- - sum(fcomp.tmp$stepps.mean.lbda * log(fcomp.tmp$stepps.mean.lbda)) # calculate shannon-weiner index
    dat.sites.stepps[i, "dom.pft.lbda"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.lbda])
    dat.sites.stepps[i, "dom.mean.lbda"   ] <- fcomp.tmp$stepps.mean.lbda[ind.fcomp.lbda]
    
    dat.sites.stepps[i, "diff.abs.lbda"   ] <- fcomp.tmp$stepps.diff.abs.lbda[ind.fcomp.lbda]
    # dat.sites.stepps[i, "var.stepps.lbda"] <- fcomp.tmp$stability.lbda[ind.fcomp.lbda]
    dat.sites.stepps[i, "var.lbda.all"] <- fcomp.tmp$variability.lbda[ind.fcomp.lbda]
  }
  # ------------
  
  # ------------
  # Drought Atlas
  # ------------
  # Find the closest stepps site
  lbda.dist <- sqrt((lbda$lon - dat.sites.stepps$lon[i])^2 + (lbda$lat - dat.sites.stepps$lat[i])^2)
  lbda.ind <- which(lbda.dist==min(lbda.dist, na.rm=T) ) 
  
  if(length(lbda.ind)>0){
    dat.sites.stepps[i, "nyrs.lbda"] <- mean(lbda$n.yrs[lbda.ind]) # If we have 2 equidistant sites, use the mean
    # dat.sites.stepps[i, "var.lbda"] <- mean(lbda$stability.lbda[lbda.ind]) # If we have 2 equidistant sites, use the mean
    dat.sites.stepps[i, "var.lbda"] <- mean(lbda$variability.lbda[lbda.ind]) # If we have 2 equidistant sites, use the mean
  }
  # ------------
  
}
dat.sites.stepps$dataset2 <- as.factor(ifelse(dat.sites.stepps$lon>-80, "STEPPS-NEUS", "STEPPS-UMW"))
dat.sites.stepps$dom.pft.1k <- as.factor(dat.sites.stepps$dom.pft.1k)
dat.sites.stepps$dom.pft.lbda <- as.factor(dat.sites.stepps$dom.pft.lbda)
dat.sites.stepps$var.stepps.1k <- dat.sites.stepps$diff.abs.1k/abs(mean(dat.sites.stepps$dom.mean.1k, na.rm=T))
dat.sites.stepps$var.stepps.lbda <- dat.sites.stepps$diff.abs.lbda/abs(mean(dat.sites.stepps$dom.mean.lbda, na.rm=T))

for(REG in unique(dat.sites.stepps$dataset2)){
  ind.reg <- which(dat.sites.stepps$dataset2==REG)
  dat.sites.stepps[ind.reg, "var.stepps.1k.region"] <- dat.sites.stepps$diff.abs.1k[ind.reg]/abs(mean(dat.sites.stepps$dom.mean.1k[ind.reg], na.rm=T))
  dat.sites.stepps[ind.reg, "var.stepps.lbda.region"] <- dat.sites.stepps$diff.abs.lbda[ind.reg]/abs(mean(dat.sites.stepps$dom.mean.lbda[ind.reg], na.rm=T))
  
}
# dat.sites.stepps$var.stepps.lbda <- dat.sites.stepps$diff.abs.lbda/abs(mean(dat.sites.stepps$dom.mean.lbda, na.rm=T))

summary(dat.sites.stepps)

plot(var.1k.all ~ var.stepps.1k, data=dat.sites.stepps)
# par(mfrow=c(2,1)); hist(log(dat.sites.stepps$var.1k.all)); hist(log(dat.sites.stepps$var.stepps.1k)); par(mfrow=c(1,1))
dim(dat.sites.stepps); dim(stepps)

ggplot(data=dat.sites.stepps) +
  geom_histogram(aes(x=log(var.stepps.1k), fill=dataset2))
ggplot(data=dat.sites.stepps) +
  geom_histogram(aes(x=log(var.stepps.1k.region), fill=dataset2))

ggplot(dat=dat.sites.stepps) +
  coord_equal() +
  geom_point(aes(x=lon, y=lat, color=log(var.1k.all))) +
  scale_color_gradient2(name="Log\nRelative\nvariability", low="#27647B", high="#CA3542", limits=range(log(dat.sites.stepps$var.1k.all), na.rm=T), midpoint=mean(log(dat.sites.stepps$var.1k.all), na.rm=T)) 
  
ggplot(dat=dat.sites.stepps) +
  coord_equal() +
  geom_point(aes(x=lon, y=lat, color=dom.mean.1k)) +
  scale_color_gradient2(name="Mean\nCover of\nDominant\nPFT", low="#27647B", high="#CA3542", limits=range(dat.sites.stepps$dom.mean.1k, na.rm=T), midpoint=mean(dat.sites.stepps$dom.mean.1k, na.rm=T)) 
  

write.csv(dat.sites.stepps, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_STEPPS.csv"), row.names=F)
# -------------------------------------------

# -------------------------------------------
# ReFAB (empirical biomass, centennially-resolved)
# -------------------------------------------
refab <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_ReFAB.csv"))
refab <- refab[!is.na(refab$lat),]
# names(refab) <- c("X", "lat", "lon", "diff.mean", "diff.abs", "n.sig")
summary(refab)

# Doing a unit correction
cols.bm <- c("refab.mean.1k", "refab.mean.slope.1k", "refab.mean.slope.abs.1k", "refab.mean.lbda", "refab.mean.slope.lbda", "refab.mean.slope.abs.lbda")
refab[,cols.bm] <- refab[,cols.bm]*0.1/2 # Native units = Mg/Ha biomass

# refab$fract.sig <- refab$n.sig/10
refab$model <- as.factor("ReFAB")
refab$class <- as.factor("biomass")
refab$var <- as.factor("biomass")
refab$type <- as.factor("empirical")
refab$resolution <- as.factor("centennial")
# refab$stability.1k <- -log(refab$refab.mean.slope.abs.1k/abs(mean(refab$refab.mean.1k, na.rm=T)))
# refab$stability.lbda <- -log(refab$refab.mean.slope.abs.lbda/abs(mean(refab$refab.mean.lbda, na.rm=T)))
refab$variability.1k <- refab$refab.mean.slope.abs.1k/abs(mean(refab$refab.mean.1k, na.rm=T))
refab$variability.lbda <- refab$refab.mean.slope.abs.lbda/abs(mean(refab$refab.mean.lbda, na.rm=T))
summary(refab)
# -----------
# -------------------------------------------

# -------------------------------------------
# 3. Compare Biomass to Composition (Full time?; are they correlated at all?)
# -------------------------------------------
# -----------
# Extracting STEPPS & LBDA info for refab
# -----------

dat.sites.refab <- data.frame(refab[,c("lat", "lon")],
                              # var.refab.1k=refab$stability.1k,
                              # var.refab.lbda=refab$stability.lbda,
                              var.refab.1k=refab$variability.1k,
                              var.refab.lbda=refab$variability.lbda)
dat.sites.refab <- dat.sites.refab[!is.na(dat.sites.refab$lat),]

# Just Refab sites
for(i in 1:nrow(dat.sites.refab)){
  # ------------
  # STEPPS Composition
  # ------------
  # Find the closest stepps site
  stepps.dist <- sqrt((stepps$lon - dat.sites.refab$lon[i])^2 + (stepps$lat - dat.sites.refab$lat[i])^2)
  stepps.ind <- which(stepps.dist==min(stepps.dist, na.rm=T) & !stepps$taxon %in% c("Deciduous", "Evergreen")) # Don't include generic Decid/Evergreen

  # # Just pull the stability data from the dominant PFT
  stepps.dist2 <- sqrt((dat.sites.stepps$lon - dat.sites.refab$lon[i])^2 + (dat.sites.stepps$lat - dat.sites.refab$lat[i])^2)
  stepps.ind2 <- which(stepps.dist2==min(stepps.dist2, na.rm=T))
  
  dist.stepps <- mean(stepps.dist[stepps.ind])
  dat.sites.refab[i, "dist.stepps"] <- dist.stepps
  
  if(dist.stepps>0.5) next 
  
  if(length(stepps.ind)>0){
    # Subset to just that site, considering only sites with abundance >=0.1% (1e-3)
    fcomp.tmp.all <- stepps[stepps.ind,]
    fcomp.tmp <- fcomp.tmp.all[fcomp.tmp.all$stepps.mean.1k>1e-3,]
    
    # Find the dominant PFT 
    ind.fcomp.1k <- which(fcomp.tmp$stepps.mean.1k==max(fcomp.tmp$stepps.mean.1k))
    
    dat.sites.refab[i, "richness.1k"   ] <- nrow(fcomp.tmp)
    dat.sites.refab[i, "H.prime.1k"    ] <- - sum(fcomp.tmp$stepps.mean.1k * log(fcomp.tmp$stepps.mean.1k)) # calculate shannon-weiner index
    dat.sites.refab[i, "dom.pft.1k"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.1k])
    dat.sites.refab[i, "dom.mean.1k"   ] <- fcomp.tmp$stepps.mean.1k[ind.fcomp.1k]
    # dat.sites.refab[i, "var.stepps.1k"] <- fcomp.tmp$stability.1k[ind.fcomp.1k]
    dat.sites.stepps[i, "var.1k.all"] <- fcomp.tmp$variability.1k[ind.fcomp.1k]
    dat.sites.refab[i, "var.stepps.1k"] <- dat.sites.stepps$var.stepps.1k[stepps.ind2]
    
    fcomp.tmp <- fcomp.tmp.all[fcomp.tmp.all$stepps.mean.lbda>1e-3,]
    ind.fcomp.lbda <- which(fcomp.tmp$stepps.mean.lbda==max(fcomp.tmp$stepps.mean.lbda))
    if(length(ind.fcomp.lbda)>0){
      dat.sites.refab[i, "richness.lbda"   ] <- nrow(fcomp.tmp)
      dat.sites.refab[i, "H.prime.lbda"    ] <- - sum(fcomp.tmp$stepps.mean.lbda * log(fcomp.tmp$stepps.mean.lbda)) # calculate shannon-weiner index
      dat.sites.refab[i, "dom.pft.lbda"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.lbda])
      dat.sites.refab[i, "dom.mean.lbda"   ] <- fcomp.tmp$stepps.mean.lbda[ind.fcomp.lbda]
      # dat.sites.refab[i, "var.stepps.lbda"] <- fcomp.tmp$stability.lbda[ind.fcomp.lbda]
      dat.sites.refab[i, "var.stepps.lbda"] <- fcomp.tmp$variability.lbda[ind.fcomp.lbda]
      dat.sites.refab[i, "var.stepps.lbda"] <- dat.sites.stepps$var.stepps.lbda[stepps.ind2]
    }
  }
  # ------------
  
  # ------------
  # Drought Atlas
  # ------------
  # Find the closest stepps site
  lbda.dist <- sqrt((lbda$lon - dat.sites.refab$lon[i])^2 + (lbda$lat - dat.sites.refab$lat[i])^2)
  lbda.ind <- which(lbda.dist==min(lbda.dist, na.rm=T) ) 
  
  if(length(lbda.ind)>0){
    dat.sites.refab[i, "nyrs.lbda"] <- mean(lbda$n.yrs[lbda.ind]) # If we have 2 equidistant sites, use the mean
    # dat.sites.refab[i, "var.lbda"] <- mean(lbda$stability[lbda.ind]) # If we have 2 equidistant sites, use
    dat.sites.refab[i, "var.lbda"] <- mean(lbda$variability[lbda.ind]) # If we have 2 equidistant sites, use the mean
  }
  # ------------
  
}
dat.sites.refab$dom.pft.1k <- as.factor(dat.sites.refab$dom.pft.1k)
dat.sites.refab$dom.pft.lbda <- as.factor(dat.sites.refab$dom.pft.lbda)
summary(dat.sites.refab)

write.csv(dat.sites.refab , file.path(path.google, "Current Data/Stability_Synthesis", "Stability_ReFAB.csv"), row.names=F)
# -----------

# -----------
# Comparing distribution of stability
# -----------
lbda.stepps <- lbda[lbda$lon >= min(stepps$lon, refab$lon, na.rm=T) & lbda$lon <= max(stepps$lon, refab$lon, na.rm=T) & 
                      lbda$lat >= min(stepps$lat, refab$lat, na.rm=T) & lbda$lat <= max(stepps$lat, refab$lat, na.rm=T), ]
summary(lbda.stepps)
var.comparison <- data.frame(lat=c(dat.sites.refab$lat, dat.sites.stepps$lat, lbda.stepps$lat),
                              lon=c(dat.sites.refab$lon, dat.sites.stepps$lon, lbda.stepps$lon),
                              dataset = c(rep("ReFAB", nrow(dat.sites.refab)), rep("STEPPS", nrow(dat.sites.stepps)), rep("LBDA", nrow(lbda.stepps))),
                              dataset2 = c(rep("ReFAB", nrow(dat.sites.refab)), dat.sites.stepps$dataset2, rep("LBDA", nrow(lbda.stepps))),
                              variability = c(dat.sites.refab$var.refab.lbda, dat.sites.stepps$var.stepps.lbda, lbda.stepps$variability.lbda)
                              )
var.comparison$dataset2 <- as.factor(ifelse(var.comparison$dataset=="STEPPS" & var.comparison$lon>-80, "STEPPS-NEUS", paste(var.comparison$dataset)))
var.comparison$dataset <- factor(var.comparison$dataset, levels=c("LBDA", "STEPPS", "ReFAB"))
var.comparison$dataset2 <- factor(var.comparison$dataset2, levels=c("LBDA", "STEPPS", "STEPPS-NEUS", "ReFAB"))
summary(var.comparison)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Comparisons_Histograms.png"), height=6, width=5, units="in", res=320)
ggplot(data=var.comparison) +
  facet_grid(dataset~., scales="free_y") +
  geom_histogram(aes(x=log(variability), fill=dataset2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Comparisons_Density.png"), height=6, width=5, units="in", res=320)
ggplot(data=var.comparison) +
  facet_grid(dataset~.) +
  geom_density(aes(x=log(variability), fill=dataset2), alpha=0.8) +
  theme_bw()
dev.off()


# Making a figure showing correlation of STEPPS & ReFAB with climate (or lack thereof)
climate.comparison <- data.frame(lat=c(dat.sites.refab$lat, dat.sites.stepps$lat),
                                 lon=c(dat.sites.refab$lon, dat.sites.stepps$lon),
                                 dataset = c(rep("ReFAB", nrow(dat.sites.refab)), rep("STEPPS", nrow(dat.sites.stepps))),
                                 # mean.ecosys = c(dat.sites.refab$m)
                                 # var.ecosys = c(dat.sites.refab$var.refab.lbda, dat.sites.stepps$var.stepps.lbda),
                                 # var.pdsi   = c(dat.sites.refab$var.lbda, dat.sites.stepps$var.lbda),
                                 var.ecosys = c(dat.sites.refab$var.refab.lbda, dat.sites.stepps$var.stepps.lbda),
                                 var.pdsi   = c(dat.sites.refab$var.lbda, dat.sites.stepps$var.lbda)
                                 )

climate.comparison.sp <- data.frame(lat=c(dat.sites.refab$lat, dat.sites.stepps$lat, lbda$lat),
                                    lon=c(dat.sites.refab$lon, dat.sites.stepps$lon, lbda$lon),
                                    dataset = c(rep("ReFAB", nrow(dat.sites.refab)), rep("STEPPS", nrow(dat.sites.stepps)), rep("LBDA", nrow(lbda))),
                                    # stability = c(refab$stability.lbda, stepps$stability.lbda, lbda$stability.lbda),
                                    variability = c(dat.sites.refab$var.refab.lbda, dat.sites.stepps$var.stepps.1k, lbda$variability.lbda)
                                    )
climate.comparison$dataset2 <- as.factor(ifelse(climate.comparison$dataset=="STEPPS" & climate.comparison$lon>-81, "STEPPS-NEUS", paste(climate.comparison$dataset)))
climate.comparison$dataset <- factor(climate.comparison$dataset, levels=c("LBDA", "STEPPS", "ReFAB"))
climate.comparison$dataset2 <- factor(climate.comparison$dataset2, levels=c("LBDA", "STEPPS", "STEPPS-NEUS", "ReFAB"))

climate.comparison.sp$dataset2 <- as.factor(ifelse(climate.comparison.sp$dataset=="STEPPS" & climate.comparison.sp$lon>-81, "STEPPS-NEUS", paste(climate.comparison.sp$dataset)))
climate.comparison.sp$dataset <- factor(climate.comparison.sp$dataset, levels=c("LBDA", "STEPPS", "ReFAB"))
climate.comparison.sp$dataset2 <- factor(climate.comparison.sp$dataset2, levels=c("LBDA", "STEPPS", "STEPPS-NEUS", "ReFAB"))

write.csv(climate.comparison, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"), row.names=F)
write.csv(climate.comparison.sp, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data_Spatial.csv"), row.names=F)
# -----------------------

# -----------------------------------------------------------------------
# Start here if you're just re-running analyses but don't need to do all the data munging
# -----------------------------------------------------------------------
stepps <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_STEPPS_aggregated.csv"))
dat.sites.refab <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_ReFAB.csv"))
dat.sites.stepps <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_STEPPS.csv"))
climate.comparison <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"))
climate.comparison.sp <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data_Spatial.csv"))

climate.comparison.sp$dataset <- car::recode(climate.comparison.sp$dataset, "'LBDA'='Drought'; 'STEPPS'='Composition'; 'ReFAB'='Biomass'")
climate.comparison.sp$dataset <- factor(climate.comparison.sp$dataset, levels=c("Drought", "Composition", "Biomass"))
summary(climate.comparison.sp)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Ecosystem_v_Climate_Data_Map.png"), height=5.5, width=5, units="in", res=320)
ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability) & climate.comparison.sp$dataset!="Drought",]) +
  facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray90") +
  geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_tile(data=climate.comparison.sp[climate.comparison.sp$dataset=="Drought" & !is.na(climate.comparison.sp$variability),], aes(x=lon, y=lat, fill=log(variability))) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  scale_fill_gradient2(name="Log\nRelative\nvariability", low="#27647B", high="#CA3542", limits=range(log(climate.comparison.sp$variability), na.rm=T), midpoint=mean(log(climate.comparison.sp$variability), na.rm=T)) +
  scale_color_gradient2(name="Log\nRelative\nvariability", low="#27647B", high="#CA3542", limits=range(log(climate.comparison.sp$variability), na.rm=T), midpoint=mean(log(climate.comparison.sp$variability), na.rm=T)) +
  theme_bw() +
  theme(panel.background=element_rect(fill="gray25"),
        panel.grid = element_blank(),
        legend.position = "top")
dev.off()  

library(cowplot)
plot.lbda <-  ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability) & climate.comparison.sp$dataset=="Drought",]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  # geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_tile(data=climate.comparison.sp[climate.comparison.sp$dataset=="Drought" & !is.na(climate.comparison.sp$variability),], aes(x=lon, y=lat, fill=log(variability))) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  scale_fill_gradientn(name="Drought\nVariability", colors=c("#018571","#80cdc1", "#f5f5f5", "#dfc27d", "#a6611a"), limits=range(log(climate.comparison.sp$variability[climate.comparison.sp$dataset=="Drought"]), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

plot.stepps <- ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability) & climate.comparison.sp$dataset=="Composition",]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  scale_color_gradientn(name="Compos.\nVariability", colors=c("#008837","#a6dba0", "#f7f7f7", "#c2a5cf", "#7b3294"), limits=range(log(climate.comparison.sp$variability[climate.comparison.sp$dataset=="Composition"]), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

plot.refab <- ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability) & climate.comparison.sp$dataset=="Biomass",]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  scale_color_gradientn(name="Biomass\nVariability", colors=c("#4dac26","#b8e186", "#f7f7f7", "#f1b6da", "#d01c8b"), limits=range(log(climate.comparison.sp$variability[climate.comparison.sp$dataset=="Biomass"]), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Ecosystem_v_Climate_Data_Map_FreeColor.png"), height=6, width=8, units="in", res=220)
plot_grid(plot.lbda, plot.stepps, plot.refab, ncol=1)
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Ecosystem_v_Climate_Data.png"), height=6, width=6, units="in", res=320)
ggplot(data=climate.comparison) +
  geom_point(aes(x=log(var.pdsi), y=log(var.ecosys), color=dataset2), size=0.5) +
  stat_smooth(aes(x=log(var.pdsi), y=log(var.ecosys), color=dataset2, fill=dataset2), method="lm") +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "LBDA_depth_STEPPS_area.png"), height=6, width=8, units="in", res=320)
ggplot(data=lbda[,]) +
  geom_tile(aes(x=lon, y=lat, fill=n.yrs)) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  # scale_fill_continuous(limits=range(climate.comparison.sp$stability, na.rm=T)) +
  theme_bw() 
dev.off()
# -----------

# -----------
# Comparing mean relative variability
# -----------
summary(climate.comparison.sp)

# STEPPS Composition; UMW vs NEUS
mean(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="STEPPS-NEUS"]), na.rm=T); sd(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="STEPPS-NEUS"]), na.rm=T)
mean(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="STEPPS"]), na.rm=T); sd(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="STEPPS"]), na.rm=T)

t.test(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="STEPPS"]), log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="STEPPS-NEUS"]), na.rm=T)

# Drought Variability
summary(climate.comparison.sp)
mean(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="LBDA"]), na.rm=T); sd(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="LBDA"]), na.rm=T)

# Biomass Variability
summary(climate.comparison.sp)
mean(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="ReFAB"]), na.rm=T); sd(log(climate.comparison.sp$variability[climate.comparison.sp$dataset2=="ReFAB"]), na.rm=T)

# -----------

# -----------
# Quantitative comparisons with climate
# -----------
# Is biomass variability correlated with climate?
bm.pdsi <- lm(log(var.refab.lbda) ~ log(var.lbda), data=dat.sites.refab)
bm.pdsi2 <- lm(log(var.refab.lbda) ~ log(var.lbda), data=dat.sites.refab[dat.sites.refab$nyrs.lbda>=900,])
summary(bm.pdsi)
summary(bm.pdsi2)
par(mfrow=c(2,2)); plot(bm.pdsi); par(mfrow=c(1,1))
par(mfrow=c(2,2)); plot(bm.pdsi2); par(mfrow=c(1,1))

# Is composition varility correlated with climate?
# fcomp.pdsi <- lm(var.stepps.lbda ~ var.lbda, data=dat.sites.refab)
# summary(fcomp.pdsi)
summary(dat.sites.stepps)
fcomp.pdsi <- lm(log(var.stepps.lbda) ~ log(var.lbda)*dataset2-1, data=dat.sites.stepps)
fcomp.pdsi2 <- lm(log(var.stepps.lbda) ~ log(var.lbda)*dataset2-1-log(var.lbda), data=dat.sites.stepps[!is.na(dat.sites.stepps$var.lbda),]) # dat.sites.stepps$nyrs.lbda>=900 & 
fcomp.pdsi3 <- lm(log(var.stepps.lbda) ~ log(var.lbda)+dataset2-1, data=dat.sites.stepps[!is.na(dat.sites.stepps$var.lbda),]) # dat.sites.stepps$nyrs.lbda>=900 & 
summary(fcomp.pdsi); anova(fcomp.pdsi)
summary(fcomp.pdsi2); anova(fcomp.pdsi2)
summary(fcomp.pdsi3); anova(fcomp.pdsi3)
par(mfrow=c(2,2)); plot(fcomp.pdsi); par(mfrow=c(1,1))
par(mfrow=c(2,2)); plot(fcomp.pdsi2); par(mfrow=c(1,1))
nrow(dat.sites.stepps[dat.sites.stepps$nyrs.lbda>=900 & !is.na(dat.sites.stepps$var.lbda),])
nrow(lbda.stepps[lbda.stepps$n.yrs>=900 & !is.na(lbda.stepps$varility.lbda),])

# Is biomass and/or composition more/less varle than overall climate?
t.test(dat.sites.refab$var.refab.lbda, dat.sites.refab$var.lbda, paired=T)
t.test(log(dat.sites.refab$var.refab.lbda), log(dat.sites.refab$var.lbda), paired=T)
# t.test(dat.sites.refab$var.stepps.lbda, dat.sites.refab$var.lbda, paired=T)
t.test(dat.sites.stepps$var.stepps.lbda, dat.sites.stepps$var.lbda, paired=T)

mean(dat.sites.refab$var.refab.lbda - dat.sites.refab$var.lbda, na.rm=T)
sd(dat.sites.refab$var.refab.lbda - dat.sites.refab$var.lbda, na.rm=T)

mean(dat.sites.stepps$var.stepps.lbda - dat.sites.stepps$var.lbda, na.rm=T)
sd(dat.sites.stepps$var.stepps.lbda - dat.sites.stepps$var.lbda, na.rm=T)

# -----------


# -----------
# Quantitative comparisons between Biomass & Composition Stability
# -----------
dat.sites <- data.frame(lat=c(dat.sites.refab$lat, dat.sites.stepps$lat), 
                        lon=c(dat.sites.refab$lon, dat.sites.stepps$lon), 
                        dataset=c(rep("ReFAB", nrow(dat.sites.refab)), rep("STEPPS", nrow(dat.sites.stepps))),
                        # stability=c(dat.sites.refab$var.refab.1k, dat.sites.stepps$var.stepps.1k),
                        variability=c(dat.sites.refab$var.refab.1k, dat.sites.stepps$var.stepps.1k),
                        H.prime = c(dat.sites.refab$H.prime.1k, dat.sites.stepps$H.prime.1k),
                        richness = c(dat.sites.refab$richness.1k, dat.sites.stepps$richness.1k),
                        dom.pft = c(paste(dat.sites.refab$dom.pft.1k), paste(dat.sites.stepps$dom.pft.1k))
                        )
dat.sites$dataset2 <- as.factor(ifelse(dat.sites$dataset=="STEPPS" & dat.sites$lon >- 80, "STEPPS-NEUS", paste(dat.sites$dataset)))
dat.sites[dat.sites$dom.pft=="NA","dom.pft"] <- NA
summary(dat.sites)

write.csv(dat.sites, file.path(path.google, "Current Data/Stability_Synthesis", "Variability_Ecosystem_v_Diversity_Data.csv"), row.names=F)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Map_STEPPS_Hprime.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites[,]) +
  facet_grid(dataset~.) +
  geom_point(aes(x=lon, y=lat, color=H.prime), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  # scale_fill_continuous(limits=range(dat.sites$H.prime, na.rm=T)) +
  theme_bw() 
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Map_STEPPS_Richness.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites[,]) +
  facet_grid(dataset~.) +
  geom_point(aes(x=lon, y=lat, color=richness), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  # scale_fill_continuous(limits=range(dat.sites$H.prime, na.rm=T)) +
  theme_bw() 
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Map_STEPPS_DominantPFT.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites[,]) +
  facet_grid(dataset~.) +
  geom_point(aes(x=lon, y=lat, color=dom.pft), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group)) +
  coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
  # scale_fill_continuous(limits=range(dat.sites$H.prime, na.rm=T)) +
  theme_bw() +
  theme(legend.position="top")
dev.off()
# ggplot(data=dat.sites[,]) +
#   facet_grid(dataset~.) +
#   geom_point(aes(x=lon, y=lat, color=richness), size=2) +
#   geom_path(data=us, aes(x=long, y=lat, group=group)) +
#   coord_equal(xlim=range(stepps$lon, na.rm=T), ylim=range(stepps$lat, na.rm=T)) +
#   # scale_fill_continuous(limits=range(dat.sites$H.prime, na.rm=T)) +
#   theme_bw() 


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Data_Variability_v_Diversity_NoTrendline.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites) +
  facet_grid(dataset2~., scales="free_y") +
  geom_point(aes(x=H.prime, y=log(variability), color=dom.pft)) +
  # stat_smooth(data=dat.sites[dat.sites$dataset=="ReFAB",], aes(x=H.prime, y=stability, color=dataset, fill=dataset), method="lm") +
  # stat_smooth(data=dat.sites[dat.sites$dataset=="STEPPS",], aes(x=H.prime, y=stability, color=dataset, fill=dataset), method="lm") +
  theme_bw() +
  theme(legend.position="top")
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Data_Variability_v_Richness_NoTrendline.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites) +
  facet_grid(dataset2~., scales="free_y") +
  geom_point(aes(x=richness, y=log(variability), color=dom.pft)) +
  # geom_boxplot(aes(x=richness, y=log(variability), fill=dom.pft), position="dodge", alpha=0.5) +
  # geom_violin(aes(x=richness, y=log(variability), fill=dom.pft), position="dodge", alpha=0.5) +
  # stat_smooth(data=dat.sites[dat.sites$dataset=="ReFAB",], aes(x=H.prime, y=stability, color=dataset, fill=dataset), method="lm") +
  # stat_smooth(data=dat.sites[dat.sites$dataset=="STEPPS",], aes(x=H.prime, y=stability, color=dataset, fill=dataset), method="lm") +
  theme_bw() +
  theme(legend.position="top")
dev.off()

# Comparing means -- is one more stable than the other?
t.test(log(dat.sites.refab$var.refab.1k), log(dat.sites.refab$var.stepps.1k), paired=T)
mean(log(dat.sites.refab$var.refab.1k) - log(dat.sites.refab$var.stepps.1k), na.rm=T)
sd(log(dat.sites.refab$var.refab.1k) - log(dat.sites.refab$var.stepps.1k), na.rm=T)

# Is biomass stability correlated with composition stability?
lm.comp <- lm(log(var.refab.1k) ~ log(var.stepps.1k), data=dat.sites.refab)
summary(lm.comp)

# Is biomass stability correlated with diversity?
lm.diversity1 <- lm(log(var.refab.1k) ~ H.prime.1k, data=dat.sites.refab)
summary(lm.diversity1)

lm.diversity1b <- lm(log(var.refab.1k) ~ H.prime.1k + I(H.prime.1k^2), data=dat.sites.refab)
summary(lm.diversity1b)
# dat.sites[dat.sites$dataset=="ReFAB" & !is.na(dat.sites$H.prime), "H.prime.quad"] <- predict(lm.diversity1b)
AIC(lm.diversity1, lm.diversity1b)  # Linear model best

# Is composition stability correlated with Diversity
lm.diversity2 <- lm(log(var.stepps.1k) ~ H.prime.1k*dataset2-1, data=dat.sites.stepps)
summary(lm.diversity2)

lm.diversity2b <- lm(log(var.stepps.1k) ~ (H.prime.1k + I(H.prime.1k^2))*dataset2-1, data=dat.sites.stepps)
# dat.sites[dat.sites$dataset=="STEPPS", "H.prime.quad"] <- predict(lm.diversity2b)
summary(lm.diversity2b)
AIC(lm.diversity2, lm.diversity2b)

# Testing residuals for normalcy
par(mfrow=c(2,2)); plot(lm.diversity2b); par(mfrow=c(1,1))
plot(predict(lm.diversity2b)~dat.sites.stepps$var.stepps.1k[!is.na(dat.sites.stepps$H.prime.1k)]); abline(a=0, b=1, col="red")
plot(resid(lm.diversity2b)~predict(lm.diversity2b))
hist(resid(lm.diversity2b))

# Some summary figures
png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Data_Biomass_v_Composition.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites.refab) +
  geom_point(aes(x=var.stepps.1k, y=var.refab.1k)) +
  stat_smooth(aes(x=var.stepps.1k, y=var.refab.1k), method="lm") +
  theme_bw()
dev.off()


# png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Data_Stability_v_Diversity_ModFit_Quad.png"), height=6, width=6, units="in", res=320)
# ggplot(data=dat.sites) +
#   facet_wrap(~dataset, scales="free", ncol=1) +
#   geom_point(aes(x=H.prime.quad, y=stability)) +
#   stat_smooth(aes(x=H.prime.quad, y=stability), method="lm") +
#   theme_bw()
# dev.off()

# -----------

# -------------------------------------------
# -------------------------------------------