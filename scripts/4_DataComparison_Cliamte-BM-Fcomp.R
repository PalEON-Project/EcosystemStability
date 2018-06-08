# -------------------------------------------
# Data-Climate Comparisons: 
# 
# H0: Climate change causes ecosystem change, as seen through shifts in 
#     biomass and composition (consistent & proportional responses; data)
#        Prediction: Locations that have experienced the most changes in 
#                    climate show the most changes in biomass and composition
#        Prediction: Places with the greatest change in composition also 
#                    have the greatest change in biomass (Zhang 2018 Nature Climate Change)
# 
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
lbda$stability.lbda <- -log(lbda$diff.abs/mean(lbda$diff.abs, na.rm=T))  # Note: Positive numbers mean MORE stable
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

# -----------
# STEPPS (empirical composition, centennially-resolved)
# -----------
stepps <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_STEPPS2.csv"))
summary(stepps)

# Change names to match up with drivers
stepps$model <- as.factor("STEPPS")
stepps$class <- as.factor("composition")
stepps$var <- stepps$taxon
stepps$type <- as.factor("empirical")
stepps$resolution <- as.factor("centennial")
stepps$stability.1k <- -log(stepps$stepps.diff.abs.1k/mean(stepps$stepps.diff.abs.1k, na.rm=T))
stepps$stability.lbda <- -log(stepps$stepps.diff.abs.lbda/mean(stepps$stepps.diff.abs.lbda, na.rm=T))
summary(stepps)
# -----------

# -----------
# ReFAB (empirical biomass, centennially-resolved)
# -----------
refab <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_ReFAB.csv"))
# names(refab) <- c("X", "lat", "lon", "diff.mean", "diff.abs", "n.sig")
summary(refab)

# Doing a unit correction
cols.bm <- c("refab.mean.1k", "refab.mean.slope.1k", "refab.mean.slope.abs.1k", "refab.mean.lbda", "refab.mean.slope.lbda", "refab.mean.slope.abs.lbda")
refab[,cols.bm] <- refab[,cols.bm]*0.1/2 # Native units = Mg/Ha biomass

refab$fract.sig <- refab$n.sig/10
refab$model <- as.factor("ReFAB")
refab$class <- as.factor("biomass")
refab$var <- as.factor("biomass")
refab$type <- as.factor("empirical")
refab$resolution <- as.factor("centennial")
refab$stability.1k <- -log(refab$refab.mean.slope.abs.1k/mean(refab$refab.mean.slope.abs.1k, na.rm=T))
refab$stability.lbda <- -log(refab$refab.mean.slope.abs.lbda/mean(refab$refab.mean.slope.abs.lbda, na.rm=T))
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
                        stab.refab.1k=refab$stability.1k,
                        stab.refab.lbda=refab$stability.lbda)
dat.sites.refab <- dat.sites.refab[!is.na(dat.sites.refab$lat),]

# Just Refab sites
for(i in 1:nrow(dat.sites.refab)){
  # ------------
  # STEPPS Composition
  # ------------
  # Find the closest stepps site
  stepps.dist <- sqrt((stepps$lon - dat.sites.refab$lon[i])^2 + (stepps$lat - dat.sites.refab$lat[i])^2)
  stepps.ind <- which(stepps.dist==min(stepps.dist, na.rm=T) & !stepps$taxon %in% c("Deciduous", "Evergreen")) # Don't include generic Decid/Evergreen
  
  dist.stepps <- mean(stepps.dist[stepps.ind])
  dat.sites.refab[i, "dist.stepps"] <- dist.stepps
  
  if(dist.stepps>0.5) next 
  
  if(length(stepps.ind)>0){
    # Subset to just that site
    fcomp.tmp <- stepps[stepps.ind,]
    
    # Find the dominant PFT 
    ind.fcomp.1k <- which(fcomp.tmp$stepps.mean.1k==max(fcomp.tmp$stepps.mean.1k))
    dat.sites.refab[i, "H.prime.1k"    ] <- - sum(fcomp.tmp$stepps.mean.1k * log(fcomp.tmp$stepps.mean.1k)) # calculate shannon-weiner index
    dat.sites.refab[i, "dom.pft.1k"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.1k])
    dat.sites.refab[i, "stab.stepps.1k"] <- fcomp.tmp$stability.1k[ind.fcomp.1k]
    
    ind.fcomp.lbda <- which(fcomp.tmp$stepps.mean.lbda==max(fcomp.tmp$stepps.mean.lbda))
    if(length(ind.fcomp.lbda)>0){
      dat.sites.refab[i, "H.prime.lbda"    ] <- - sum(fcomp.tmp$stepps.mean.lbda * log(fcomp.tmp$stepps.mean.lbda)) # calculate shannon-weiner index
      dat.sites.refab[i, "dom.pft.lbda"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.lbda])
      dat.sites.refab[i, "stab.stepps.lbda"] <- fcomp.tmp$stability.lbda[ind.fcomp.lbda]
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
    dat.sites.refab[i, "stab.lbda"] <- mean(lbda$stability[lbda.ind]) # If we have 2 equidistant sites, use the mean
  }
  # ------------
  
}
dat.sites.refab$dom.pft.1k <- as.factor(dat.sites.refab$dom.pft.1k)
dat.sites.refab$dom.pft.lbda <- as.factor(dat.sites.refab$dom.pft.lbda)
summary(dat.sites.refab)

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
  fcomp.tmp <- stepps[stepps.ind,]
  
  # Find the dominant PFT 
  ind.fcomp.1k <- which(fcomp.tmp$stepps.mean.1k==max(fcomp.tmp$stepps.mean.1k))
  dat.sites.stepps[i, "H.prime.1k"    ] <- - sum(fcomp.tmp$stepps.mean.1k * log(fcomp.tmp$stepps.mean.1k)) # calculate shannon-weiner index
  dat.sites.stepps[i, "dom.pft.1k"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.1k])
  dat.sites.stepps[i, "stab.stepps.1k"] <- fcomp.tmp$stability.1k[ind.fcomp.1k]
  
  ind.fcomp.lbda <- which(fcomp.tmp$stepps.mean.lbda==max(fcomp.tmp$stepps.mean.lbda))
  if(length(ind.fcomp.lbda)>0){
    dat.sites.stepps[i, "H.prime.lbda"    ] <- - sum(fcomp.tmp$stepps.mean.lbda * log(fcomp.tmp$stepps.mean.lbda)) # calculate shannon-weiner index
    dat.sites.stepps[i, "dom.pft.lbda"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.lbda])
    dat.sites.stepps[i, "stab.stepps.lbda"] <- fcomp.tmp$stability.lbda[ind.fcomp.lbda]
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
    dat.sites.stepps[i, "stab.lbda"] <- mean(lbda$stability[lbda.ind]) # If we have 2 equidistant sites, use the mean
  }
  # ------------
  
}
dat.sites.stepps$dom.pft.1k <- as.factor(dat.sites.stepps$dom.pft.1k)
dat.sites.stepps$dom.pft.lbda <- as.factor(dat.sites.stepps$dom.pft.lbda)
summary(dat.sites.stepps)

dim(dat.sites.stepps); dim(stepps)
# -----------

# -----------
# Evalutating correlation with Climate
# -----------
bm.pdsi <- lm(stab.refab.lbda ~ stab.lbda, data=dat.sites.refab)
bm.pdsi2 <- lm(stab.refab.lbda ~ stab.lbda, data=dat.sites.refab[dat.sites.refab$nyrs.lbda>=900,])
summary(bm.pdsi)
summary(bm.pdsi2)

# fcomp.pdsi <- lm(stab.stepps.lbda ~ stab.lbda, data=dat.sites.refab)
# summary(fcomp.pdsi)
fcomp.pdsi <- lm(stab.stepps.lbda ~ stab.lbda, data=dat.sites.stepps)
fcomp.pdsi2 <- lm(stab.stepps.lbda ~ stab.lbda, data=dat.sites.stepps[dat.sites.stepps$nyrs.lbda>=900 & !is.na(dat.sites.stepps$stab.lbda),])
summary(fcomp.pdsi)

summary(fcomp.pdsi2); 
nrow(dat.sites.stepps[dat.sites.stepps$nyrs.lbda>=900 & !is.na(dat.sites.stepps$stab.lbda),])


t.test(dat.sites.refab$stab.refab.lbda, dat.sites.refab$stab.lbda, paired=T)
t.test(dat.sites.refab$stab.stepps.lbda, dat.sites.refab$stab.lbda, paired=T)
t.test(dat.sites.stepps$stab.stepps.lbda, dat.sites.stepps$stab.lbda, paired=T)

mean(dat.sites.refab$stab.refab.lbda - dat.sites.refab$stab.lbda, na.rm=T)
sd(dat.sites.refab$stab.refab.lbda - dat.sites.refab$stab.lbda, na.rm=T)

mean(dat.sites.stepps$stab.stepps.lbda - dat.sites.stepps$stab.lbda, na.rm=T)
sd(dat.sites.stepps$stab.stepps.lbda - dat.sites.stepps$stab.lbda, na.rm=T)
# -----------


# -----------
# Quantitative comparisons between Biomass & Composition Stability
# -----------
# Comparing means -- is one more stable than the other?
t.test(dat.sites.refab$stab.refab.1k, dat.sites.refab$stab.stepps.1k, paired=T)
mean(dat.sites.refab$stab.refab.1k - dat.sites.refab$stab.stepps.1k)
sd(dat.sites.refab$stab.refab.1k - dat.sites.refab$stab.stepps.1k)

# Is biomass stability correlated with composition stability?
lm.comp <- lm(stab.refab.1k ~ stab.stepps.1k, data=dat.sites.refab)
summary(lm.comp)

# Is biomass stability correlated with diversity?
lm.diversity <- lm(stab.refab.1k ~ H.prime.1k, data=dat.sites.refab)
summary(lm.diversity)

# Is composition stability correlated with Diversity
lm.diversity2 <- lm(stab.stepps.1k ~ H.prime.1k, data=dat.sites.refab)
summary(lm.diversity2)

lm.diversity2b <- lm(stab.stepps.1k ~ H.prime.1k, data=dat.sites.stepps)
summary(lm.diversity2b)

# Some summary figures
png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Data_Biomass_v_Composition.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites.refab) +
  geom_point(aes(x=stab.stepps.1k, y=stab.refab.1k)) +
  stat_smooth(aes(x=stab.stepps.1k, y=stab.refab.1k), method="lm") +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Data_Biomass_v_Diversity.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites.refab) +
  geom_point(aes(x=H.prime.1k, y=stab.refab.1k)) +
  stat_smooth(aes(x=H.prime.1k, y=stab.refab.1k), method="lm") +
  theme_bw()
dev.off()
# -----------

# -------------------------------------------
# -------------------------------------------
