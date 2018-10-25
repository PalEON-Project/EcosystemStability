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
lbda$stability.lbda <- -log(lbda$diff.abs/abs(mean(lbda$lbda.mean, na.rm=T)))  # Note: Positive numbers mean MORE stable
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
stepps$stability.1k <- -log(stepps$stepps.diff.abs.1k/abs(mean(stepps$stepps.mean.1k, na.rm=T)))
stepps$stability.lbda <- -log(stepps$stepps.diff.abs.lbda/abs(mean(stepps$stepps.mean.lbda, na.rm=T)))
stepps$variability.1k <- stepps$stepps.diff.abs.1k/abs(mean(stepps$stepps.mean.1k, na.rm=T))
stepps$variability.lbda <- stepps$stepps.diff.abs.lbda/abs(mean(stepps$stepps.mean.lbda, na.rm=T))
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
refab$stability.1k <- -log(refab$refab.mean.slope.abs.1k/abs(mean(refab$refab.mean.1k, na.rm=T)))
refab$stability.lbda <- -log(refab$refab.mean.slope.abs.lbda/abs(mean(refab$refab.mean.lbda, na.rm=T)))
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
                              stab.refab.1k=refab$stability.1k,
                              stab.refab.lbda=refab$stability.lbda,
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
    dat.sites.refab[i, "stab.stepps.1k"] <- fcomp.tmp$stability.1k[ind.fcomp.1k]
    dat.sites.refab[i, "var.stepps.1k"] <- fcomp.tmp$variability.1k[ind.fcomp.1k]
    
    fcomp.tmp <- fcomp.tmp.all[fcomp.tmp.all$stepps.mean.lbda>1e-3,]
    ind.fcomp.lbda <- which(fcomp.tmp$stepps.mean.lbda==max(fcomp.tmp$stepps.mean.lbda))
    if(length(ind.fcomp.lbda)>0){
      dat.sites.refab[i, "richness.lbda"   ] <- nrow(fcomp.tmp)
      dat.sites.refab[i, "H.prime.lbda"    ] <- - sum(fcomp.tmp$stepps.mean.lbda * log(fcomp.tmp$stepps.mean.lbda)) # calculate shannon-weiner index
      dat.sites.refab[i, "dom.pft.lbda"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.lbda])
      dat.sites.refab[i, "stab.stepps.lbda"] <- fcomp.tmp$stability.lbda[ind.fcomp.lbda]
      dat.sites.refab[i, "var.stepps.lbda"] <- fcomp.tmp$variability.lbda[ind.fcomp.lbda]
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
    dat.sites.refab[i, "stab.lbda"] <- mean(lbda$stability[lbda.ind]) # If we have 2 equidistant sites, use
    dat.sites.refab[i, "var.lbda"] <- mean(lbda$variability[lbda.ind]) # If we have 2 equidistant sites, use the mean
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
  fcomp.tmp.all <- stepps[stepps.ind,]
  fcomp.tmp <- fcomp.tmp.all[fcomp.tmp.all$stepps.mean.1k>1e-3,]
  
  # Find the dominant PFT 
  ind.fcomp.1k <- which(fcomp.tmp$stepps.mean.1k==max(fcomp.tmp$stepps.mean.1k))
  
  dat.sites.stepps[i, "richness.1k"   ] <- nrow(fcomp.tmp)
  dat.sites.stepps[i, "H.prime.1k"    ] <- - sum(fcomp.tmp$stepps.mean.1k * log(fcomp.tmp$stepps.mean.1k)) # calculate shannon-weiner index
  dat.sites.stepps[i, "dom.pft.1k"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.1k])
  dat.sites.stepps[i, "stab.stepps.1k"] <- fcomp.tmp$stability.1k[ind.fcomp.1k]
  dat.sites.stepps[i, "var.stepps.1k"] <- fcomp.tmp$variability.1k[ind.fcomp.1k]
  
  fcomp.tmp <- fcomp.tmp.all[fcomp.tmp.all$stepps.mean.lbda>1e-3,]
  ind.fcomp.lbda <- which(fcomp.tmp$stepps.mean.lbda==max(fcomp.tmp$stepps.mean.lbda))
  if(length(ind.fcomp.lbda)>0){
    dat.sites.stepps[i, "richness.lbda"   ] <- nrow(fcomp.tmp)
    dat.sites.stepps[i, "H.prime.lbda"    ] <- - sum(fcomp.tmp$stepps.mean.lbda * log(fcomp.tmp$stepps.mean.lbda)) # calculate shannon-weiner index
    dat.sites.stepps[i, "dom.pft.lbda"    ] <- paste(fcomp.tmp$taxon[ind.fcomp.lbda])
    dat.sites.stepps[i, "stab.stepps.lbda"] <- fcomp.tmp$stability.lbda[ind.fcomp.lbda]
    dat.sites.stepps[i, "var.stepps.lbda"] <- fcomp.tmp$variability.lbda[ind.fcomp.lbda]
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
    dat.sites.stepps[i, "stab.lbda"] <- mean(lbda$stability.lbda[lbda.ind]) # If we have 2 equidistant sites, use the mean
    dat.sites.stepps[i, "var.lbda"] <- mean(lbda$variability.lbda[lbda.ind]) # If we have 2 equidistant sites, use the mean
  }
  # ------------
  
}
dat.sites.stepps$dom.pft.1k <- as.factor(dat.sites.stepps$dom.pft.1k)
dat.sites.stepps$dom.pft.lbda <- as.factor(dat.sites.stepps$dom.pft.lbda)
summary(dat.sites.stepps)

dim(dat.sites.stepps); dim(stepps)

write.csv(dat.sites.stepps, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_STEPPS.csv"), row.names=F)
write.csv(dat.sites.refab , file.path(path.google, "Current Data/Stability_Synthesis", "Stability_ReFAB.csv"), row.names=F)
# -----------

# -----------
# Comparing distribution of stability
# -----------
lbda.stepps <- lbda[lbda$lon >= min(stepps$lon, refab$lon, na.rm=T) & lbda$lon <= max(stepps$lon, refab$lon, na.rm=T) & 
                      lbda$lat >= min(stepps$lat, refab$lat, na.rm=T) & lbda$lat <= max(stepps$lat, refab$lat, na.rm=T), ]
summary(lbda.stepps)
stab.comparison <- data.frame(lat=c(dat.sites.refab$lat, dat.sites.stepps$lat, lbda.stepps$lat),
                              lon=c(dat.sites.refab$lon, dat.sites.stepps$lon, lbda.stepps$lon),
                              dataset = c(rep("ReFAB", nrow(dat.sites.refab)), rep("STEPPS", nrow(dat.sites.stepps)), rep("LBDA", nrow(lbda.stepps))),
                              stability = c(dat.sites.refab$stab.refab.lbda, dat.sites.stepps$stab.stepps.lbda, lbda.stepps$stability.lbda),
                              variability = c(dat.sites.refab$var.refab.lbda, dat.sites.stepps$var.stepps.lbda, lbda.stepps$variability.lbda)
                              
                              )
stab.comparison$dataset <- factor(stab.comparison$dataset, levels=c("LBDA", "STEPPS", "ReFAB"))

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Comparisons_Histograms.png"), height=6, width=5, units="in", res=320)
ggplot(data=stab.comparison) +
  facet_grid(dataset~.) +
  geom_histogram(aes(x=log(variability))) +
  theme_bw()
dev.off()


# Making a figure showing correlation of STEPPS & ReFAB with climate (or lack thereof)
climate.comparison <- data.frame(lat=c(dat.sites.refab$lat, dat.sites.stepps$lat),
                                 lon=c(dat.sites.refab$lon, dat.sites.stepps$lon),
                                 dataset = c(rep("ReFAB", nrow(dat.sites.refab)), rep("STEPPS", nrow(dat.sites.stepps))),
                                 stab.ecosys = c(dat.sites.refab$stab.refab.lbda, dat.sites.stepps$stab.stepps.lbda),
                                 stab.pdsi   = c(dat.sites.refab$stab.lbda, dat.sites.stepps$stab.lbda),
                                 var.ecosys = c(dat.sites.refab$var.refab.lbda, dat.sites.stepps$var.stepps.lbda),
                                 var.pdsi   = c(dat.sites.refab$var.lbda, dat.sites.stepps$var.lbda)
                                 )

climate.comparison.sp <- data.frame(lat=c(refab$lat, stepps$lat, lbda$lat),
                                    lon=c(refab$lon, stepps$lon, lbda$lon),
                                    dataset = c(rep("ReFAB", nrow(refab)), rep("STEPPS", nrow(stepps)), rep("LBDA", nrow(lbda))),
                                    stability = c(refab$stability.lbda, stepps$stability.lbda, lbda$stability.lbda),
                                    variability = c(refab$variability.lbda, stepps$variability.lbda, lbda$variability.lbda)
                                    )
climate.comparison$dataset <- factor(climate.comparison$dataset, levels=c("LBDA", "STEPPS", "ReFAB"))
climate.comparison.sp$dataset <- factor(climate.comparison.sp$dataset, levels=c("LBDA", "STEPPS", "ReFAB"))

write.csv(climate.comparison, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"), row.names=F)
write.csv(climate.comparison.sp, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data_Spatial.csv"), row.names=F)

climate.comparison <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"))
climate.comparison.sp <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data_Spatial.csv"))

climate.comparison.sp$dataset <- car::recode(climate.comparison.sp$dataset, "'LBDA'='Drought'; 'STEPPS'='Composition'; 'ReFAB'='Biomass'")
climate.comparison.sp$dataset <- factor(climate.comparison.sp$dataset, levels=c("Drought", "Composition", "Biomass"))

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Ecosystem_v_Climate_Data_Map.png"), height=5.5, width=5, units="in", res=320)
ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability),]) +
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
        panel.grid = element_blank())
dev.off()  

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Ecosystem_v_Climate_Data.png"), height=6, width=6, units="in", res=320)
ggplot(data=climate.comparison) +
  geom_point(aes(x=log(var.pdsi), y=log(var.ecosys), color=dataset), size=0.5) +
  stat_smooth(aes(x=log(var.pdsi), y=log(var.ecosys), color=dataset, fill=dataset), method="lm") +
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
# Quantitative comparisons with climate
# -----------
# Is biomass stability correlated with climate?
bm.pdsi <- lm(stab.refab.lbda ~ stab.lbda, data=dat.sites.refab)
bm.pdsi2 <- lm(stab.refab.lbda ~ stab.lbda, data=dat.sites.refab[dat.sites.refab$nyrs.lbda>=900,])
summary(bm.pdsi)
summary(bm.pdsi2)

# Is composition stability correlated with climate?
# fcomp.pdsi <- lm(stab.stepps.lbda ~ stab.lbda, data=dat.sites.refab)
# summary(fcomp.pdsi)
fcomp.pdsi <- lm(stab.stepps.lbda ~ stab.lbda, data=dat.sites.stepps)
fcomp.pdsi2 <- lm(stab.stepps.lbda ~ stab.lbda, data=dat.sites.stepps[dat.sites.stepps$nyrs.lbda>=900 & !is.na(dat.sites.stepps$stab.lbda),])
summary(fcomp.pdsi)
summary(fcomp.pdsi2); 
nrow(dat.sites.stepps[dat.sites.stepps$nyrs.lbda>=900 & !is.na(dat.sites.stepps$stab.lbda),])
nrow(lbda.stepps[lbda.stepps$n.yrs>=900 & !is.na(lbda.stepps$stability.lbda),])

# Is biomass and/or composition more/less stable than overall climate?
t.test(dat.sites.refab$stab.refab.lbda, dat.sites.refab$stab.lbda, paired=T)
# t.test(dat.sites.refab$stab.stepps.lbda, dat.sites.refab$stab.lbda, paired=T)
t.test(dat.sites.stepps$stab.stepps.lbda, dat.sites.stepps$stab.lbda, paired=T)

mean(dat.sites.refab$stab.refab.lbda - dat.sites.refab$stab.lbda, na.rm=T)
sd(dat.sites.refab$stab.refab.lbda - dat.sites.refab$stab.lbda, na.rm=T)

mean(dat.sites.stepps$stab.stepps.lbda - dat.sites.stepps$stab.lbda, na.rm=T)
sd(dat.sites.stepps$stab.stepps.lbda - dat.sites.stepps$stab.lbda, na.rm=T)

# -----------


# -----------
# Quantitative comparisons between Biomass & Composition Stability
# -----------
dat.sites <- data.frame(lat=c(dat.sites.refab$lat, dat.sites.stepps$lat), 
                        lon=c(dat.sites.refab$lon, dat.sites.stepps$lon), 
                        dataset=c(rep("ReFAB", nrow(dat.sites.refab)), rep("STEPPS", nrow(dat.sites.stepps))),
                        stability=c(dat.sites.refab$stab.refab.1k, dat.sites.stepps$stab.stepps.1k),
                        variability=c(dat.sites.refab$var.refab.1k, dat.sites.stepps$var.stepps.1k),
                        H.prime = c(dat.sites.refab$H.prime.1k, dat.sites.stepps$H.prime.1k),
                        richness = c(dat.sites.refab$richness.1k, dat.sites.stepps$richness.1k),
                        dom.pft = c(paste(dat.sites.refab$dom.pft.1k), paste(dat.sites.stepps$dom.pft.1k))
                        )
dat.sites[dat.sites$dom.pft=="NA","dom.pft"] <- NA
summary(dat.sites)

write.csv(dat.sites, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Diversity_Data.csv"), row.names=F)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Map_STEPPS_Hprime.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites[,]) +
  facet_grid(dataset~.) +
  geom_point(aes(x=lon, y=lat, color=H.prime), size=2) +
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
  facet_grid(dataset~., scales="free_y") +
  geom_point(aes(x=H.prime, y=log(variability), color=dom.pft)) +
  # stat_smooth(data=dat.sites[dat.sites$dataset=="ReFAB",], aes(x=H.prime, y=stability, color=dataset, fill=dataset), method="lm") +
  # stat_smooth(data=dat.sites[dat.sites$dataset=="STEPPS",], aes(x=H.prime, y=stability, color=dataset, fill=dataset), method="lm") +
  theme_bw() +
  theme(legend.position="top")
dev.off()

# Comparing means -- is one more stable than the other?
t.test(dat.sites.refab$stab.refab.1k, dat.sites.refab$stab.stepps.1k, paired=T)
mean(dat.sites.refab$stab.refab.1k - dat.sites.refab$stab.stepps.1k, na.rm=T)
sd(dat.sites.refab$stab.refab.1k - dat.sites.refab$stab.stepps.1k, na.rm=T)

# Is biomass stability correlated with composition stability?
lm.comp <- lm(stab.refab.1k ~ stab.stepps.1k, data=dat.sites.refab)
summary(lm.comp)

# Is biomass stability correlated with diversity?
lm.diversity1 <- lm(stab.refab.1k ~ H.prime.1k, data=dat.sites.refab)
summary(lm.diversity1)

lm.diversity1b <- lm(stab.refab.1k ~ H.prime.1k + I(H.prime.1k^2), data=dat.sites.refab)
summary(lm.diversity1b)
dat.sites[dat.sites$dataset=="ReFAB" & !is.na(dat.sites$H.prime), "H.prime.quad"] <- predict(lm.diversity1b)

# Is composition stability correlated with Diversity
lm.diversity2 <- lm(stab.stepps.1k ~ H.prime.1k + I(H.prime.1k^2), data=dat.sites.stepps)
summary(lm.diversity2)

lm.diversity2b <- lm(stab.stepps.1k ~ H.prime.1k + I(H.prime.1k^2), data=dat.sites.stepps)
dat.sites[dat.sites$dataset=="STEPPS", "H.prime.quad"] <- predict(lm.diversity2b)
summary(lm.diversity2b)

# Testing residuals for normalcy
plot(predict(lm.diversity2b)~dat.sites.stepps$stab.stepps.1k[!is.na(dat.sites.stepps$H.prime.1k)]); abline(a=0, b=1, col="red")
plot(resid(lm.diversity2b)~predict(lm.diversity2b))
hist(resid(lm.diversity2b))

# Some summary figures
png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Data_Biomass_v_Composition.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites.refab) +
  geom_point(aes(x=stab.stepps.1k, y=var.refab.1k)) +
  stat_smooth(aes(x=stab.stepps.1k, y=var.refab.1k), method="lm") +
  theme_bw()
dev.off()


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Data_Stability_v_Diversity_ModFit_Quad.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites) +
  facet_wrap(~dataset, scales="free", ncol=1) +
  geom_point(aes(x=H.prime.quad, y=stability)) +
  stat_smooth(aes(x=H.prime.quad, y=stability), method="lm") +
  theme_bw()
dev.off()

# -----------

# -------------------------------------------
# -------------------------------------------