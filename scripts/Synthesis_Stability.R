# --------------------------------------------
# Script Description
# --------------------------------------------
# Purpose: 
#
# Analysis Details:
#
#
# Workflow
# 0. Define file paths etc
# 1. Read in & format datasets
# --------------------------------------------

rm(list=ls())

# --------------------------------------------
# 0. Define file paths etc
# --------------------------------------------
# Path to github repository/working directory
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

# Path to where data are; lets just pull straight from the Google Drive folder
path.data <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/Current Data/Stability_Index/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"

# --------------------------------------------


# --------------------------------------------
# 1. Read in and merge datasets
# --------------------------------------------
# load in the paleon domain info;
# This got generated using domain_environment_extraction.R
paleon <- read.csv(file.path(path.repo, "data/paleon_models_environment_master.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
summary(paleon)

stab.clim <- read.csv(file.path(path.data, "PalEON_ClimateStability.csv"))
stab.nada <- read.csv(file.path(path.data, "NADA_Stability.csv"))
stab.lbda <- read.csv(file.path(path.data, "LBDA_Stability.csv"))
stab.fcomp <- read.csv("data/si_fcomp_multi.csv")
stab.refab <- read.csv(file.path(path.data, "refab.stability.df.csv"))

# biomass models all saved by species
stab.bm.ed <- read.csv(file.path(path.data, "ed.stability.csv"))
stab.bm.link <- read.csv(file.path(path.data, "linkages.stability.csv"))
stab.bm.lpjg <- read.csv(file.path(path.data, "lpj.guess.stability.csv"))
stab.bm.lpjw <- read.csv(file.path(path.data, "lpj.wsl.stability.csv"))
stab.bm.trif <- read.csv(file.path(path.data, "triffid.stability.csv"))

summary(stab.fcomp)
summary(stab.clim)
summary(stab.nada)
summary(stab.lbda)
summary(stab.refab)

# Rename some columns so we can do some merging:
names(paleon)[3] <- "Site"

stab.clim <- stab.clim[,c(1:4, 26:ncol(stab.clim))]
names(stab.clim)[3] <- "Site"
names(stab.clim)[5:ncol(stab.clim)] <- substr(names(stab.clim)[5:ncol(stab.clim)], 6, nchar(names(stab.clim)[5:ncol(stab.clim)]))
summary(stab.clim)

names(stab.nada)[3:5] <- c("nada.end", "nada.start", "nada")
names(stab.lbda)[3:5] <- c("lbda.start", "lbda.end", "lbda")
# stab.drought[,c("lat", "lon")] <- stab.drought[,c("lat", "lon")]-0.25
summary(stab.nada)
summary(stab.lbda)

names(stab.refab) <- c("Site", "biomass", "lat", "lon")
summary(stab.refab)

summary(stab.fcomp)

# --------------------------------------------
# 2. Exploratory Comparisions
#    A. Climate
#    B. Composition
#    C. Biomass
# --------------------------------------------
# Checking out spatial alignment of products
library(ggplot2)
us <- map_data("state")

ggplot() + 
  geom_tile(data=stab.lbda, aes(x=lon, y=lat, fill=lbda))+
  geom_raster(data=stab.clim, aes(x=lon, y=lat, fill=pdsi.ann), fill="red", alpha=0.5) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(stab.lbda$lon), ylim=range(stab.lbda$lat))


# The PalEON grid is at the same resolution as the Living, Blended Atlas, so those comparisons are easy
summary(stab.lbda)
summary(stab.clim)

stab.lbda2 <- merge(stab.lbda, stab.clim)
summary(stab.lbda2)

ggplot() + 
  geom_raster(data=stab.clim[stab.clim$pdsi.ann<quantile(stab.clim$pdsi.ann, 1-3/nrow(stab.clim)),], aes(x=lon, y=lat, fill=pdsi.ann)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(stab.lbda$lon), ylim=range(stab.lbda$lat))

# Remove our 3 weirdos
stab.lbda2[stab.lbda2$pdsi.ann>quantile(stab.lbda2$pdsi.ann, 1-3/nrow(stab.lbda2), na.rm=T), c("pdsi.ann", "pdsi.ann.cent","pdsi.lbda", "pdsi.lbda.nyr")] <- NA


# stab.lbda2$lbda.std <- stab.lbda2$lbda/(stab.lbda2$n.yrs)
# stab.lbda2$pdsi.std <- stab.lbda2$pdsi.lbda/(stab.lbda2$n.yrs)

stab.lbda2$lbda.std <- (stab.lbda2$lbda/mean(stab.lbda2$lbda, na.rm=T))
stab.lbda2$pdsi.lbda.std <- (stab.lbda2$pdsi.lbda/mean(stab.lbda2$pdsi.lbda, na.rm=T))

stab.lbda2$tair.std <- (stab.lbda2$tair.ann/mean(stab.lbda2$tair.ann, na.rm=T))
stab.lbda2$precip.std <- (stab.lbda2$precip.ann/mean(stab.lbda2$precip.ann, na.rm=T))
stab.lbda2$pdsi.std <- (stab.lbda2$pdsi.ann/mean(stab.lbda2$pdsi.ann, na.rm=T))

# paleon$tair.std <- (paleon$tair.yr.set/mean(paleon$tair.yr.set, na.rm=T))
# paleon$precip.std <- (paleon$precip.yr.set/mean(paleon$precip.yr.set, na.rm=T))
# 
# met.comparison.raw <- stack(paleon[,c("tair.std", "precip.std")])
# met.comparison.raw[,c("Site", "lon", "lat", "umw")] <- paleon[,c("Site", "lon", "lat", "umw")]
# 
# ggplot(data=met.comparison.raw[,]) + 
#   facet_wrap(~ind, ncol=1) +
#   geom_tile(aes(x=lon, y=lat, fill=values))+
#   geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1) + 
#   coord_equal(xlim=range(stab.lbda$lon), ylim=range(stab.lbda$lat), expand=0)
# 

# Stacking stuff to make a handy figure
library(car)
stab.comparison <- stack(stab.lbda2[,c("lbda.std", "pdsi.lbda.std", "tair.std", "precip.std", "pdsi.std")])
stab.comparison[,c("Site", "lon", "lat", "umw")] <- stab.lbda2[,c("Site", "lon", "lat", "umw")]
stab.comparison$ind <- recode(stab.comparison$ind, "'lbda.std'='LBDA'; 'pdsi.lbda.std'='PDSI (LBDA Overlap)'; 'tair.std'='Tair'; 'precip.std'='Precip'; 'pdsi.std'='PDSI (all)'")
stab.comparison$ind <- factor(stab.comparison$ind, levels=c("Tair", "Precip", "PDSI (all)", "PDSI (LBDA Overlap)", "LBDA"))
summary(stab.comparison)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "MetStability_Relative.png"), height=3, width=6, unit="in", res=320)
ggplot(data=stab.comparison[,]) + 
  facet_wrap(~ind, ncol=3) +
  geom_tile(aes(x=lon, y=lat, fill=values))+
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1) + 
  coord_equal(xlim=range(stab.lbda$lon), ylim=range(stab.lbda$lat), expand=0)
dev.off()



ggplot(data=stab.lbda2) +
  geom_point(aes(x=lbda.std, y=pdsi.std)) +
  stat_smooth(aes(x=lbda.std, y=pdsi.std), method="lm")

ggplot(data=stab.lbda2[stab.lbda2$pdsi.std2<5,]) +
  geom_point(aes(x=lbda.std2, y=pdsi.std2)) +
  stat_smooth(aes(x=lbda.std2, y=pdsi.std2), method="lm")

ggplot(data=stab.lbda2) +
  geom_point(aes(x=lbda.std, y=tair.jja/1161)) +
  stat_smooth(aes(x=lbda.std, y=tair.jja/1161), method="lm")

ggplot(data=stab.lbda2) +
  geom_point(aes(x=lbda.std, y=precip.ann/1161)) +
  stat_smooth(aes(x=lbda.std, y=precip.ann/1161), method="lm")

lbda.pdsi.lm <- lm(pdsi.std ~ lbda.std, data=stab.lbda2 )
summary(lbda.pdsi.lm)

lbda.tair.lm <- lm(tair.jja/1161 ~ lbda.std, data=stab.lbda2 )
summary(lbda.tair.lm)

lbda.ppt.lm <- lm(precip.jja/1161 ~ lbda.std, data=stab.lbda2 )
summary(lbda.ppt.lm)

par(mfrow=c(3,1))
hist(stab.lbda2$lbda)
hist(stab.lbda2$lbda.std)
hist(stab.lbda2$pdsi.ann/1161)


# Finding out how to re-center points to line up with NADA
ggplot() + 
  geom_tile(data=stab.nada, aes(x=lon, y=lat, fill=nada))+
  geom_raster(data=stab.clim, aes(x=lon, y=lat, fill=pdsi.ann)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(stab.nada$lon), ylim=range(stab.nada$lat))


# --------------------------------------------
