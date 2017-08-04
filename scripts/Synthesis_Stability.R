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

library(ggplot2)
us <- map_data("state")
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

names(stab.refab) <- c("Site", "biomass", "lon", "lat")
summary(stab.refab)

# Merge models into 1
stab.bm.ed$Model   <- "ed2"
stab.bm.link$Model <- "linkages"
stab.bm.lpjg$Model <- "lpj-guess"
stab.bm.lpjw$Model <- "lpj-wsl"
stab.bm.trif$Model <- "triffid"

stab.bm.model <- rbind(stab.bm.ed, stab.bm.link, stab.bm.lpjg, stab.bm.lpjw, stab.bm.trif)
stab.bm.model$Model <- as.factor(stab.bm.model$Model)
stab.bm.model$Site  <- as.factor(paste0("lat", stab.bm.model$lat, "lon", stab.bm.model$lon))
names(stab.bm.model)[1] <- "biomass"
summary(stab.bm.model)


summary(stab.fcomp)

# --------------------------------------------
# 2. Exploratory Comparisions
#    A. Climate
#    B. Composition
#    C. Biomass
# --------------------------------------------
# Checking out spatial alignment of products

# ggplot() + 
#   geom_tile(data=stab.lbda, aes(x=lon, y=lat, fill=lbda))+
#   geom_raster(data=stab.clim, aes(x=lon, y=lat, fill=pdsi.ann), fill="red", alpha=0.5) +
#   geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
#   coord_equal(xlim=range(stab.lbda$lon), ylim=range(stab.lbda$lat))
# 

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

png(file.path(path.google, "Current Figures/Stability_Synthesis", "MetStability_Relative.png"), height=4, width=10, unit="in", res=320)
ggplot(data=stab.comparison[,]) + 
  facet_wrap(~ind, ncol=3) +
  geom_tile(aes(x=lon, y=lat, fill=values))+
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1) + 
  ggtitle("PDSI Stability Index: Reconstruction vs. Drivers") +
  coord_equal(xlim=range(stab.lbda$lon), ylim=range(stab.lbda$lat), expand=0)
dev.off()

# Doing some exploratory graphing of stability in met variables
ggplot(data=stab.lbda2, aes(x=lbda, y=pdsi.lbda)) +
  geom_point() +
  stat_smooth(method="lm")

ggplot(data=stab.lbda2, aes(x=lbda.std, y=pdsi.lbda.std)) +
  geom_point() +
  stat_smooth(method="lm")

ggplot(data=stab.lbda2, aes(x=pdsi.std, y=pdsi.lbda.std)) +
  geom_point() +
  stat_smooth(method="lm")

ggplot(data=stab.lbda2, aes(x=precip.std, y=pdsi.std)) +
  geom_point() +
  stat_smooth(method="lm")

ggplot(data=stab.lbda2, aes(x=tair.std, y=pdsi.std)) +
  geom_point() +
  stat_smooth(method="lm")


# Looking at the correlations between empirical & model
lbda.pdsi.lm <- lm(pdsi.lbda.std ~ lbda.std, data=stab.lbda2 )
summary(lbda.pdsi.lm) # No correlation (at all!) between model & driver instability

# Making sure that the poor coorealtion with osberved is likely to hold true for logner time periods as well
pdsi.lm <- lm(pdsi.lbda.std ~ pdsi.std, data=stab.lbda2 )
summary(pdsi.lm) # Not perfect, but tight correlation between the LBDA overlap stability & with the full time
# --------------------------------------------


# --------------------------------------------
# Now to compare ecosystem stability with models & data against each other & with met 
# --------------------------------------------
# Lets start with biomass since that's only 1-dimension
# biomass models all saved by species
summary(stab.bm.model)
summary(stab.refab)

# Standardizing the biomass stability to itself
stab.refab$bm.std <- stab.refab$biomass/mean(stab.refab$biomass)

for(mod in unique(stab.bm.model$Model)){
  stab.bm.model[stab.bm.model$Model==mod, "bm.std"] <- stab.bm.model[stab.bm.model$Model==mod, "biomass"]/mean(stab.bm.model[stab.bm.model$Model==mod, "biomass"], na.rm=T)
}
summary(stab.bm.model)

summary(stab.lbda2)

# Making a pdsi data frame that will correspond to the biomass stability one
pdsi.stab <- data.frame(stab.lbda2[,c("Site", "lon", "lat", "umw")], 
                        Model=rep(unique(stab.bm.model$Model), each=nrow(stab.lbda2)),
                        pdsi.stab=stab.lbda2$pdsi.lbda.std)
pdsi.stab <- rbind(pdsi.stab, data.frame(stab.lbda2[,c("Site", "lon", "lat", "umw")], 
                                         Model=rep("refab", each=nrow(stab.lbda2)),
                                         pdsi.stab=stab.lbda2$lbda.std))

pdsi.stab <- pdsi.stab[pdsi.stab$lon>=min(stab.lbda2[stab.lbda2$umw=="y", "lon"], na.rm=T)-1 & 
                         pdsi.stab$lon<=max(stab.lbda2[stab.lbda2$umw=="y", "lon"], na.rm=T)+1 & 
                         pdsi.stab$lat>=min(stab.lbda2[stab.lbda2$umw=="y", "lat"], na.rm=T)-1 & 
                         pdsi.stab$lat<=max(stab.lbda2[stab.lbda2$umw=="y", "lat"], na.rm=T)+1   , ]

summary(pdsi.stab)


# Merging empirical and model biomass stability
stab.refab$Model <- as.factor("refab")
stab.refab$umw <- as.factor("y")

stab.bm.model <- merge(stab.bm.model, paleon[,c("Site", "lon", "lat", "umw")])

stab.bm <- merge(stab.bm.model, stab.refab, all=T)
stab.bm <- stab.bm[stab.bm$umw=="y",]

# Since we're log-transforming thing, make 0 something really tiny
stab.bm[stab.bm$bm.std==0 & !is.na(stab.bm$bm.std), "bm.std"] <- 1e-4 

# put refab first so it's always our reference 
stab.bm$Model <- factor(stab.bm$Model, levels=c("refab", "ed2", "linkages", "lpj-guess", "lpj-wsl", "triffid"))
summary(stab.bm)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "BiomassStability_Maps.png"), height=6, width=8, units = "in", res=320)
ggplot(data=stab.bm) +
  facet_wrap(~Model, ncol=2) +
  geom_tile(data=pdsi.stab, aes(x=lon, y=lat, fill=pdsi.stab)) +
  geom_point(data=stab.bm, aes(x=lon, y=lat, color=log(bm.std)), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group), color="black") + 
  scale_fill_gradient2(low = "blue2", high = "red2", mid = "gray80", midpoint = 1) + 
  scale_color_gradient(low = "darkgreen", high="lightgreen") + 
  coord_equal(xlim=range(pdsi.stab$lon), 
              ylim=c(min(pdsi.stab$lat), max(pdsi.stab$lat)+0.25),
              expand=0) +
  ggtitle("Stability Index: PDSI vs. Biomass") +
  # guides(fill=guide_legend(aes.overide=list(alpha=0.5))) +
  theme(panel.background=element_blank())
dev.off()


# Extracting point information for each grid cell
for(i in 1:nrow(stab.bm)){
  pdsi.ind <- which(pdsi.stab$Model==stab.bm$Model[i] & 
                      pdsi.stab$lat-0.25<stab.bm$lat[i] & pdsi.stab$lat+0.25>=stab.bm$lat[i] &
                      pdsi.stab$lon-0.25<stab.bm$lon[i] & pdsi.stab$lon+0.25>=stab.bm$lon[i])
  
  if(length(pdsi.ind)==0) next
  
  stab.bm[i,"pdsi.stab"] <- pdsi.stab[pdsi.ind, "pdsi.stab"]
}
summary(stab.bm)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "BiomassStability_Maps.png"), height=6, width=8, units = "in", res=320)
ggplot(data=stab.bm) + 
  geom_point(aes(x=pdsi.stab, y=log(bm.std), color=Model)) +
  stat_smooth(aes(x=pdsi.stab, y=log(bm.std), color=Model, fill=Model), method="lm")
dev.off()


summary(stab.bm)

# Testing things with a linear model
bm.stab.lm <- lm(log(bm.std) ~ pdsi.stab*Model, data=stab.bm[complete.cases(stab.bm),])
summary(bm.stab.lm) # No relationships between environmental stability & biomass; this 
anova(bm.stab.lm)

# --------------------------------------------
