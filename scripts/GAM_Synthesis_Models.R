# -------------------------------------------
# Assessing the scales and magnitude of climate change variability over the past millennium in the MIP drivers
# Author: Christy Rollinson, crollinson@gmail.com

# 1. Stastical detection of significant change
#    - Use loess or TP regression splines
# 2. Comparison with paleoclimate reconstructions
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra); library(scales); library(grid)
library(mgcv)
library(plyr); library(parallel)

# Setting the path to this repository
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"

# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
path.gamm.func <- "~/Desktop/Research/R_Functions/"
# mip.utils <- "~/Dropbox/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

path.co2 <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/co2/paleon_annual_co2.nc"
# -------------------------------------------

# -------------------------------------------
# Define some useful variables
# -------------------------------------------
# Set up some time variables just to help with indexing
yrs <- 850:2010
mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)

us <- map_data("state")
# -------------------------------------------

# -------------------------------------------
# Comparing met & Ecosystem stability
# -------------------------------------------
stab.met   <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers.csv"))
stab.models <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models.csv"))

# Standardizing the derivatives relative to the mean
stab.met$tair.deriv.std <- stab.met$tair.deriv/mean(stab.met$tair.deriv)
stab.met$precip.deriv.std <- stab.met$precip.deriv/mean(stab.met$precip.deriv)
stab.met$pdsi.deriv.std <- stab.met$pdsi.deriv/mean(stab.met$pdsi.deriv)

for(mod in unique(stab.models$Model)){
  stab.models[stab.models$Model==mod,"bm.deriv.std"] <- stab.models[stab.models$Model==mod,"deriv.bm"]/mean(stab.models[stab.models$Model==mod,"deriv.bm"], na.rm=T)
  stab.models[stab.models$Model==mod,"npp.deriv.std"] <- stab.models[stab.models$Model==mod,"deriv.npp"]/mean(stab.models[stab.models$Model==mod,"deriv.npp"], na.rm=T)
}

stab.syn <- merge(stab.models, stab.met)
summary(stab.syn)

summary(stab.met)

plot.precip <- ggplot(data=stab.met) +
  # facet_wrap(~Model) +
  geom_tile(aes(x=lon, y=lat, fill=log(precip.deriv.std))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  # geom_point(data=stab.models, aes(x=lon, y=lat, color=bm.deriv.std)) +
  scale_fill_gradient2(name="stability", low = "blue2", high = "red2", mid = "gray80", midpoint = 0) + 
  # scale_color_gradient(low = "darkgreen", high="lightgreen") + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("Precipitation") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face=2))
  
plot.temp <- ggplot(data=stab.met) +
  # facet_wrap(~Model) +
  geom_tile(aes(x=lon, y=lat, fill=log(tair.deriv.std))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  # geom_point(data=stab.models, aes(x=lon, y=lat, color=bm.deriv.std)) +
  scale_fill_gradient2(name="stability", low = "blue2", high = "red2", mid = "gray80", midpoint = 0) + 
  # scale_color_gradient(low = "darkgreen", high="lightgreen") + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("Temperature") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face=2))



plot.ed.bm <- ggplot(data=stab.models[stab.models$Model=="ED2",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=log(bm.deriv.std)), size=2) +
  scale_color_gradient2(name="stability", low = "blue2", high = "red2", mid = "gray80", midpoint = 0) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("ED2") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.lpjg.bm <- ggplot(data=stab.models[stab.models$Model=="LPJ-GUESS",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=log(bm.deriv.std)), size=2) +
  scale_color_gradient2(name="stability", low = "blue2", high = "red2", mid = "gray80", midpoint = 0) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LPJ-GUESS") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.lpjw.bm <- ggplot(data=stab.models[stab.models$Model=="LPJ-WSL",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=log(bm.deriv.std)), size=2) +
  scale_color_gradient2(name="stability", low = "blue2", high = "red2", mid = "gray80", midpoint = 0) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LPJ-WSL") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.link.bm <- ggplot(data=stab.models[stab.models$Model=="LINKAGES",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=log(bm.deriv.std)), size=2) +
  scale_color_gradient2(name="stability", low = "blue2", high = "red2", mid = "gray80", midpoint = 0) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LINKAGES") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))



png(file.path(path.google, "Current Figures/Stability_GAMs", "Stability_Biomass_Region.png"), height=6, width=8, units="in", res=320)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
print(plot.temp   , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot.precip , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(plot.ed.bm  , vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(plot.link.bm, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(plot.lpjg.bm, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(plot.lpjw.bm, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()


plot.precip.nyr <- ggplot(data=stab.met) +
  # facet_wrap(~Model) +
  geom_tile(aes(x=lon, y=lat, fill=precip.nyr)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  # geom_point(data=stab.models, aes(x=lon, y=lat, color=bm.deriv.std)) +
  scale_fill_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80") + 
  # scale_color_gradient(low = "darkgreen", high="lightgreen") + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("Precipitation") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face=2))

plot.temp.nyr <- ggplot(data=stab.met) +
  # facet_wrap(~Model) +
  geom_tile(aes(x=lon, y=lat, fill=tair.nyr)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  # geom_point(data=stab.models, aes(x=lon, y=lat, color=bm.deriv.std)) +
  scale_fill_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80") + 
  # scale_color_gradient(low = "darkgreen", high="lightgreen") + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("Temperature") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face=2))



plot.ed.bm.nyr <- ggplot(data=stab.models[stab.models$Model=="ED2",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=bm.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("ED2") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.lpjg.bm.nyr <- ggplot(data=stab.models[stab.models$Model=="LPJ-GUESS",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=bm.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LPJ-GUESS") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.lpjw.bm.nyr <- ggplot(data=stab.models[stab.models$Model=="LPJ-WSL",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=bm.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LPJ-WSL") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.link.bm.nyr <- ggplot(data=stab.models[stab.models$Model=="LINKAGES",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=bm.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LINKAGES") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))



png(file.path(path.google, "Current Figures/Stability_GAMs", "Stability_Biomass_Region_Nyrs.png"), height=6, width=8, units="in", res=320)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
print(plot.temp.nyr   , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot.precip.nyr , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(plot.ed.bm.nyr  , vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(plot.link.bm.nyr, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(plot.lpjg.bm.nyr, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(plot.lpjw.bm.nyr, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()


plot.ed.npp.nyr <- ggplot(data=stab.models[stab.models$Model=="ED2",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=npp.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("ED2") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.lpjg.npp.nyr <- ggplot(data=stab.models[stab.models$Model=="LPJ-GUESS",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=npp.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LPJ-GUESS") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.lpjw.npp.nyr <- ggplot(data=stab.models[stab.models$Model=="LPJ-WSL",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=npp.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LPJ-WSL") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))

plot.link.npp.nyr <- ggplot(data=stab.models[stab.models$Model=="LINKAGES",]) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") +
  geom_point(aes(x=lon, y=lat, color=npp.nyr), size=2) +
  scale_color_gradient2(name="# Yrs", low = "blue2", high = "red2", mid = "gray80", limits=range(0,1000)) + 
  coord_equal(xlim=c(min(stab.met$lon)-0.25, max(stab.met$lon)+0.25), 
              ylim=c(min(stab.met$lat)-0.25, max(stab.met$lat)+0.25),
              expand=0) +
  ggtitle("LINKAGES") +
  theme_bw() + theme(plot.title=element_text(hjust=0.5))


png(file.path(path.google, "Current Figures/Stability_GAMs", "Stability_NPP_Region_Nyrs.png"), height=6, width=8, units="in", res=320)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
print(plot.temp.nyr   , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot.precip.nyr , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(plot.ed.npp.nyr  , vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(plot.link.npp.nyr, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(plot.lpjg.npp.nyr, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(plot.lpjw.npp.nyr, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()


ggplot(data=stab.models) +
  facet_wrap(~Model) +
  geom_histogram(aes(x=bm.deriv.std)) +
  geom_vline(xintercept=1, color="red") +
  theme_bw()

ggplot(data=stab.models) +
  facet_wrap(~Model) +
  geom_histogram(aes(x=log(npp.deriv.std))) +
  geom_vline(xintercept=0, color="red") +
  theme_bw()

ggplot(data=stab.met) +
  # facet_wrap(~Model) +
  geom_histogram(aes(x=log(pdsi.deriv.std))) +
  geom_vline(xintercept=0, color="red") +
  theme_bw()

ggplot(data=stab.met) +
  # facet_wrap(~Model) +
  geom_histogram(aes(x=log(tair.deriv.std))) +
  geom_vline(xintercept=0, color="red") +
  theme_bw()

ggplot(data=stab.met) +
  # facet_wrap(~Model) +
  geom_histogram(aes(x=log(precip.deriv.std))) +
  geom_vline(xintercept=0, color="red") +
  theme_bw()



ggplot(data=stab.syn, aes(x=log(pdsi.deriv.std), y=log(bm.deriv.std), color=Model, fill=Model)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()

pdsi.bm <- lm(log(bm.deriv.std) ~ log(pdsi.deriv.std)*Model - log(pdsi.deriv.std), data=stab.syn)
summary(pdsi.bm)


ggplot(data=stab.syn, aes(x=log(tair.deriv.std), y=log(bm.deriv.std), color=Model, fill=Model)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()

tair.bm <- lm(log(bm.deriv.std) ~ log(tair.deriv.std)*Model - log(tair.deriv.std), data=stab.syn)
summary(tair.bm)

ggplot(data=stab.syn, aes(x=log(tair.deriv.std), y=log(npp.deriv.std), color=Model, fill=Model)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()

tair.npp <- lm(log(npp.deriv.std) ~ log(tair.deriv.std)*Model - log(tair.deriv.std), data=stab.syn)
summary(tair.npp)

ggplot(data=stab.syn, aes(x=log(npp.deriv.std), y=log(bm.deriv.std), color=Model, fill=Model)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()

bm.npp <- lm(log(bm.deriv.std) ~ log(npp.deriv.std)*Model - log(npp.deriv.std), data=stab.syn)
summary(bm.npp)


ggplot(data=stab.syn, aes(x=log(pdsi.deriv), y=log(deriv.npp), color=Model, fill=Model)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()


ggplot(data=stab.syn, aes(x=log(precip.deriv), y=log(deriv.bm), color=Model, fill=Model)) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()


ggplot(data=stab.syn, aes(x=log(precip.deriv), y=log(pdsi.deriv))) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()
ggplot(data=stab.syn, aes(x=log(tair.deriv), y=log(pdsi.deriv))) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()
ggplot(data=stab.syn, aes(x=log(tair.deriv), y=log(precip.deriv))) +
  geom_point() +
  stat_smooth(method="lm") +
  theme_bw()




# -------------------------------------------

