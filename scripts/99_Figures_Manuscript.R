# -------------------------------------------
# Manuscript Figures in one place to make reproduction and manipulation easy

# Manuscript Figure List:
# 1. Methods: Dataset Maps -- Model region + STEPPS + ReFAB
# 2. Methods: GAM illustrative example
# 3. Results: Map of empirical vs. driver PDSI stability
# 4. Results: Scatterplot of ecosystem vs. climate stability (models + data)
#    4.a. Composition
#    4.b. Biomass (separate panels for Linkages & other data)
# 5. Results: Scatterplots by models (group by complexity/responses)
# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# Load Libaries & Set file paths
# -------------------------------------------
library(ggplot2); library(gridExtra); library(scales)
library(sp); library(raster)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"


# Set up a table with colors for models & data
# Using a CB palette from here: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
dat.colors <- data.frame(model=c("LBDA", "STEPPS", "ReFAB", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"),
                         type=c(rep("emipirical", 3), rep("model", 7)),
                         color=c(rep("#000000", 3), rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
dat.colors
# -------------------------------------------


# -------------------------------------------
# 1. Methods: Dataset Maps -- Model region + STEPPS + ReFAB
# -------------------------------------------
us <- map_data("state")

drivers.all <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
stepps <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_STEPPS.csv"))
refab <- read.csv(file.path(path.google, "Current Data", "refab.mean.slope.csv"))
models1 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))


png(file.path(path.google, "Manuscript/Figures", "1_Model_v_Data_Spatial_Extent.png"), height=5, width=10, units="in", res=320)
ggplot() +
  geom_tile(data=drivers.all, aes(x=lon, y=lat, fill=tair.yr.set-273.15)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray20") + 
  geom_point(data=stepps[stepps$pft=="OTHER.HARDWOOD",], aes(x=lon, y=lat, color="STEPPS"), shape=19, size=0.75, alpha=0.75) +
  geom_point(data=models1[models1$Model=="ED2",], aes(x=lon, y=lat, color="Models"), shape=19, size=1.25) +
  geom_point(data=refab, aes(x=lon, y=lat, color="ReFAB"), size=1.75) +
  scale_color_manual(name="Datasets", values=c("gray40", "black", "darkgoldenrod3")) +
  scale_fill_gradient2(name="Mean Temp \n(1800-1850)", low = "blue", high = "red", mid = "white", midpoint = mean(drivers.all$tair.yr.set-273.15, na.rm=T)) +
  labs(x="Longitude", y="Latitude") +
  # guides(color=guide_legend(aes.overide=list(size=10)))+
  # guides(shape=guide_legend(aes.override=list(size=20))) +
  coord_equal(xlim=range(drivers.all$lon), ylim=range(drivers.all$lat), expand=0) +
  theme(legend.position="top",
        legend.key=element_blank(),
        panel.background = element_rect(fill="gray80"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
dev.off()

rm(drivers.all, stepps, refab, models1)
# -------------------------------------------

# -------------------------------------------
# 2. Methods: GAM illustrative example
# -------------------------------------------
load(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Example_ED2_HarvardForest.RData"))

df.ed <- data.frame(Year=850:1850, Biomass = ed.list$bm, 
                    gam.mean = ed.out$bm$gam.post$ci$mean, 
                    gam.lwr = ed.out$bm$gam.post$ci$lwr,
                    gam.upr = ed.out$bm$gam.post$ci$upr)
summary(df.ed)


png(file.path(path.google, "Manuscript/Figures", "2_HF_TimeSeries_Biomass_ED2_GAM.png"), height=9, width=10, units="in", res=220)
ggplot(data=df.ed) +
  geom_line(aes(x=Year, y=Biomass), size=0.5, color="black") +
  geom_ribbon(aes(x=Year, ymin=gam.lwr, ymax=gam.upr), alpha=0.5, fill="blue") +
  geom_line(aes(x=Year, y=gam.mean), size=5, color="blue") +
  geom_point(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, y=gam.mean), size=10, color="blue4") +
  geom_linerange(data=df.ed[df.ed$Year %in% seq(850, 1850, 100),], aes(x=Year, ymin=gam.lwr, ymax=gam.upr), size=5, color="blue4") +
  scale_y_continuous(limits=range(df.ed$Biomass), name="Biomass") +
  theme_bw() +
  theme(axis.text=element_text(size=rel(2)),
        axis.title=element_text(size=rel(2.5)))
dev.off()
# -------------------------------------------

# -------------------------------------------
# 3. Results: Map of empirical vs. driver PDSI stability
# -------------------------------------------

# png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_PDSI_diff_abs.png"), height=6, width=8, units="in", res=320)
ggplot(data=dat.all[dat.all$class=="climate" & dat.all$var=="pdsi",]) +
  facet_wrap(~model, ncol=1) +
  geom_tile(aes(x=lon, y=lat, fill=log(diff.abs))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(dat.all$lon), ylim=range(dat.all$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="climate" & dat.all$var=="pdsi","diff.abs"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray50"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("PDSI: empirical v. model; NOTE: LOG SCALE!!")
# dev.off()
# -------------------------------------------


# -------------------------------------------
# 4. Results: Scatterplot of ecosystem vs. climate stability (models + data)
#    4.a. Composition
#    4.b. Biomass (separate panels for Linkages & other data)
# -------------------------------------------
# -------------------------------------------


# -------------------------------------------
# 5. Results: Scatterplots by models (group by complexity/responses)
# -------------------------------------------
# -------------------------------------------
