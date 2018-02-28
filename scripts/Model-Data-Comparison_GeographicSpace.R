# -------------------------------------------
# Model-Data comparisons: 1. Comparisons of Model & Data Stability in Geographic Space

# Here we will:
# 1. Look at patterns of stability in composition, biomass & climate in geographic space
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


# Set up a table with colors for models & data
# Using a CB palette from here: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
dat.colors <- data.frame(model=c("LBDA", "STEPPS", "ReFAB", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"),
                         type=c(rep("emipirical", 3), rep("model", 7)),
                         color=c(rep("#000000", 3), rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
dat.colors
# -------------------------------------------

# -------------------------------------------
# Read in & align data
# -------------------------------------------
# -----------
# Empirical Drought Reconstruction (LBDA, annually-resolved)
# -----------
lbda <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_LBDA_100.csv"))
lbda$model <- as.factor("LBDA")
lbda$class <- as.factor("climate")
lbda$var <- as.factor("pdsi")
lbda$type <- as.factor("empirical")
lbda$resolution <- as.factor("annual")
names(lbda)[which(names(lbda)=="n.yrs.sig")] <- "n.sig"
summary(lbda)

us <- map_data("state")

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
# Model Drivers
# -----------
drivers.all <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
drivers <- stack(drivers.all[,c("pdsi.deriv", "tair.deriv", "precip.deriv")])
names(drivers) <- c("deriv.abs", "var")
drivers[,"diff.abs"] <- stack(drivers.all[,c("pdsi.diff", "tair.diff", "precip.diff")])[,1]
drivers$var <- as.factor(unlist(lapply(stringr::str_split(drivers$var, "[.]"), function(x) {x[1]})))
drivers$n.sig <- stack(drivers.all[,c("pdsi.nyr", "tair.nyr", "precip.nyr")])[,1]
drivers$fract.sig <- drivers$n.sig/1000
drivers[,c("lon", "lat")] <- drivers.all[,c("lon", "lat")]
drivers$type <- as.factor("model")
drivers$resolution <- as.factor("annual")
drivers$var <- car::recode(drivers$var, "'tair'='temp'") # Since people have a hard time figuring out "Tair"
drivers$class <- as.factor("climate")
drivers$model <- as.factor("drivers")
drivers <- drivers[,c("lon", "lat", "model", "class", "var", "type", "resolution", "diff.abs", "deriv.abs", "n.sig", "fract.sig")]
summary(drivers)

# Load in the PDSI calc that's match in temporal extent to LBDA
pdsi2 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_PDSI_Drivers_LBDA_time_100.csv"))
pdsi2$var <- as.factor("pdsi")
pdsi2$type <- as.factor("model")
pdsi2$resolution <- as.factor("annual")
pdsi2$class <- as.factor("climate")
pdsi2$model <- as.factor("drivers-modified")
names(pdsi2)[which(names(pdsi2)=="n.yrs.sig")] <- "n.sig"
summary(pdsi2)

drivers <- merge(drivers, pdsi2, all=T)
summary(drivers)
# -----------

# -----------
# STEPPS (empirical composition, centennially-resolved)
# -----------
stepps <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_STEPPS.csv"))

# Need to convert this to latlon
# albers <- sp::CRS("+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") # Define the spatial projection
# stepps <- sp::SpatialPointsDataFrame(coords=stepps[,c("x", "y")], data=stepps, proj4string = albers) # Create a spatial object
# stepps <- sp::spTransform(stepps, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # spatial projection transformation
# stepps <- data.frame(stepps) # Pull out the data frame
# names(stepps)[which(names(stepps) %in% c("x.1", "y.1"))] <- c("lon", "lat")

# Change names to match up with drivers
names(stepps)[which(names(stepps)=="sig")] <- c("n.sig")
names(stepps)[which(names(stepps)=="deriv.abs")] <- c("diff.abs")
stepps$deriv.abs <- NA

stepps$fract.sig <- stepps$n.sig/10
stepps$model <- as.factor("STEPPS")
stepps$class <- as.factor("composition")
stepps$var <- stepps$pft
stepps$type <- as.factor("empirical")
stepps$resolution <- as.factor("centennial")
summary(stepps)
# -----------

# -----------
# ReFAB (empirical biomass, centennially-resolved)
# -----------
refab.means <- read.csv(file.path(path.google, "Current Data", "biomass.means.csv"))
names(refab.means)[4] <- "value" 
summary(refab.means)

refab <- read.csv(file.path(path.google, "Current Data", "refab.mean.slope.csv"))
names(refab) <- c("X", "lat", "lon", "diff.mean", "diff.abs", "n.sig")
refab$deriv.abs <- NA

refab <- merge(refab, refab.means[,c("lon", "lat", "value")])
summary(refab)

# Doing a unit correction
refab[,c("value", "diff.mean", "diff.abs")] <- refab[,c("value", "diff.mean", "diff.abs")]*0.1/2 # Native units = Mg/Ha biomass

refab$fract.sig <- refab$n.sig/10
refab$model <- as.factor("ReFAB")
refab$class <- as.factor("biomass")
refab$var <- as.factor("biomass")
refab$type <- as.factor("empirical")
refab$resolution <- as.factor("centennial")
summary(refab)
# -----------

# -----------
# Model Output - Biomass (annually-resovled)
# -----------
models1 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))
models1 <- models1[,c("lon", "lat", "Model", "mean.bm", "diff.bm", "deriv.bm", "nyr.bm")] # Add in composition once you do it
names(models1) <- c("lon", "lat", "model", "value", "diff.abs", "deriv.abs", "n.sig")
models1$fract.sig <- models1$n.sig/1000
models1[models1$value<=0 & !is.na(models1$value),c("value", "diff.abs", "deriv.abs", "n.sig", "fract.sig")] <- NA

models1$class <- as.factor("biomass")
models1$var <- as.factor("biomass")
models1$type <- as.factor("model")
models1$resolution <- as.factor("annual")

summary(models1)
# -----------

# -----------
# Model Ouput - Fcomp (annually-resolved)
# -----------
load(file.path(path.data, "PalEON_siteInfo_all.RData"))
paleon$cell <- 1:nrow(paleon)
summary(paleon)

fcomp.diff  <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_FCOMP_DIFF.csv"))
fcomp.deriv <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_FCOMP_DERIV.csv"))

# Fixing some column names
names(fcomp.diff)[which(names(fcomp.diff) %in% c("deriv.abs", "deriv.mean"))] <- c("diff.abs", "diff.mean")
summary(fcomp.diff)

names(fcomp.deriv)[which(names(fcomp.deriv)=="n.yrs.sig")] <- "n.sig"
summary(fcomp.deriv)

# Adding lat/lon into the fcomp data frames
fcomp.diff  <- merge(fcomp.diff , paleon[,c("cell", "lon", "lat")], all.x=T)
fcomp.deriv <- merge(fcomp.deriv, paleon[,c("cell", "lon", "lat")], all.x=T)

# Merging the two fcomp data frames together
fcomp <- merge(fcomp.diff [,c("lon", "lat", "pft", "model", "value", "diff.abs" , "diff.mean" )], 
               fcomp.deriv[,c("lon", "lat", "pft", "model", "value", "deriv.abs", "deriv.mean", "n.sig", "fract.sig")])
summary(fcomp)

fcomp$var <- fcomp$pft
fcomp$class <- as.factor("composition")
fcomp$type <- as.factor("model")
fcomp$resolution <- as.factor("annual")
summary(fcomp)
# -----------

# -----------
# Model Output (centennially-resolved)
# # NOTE: Derivative needs to be the mean ABSOLUTE VALUE!  NEEDS TO BE RE-RUN!!
# # NOTE: This isn't accurate and I think we'll just roll with what I'd already done for biomass
# -----------
# models <- c("ed", "lpj", "lpj.wsl", "triffid", "linkages")
# models.bm.100 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_Biomass_Centennial.csv"))
# 
# models.bm <- stack(models.bm.100[,c(paste0(models, ".mean.deriv"))])
# names(models.bm) <- c("deriv.abs", "model")
# models.bm$deriv.abs <- abs(models.bm$deriv.abs) # NOTE: SHOULDN"T HAVE TO DO THIS! FIX IN FUTURE VERSIONS
# models.bm$model <- as.factor(unlist(lapply(stringr::str_split(models.bm$model, "[.]"), function(x) {paste(x[1:(length(x)-2)], collapse=".")})))
# models.bm$n.sig <- stack(models.bm.100[,c(paste0(models, ".nsig"))])[,1]
# models.bm$fract.sig <- models.bm$n.sig/10
# 
# models.bm[,c("lat", "lon")] <- models.bm.100[,c("X.1", "X.2")]
# 
# models.bm$model <- car::recode(models.bm$model, "'ed'='ED2'; 'lpj'='LPJ-GUESS'; 'lpj.wsl'='LPJ-WSL'; 'linkages'='LINKAGES'; 'triffid'='TRIFFID'")
# 
# models.bm$class <- as.factor("biomass")
# models.bm$var <- as.factor("biomass")
# models.bm$type <- as.factor("model")
# models.bm$resolution <- as.factor("centennial")
# 
# summary(models.bm)
# -----------


# Putting everything into 1 data frame to see how this goes
# This probably makes the most sense for mapping at least
cols.bind <- c("lon", "lat", "model", "class", "var", "type", "resolution", "diff.abs", "deriv.abs", "n.sig", "fract.sig")

dat.all <- rbind(lbda     [,cols.bind], 
                 drivers  [,cols.bind], 
                 stepps   [,cols.bind], 
                 refab    [,cols.bind], 
                 models1  [,cols.bind], 
                 # models.bm[,cols.bind],
                 fcomp    [,cols.bind])
dat.all <- dat.all[!is.na(dat.all$lon),]
summary(dat.all)

# Getting a relative stability so that we can ignore actual value and just look at spatial patterning
for(mod in unique(dat.all$model)){
  for(var in unique(dat.all$var[dat.all$model==mod])){
    for(res in unique(dat.all$resolution[dat.all$model==mod & dat.all$var==var])){
      deriv.subset <- dat.all[dat.all$model==mod & dat.all$var==var & dat.all$resolution==res, "deriv.abs"]
      dat.all[dat.all$model==mod & dat.all$var==var & dat.all$resolution==res, "deriv.rel"] <- deriv.subset/mean(deriv.subset, na.rm=T)
      
      diff.subset <- dat.all[dat.all$model==mod & dat.all$var==var & dat.all$resolution==res, "diff.abs"]
      dat.all[dat.all$model==mod & dat.all$var==var & dat.all$resolution==res, "diff.rel"] <- diff.subset/mean(diff.subset, na.rm=T)
      
    }
  }
}

dat.colors$model <- factor(dat.colors$model, levels(dat.all$model))
dat.colors[,] <- dat.colors[order(dat.colors$model),]
# quick comparison of deriv vs. diff for all variables & dataset
diff.deriv <- lm(diff.abs ~ deriv.abs, data=dat.all)
summary(diff.deriv)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Diff_v_Deriv_AllVars.png"), height=8, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$class %in% c("climate", "biomass"),]) +
  facet_wrap(~var, scales="free") +
  geom_point(aes(x=diff.abs, y=deriv.abs, color=model), alpha=0.2, size=0.5) +
  geom_abline(intercept=0, slope=1, color="black") +
  stat_smooth(aes(x=diff.abs, y=deriv.abs, color=model), method="lm", se=F) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(dat.all[dat.all$class %in% c("climate", "biomass"),"model"]),"color"])) +
  theme_bw()
dev.off()

# -------------------------------------------


# -------------------------------------------
# Exploratory graphing of datasets
# -------------------------------------------
us <- map_data("state")


# -----------
# Maping Datasets for comparison
# -----------
png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Spatial_Extent.png"), height=5, width=10, units="in", res=320)
ggplot() +
  geom_tile(data=drivers.all, aes(x=lon, y=lat, fill=tair.yr.set-273.15)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  geom_point(data=stepps[stepps$pft=="OTHER.HARDWOOD",], aes(x=lon, y=lat, color="STEPPS"), shape=19, size=0.75, alpha=0.75) +
  geom_point(data=models1[models1$model=="ED2",], aes(x=lon, y=lat, color="Models"), shape=19, size=1.25) +
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
# -----------


# -----------
# Mapping Met
# -----------
# Looking explicitly at model & driver PDSI -- Note the Log scale indiciating the orders of magnitude difference
png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_PDSI_diff_abs.png"), height=6, width=8, units="in", res=320)
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
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_PDSI_deriv_abs.png"), height=6, width=8, units="in", res=320)
ggplot(data=dat.all[dat.all$class=="climate" & dat.all$var=="pdsi",]) +
  facet_wrap(~model, ncol=1) +
  geom_tile(aes(x=lon, y=lat, fill=log(deriv.abs))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(dat.all$lon), ylim=range(dat.all$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="climate" & dat.all$var=="pdsi","deriv.abs"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray50"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("PDSI: empirical v. model; NOTE: LOG SCALE!!")
dev.off()


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_PDSI_sig_fract.png"), height=6, width=8, units="in", res=320)
ggplot(data=dat.all[dat.all$class=="climate" & dat.all$var=="pdsi",]) +
  facet_wrap(~model, ncol=2) +
  geom_tile(aes(x=lon, y=lat, fill=fract.sig)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(dat.all$lon), ylim=range(dat.all$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(dat.all[dat.all$class=="climate" & dat.all$var=="pdsi","fract.sig"], na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray50"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("PDSI: empirical v. model")
dev.off()


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Climate_diff_rel.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$class=="climate" ,]) +
  facet_grid(type~var) +
  geom_tile(aes(x=lon, y=lat, fill=log(diff.rel))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(dat.all$lon), ylim=range(dat.all$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="climate","diff.rel"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.background = element_rect(fill="gray50"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Relative Climate Stability: empirical v. model (850 - 1850 A.D.); NOTE: LOG SCALE!!")
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Climate_deriv_rel.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$class=="climate" & !dat.all$model=="drivers-modified",]) +
  facet_grid(type~var) +
  geom_tile(aes(x=lon, y=lat, fill=log(deriv.rel))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(dat.all$lon), ylim=range(dat.all$lat), expand=0) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="climate","deriv.rel"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.background = element_rect(fill="gray50"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Relative Climate Stability: empirical v. model (850 - 1850 A.D.); NOTE: LOG SCALE!!")
dev.off()
# -----------


# -----------
# Mapping Biomass
# -----------
lon.range <- c(min(dat.all[dat.all$model=="ReFAB","lon"])-1, max(dat.all[dat.all$model=="ReFAB","lon"])+1)
lat.range <- c(min(dat.all[dat.all$model=="ReFAB","lat"])-1, max(dat.all[dat.all$model=="ReFAB","lat"])+1.5)
lat.lon.ind <- dat.all$lat>=lat.range[1] & dat.all$lat<=lat.range[2] & dat.all$lon>=lon.range[1] & dat.all$lon<=lon.range[2]

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_diff_abs_century_full.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$var=="biomass" & !is.na(dat.all$diff.abs),]) +
  facet_wrap(~model) +
  geom_point(aes(x=lon, y=lat, color=log(diff.abs))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
  coord_equal(xlim=range(dat.all$lon), 
              ylim=range(dat.all$lat), 
              expand=0) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="biomass" & dat.all$resolution=="centennial", "diff.abs"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray80"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Centennial-Scale Climate Stability: empirical v. model")
dev.off()


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_diff_abs_century.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$var=="biomass" & lat.lon.ind & !is.na(dat.all$diff.abs),]) +
  facet_wrap(~model) +
  geom_point(aes(x=lon, y=lat, color=log(diff.abs))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
  coord_equal(xlim=lon.range, 
              ylim=lat.range, 
              expand=0) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="biomass" & lat.lon.ind, "diff.abs"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray80"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Centennial-Scale Climate Stability: empirical v. model")
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_deriv_abs_century.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$var=="biomass" & lat.lon.ind,]) +
  facet_wrap(~model) +
  geom_point(aes(x=lon, y=lat, color=log(deriv.abs))) +
  geom_point(data=dat.all[dat.all$model=="ReFAB",], aes(x=lon, y=lat, color=log(diff.abs))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
  coord_equal(xlim=lon.range, 
              ylim=lat.range, 
              expand=0) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="biomass" & lat.lon.ind, "deriv.abs"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray80"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Centennial-Scale Climate Stability: empirical v. model")
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_diff_rel.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$var=="biomass" & lat.lon.ind & !is.na(dat.all$diff.rel),]) +
  facet_wrap(~model) +
  geom_point(aes(x=lon, y=lat, color=log(diff.rel))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
  coord_equal(xlim=lon.range, 
              ylim=lat.range, 
              expand=0) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="biomass" & lat.lon.ind, "deriv.rel"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray80"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Relative Centennial-Scale Climate Stability: empirical v. model")
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_diff_rel.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$var=="biomass" & lat.lon.ind,]) +
  facet_wrap(~model) +
  geom_point(aes(x=lon, y=lat, color=log(diff.rel))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
  coord_equal(xlim=lon.range, 
              ylim=lat.range, 
              expand=0) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="biomass" & lat.lon.ind, "diff.rel"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray80"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Relative Centennial-Scale Climate Stability: empirical v. model")
dev.off()
# -----------


# -----------
# Mapping Composition
# -----------
summary(dat.all[dat.all$class=="composition",])

pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Composition_diff_abs_extent_full.pdf"))
for(mod in unique(dat.all[dat.all$class=="composition","model"])){
  cols <- ifelse(mod %in% c("STEPPS", "LINKAGES"), 3, 2)
  print(
    ggplot(data=dat.all[dat.all$class=="composition" & dat.all$model==mod & dat.all$var!="Total",]) +
      facet_wrap(~var, ncol=cols) +
      geom_point(aes(x=lon, y=lat, color=log(diff.abs))) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
      coord_equal(xlim=range(dat.all$lon), 
                  ylim=range(dat.all$lat), 
                  expand=0) +
      scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="composition", "diff.abs"]), na.rm=T)) +
      theme_bw() +
      theme(legend.position="right",
            panel.background = element_rect(fill="gray80"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      ggtitle(paste0("Centennial-Scale Climate Stability: ", mod))
  )
}
dev.off()


lon.range <- c(min(dat.all[dat.all$model=="STEPPS","lon"])-1, max(dat.all[dat.all$model=="STEPPS","lon"])+1)
lat.range <- c(min(dat.all[dat.all$model=="STEPPS","lat"])-1, max(dat.all[dat.all$model=="STEPPS","lat"])+1.5)
lat.lon.ind <- dat.all$lat>=lat.range[1] & dat.all$lat<=lat.range[2] & dat.all$lon>=lon.range[1] & dat.all$lon<=lon.range[2]

pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Composition_diff_abs_extent_stepps.pdf"))
for(mod in unique(dat.all[dat.all$class=="composition","model"])){
  cols <- ifelse(mod %in% c("STEPPS", "LINKAGES"), 3, 2)
  print(
    ggplot(data=dat.all[dat.all$class=="composition" & dat.all$model==mod & dat.all$var!="Total" & lat.lon.ind,]) +
      facet_wrap(~var, ncol=cols) +
      geom_point(aes(x=lon, y=lat, color=log(diff.abs))) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
      coord_equal(xlim=lon.range, 
                  ylim=lat.range, 
                  expand=0) +
      scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="composition" & lat.lon.ind, "diff.abs"]), na.rm=T)) +
      theme_bw() +
      theme(legend.position="right",
            panel.background = element_rect(fill="gray80"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      ggtitle(paste0("Centennial-Scale Climate Stability: ", mod))
  )
}
dev.off()
# -----------

# -------------------------------------------


# -------------------------------------------
# Model-Data Comparisons: Model ~ Data
# -------------------------------------------
# Coming up with a table of the coordinates for the models
coords.models <- drivers.all[!is.na(drivers.all$umw), c("lon", "lat")]
summary(coords.models)

# -----------
# Looking at PDSI
# -----------
summary(drivers.all)
summary(pdsi2)
summary(lbda)

dim(drivers[drivers$model=="drivers-modified" & drivers$var=="pdsi",])
dim(lbda)
lbda2 <- lbda[,c("lon", "lat", "diff.abs")]
names(lbda2)[which(names(lbda2)=="diff.abs")] <- "lbda"

dat.pdsi <- merge(drivers.all[,c("lon", "lat", "umw")], pdsi2[,c("lon", "lat", "n.yrs", "diff.abs")])
dat.pdsi <- merge(dat.pdsi, lbda2)
summary(dat.pdsi)

ggplot(data=dat.pdsi) +
  geom_point(aes(x=log(lbda), y=log(diff.abs)), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(lbda), y=log(diff.abs)), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,quantile(models1$deriv.abs, 0.97, na.rm=T))) +
  theme_bw()

ggplot(data=dat.pdsi[dat.pdsi$umw=="y",]) +
  geom_point(aes(x=log(lbda), y=log(diff.abs)), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(lbda), y=log(diff.abs)), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,quantile(models1$deriv.abs, 0.97, na.rm=T))) +
  theme_bw()

# for the whole region
lm.pdsi <- lm(log(diff.abs) ~ log(lbda), data=dat.pdsi)
summary(lm.pdsi)

# For where we have good data
lm.pdsi2 <- lm(log(diff.abs) ~ log(lbda), data=dat.pdsi[dat.pdsi$umw=="y",])
summary(lm.pdsi2 )

hist(resid(lm.pdsi))
# -----------

# -----------
# Biomass first (because it's easiest)
# -----------
for(i in 1:nrow(refab)){
  # Find out which model grid point is closest to the refab site
  lon.dist <- refab$lon[i] - coords.models$lon
  lat.dist <- refab$lat[i] - coords.models$lat
  ind.models <- which(sqrt(lon.dist^2 + lat.dist^2)==min(sqrt(lon.dist^2 + lat.dist^2)))

  models1[models1$lon==coords.models$lon[ind.models] & models1$lat==coords.models$lat[ind.models], "refab"] <- refab$diff.abs[i]
}
summary(models1)


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_DerivAbs_Biomass.png"), height=6, width=6, units="in", res=320)
ggplot(data=models1) +
  geom_point(aes(x=refab, y=deriv.abs, color=model), size=2) +
  stat_smooth(aes(x=refab, y=deriv.abs, color=model, fill=model), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model),"color"])) +
  coord_cartesian(ylim=c(0,quantile(models1$deriv.abs, 0.97, na.rm=T))) +
  theme_bw()
dev.off()

# png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_DerivAbs_Biomass.png"), height=6, width=6, units="in", res=320)
ggplot(data=models1[models1$model!="LINKAGES",]) +
  geom_point(aes(x=refab, y=deriv.abs, color=model), size=2) +
  stat_smooth(aes(x=refab, y=deriv.abs, color=model, fill=model), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model[models1$model!="LINKAGES"]),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model[models1$model!="LINKAGES"]),"color"])) +
  coord_cartesian(ylim=c(0,quantile(models1$deriv.abs[models1$model!="LINKAGES"], 0.97, na.rm=T))) +
  theme_bw()
# dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_DiffAbs_Biomass.png"), height=6, width=6, units="in", res=320)
ggplot(data=models1) +
  geom_point(aes(x=log(refab), y=log(diff.abs), color=model), size=2) +
  stat_smooth(aes(x=log(refab), y=log(diff.abs), color=model, fill=model), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model),"color"])) +
  # coord_cartesian(ylim=c(0,quantile(models1$diff.abs, 0.97, na.rm=T))) +
  theme_bw()
dev.off()

ggplot(data=models1[models1$model!="LINKAGES",]) +
  geom_point(aes(x=log(refab), y=log(diff.abs), color=model), size=2) +
  stat_smooth(aes(x=log(refab), y=log(diff.abs), color=model, fill=model), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model[models1$model!="LINKAGES"]),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models1$model[models1$model!="LINKAGES"]),"color"])) +
  # coord_cartesian(ylim=c(0,quantile(models1$deriv.abs[models1$model!="LINKAGES"], 0.97, na.rm=T))) +
  theme_bw()

lm.refab <- lm(log(refab) ~ model*log(diff.abs) - log(diff.abs) - 1, data=models1)
bm.summary <- round(summary(lm.refab)$coefficients, 3)
bm.summary
write.csv(bm.summary, file.path(path.google, "Current Data/Stability_Synthesis", "Model_v_Data_Biomass.csv"), row.names=T)

mean(models1[models1$model!="LINKAGES", "diff.abs"], na.rm=T); 
mean(models1[models1$model!="LINKAGES", "refab"], na.rm=T)
mean(models1[models1$model=="LINKAGES", "diff.abs"], na.rm=T); 

refab.ed <- lm(log(diff.abs) ~ log(refab), data=models1[models1$model=="ED2",])
summary(refab.ed)
hist(resid(refab.ed))

refab.lpjg <- lm(log(diff.abs) ~ log(refab), data=models1[models1$model=="LPJ-GUESS",])
summary(refab.lpjg)

refab.lpjw <- lm(log(diff.abs) ~ log(refab), data=models1[models1$model=="LPJ-WSL",])
summary(refab.lpjw)

refab.link <- lm(log(diff.abs) ~ log(refab), data=models1[models1$model=="LINKAGES",])
summary(refab.link)

refab.triff <- lm(log(diff.abs) ~ log(refab), data=models1[models1$model=="TRIFFID",])
summary(refab.triff)

# -----------


# -----------
# STEPPS Composition: 
# - look at only only compare dominant PFTs
# -----------
summary(stepps)
coords.stepps <- stepps[stepps$pft=="OTHER.HARDWOOD", c("lon", "lat")]

# create a dataframe where we'll store all of the fcomp
fcomp.stab <- coords.models
fcomp.stab$model <- paste(unique(fcomp$model)[1])
for(mod in unique(fcomp$model)){
  # If we need to add the model, do so first
  if(length(which(fcomp.stab$model==mod))==0) fcomp.stab <- rbind(fcomp.stab, data.frame(coords.models, model=mod, pft.abs=NA, deriv.abs=NA, diff.abs=NA))
  
  for(lon in unique(coords.models$lon)){
    for(lat in unique(coords.models[coords.models$lon==lon, "lat"])){
      df.fcomp <- fcomp[fcomp$model==mod & fcomp$lon==lon & fcomp$lat==lat,]
      if(nrow(df.fcomp)==0) next
      
      # Make sure to get rid of the "total" column
      df.fcomp <- df.fcomp[df.fcomp$pft!="Total",]
      
      ind.fcomp <- which(df.fcomp$value==max(df.fcomp$value))[1]
      ind.stab <- which(fcomp.stab$model==mod & fcomp.stab$lon==lon & fcomp.stab$lat==lat)
      
      fcomp.stab[ind.stab, "pft.abs"] <- df.fcomp$pft[ind.fcomp]
      fcomp.stab[ind.stab, "deriv.abs"] <- df.fcomp$deriv.abs[ind.fcomp]
      fcomp.stab[ind.stab, "diff.abs"] <- df.fcomp$diff.abs[ind.fcomp]
    }
  }
}
fcomp.stab$model <- as.factor(fcomp.stab$model)
fcomp.stab$pft.abs <- as.factor(fcomp.stab$pft.abs)
summary(fcomp.stab)

for(i in 1:nrow(coords.stepps)){
  # Find out which model grid point is closest to the STEPPS site
  lon.stepps <- coords.stepps$lon[i]
  lat.stepps <- coords.stepps$lat[i]

  lon.dist <- lon.stepps - coords.models$lon
  lat.dist <- lat.stepps - coords.models$lat
  ind.models <- which(sqrt(lon.dist^2 + lat.dist^2)==min(sqrt(lon.dist^2 + lat.dist^2)))
  
  lon.models <- coords.models$lon[ind.models]
  lat.models <- coords.models$lat[ind.models]
  
  # Figure out which STEPPS PFT to pull
  stepps.now <- stepps[stepps$lon==lon.stepps & stepps$lat==lat.stepps,]
  stepps.now <- stepps.now[!stepps.now$pft %in% c("Deciduous", "Evergreen"),]
  stepps.ind <- which(stepps.now$value==max(stepps.now$value, na.rm=T))[1]
  # stepps.ind <- which(stepps.now$value==max(stepps.now$value, na.rm=T))[1]
  
  fcomp.stab[fcomp.stab$lon==lon.models & fcomp.stab$lat==lat.models, "pft.stepps"] <- stepps.now$pft[stepps.ind]
  fcomp.stab[fcomp.stab$lon==lon.models & fcomp.stab$lat==lat.models, "stepps"] <- stepps.now$diff.abs[stepps.ind]
}
summary(fcomp.stab)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_DerivAbs_Fcomp.png"), height=6, width=6, units="in", res=320)
ggplot(data=fcomp.stab) +
  geom_point(aes(x=stepps, y=deriv.abs, color=model), size=2) +
  stat_smooth(aes(x=stepps, y=deriv.abs, color=model, fill=model), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model),"color"])) +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_DiffAbs_Fcomp.png"), height=6, width=6, units="in", res=320)
ggplot(data=fcomp.stab) +
  geom_point(aes(x=stepps, y=diff.abs, color=model), size=2) +
  stat_smooth(aes(x=stepps, y=diff.abs, color=model, fill=model), method="lm") +
  geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model),"color"])) +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

lm.stepps <- lm(log(stepps) ~ model*log(diff.abs) - log(diff.abs) - 1, data=fcomp.stab)
comp.summary <- round(summary(lm.stepps)$coefficients, 3)
comp.summary
write.csv(comp.summary, file.path(path.google, "Current Data/Stability_Synthesis", "Model_v_Data_Composition.csv"), row.names=T)

stepps.ed0 <- lm(log(stepps) ~ log(diff.abs), data=fcomp.stab[fcomp.stab$model=="ED2",])
summary(stepps.ed0)

stepps.ed <- lm(log(stepps) ~ log(diff.abs), data=fcomp.stab[fcomp.stab$model=="ED2",])
summary(stepps.ed)

stepps.lpjg <- lm(log(stepps) ~ log(diff.abs), data=fcomp.stab[fcomp.stab$model=="LPJ-GUESS",])
summary(stepps.lpjg)

stepps.lpjw <- lm(log(stepps) ~ log(diff.abs), data=fcomp.stab[fcomp.stab$model=="LPJ-WSL",])
summary(stepps.lpjw)

stepps.trif <- lm(log(stepps) ~ log(diff.abs), data=fcomp.stab[fcomp.stab$model=="TRIFFID",])
summary(stepps.trif)


# -----------

# -------------------------------------------
