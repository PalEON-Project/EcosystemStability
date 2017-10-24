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
library(sp)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"
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
# -----------

# -----------
# Model Drivers
# -----------
drivers.all <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
drivers <- stack(drivers.all[,c("pdsi.deriv", "tair.deriv", "precip.deriv")])
names(drivers) <- c("deriv.abs", "var")
drivers$var <- as.factor(unlist(lapply(stringr::str_split(drivers$var, "[.]"), function(x) {x[1]})))
drivers$n.sig <- stack(drivers.all[,c("pdsi.nyr", "tair.nyr", "precip.nyr")])[,1]
drivers$fract.sig <- drivers$n.sig/1000
drivers[,c("lon", "lat")] <- drivers.all[,c("lon", "lat")]
drivers$type <- as.factor("model")
drivers$resolution <- as.factor("annual")
drivers$var <- car::recode(drivers$var, "'tair'='temp'") # Since people have a hard time figuring out "Tair"
drivers$class <- as.factor("climate")
drivers$model <- as.factor("drivers")
drivers <- drivers[,c("lon", "lat", "model", "class", "var", "type", "resolution", "deriv.abs", "n.sig", "fract.sig")]
summary(drivers)
# -----------

# -----------
# STEPPS (empirical composition, centennially-resolved)
# -----------
stepps <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_STEPPS2.csv"))

# Need to convert this to latlon
albers <- sp::CRS("+proj=aea +lat_1=42.122774 +lat_2=49.01518 +lat_0=45.568977 +lon_0=-83.248627 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") # Define the spatial projection
stepps <- sp::SpatialPointsDataFrame(coords=stepps[,c("x", "y")], data=stepps, proj4string = albers) # Create a spatial object
stepps <- sp::spTransform(stepps, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # spatial projection transformation
stepps <- data.frame(stepps) # Pull out the data frame
names(stepps)[which(names(stepps) %in% c("x.1", "y.1"))] <- c("lon", "lat")
names(stepps)[which(names(stepps)=="sig")] <- c("n.sig")

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
# # NOTE: Derivative needs to be the mean ABSOLUTE VALUE!  NEEDS TO BE RE-RUN!!
# -----------
refab <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_ReFAB.csv"))
refab$refab.mean.slope <- abs(refab$refab.mean.slope) # NOTE: SHOULDN"T HAVE TO DO THIS! FIX IN FUTURE VERSIONS
names(refab) <- c("X", "lat", "lon", "deriv.abs", "n.sig")

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
models1 <- models1[,c("lon", "lat", "Model", "deriv.bm", "bm.nyr")] # Add in composition once you do it
names(models1) <- c("lon", "lat", "model", "deriv.abs", "n.sig")
models1$fract.sig <- models1$n.sig/1000

models1$class <- as.factor("biomass")
models1$var <- as.factor("biomass")
models1$type <- as.factor("model")
models1$resolution <- as.factor("annual")

summary(models1)
# -----------

# -----------
# Model Ouput - Fcomp (annually-resolved)
# -----------
fcomp <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_FCOMP_annual.csv"))
names(fcomp)[which(names(fcomp)=="n.yrs.sig")] <- "n.sig"
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
cols.bind <- c("lon", "lat", "model", "class", "var", "type", "resolution", "deriv.abs", "n.sig", "fract.sig")

dat.all <- rbind(lbda     [,cols.bind], 
                 drivers  [,cols.bind], 
                 stepps   [,cols.bind], 
                 refab    [,cols.bind], 
                 models1  [,cols.bind], 
                 # models.bm[,cols.bind],
                 fcomp    [,cols.bind])
summary(dat.all)

# Getting a relative stability so that we can ignore actual value and just look at spatial patterning
for(mod in unique(dat.all$model)){
  for(var in unique(dat.all$var[dat.all$model==mod])){
    for(res in unique(dat.all$resolution[dat.all$model==mod & dat.all$var==var])){
      deriv.subset <- dat.all[dat.all$model==mod & dat.all$var==var & dat.all$resolution==res, "deriv.abs"]
      dat.all[dat.all$model==mod & dat.all$var==var & dat.all$resolution==res, "deriv.rel"] <- deriv.subset/mean(deriv.subset, na.rm=T)
    }
  }
}

# -------------------------------------------


# -------------------------------------------
# Exploratory graphing of datasets
# -------------------------------------------
us <- map_data("state")


# -----------
# Mapping Met
# -----------
# Looking explicitly at model & driver PDSI -- Note the Log scale indiciating the orders of magnitude difference
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
  facet_wrap(~model, ncol=1) +
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


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Climate_deriv_rel.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$class=="climate" ,]) +
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


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_deriv_abs_century.png"), height=5, width=10, units="in", res=320)
ggplot(data=dat.all[dat.all$var=="biomass" & lat.lon.ind,]) +
  facet_wrap(~model) +
  geom_point(aes(x=lon, y=lat, color=log(deriv.abs))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
  coord_equal(xlim=lon.range, 
              ylim=lat.range, 
              expand=0) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(dat.all[dat.all$class=="biomass" & dat.all$resolution=="centennial" & lat.lon.ind, "deriv.abs"]), na.rm=T)) +
  theme_bw() +
  theme(legend.position="right",
        panel.background = element_rect(fill="gray80"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  ggtitle("Centennial-Scale Climate Stability: empirical v. model")
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_deriv_rel.png"), height=10, width=6, units="in", res=320)
ggplot(data=dat.all[dat.all$var=="biomass" & lat.lon.ind,]) +
  facet_wrap(~model) +
  geom_point(aes(x=lon, y=lat, color=deriv.rel)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
  coord_equal(xlim=lon.range, 
              ylim=lat.range, 
              expand=0) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(dat.all[dat.all$class=="biomass" & lat.lon.ind, "deriv.rel"], na.rm=T)) +
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

pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Composition_deriv_abs_extent_full.pdf"))
for(mod in unique(dat.all[dat.all$class=="composition","model"])){
  cols <- ifelse(mod=="STEPPS", 3, 2)
  print(
    ggplot(data=dat.all[dat.all$class=="composition" & dat.all$model==mod,]) +
      facet_wrap(~var, ncol=cols) +
      geom_point(aes(x=lon, y=lat, color=deriv.abs)) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
      coord_equal(xlim=range(dat.all$lon), 
                  ylim=range(dat.all$lat), 
                  expand=0) +
      scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(dat.all[dat.all$class=="composition", "deriv.abs"], na.rm=T)) +
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

pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Composition_deriv_abs_extent_stepps.pdf"))
for(mod in unique(dat.all[dat.all$class=="composition","model"])){
  cols <- ifelse(mod=="STEPPS", 3, 2)
  print(
    ggplot(data=dat.all[dat.all$class=="composition" & dat.all$model==mod  & lat.lon.ind,]) +
      facet_wrap(~var, ncol=cols) +
      geom_point(aes(x=lon, y=lat, color=deriv.abs)) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray25") +
      coord_equal(xlim=lon.range, 
                  ylim=lat.range, 
                  expand=0) +
      scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(dat.all[dat.all$class=="composition" & lat.lon.ind, "deriv.abs"], na.rm=T)) +
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
