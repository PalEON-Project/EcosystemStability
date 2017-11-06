# -------------------------------------------
# Model-Data comparisons: Comparisons of Model & Data Stability in Climate & Environmental Space

# Here we will:
# 2. Look at patterns of stability in environmental space
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
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"
# -------------------------------------------

# -------------------------------------------
# Read in & align predictor layers
# -------------------------------------------
# -----------
# Empirical Drought Reconstruction (LBDA, annually-resolved)
# -----------
lbda <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_LBDA_100.csv"))
coordinates(lbda) <- lbda[,c("lon", "lat")]

# Set up a blank raster
lbda.rast <- raster(ncol=length(unique(lbda$lon)), nrow=length(unique(lbda$lat)))
extent(lbda.rast) <- extent(lbda)
projection(lbda.rast) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Rasterize the layer we want
lbda.deriv <- rasterize(lbda, lbda.rast, lbda$deriv.abs, fun=mean)
# projection(lbda.deriv) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
names(lbda.deriv) <- "deriv.pdsi"
lbda.deriv
# -----------

# -----------
# Model Drivers
# -----------
drivers.all <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
coordinates(drivers.all) <- drivers.all[,c("lon", "lat")]
drivers.all$whc.tot[drivers.all$whc.tot>1e3] <- NA
summary(drivers.all)

# Set up a blank raster
domain <- raster(ncol=length(unique(drivers.all$lon)), nrow=length(unique(drivers.all$lat)))
extent(domain) <- extent(drivers.all)
projection(domain) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Rasterize the layer we want
drivers <- rasterize(drivers.all, domain, drivers.all$depth, fun=mean)
# projection(drivers) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
names(drivers) <- "soil.depth"
drivers <- stack(drivers)
drivers

drivers$whc <- rasterize(drivers.all, domain, drivers.all$whc.tot, fun=mean)
drivers$tair.sett <- rasterize(drivers.all, domain, drivers.all$tair.yr.set, fun=mean)
drivers$precip.sett <- rasterize(drivers.all, domain, drivers.all$precip.yr.set, fun=mean)
drivers$deriv.pdsi <- rasterize(drivers.all, domain, drivers.all$pdsi.deriv, fun=mean)
drivers$deriv.tair <- rasterize(drivers.all, domain, drivers.all$tair.deriv, fun=mean)
drivers$deriv.precip <- rasterize(drivers.all, domain, drivers.all$precip.deriv, fun=mean)
drivers
# plot(drivers)
# -----------
# -------------------------------------------

# -------------------------------------------
# Read in, format, & extract the point data
# -------------------------------------------
# -----------
# STEPPS (empirical composition, centennially-resolved)
# -----------
stepps <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_STEPPS.csv"))
coordinates(stepps) <- stepps[,c("lon", "lat")]
projection(stepps) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(stepps)

# plot(drivers$whc)
# plot(stepps, add=T, col="black", cex=0.2)
# -----------

# -----------
# ReFAB (empirical biomass, centennially-resolved)
# # NOTE: Derivative needs to be the mean ABSOLUTE VALUE!  NEEDS TO BE RE-RUN!!
# -----------
refab <- read.csv(file.path(path.google, "Current Data", "refab.mean.slope.csv"))
names(refab) <- c("X", "lat", "lon", "deriv.mean", "deriv.abs", "n.sig")
refab <- refab[!is.na(refab$lat),]
coordinates(refab) <- refab[,c("lon", "lat")]
projection(refab) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(refab)
# -----------

# -----------
# Model Output
# -----------
models1 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))
models1 <- models1[,c("lon", "lat", "Model", "deriv.bm", "bm.nyr")] # Add in composition once you do it
names(models1) <- c("lon", "lat", "model", "deriv.abs", "n.sig")
coordinates(models1) <- models1[,c("lon", "lat")]
projection(models1) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(models1)
# -----------


# -----------
# Model Ouput - Fcomp (annually-resolved)
# -----------
fcomp <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_FCOMP_annual.csv"))
coordinates(fcomp) <- fcomp[,c("lon", "lat")]
projection(fcomp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(fcomp)
# -----------
# -------------------------------------------


# -------------------------------------------
# Extract environment data for the models & data
# Note: Using simple extractions rather than nearest point because the gridded predictors are supposed to be the average of the grid cell
# -------------------------------------------
# -----------
# Empirical: STEPPS, REFAB
# -----------
# Drought Stability
stepps$deriv.pdsi <- extract(lbda.deriv, stepps)
refab$deriv.pdsi <- extract(lbda.deriv, refab)
summary(stepps)
summary(refab)

# Drivers
drivers.stepps <- data.frame(extract(drivers, stepps))
drivers.stepps <- drivers.stepps[,names(drivers.stepps)[!names(drivers.stepps)=="deriv.pdsi"]]
stepps <- data.frame(stepps)
stepps[,names(drivers.stepps)] <- drivers.stepps
stepps[,c("deriv.tair", "deriv.precip")] <- NA
summary(stepps)

drivers.refab <- data.frame(extract(drivers, refab))
drivers.refab <- drivers.refab[,names(drivers.refab)[!names(drivers.refab)=="deriv.pdsi"]]
refab <- data.frame(refab)
refab[,names(drivers.refab)] <- drivers.refab
refab[,c("deriv.tair", "deriv.precip")] <- NA
summary(refab)
# -----------

# -----------
# Models: models1, fcomp
# -----------
# Drivers
drivers.models1 <- extract(drivers, models1)
models1 <- data.frame(models1)
models1[,names(drivers)] <- drivers.models1
summary(models1)

drivers.fcomp <- extract(drivers, fcomp)
fcomp <- data.frame(fcomp)
fcomp[,names(drivers)] <- drivers.fcomp
summary(fcomp)
# -----------
# -------------------------------------------


# -------------------------------------------
# Model-Data Comparisons: Biomass
# -------------------------------------------
# -----------
# Data fromatting and exploration
# -----------
refab$model <- as.factor("ReFAB")
refab$class <- as.factor("biomass")
refab$var <- as.factor("biomass")
refab$type <- as.factor("empirical")
refab$resolution <- as.factor("centennial")
summary(refab)

models1$class <- as.factor("biomass")
models1$var <- as.factor("biomass")
models1$type <- as.factor("model")
models1$resolution <- as.factor("annual")
summary(models1)

cols.bind <- c("lon", "lat", "type", "model", "deriv.abs", names(drivers))
stab.bm <- rbind(refab[,cols.bind], models1[,cols.bind])
summary(stab.bm)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_PDSI.png"), height=6, width=6, units="in", res=320)
ggplot(data=stab.bm) +
  # facet_wrap(~model, scales="free") +
  geom_point(aes(x=log(deriv.pdsi), y=log(deriv.abs), color=model), size=2) +
  stat_smooth(aes(x=log(deriv.pdsi), y=log(deriv.abs), color=model, fill=model), method="lm") +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_Biomass_Temp.png"), height=6, width=6, units="in", res=320)
ggplot(data=stab.bm[stab.bm$type=="model",]) +
  # facet_wrap(~model, scales="free") +
  geom_point(aes(x=log(deriv.tair), y=log(deriv.abs), color=model), size=2) +
  stat_smooth(aes(x=log(deriv.tair), y=log(deriv.abs), color=model, fill=model), method="lm") +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_Biomass_Precip.png"), height=6, width=6, units="in", res=320)
ggplot(data=stab.bm[stab.bm$type=="model",]) +
  # facet_wrap(~model, scales="free") +
  geom_point(aes(x=log(deriv.precip), y=log(deriv.abs), color=model), size=2) +
  stat_smooth(aes(x=log(deriv.precip), y=log(deriv.abs), color=model, fill=model), method="lm") +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()
# -----------

# -----------
# Linear models
# -----------
# Simple comparison with climate
bm.pdsi <- lm(log(deriv.abs) ~ model*log(deriv.pdsi), data=stab.bm)
summary(bm.pdsi)
anova(bm.pdsi)

bm.pdsi.link <- lm(log(deriv.abs) ~ log(deriv.pdsi), data=stab.bm[stab.bm$model=="LINKAGES",])
summary(bm.pdsi.link)

# Checking our residuals
hist(residuals(bm.pdsi))
plot(residuals(bm.pdsi) ~ predict(bm.pdsi)); abline(h=0, col="red")

# Multi-variate analysis
bm.env <- lm(log(deriv.abs) ~ model*log(deriv.pdsi)*tair.sett*precip.sett*whc, data=stab.bm)
summary(bm.env)
anova(bm.env)

hist(residuals(bm.env))
plot(residuals(bm.env) ~ predict(bm.env)); abline(h=0, col="red")

bm.env.refab <- lm(log(deriv.abs) ~ log(deriv.pdsi)*tair.sett*precip.sett*whc, data=stab.bm[stab.bm$model=="ReFAB",])
anova(bm.env.refab)

bm.env.ed2 <- lm(log(deriv.abs) ~ log(deriv.pdsi)*tair.sett*precip.sett*whc, data=stab.bm[stab.bm$model=="ED2",])
anova(bm.env.ed2)

bm.env.refab <- lm(log(deriv.abs) ~ log(deriv.pdsi)*tair.sett*precip.sett*whc, data=stab.bm[stab.bm$model=="ReFAB",])
anova(bm.env.refab)

bm.env.lpjg <- lm(log(deriv.abs) ~ log(deriv.pdsi)*tair.sett*precip.sett*whc, data=stab.bm[stab.bm$model=="LPJ-GUESS",])
anova(bm.env.lpjg)

bm.env.lpjw <- lm(log(deriv.abs) ~ log(deriv.pdsi)*tair.sett*precip.sett*whc, data=stab.bm[stab.bm$model=="LPJ-WSL",])
anova(bm.env.lpjw)

bm.env.link <- lm(log(deriv.abs) ~ log(deriv.pdsi)*tair.sett*precip.sett*whc, data=stab.bm[stab.bm$model=="LINKAGES",])
anova(bm.env.link)

# -----------


# -------------------------------------------
