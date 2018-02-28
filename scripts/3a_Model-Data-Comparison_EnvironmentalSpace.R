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
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"

dat.colors <- data.frame(model=c("LBDA", "STEPPS", "ReFAB", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"),
                         type=c(rep("emipirical", 3), rep("model", 7)),
                         color=c(rep("#000000", 3), rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
dat.colors
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

lbda.diff <- rasterize(lbda, lbda.rast, lbda$diff.abs, fun=mean)
# projection(lbda.deriv) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
names(lbda.diff) <- "diff.pdsi"
lbda.diff
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
drivers$diff.pdsi <- rasterize(drivers.all, domain, drivers.all$pdsi.diff, fun=mean)
drivers$diff.tair <- rasterize(drivers.all, domain, drivers.all$tair.diff, fun=mean)
drivers$diff.precip <- rasterize(drivers.all, domain, drivers.all$precip.diff, fun=mean)
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
names(stepps)[which(names(stepps)=="sig")] <- c("n.sig")
names(stepps)[which(names(stepps)=="deriv.abs")] <- c("diff.abs")
stepps$deriv.abs <- NA

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
refab.means <- read.csv(file.path(path.google, "Current Data", "biomass.means.csv"))
names(refab.means)[4] <- "value" 
summary(refab.means)

refab <- read.csv(file.path(path.google, "Current Data", "refab.mean.slope.csv"))
names(refab) <- c("X", "lat", "lon", "diff.mean", "diff.abs", "n.sig")
refab$deriv.abs <- NA

refab <- merge(refab, refab.means[,c("lon", "lat", "value")])
summary(refab)

# Doing a unit correction
refab[,c("value", "diff.mean", "diff.abs")] <- refab[,c("value", "diff.mean", "diff.abs")]*0.1

refab <- refab[!is.na(refab$lat),]
coordinates(refab) <- refab[,c("lon", "lat")]
projection(refab) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
summary(refab)
# -----------

# -----------
# Model Output
# -----------
models1 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))
models1 <- models1[,c("lon", "lat", "Model", "mean.bm", "diff.bm", "deriv.bm", "nyr.bm")] # Add in composition once you do it
names(models1) <- c("lon", "lat", "model", "value", "diff.abs", "deriv.abs", "n.sig")
models1$fract.sig <- models1$n.sig/1000

models1[models1$value<=0 & !is.na(models1$value),c("value", "diff.abs", "deriv.abs", "n.sig", "fract.sig")] <- NA

coordinates(models1) <- models1[,c("lon", "lat")]
projection(models1) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
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

# Making a single spatial fcomp
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
stepps$diff.pdsi <- extract(lbda.diff, stepps)
refab$diff.pdsi <- extract(lbda.diff, refab)
summary(stepps)
summary(refab)

# Drivers
drivers.stepps <- data.frame(extract(drivers, stepps))
drivers.stepps <- drivers.stepps[,names(drivers.stepps)[!names(drivers.stepps) %in% c("deriv.pdsi", "diff.pdsi")]]
stepps <- data.frame(stepps)
stepps[,names(drivers.stepps)] <- drivers.stepps
stepps[,c("deriv.tair", "deriv.precip", "diff.tair", "diff.precip")] <- NA
summary(stepps)

drivers.refab <- data.frame(extract(drivers, refab))
drivers.refab <- drivers.refab[,names(drivers.refab)[!names(drivers.refab) %in% c("deriv.pdsi", "diff.pdsi")]]
refab <- data.frame(refab)
refab[,names(drivers.refab)] <- drivers.refab
refab[,c("deriv.tair", "deriv.precip", "diff.tair", "diff.precip")] <- NA
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
# Model-Data Comparisons: Composition
# -------------------------------------------
# -----------
# Data fromatting: Condense down to 1 value per grid cell
# -----------
coords.models <- data.frame(drivers.all[!is.na(drivers.all$umw), c("lon", "lat")])[,1:2]
summary(coords.models)

fcomp.stab <- coords.models
fcomp.stab$model <- paste(unique(fcomp$model)[1])
fcomp.stab$type <- "model"
for(mod in unique(fcomp$model)){
  # If we need to add the model, do so first
  if(length(which(fcomp.stab$model==mod))==0) fcomp.stab <- rbind(fcomp.stab, 
                                                                  data.frame(coords.models, model=mod, type="model", pft.abs=NA, diff.abs=NA, deriv.abs=NA,
                                                                             whc=NA, tair.sett=NA, precip.sett=NA, 
                                                                             deriv.pdsi=NA, deriv.tair=NA, deriv.precip=NA, 
                                                                             diff.pdsi=NA, diff.tair=NA, diff.precip=NA))
  
  for(lon in unique(coords.models$lon)){
    for(lat in unique(coords.models[coords.models$lon==lon, "lat"])){
      df.fcomp <- fcomp[fcomp$model==mod & fcomp$lon==lon & fcomp$lat==lat,]
      if(nrow(df.fcomp)==0) next
      
      # Make sure to get rid of the "total" column
      df.fcomp <- df.fcomp[df.fcomp$pft!="Total",]
      
      ind.fcomp <- which(df.fcomp$value==max(df.fcomp$value))[1]
      ind.stab <- which(fcomp.stab$model==mod & fcomp.stab$lon==lon & fcomp.stab$lat==lat)
      
      fcomp.stab[ind.stab, "pft.abs"     ] <- paste(df.fcomp$pft[ind.fcomp])
      fcomp.stab[ind.stab, "diff.abs"    ] <- df.fcomp$diff.abs[ind.fcomp]
      fcomp.stab[ind.stab, "deriv.abs"   ] <- df.fcomp$deriv.abs[ind.fcomp]
      fcomp.stab[ind.stab, "whc"         ] <- df.fcomp$whc[ind.fcomp]
      fcomp.stab[ind.stab, "tair.sett"   ] <- df.fcomp$tair.sett[ind.fcomp]
      fcomp.stab[ind.stab, "precip.sett" ] <- df.fcomp$precip.sett[ind.fcomp]
      fcomp.stab[ind.stab, "deriv.pdsi"  ] <- df.fcomp$deriv.pdsi[ind.fcomp]
      fcomp.stab[ind.stab, "deriv.tair"  ] <- df.fcomp$deriv.tair[ind.fcomp]
      fcomp.stab[ind.stab, "deriv.precip"] <- df.fcomp$deriv.precip[ind.fcomp]
      fcomp.stab[ind.stab, "diff.pdsi"   ] <- df.fcomp$diff.pdsi[ind.fcomp]
      fcomp.stab[ind.stab, "diff.tair"   ] <- df.fcomp$diff.tair[ind.fcomp]
      fcomp.stab[ind.stab, "diff.precip" ] <- df.fcomp$diff.precip[ind.fcomp]
      
      # Adding the dominant PFT into the model output
      models1[models1$model==mod & models1$lon==lon & models1$lat==lat,"pft.abs"] <- paste(df.fcomp$pft[ind.fcomp])
    }
  }
}

coords.stepps <- stepps[stepps$pft=="OTHER.HARDWOOD", c("lon", "lat")]
fcomp.stab <- rbind(fcomp.stab, 
                    data.frame(coords.stepps, model="STEPPS", type="empirical", pft.abs=NA, diff.abs=NA, deriv.abs=NA,
                               whc=NA, tair.sett=NA, precip.sett=NA, 
                               deriv.pdsi=NA, deriv.tair=NA, deriv.precip=NA,
                               diff.pdsi=NA , diff.tair=NA , diff.precip=NA))
for(lon in unique(coords.stepps$lon)){
  for(lat in unique(coords.stepps[coords.stepps$lon==lon, "lat"])){
    ind.stab <- which(fcomp.stab$model=="STEPPS" & fcomp.stab$lon==lon & fcomp.stab$lat==lat)
    
    df.fcomp <- stepps[stepps$lon==lon & stepps$lat==lat,]
    # if(nrow(df.fcomp)==0) next
    
    # Make sure to get rid of the total Decid/Evergreen column
    df.fcomp <- df.fcomp[!df.fcomp$pft %in% c("Deciduous", "Evergreen"),]
    
    ind.fcomp <- which(df.fcomp$value==max(df.fcomp$value))[1]

    fcomp.stab[ind.stab, "pft.abs"     ] <- paste(df.fcomp$pft[ind.fcomp])
    fcomp.stab[ind.stab, "diff.abs"    ] <- df.fcomp$diff.abs[ind.fcomp]
    fcomp.stab[ind.stab, "deriv.abs"   ] <- df.fcomp$deriv.abs[ind.fcomp]
    fcomp.stab[ind.stab, "whc"         ] <- df.fcomp$whc[ind.fcomp]
    fcomp.stab[ind.stab, "tair.sett"   ] <- df.fcomp$tair.sett[ind.fcomp]
    fcomp.stab[ind.stab, "precip.sett" ] <- df.fcomp$precip.sett[ind.fcomp]
    fcomp.stab[ind.stab, "deriv.pdsi"  ] <- df.fcomp$deriv.pdsi[ind.fcomp]
    fcomp.stab[ind.stab, "deriv.tair"  ] <- df.fcomp$deriv.tair[ind.fcomp]
    fcomp.stab[ind.stab, "deriv.precip"] <- df.fcomp$deriv.precip[ind.fcomp]
    fcomp.stab[ind.stab, "diff.pdsi"  ] <- df.fcomp$diff.pdsi[ind.fcomp]
    fcomp.stab[ind.stab, "diff.tair"  ] <- df.fcomp$diff.tair[ind.fcomp]
    fcomp.stab[ind.stab, "diff.precip"] <- df.fcomp$diff.precip[ind.fcomp]
  }
}

fcomp.stab$model <- as.factor(fcomp.stab$model)
fcomp.stab$type <- as.factor(fcomp.stab$type)
fcomp.stab$pft.abs <- as.factor(fcomp.stab$pft.abs)
summary(fcomp.stab)

fcomp.stab$model <- factor(fcomp.stab$model, levels=c("STEPPS", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Composition_PDSI.png"), height=6, width=6, units="in", res=320)
ggplot(data=fcomp.stab) +
  # facet_wrap(~model, scales="free") +
  geom_point(aes(x=log(diff.pdsi), y=log(diff.abs), color=model), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(diff.pdsi), y=log(diff.abs), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw() +
  theme(legend.position=c(0.8, 0.2))
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_Composition_Temp.png"), height=6, width=6, units="in", res=320)
ggplot(data=fcomp.stab[fcomp.stab$type=="model",]) +
  # facet_wrap(~model, scales="free") +
  geom_point(aes(x=log(diff.tair), y=log(diff.abs), color=model), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(diff.tair), y=log(diff.abs), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model[fcomp.stab$type=="model"]),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model[fcomp.stab$type=="model"]),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_Composition_Precip.png"), height=6, width=6, units="in", res=320)
ggplot(data=fcomp.stab[fcomp.stab$type=="model",]) +
  # facet_wrap(~model, scales="free") +
  geom_point(aes(x=log(diff.precip), y=log(diff.abs), color=model), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(diff.precip), y=log(diff.abs), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model[fcomp.stab$type=="model"]),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(fcomp.stab$model[fcomp.stab$type=="model"]),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()
# -----------

# -----------
# Linear models
# -----------
# Simple comparison with climate
comp.pdsi.pft <- lm(log(diff.abs) ~ model*log(diff.pdsi)*pft.abs, data=fcomp.stab)
anova(comp.pdsi.pft)

# comp.summary <- round(summary(comp.pdsi.pft)$coefficients, 3)
comp.summary <- round(data.frame(anova(comp.pdsi.pft)), 3)
comp.summary

write.csv(comp.summary, file.path(path.google, "Current Data/Stability_Synthesis", "ANOVA_Composition_PFT.csv"), row.names=T)

comp.pdsi.pft2 <- lm(log(diff.abs) ~ model*log(diff.pdsi), data=fcomp.stab)
anova(comp.pdsi.pft2)
comp.summary2 <- round(summary(comp.pdsi.pft2)$coefficients, 3)
comp.summary2

write.csv(comp.summary2, file.path(path.google, "Current Data/Stability_Synthesis", "ANOVA_Composition_RelEffects.csv"), row.names=T)

comp.pdsi.pft3 <- lm(log(diff.abs) ~ model*log(diff.pdsi) - log(diff.pdsi) -1, data=fcomp.stab)
anova(comp.pdsi.pft3)
comp.summary3 <- round(summary(comp.pdsi.pft3)$coefficients, 3)
comp.summary3

write.csv(comp.summary3, file.path(path.google, "Current Data/Stability_Synthesis", "ANOVA_Composition_AbsoltueEffects.csv"), row.names=T)


fcomp.pdsi0 <- lm(log(diff.abs) ~ model*log(diff.pdsi), data=fcomp.stab)
summary(fcomp.pdsi0)

fcomp.pdsi <- lm(log(diff.abs) ~ model*log(diff.pdsi), data=fcomp.stab)
summary(fcomp.pdsi)

anova(fcomp.pdsi)


pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Composition_v_PDSI.pdf"))
for(mod in unique(fcomp.stab$model)){
  # Only show pfts with >1 data point
  pfts.use <- summary(fcomp.stab[fcomp.stab$model==mod & !is.na(fcomp.stab$pft.abs), "pft.abs"])
  pfts.use <- names(pfts.use)[which(pfts.use>1)]
  print(
    ggplot(data=fcomp.stab[fcomp.stab$model==mod & fcomp.stab$pft.abs %in% pfts.use,]) +
      # facet_wrap(~model, scales="free") +
      geom_point(aes(x=log(diff.pdsi), y=log(diff.abs), color=pft.abs), size=2) +
      stat_smooth(aes(x=log(diff.pdsi), y=log(diff.abs), color=pft.abs, fill=pft.abs), method="lm") +
      # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
      # coord_cartesian(ylim=c(0,0.2)) +
      ggtitle(mod) +
      theme_bw()
    
  )
}
dev.off()

fcomp.pdsi.stepps <- lm(log(diff.abs) ~ log(diff.pdsi), data=fcomp.stab[fcomp.stab$model=="STEPPS",])
summary(fcomp.pdsi.stepps)

fcomp.pdsi.stepps2 <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=fcomp.stab[fcomp.stab$model=="STEPPS",])
summary(fcomp.pdsi.stepps2)
anova(fcomp.pdsi.stepps2)

fcomp.pdsi.ed2 <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs-log(diff.pdsi), data=fcomp.stab[fcomp.stab$model=="ED2",])
summary(fcomp.pdsi.ed2)
anova(fcomp.pdsi.ed2)

fcomp.pdsi.link0 <- lm(log(diff.abs) ~ log(diff.pdsi) , data=fcomp.stab[fcomp.stab$model=="LINKAGES",])
summary(fcomp.pdsi.link0)

fcomp.pdsi.link <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=fcomp.stab[fcomp.stab$model=="LINKAGES",])
summary(fcomp.pdsi.link)
anova(fcomp.pdsi.link)

fcomp.pdsi.lpjw <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=fcomp.stab[fcomp.stab$model=="LPJ-WSL",])
anova(fcomp.pdsi.lpjw)
summary(fcomp.pdsi.lpjw)

fcomp.pdsi.lpjg <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=fcomp.stab[fcomp.stab$model=="LPJ-GUESS",])
summary(fcomp.pdsi.lpjg)
anova(fcomp.pdsi.lpjg)

fcomp.pdsi.trif <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=fcomp.stab[fcomp.stab$model=="TRIFFID",])
summary(fcomp.pdsi.trif)
# B = Broadleaf
# N = Needleleaf
# S = shrub

# Checking our residuals
hist(residuals(fcomp.pdsi))
plot(residuals(fcomp.pdsi) ~ predict(fcomp.pdsi)); abline(h=0, col="red")

# Multi-variate analysis
fcomp.env <- lm(log(deriv.abs) ~ model*log(deriv.pdsi)*tair.sett*precip.sett, data=fcomp.stab)
# summary(fcomp.env)
anova(fcomp.env)

hist(residuals(fcomp.env))
plot(residuals(fcomp.env) ~ predict(fcomp.env)); abline(h=0, col="red")


# Looking at complicated sets of predictors (don't try explaining results by PFT)
fcomp.env.stepps <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*pft.abs, data=fcomp.stab[fcomp.stab$model=="STEPPS",])
anova.fcomp.stepps <- data.frame(anova(fcomp.env.stepps))
names(anova.fcomp.stepps)[ncol(anova.fcomp.stepps)] <- "P.STEPPS"
anova.fcomp.stepps$factor <- row.names(anova.fcomp.stepps)

fcomp.env.ed2 <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*pft.abs, data=fcomp.stab[fcomp.stab$model=="ED2",])
anova.fcomp.ed2 <- data.frame(anova(fcomp.env.ed2))
names(anova.fcomp.ed2)[ncol(anova.fcomp.ed2)] <- "P.ED2"
anova.fcomp.ed2$factor <- row.names(anova.fcomp.ed2)

fcomp.env.link <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*pft.abs, data=fcomp.stab[fcomp.stab$model=="LINKAGES",])
anova.fcomp.link <- data.frame(anova(fcomp.env.link))
names(anova.fcomp.link)[ncol(anova.fcomp.link)] <- "P.LINKAGES"
anova.fcomp.link$factor <- row.names(anova.fcomp.link)
anova(fcomp.env.link)

fcomp.env.lpjg <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*pft.abs, data=fcomp.stab[fcomp.stab$model=="LPJ-GUESS",])
anova.fcomp.lpjg <- data.frame(anova(fcomp.env.lpjg))
names(anova.fcomp.lpjg)[ncol(anova.fcomp.lpjg)] <- "P.LPJ.GUESS"
anova.fcomp.lpjg$factor <- row.names(anova.fcomp.lpjg)
anova(fcomp.env.lpjg)

fcomp.env.lpjw <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*pft.abs, data=fcomp.stab[fcomp.stab$model=="LPJ-WSL",])
anova.fcomp.lpjw <- data.frame(anova(fcomp.env.lpjw))
names(anova.fcomp.lpjw)[ncol(anova.fcomp.lpjw)] <- "P.LPJ.WSL"
anova.fcomp.lpjw$factor <- row.names(anova.fcomp.lpjw)
anova(fcomp.env.lpjw)

fcomp.env.triff <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*pft.abs, data=fcomp.stab[fcomp.stab$model=="TRIFFID",])
anova.fcomp.triff <- data.frame(anova(fcomp.env.triff))
names(anova.fcomp.triff)[ncol(anova.fcomp.triff)] <- "P.TRIFFID"
anova.fcomp.triff$factor <- row.names(anova.fcomp.triff)
anova(fcomp.env.triff)

# Creating an anova table
anova.fcomp <- merge(anova.fcomp.stepps[,c("factor", "P.STEPPS")], anova.fcomp.ed2[,c("factor", "P.ED2")], all=T)
anova.fcomp <- merge(anova.fcomp, anova.fcomp.link [,c("factor", "P.LINKAGES") ], all=T)
anova.fcomp <- merge(anova.fcomp, anova.fcomp.lpjg [,c("factor", "P.LPJ.GUESS")], all=T)
anova.fcomp <- merge(anova.fcomp, anova.fcomp.lpjw [,c("factor", "P.LPJ.WSL")  ], all=T)
anova.fcomp <- merge(anova.fcomp, anova.fcomp.triff[,c("factor", "P.TRIFFID")  ], all=T)

anova.fcomp[,2:ncol(anova.fcomp)] <- round(anova.fcomp[,2:ncol(anova.fcomp)],3)
anova.fcomp <- anova.fcomp[anova.fcomp$factor!="Residuals",]
anova.fcomp$level <- unlist(lapply(stringr::str_split(anova.fcomp$factor, ":"), length))
anova.fcomp <- anova.fcomp[order(anova.fcomp$level),]
anova.fcomp

write.csv(anova.fcomp, file.path(path.google, "Current Data/Stability_Synthesis", "ANOVA_Composition_Environment_Model.csv"), row.names=F)

# -----------
# -------------------------------------------

# -------------------------------------------
# Model-Data Comparisons: Biomass
# -------------------------------------------
# -----------
# Data fromatting and exploration
# -----------
# Add the dominant PFT from the STEPPS model to the output
for(i in 1:nrow(refab)){
  # Find out which model grid point is closest to the STEPPS site
  lon.refab <- refab$lon[i]
  lat.refab <- refab$lat[i]
  
  lon.dist <- lon.refab - stepps$lon
  lat.dist <- lat.refab - stepps$lat
  ind.stepps <- which(sqrt(lon.dist^2 + lat.dist^2)==min(sqrt(lon.dist^2 + lat.dist^2)))
  
  # Figure out which STEPPS PFT to pull
  stepps.now <- stepps[ind.stepps,]
  stepps.now <- stepps.now[!stepps.now$pft %in% c("Deciduous", "Evergreen"),]
  stepps.ind <- which(stepps.now$value==max(stepps.now$value, na.rm=T))[1]
  # stepps.ind <- which(stepps.now$value==max(stepps.now$value, na.rm=T))[1]
  
  refab[refab$lon==lon.refab & refab$lat==lat.refab, "pft.abs"] <- stepps.now$pft[stepps.ind]
}
summary(refab)

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

cols.bind <- c("lon", "lat", "type", "model", "value", "diff.abs", "deriv.abs", "pft.abs", names(drivers))
stab.bm <- rbind(refab[,cols.bind], models1[,cols.bind])
summary(stab.bm)

stab.bm$scale <- as.factor(ifelse(stab.bm$model=="LINKAGES", "LINKAGES", "other"))
summary(stab.bm)


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_PDSI.png"), height=6, width=6, units="in", res=320)
ggplot(data=stab.bm) +
  facet_wrap(~scale, scales="free") +
  geom_point(aes(x=log(diff.pdsi), y=log(diff.abs), color=model), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(diff.pdsi), y=log(diff.abs), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_Biomass_Temp.png"), height=6, width=6, units="in", res=320)
ggplot(data=stab.bm[stab.bm$type=="model",]) +
  facet_wrap(~scale, scales="free") +
  geom_point(aes(x=log(diff.tair), y=log(diff.abs), color=model), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(diff.tair), y=log(diff.abs), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model[stab.bm$type=="model"]),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model[stab.bm$type=="model"]),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_Biomass_Precip.png"), height=6, width=6, units="in", res=320)
ggplot(data=stab.bm[stab.bm$type=="model",]) +
  facet_wrap(~scale, scales="free") +
  geom_point(aes(x=log(diff.precip), y=log(diff.abs), color=model), size=1, alpha=0.5) +
  stat_smooth(aes(x=log(diff.precip), y=log(diff.abs), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model[stab.bm$type=="model"]),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model[stab.bm$type=="model"]),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_v_Data_Biomass_meanBiomass.png"), height=6, width=6, units="in", res=320)
ggplot(data=stab.bm) +
  facet_wrap(~scale, scales="free") +
  geom_point(aes(x=value, y=log(diff.abs), color=model), size=1, alpha=0.5) +
  stat_smooth(aes(x=value, y=log(diff.abs), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.bm$model),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  # coord_cartesian(ylim=c(0,0.2)) +
  theme_bw()
dev.off()

# -----------

# -----------
# Linear models
# -----------
# Simple comparison with climate
# stab.bm[stab.bm$diff.abs==0 & !is.na(stab.bm$diff.abs), "diff.abs"] <- 1e-10
bm.pdsi.bm <- lm(log(diff.abs) ~ model*log(diff.pdsi), data=stab.bm)
bm.summary <- round(summary(bm.pdsi.bm)$coefficients, 3)
bm.summary

write.csv(bm.summary, file.path(path.google, "Current Data/Stability_Synthesis", "ANOVA_Biomass_RelEffects.csv"), row.names=T)

bm.pdsi.bm2 <- lm(log(diff.abs) ~ model*log(diff.pdsi)-1 - log(diff.pdsi), data=stab.bm)
bm.summary2 <- round(summary(bm.pdsi.bm2)$coefficients, 3)
bm.summary2

write.csv(bm.summary2, file.path(path.google, "Current Data/Stability_Synthesis", "ANOVA_Biomass_AbsoltueEffects.csv"), row.names=T)


bm.pdsi <- lm(log(diff.abs) ~ model*log(diff.pdsi), data=stab.bm)
summary(bm.pdsi)
anova(bm.pdsi)

# Checking our residuals
hist(residuals(bm.pdsi))
plot(residuals(bm.pdsi) ~ predict(bm.pdsi)); abline(h=0, col="red")

bm.pdsi.link1 <- lm(log(diff.abs) ~ log(diff.pdsi), data=stab.bm[stab.bm$model=="LINKAGES",])
summary(bm.pdsi.link1)

bm.pdsi.triff1 <- lm(log(diff.abs) ~ log(diff.pdsi), data=stab.bm[stab.bm$model=="TRIFFID",])
summary(bm.pdsi.triff1)

# Looking at biomass-backgroudn biomass
bm.bm <- lm(log(diff.abs) ~ model*value, data=stab.bm)
summary(bm.bm)

hist(residuals(bm.bm))
plot(residuals(bm.bm) ~ predict(bm.bm)); abline(h=0, col="red")



# Breaking PDSI-Biomass down by PFT within a model
pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Biomass_v_PDSI_PFTs.pdf"))
for(mod in unique(stab.bm$model)){
  # Only show pfts with >1 data point
  pfts.use <- summary(stab.bm[stab.bm$model==mod & !is.na(stab.bm$pft.abs), "pft.abs"])
  pfts.use <- names(pfts.use)[which(pfts.use>1)]
  print(
    ggplot(data=stab.bm[stab.bm$model==mod & stab.bm$pft.abs %in% pfts.use,]) +
      # facet_wrap(~model, scales="free") +
      geom_point(aes(x=log(diff.pdsi), y=log(diff.abs), color=pft.abs), size=2) +
      stat_smooth(aes(x=log(diff.pdsi), y=log(diff.abs), color=pft.abs, fill=pft.abs), method="lm") +
      # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
      # coord_cartesian(ylim=c(0,0.2)) +
      ggtitle(mod) +
      theme_bw()
    
  )
}
dev.off()

pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Biomass_v_meanBiomass_PFTs.pdf"))
for(mod in unique(stab.bm$model)){
  # Only show pfts with >1 data point
  pfts.use <- summary(stab.bm[stab.bm$model==mod & !is.na(stab.bm$pft.abs), "pft.abs"])
  pfts.use <- names(pfts.use)[which(pfts.use>1)]
  print(
    ggplot(data=stab.bm[stab.bm$model==mod & stab.bm$pft.abs %in% pfts.use,]) +
      # facet_wrap(~model, scales="free") +
      geom_point(aes(x=value, y=log(diff.abs), color=pft.abs), size=2) +
      stat_smooth(aes(x=value, y=log(diff.abs), color=pft.abs, fill=pft.abs), method="lm") +
      # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
      # coord_cartesian(ylim=c(0,0.2)) +
      ggtitle(mod) +
      theme_bw()
    
  )
}
dev.off()

# Stability vs climate
bm.pdsi.refab <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=stab.bm[stab.bm$model=="ReFAB",])
summary(bm.pdsi.refab)
anova(bm.pdsi.refab)

bm.pdsi.refab2 <- lm(log(diff.abs) ~ log(diff.pdsi)*log(value) - log(diff.pdsi), data=stab.bm[stab.bm$model=="ReFAB",])
anova(bm.pdsi.refab2)
summary(bm.pdsi.refab2)

# bm.pdsi.refab3 <- lm(log(diff.abs) ~ log(value) - log(diff.pdsi), data=stab.bm[stab.bm$model=="ReFAB",])
# anova(bm.pdsi.refab3)
# summary(bm.pdsi.refab3)


bm.pdsi.ed2 <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=stab.bm[stab.bm$model=="ED2",])
anova(bm.pdsi.ed2)
summary(bm.pdsi.ed2)

bm.pdsi.ed22 <- lm(log(diff.abs) ~ log(diff.pdsi)*log(value) - log(diff.pdsi), data=stab.bm[stab.bm$model=="ED2",])
anova(bm.pdsi.ed22)
summary(bm.pdsi.ed22)

bm.pdsi.lpjg <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=stab.bm[stab.bm$model=="LPJ-GUESS",])
anova(bm.pdsi.lpjg)
summary(bm.pdsi.lpjg)

bm.pdsi.lpjg2 <- lm(log(diff.abs) ~ log(diff.pdsi)*log(value) - log(diff.pdsi), data=stab.bm[stab.bm$model=="LPJ-GUESS",])
anova(bm.pdsi.lpjg2)
summary(bm.pdsi.lpjg2)

bm.pdsi.lpjw <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=stab.bm[stab.bm$model=="LPJ-WSL",])
anova(bm.pdsi.lpjw)
summary(bm.pdsi.lpjw)

bm.pdsi.lpjw2 <- lm(log(diff.abs) ~ log(diff.pdsi)*log(value) - log(diff.pdsi), data=stab.bm[stab.bm$model=="LPJ-WSL",])
anova(bm.pdsi.lpjw2)
summary(bm.pdsi.lpjw2)

# bm.pdsi.lpjw3 <- lm(log(diff.abs) ~ log(value), data=stab.bm[stab.bm$model=="LPJ-WSL",])
# anova(bm.pdsi.lpjw3)
# summary(bm.pdsi.lpjw3)

bm.pdsi.link <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=stab.bm[stab.bm$model=="LINKAGES",])
anova(bm.pdsi.link)
summary(bm.pdsi.link)

bm.pdsi.link2 <- lm(log(diff.abs) ~ log(diff.pdsi)*log(value) - log(diff.pdsi), data=stab.bm[stab.bm$model=="LINKAGES",])
anova(bm.pdsi.link2)
summary(bm.pdsi.link2)

bm.pdsi.triff <- lm(log(diff.abs) ~ log(diff.pdsi)*pft.abs - log(diff.pdsi), data=stab.bm[stab.bm$model=="TRIFFID",])
anova(bm.pdsi.triff)
summary(bm.pdsi.triff)

bm.pdsi.triff2 <- lm(log(diff.abs) ~ log(diff.pdsi)*log(value) - log(diff.pdsi), data=stab.bm[stab.bm$model=="TRIFFID",])
anova(bm.pdsi.triff2)
summary(bm.pdsi.triff2)


# Stability vs. Biomass
bm.bm.refab <- lm(log(diff.abs) ~ value*pft.abs - value, data=stab.bm[stab.bm$model=="ReFAB",])
summary(bm.bm.refab)

bm.bm.ed2 <- lm(log(diff.abs) ~ value*pft.abs - value, data=stab.bm[stab.bm$model=="ED2",])
summary(bm.bm.ed2)

bm.bm.lpjg <- lm(log(diff.abs) ~ value*pft.abs - value, data=stab.bm[stab.bm$model=="LPJ-GUESS",])
summary(bm.bm.lpjg)

bm.bm.lpjw <- lm(log(diff.abs) ~ value*pft.abs - value, data=stab.bm[stab.bm$model=="LPJ-WSL",])
summary(bm.bm.lpjw)

bm.bm.link <- lm(log(diff.abs) ~ value*pft.abs - value, data=stab.bm[stab.bm$model=="LINKAGES",])
summary(bm.bm.link)

bm.bm.triff <- lm(log(diff.abs) ~ value*pft.abs - value, data=stab.bm[stab.bm$model=="TRIFFID",])
summary(bm.bm.triff)

# Multi-variate analysis
bm.env <- lm(log(diff.abs) ~ model*log(diff.pdsi)*tair.sett*precip.sett*value, data=stab.bm)
# summary(bm.env)
anova(bm.env)

hist(residuals(bm.env))
plot(residuals(bm.env) ~ predict(bm.env)); abline(h=0, col="red")

bm.env.refab <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*value, data=stab.bm[stab.bm$model=="ReFAB",])
anova.bm.refab <- data.frame(anova(bm.env.refab))
names(anova.bm.refab)[ncol(anova.bm.refab)] <- "P.ReFAB"
anova.bm.refab$factor <- row.names(anova.bm.refab)

bm.env.ed2 <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*value, data=stab.bm[stab.bm$model=="ED2",])
anova.bm.ed2 <- data.frame(anova(bm.env.ed2))
names(anova.bm.ed2)[ncol(anova.bm.ed2)] <- "P.ED2"
anova.bm.ed2$factor <- row.names(anova.bm.ed2)

bm.env.lpjg <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*value, data=stab.bm[stab.bm$model=="LPJ-GUESS",])
anova.bm.lpjg <- data.frame(anova(bm.env.lpjg))
names(anova.bm.lpjg)[ncol(anova.bm.lpjg)] <- "P.LPJ.GUESS"
anova.bm.lpjg$factor <- row.names(anova.bm.lpjg)

bm.env.lpjw <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*value, data=stab.bm[stab.bm$model=="LPJ-WSL",])
anova.bm.lpjw <- data.frame(anova(bm.env.lpjw))
names(anova.bm.lpjw)[ncol(anova.bm.lpjw)] <- "P.LPJ.WSL"
anova.bm.lpjw$factor <- row.names(anova.bm.lpjw)

bm.env.link <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*value, data=stab.bm[stab.bm$model=="LINKAGES",])
anova.bm.link <- data.frame(anova(bm.env.link))
names(anova.bm.link)[ncol(anova.bm.link)] <- "P.LINKAGES"
anova.bm.link$factor <- row.names(anova.bm.link)

bm.env.triff <- lm(log(diff.abs) ~ log(diff.pdsi)*tair.sett*precip.sett*value, data=stab.bm[stab.bm$model=="TRIFFID",])
anova.bm.triff <- data.frame(anova(bm.env.triff))
names(anova.bm.triff)[ncol(anova.bm.triff)] <- "P.TRIFFID"
anova.bm.triff$factor <- row.names(anova.bm.triff)

# Creating an anova table
anova.bm <- merge(anova.bm.refab[,c("factor", "P.ReFAB")], anova.bm.ed2[,c("factor", "P.ED2")], all=T)
anova.bm <- merge(anova.bm, anova.bm.link [,c("factor", "P.LINKAGES") ], all=T)
anova.bm <- merge(anova.bm, anova.bm.lpjg [,c("factor", "P.LPJ.GUESS")], all=T)
anova.bm <- merge(anova.bm, anova.bm.lpjw [,c("factor", "P.LPJ.WSL")  ], all=T)
anova.bm <- merge(anova.bm, anova.bm.triff[,c("factor", "P.TRIFFID")  ], all=T)

anova.bm[,2:ncol(anova.bm)] <- round(anova.bm[,2:ncol(anova.bm)],3)
anova.bm <- anova.bm[anova.bm$factor!="Residuals",]
anova.bm$level <- unlist(lapply(str_split(anova.bm$factor, ":"), length))
anova.bm <- anova.bm[order(anova.bm$level),]
head(anova.bm); tail(anova.bm)
 
write.csv(anova.bm, file.path(path.google, "Current Data/Stability_Synthesis", "ANOVA_Biomass_Environment_Model.csv"), row.names=F)
# -----------


# -------------------------------------------
