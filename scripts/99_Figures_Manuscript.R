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
library(ggplot2); library(gridExtra); library(grid); library(scales)
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
# 1. Methods: Dataset Maps + time series
#
# see following post for making unbalanced multiplots: 
# https://stackoverflow.com/questions/5244014/arrange-plots-in-a-layout-which-cannot-be-achieved-by-parmfrow
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
# -------------------------------------------

# png(file.path(path.google, "Manuscript/Figures", "1a_Data_SpatialExtent.png"), height=5, width=10, units="in", res=320)
# dev.off()

# -----------------
# Stacked RAW time series of STEPPS, ReFAB, & LBDA (PDSI) for one point
# -----------------
# ------
# Find the point with a high rate of change in STEPPS
# ------
# stepps[which(stepps$deriv.abs==max(stepps$deriv.abs, na.rm=T)),]

# refab.stepps <- refab[refab$lon<max(stepps[stepps$lat<45, "lon"]),]
# refab.stepps[which(refab.stepps$refab.mean.slope.abs==max(refab.stepps$refab.mean.slope.abs, na.rm=T)),]

# Load in RAW data & subset for just the point we want to graph
# NOTE: Loading from old locations & using UNDERC for now
stepps.hips <- readRDS("~/Dropbox/PalEON_CR/PalEON_MIP_Site/Analyses/Change-and-Stability/raw_data/Benchmarks/STEPPS2/r_hips_v1.0.RDS")
stepps.hips$taxon <- as.factor(stepps.hips$taxon)
stepps.hips$year <- 2000-stepps.hips$time*100
# stepps.hips <- stepps.hips[stepps.hips$year>=850,]
summary(stepps.hips)

stepps.hips2 <- aggregate(stepps.hips$prop, by=stepps.hips[,c("site", "year", "taxon")], FUN=mean)
names(stepps.hips2)[which(names(stepps.hips2)=="x")] <- "prop.mean"
stepps.hips2$prop.025 <- aggregate(stepps.hips$prop, by=stepps.hips[,c("site", "year", "taxon")], FUN=quantile, 0.025)[,"x"]
stepps.hips2$prop.975 <- aggregate(stepps.hips$prop, by=stepps.hips[,c("site", "year", "taxon")], FUN=quantile, 0.975)[,"x"]
summary(stepps.hips2)

#Lat/Lon for HIPS
hips.latlon <- data.frame(site=c("PBL", "PDL", "PUN"), lat=c(46.28, 47.17, 46.22), lon=c(-94.58, -95.17, -89.53))
stepps.hips2 <- merge(stepps.hips2, hips.latlon)
summary(stepps.hips2)

# Finding the site with the most change in any 1 taxon
hips.change <- data.frame(site=rep(unique(stepps.hips2$site), each=length(unique(stepps.hips2$taxon))),
                          taxon=unique(stepps.hips2$taxon))
summary(hips.change)
for(SITE in unique(hips.change$site)){
  tmp <- stepps.hips2[stepps.hips2$site==SITE,]
  for(TAX in unique(hips.change$taxon)){
    hips.change[hips.change$site==SITE & hips.change$taxon==TAX,"mean"] <- mean(tmp[tmp$taxon==TAX, "prop.mean"])
    
    tax.diff <- diff(tmp[tmp$taxon==TAX, "prop.mean"])
    hips.change[hips.change$site==SITE & hips.change$taxon==TAX,"diff"] <- sum(abs(tax.diff))
  }
}
summary(hips.change)
change.max <- hips.change[hips.change$diff==max(hips.change$diff),]
site.max.dat <- hips.change[hips.change$site==change.max$site,]
pft.max <- site.max.dat[site.max.dat$mean==max(site.max.dat$mean),]

pft.max
change.max
stepps.summary <- stepps.hips2[stepps.hips2$site==site.max.dat$site & stepps.hips2$taxon==pft.max$taxon & stepps.hips2$year>=850 & stepps.hips2$year<=1850, ]
stepps.summary$data <- "Composition"
site.use <- hips.latlon[hips.latlon$site==paste(change.max$site),]
# ------


load(file.path(path.google, "Current Data", "refab.sites.lat.lon.Rdata"))
refab.latlon <- lat.lon.df; rm(lat.lon.df)

load(file.path(path.google, "Current Data", "ReFAB.all.samps.list.Rdata"))
summary(refab.latlon)

load("~/Desktop/Research/PalEON_EcosystemStability/data/refab.biomass.CI13.Rdata")
refab.raw <- all.samps.list; rm(all.samps.list)
refab.raw <- refab.raw[as.numeric(summary(refab.raw)[,1])>0]
# summary(refab.raw)
dim(refab.latlon); length(refab.raw)
summary(refab.latlon)

latlon.max <- hips.latlon[hips.latlon$site==paste(change.max$site),]
refab.latlon$dist <- sqrt((refab.latlon$lon - latlon.max$lon)^2 + (refab.latlon$lat - latlon.max$lat)^2)
summary(refab.latlon)s

refab.raw2 <- refab.raw[[which(refab.latlon$dist==min(refab.latlon$dist))]]
refab.raw2 <- refab.raw2[,90:100] # following Ann's Annotation in 1_biomass-data-stability.R
summary(refab.raw2)
refab.summary <- data.frame(data="Biomass",
                            year = seq(850, 1850, by=100),
                            mean=apply(refab.raw2, 2, mean),
                            lo = apply(refab.raw2, 2, quantile, 0.025),
                            hi = apply(refab.raw2, 2, quantile, 0.975))
summary(refab.summary)
#

# path.data <- "~/Desktop/Research/PalEON_MIP_Region/NADA/"
lbda.nc <- ncdf4::nc_open("~/Desktop/Research/PalEON_MIP_Region/NADA/nada_hd2_cl.nc")

lbda.lat <- ncdf4::ncvar_get(lbda.nc, "lat")
lbda.lon <- ncdf4::ncvar_get(lbda.nc, "lon")
lbda.time <- ncdf4::ncvar_get(lbda.nc, "time")
lbda.nc$dim$time$units # Goes from past to present: 0-2005

lbda.res <- mean(diff(lbda.lon))
lbda.lon.ind <- which(lbda.lon+lbda.res/2>=site.use$lon & lbda.lon-lbda.res/2<=site.use$lon)
lbda.lat.ind <- which(lbda.lat+lbda.res/2>=site.use$lat & lbda.lat-lbda.res/2<=site.use$lat)
lbda.time.ind <- which(lbda.time>=850 & lbda.time<=1850)
lbda.lon  <- lbda.lon [lbda.lon.ind]
lbda.lat  <- lbda.lat [lbda.lat.ind]
lbda.time <- lbda.time[lbda.time.ind]

lbda.raw <- ncdf4::ncvar_get(lbda.nc, "pdsi")[lbda.time.ind, lbda.lat.ind, lbda.lon.ind]
lbda.raw[lbda.raw < -99] <- NA
summary(lbda.raw)

lbda.summary <- data.frame(data="Drought", year=850:1850, raw=lbda.raw)
summary(lbda.summary)

# Running the 100-year gam on the LBDA so we can compare it to paleo data
library(mgcv)
path.gamm.func <- "~/Desktop/Research/R_Functions/"
source(file.path(path.gamm.func, "Calculate_GAMM_Derivs.R"))
source(file.path(path.gamm.func, "Calculate_GAMM_Posteriors.R"))

calc.stability <- function(x, width=100){
  dat.tmp <- data.frame(Y=x, Year=1:length(x))
  k.use=round(length(x[!is.na(x)])/width, 0)
  mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
  
  yrs.cent <- rev(seq(length(x), 1, by=-100)) # Go backwards so everythign lines up in 1850
  
  mod.out <- list()
  mod.out$gam.post <- post.distns(mod.gam, newdata=dat.tmp[, ], vars="Year", return.sims=T)$sims # Note col1=X (index); col2=year
  mod.out$mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
  
  return(mod.out)
}

lbda.gam <- calc.stability(lbda.summary[!is.na(lbda.summary$raw), "raw"], width=100)
summary(lbda.gam)
summary(lbda.gam$gam.post[,1:15])
summary(lbda.gam$mod.deriv)
dim(lbda.gam$gam.post)

lbda.summary[!is.na(lbda.summary$raw), "mean"] <- apply(lbda.gam$gam.post[,3:ncol(lbda.gam$gam.post)], 1, mean)
lbda.summary[!is.na(lbda.summary$raw), "lo"  ] <- apply(lbda.gam$gam.post[,3:ncol(lbda.gam$gam.post)], 1, quantile, 0.025)
lbda.summary[!is.na(lbda.summary$raw), "hi"  ] <- apply(lbda.gam$gam.post[,3:ncol(lbda.gam$gam.post)], 1, quantile, 0.975)

# Stacking things together
summary(stepps.summary)
summary(refab.summary)
summary(lbda.summary)

dat.raw <- data.frame(data= c(rep("Composition", nrow(stepps.summary)), rep("Biomass", nrow(refab.summary)), rep("PDSI", nrow(lbda.summary))),
                      year= c(stepps.summary$year, refab.summary$year, lbda.summary$year),
                      raw = c(rep(NA, nrow(stepps.summary) + nrow(refab.summary)), lbda.summary$raw),
                      mean= c(stepps.summary$prop.mean, refab.summary$mean, lbda.summary$mean),
                      lo  = c(stepps.summary$prop.025, refab.summary$lo, lbda.summary$lo),
                      hi  = c(stepps.summary$prop.975, refab.summary$hi, lbda.summary$hi))
dat.raw$data <- factor(dat.raw$data, levels=c("Composition", "Biomass", "PDSI"))
summary(dat.raw)

plot.time <- ggplot(data=dat.raw) +
  facet_grid(data~., scales="free_y", switch="y") +
  geom_line(aes(x=year, y=raw), size=0.1, color="black") +
  geom_pointrange(data=dat.raw[dat.raw$data %in% c("Composition", "Biomass"),], aes(x=year, y=mean, ymin=lo, ymax=hi), size=1) +
  # geom_point(data=dat.raw[dat.raw$data %in% c("Composition", "Biomass"),], aes(x=year, y=mean), size=5) +
  geom_pointrange(data=dat.raw[dat.raw$data %in% c("PDSI") & dat.raw$year %in% seq(850, 1850, by=100),], 
                  aes(x=year, y=mean, ymin=lo, ymax=hi), size=1, color="black") +
  geom_line(aes(x=year, y=mean), color="blue", size=2) +
  theme_bw() +
  theme(strip.placement = "outside",
        axis.title.y = element_blank())

# lbda.raw <- 
# lbda.gam <- 
# -----------

# -----------
# Map of STEPPS vs. ReFAB vs. model grid 
# -----------
us <- map_data("state")

drivers.all <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
stepps <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_STEPPS.csv"))
refab <- read.csv(file.path(path.google, "Current Data", "refab.mean.slope.csv"))
# models1 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))

map.dat <- ggplot() +
  geom_tile(data=drivers.all, aes(x=lon, y=lat, fill=tair.yr.set-273.15)) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray20") + 
  geom_point(data=stepps[stepps$pft=="OTHER.HARDWOOD",], aes(x=lon, y=lat, color="STEPPS", shape="STEPPS"), size=1, alpha=0.75) +
  # geom_point(data=models1[models1$Model=="ED2",], aes(x=lon, y=lat, color="Models"), shape=19, size=1.25) +
  geom_point(data=refab, aes(x=lon, y=lat, color="ReFAB", shape="ReFAB"), size=1.75) +
  geom_point(data=hips.latlon[hips.latlon$site==paste(hips.change[hips.change$diff==max(hips.change$diff),"site"]),], aes(x=lon, y=lat), color="darkgoldenrod", size=5) +
  scale_color_manual(name="Datasets", values=c("black", "gray50")) +
  scale_shape_manual(name="Datasets", values=c(19,15)) +
  scale_fill_gradient2(name="Mean Temp \n(1800-1850)", low = "blue", high = "red", mid = "white", midpoint = mean(drivers.all$tair.yr.set-273.15, na.rm=T)) +
  labs(x="Longitude", y="Latitude") +
  guides(color=guide_legend(aes.overide=list(size=20)))+
  guides(shape=guide_legend(aes.override=list(size=20))) +
  coord_equal(xlim=range(refab$lon, na.rm=T)+c(-1.5,1.5), ylim=range(refab$lat,na.rm=T)+c(-1,2), expand=0) +
  theme(legend.position="top",
        legend.key=element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
# -----------



# ----------------------
# The Figure!
# ----------------------
png(file.path(path.google, "Manuscript/Figures", "1_Map_TimeSeries_Tall.png"), height=8, width=8, units="in", res=220)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1, heights=c(2.5, 1))))
print(map.dat  , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.time, vp = viewport(layout.pos.row = 2, layout.pos.col=1))
dev.off()

png(file.path(path.google, "Manuscript/Figures", "1_Map_TimeSeries_Wide.png"), height=4, width=8, units="in", res=220)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2, widths=c(2, 1))))
print(map.dat  , vp = viewport(layout.pos.row = 1, layout.pos.col=1))
print(plot.time, vp = viewport(layout.pos.row = 1, layout.pos.col=2))
dev.off()

# ----------------------


rm(drivers.all, stepps, refab, models1)
# -------------------------------------------








# -------------------------------------------
# 3. Results: Map of empirical vs. driver PDSI stability
# -------------------------------------------
library(ggplot2)
us <- map_data("state")

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
drivers[,"diff.abs"] <- stack(drivers.all[,c("pdsi.diff", "tair.diff", "precip.diff")])[,1]
drivers$var <- as.factor(unlist(lapply(stringr::str_split(drivers$var, "[.]"), function(x) {x[1]})))
drivers$n.sig <- stack(drivers.all[,c("pdsi.nyr", "tair.nyr", "precip.nyr")])[,1]
drivers$fract.sig <- drivers$n.sig/1000
drivers[,c("lon", "lat")] <- drivers.all[,c("lon", "lat")]
drivers$type <- as.factor("model")
drivers$resolution <- as.factor("annual")
drivers$var <- car::recode(drivers$var, "'tair'='temp'") # Since people have a hard time figuring out "Tair"
drivers$class <- as.factor("climate")
drivers$model <- as.factor("drivers (850-2010)")
drivers <- drivers[,c("lon", "lat", "model", "class", "var", "type", "resolution", "diff.abs", "deriv.abs", "n.sig", "fract.sig")]
summary(drivers)

# Load in the PDSI calc that's match in temporal extent to LBDA
pdsi2 <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_PDSI_Drivers_LBDA_time_100.csv"))
pdsi2$var <- as.factor("pdsi")
pdsi2$type <- as.factor("model")
pdsi2$resolution <- as.factor("annual")
pdsi2$class <- as.factor("climate")
pdsi2$model <- as.factor("drivers (LBDA time)")
names(pdsi2)[which(names(pdsi2)=="n.yrs.sig")] <- "n.sig"
summary(pdsi2)

drivers <- merge(drivers, pdsi2, all=T)
summary(drivers)
# -----------


cols.bind <- c("lon", "lat", "model", "class", "var", "type", "resolution", "diff.abs", "deriv.abs", "n.sig", "fract.sig")

stab.clim <- rbind(lbda   [,cols.bind], 
                   drivers[,cols.bind]
                   )

stab.clim <- stab.clim[!is.na(stab.clim$diff.abs),]
stab.clim$model <- factor(stab.clim$model, levels=c("LBDA", "drivers (LBDA time)", "drivers (850-2010)"))
summary(stab.clim)

LETTERS[1:3]
panel.labs <- data.frame(model=levels(stab.clim$model), lab=LETTERS[1:3])

png(file.path(path.google, "Manuscript/Figures", "3_PDSI_Stability_GeographicSpace.png"), height=10, width=8, units="in", res=220)
ggplot(data=stab.clim[stab.clim$var=="pdsi",]) +
  facet_wrap(~model, ncol=1) +
  geom_tile(aes(x=lon, y=lat, fill=log(diff.abs))) +
  geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
  coord_equal(xlim=range(stab.clim$lon), ylim=range(stab.clim$lat), expand=0) +
  geom_text(data=panel.labs, aes(x=-63, y=37, label=paste0(lab, ")")), fontface="bold") +
  scale_fill_gradient(low="gray90", high="red3") +
  labs(x="longitude", y="latitude") +
  # scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(log(stab.clim[stab.clim$var=="pdsi","diff.abs"]), na.rm=T)) +
  theme_bw() +
  # theme(strip.text = element_blank()) +
  theme(legend.position="top",
        panel.background = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
dev.off()
# -------------------------------------------


# -------------------------------------------
# 4. Results: Scatterplot of ecosystem vs. climate stability (models + data)
#    4.a. Composition
#    4.b. Biomass (separate panels for Linkages & other data)
# -------------------------------------------
drivers.all <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
summary(drivers.all)

fcomp.stab <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "StabilityComparision_EnvironmentalSpace_Composition.csv"))
biom.stab  <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "StabilityComparision_EnvironmentalSpace_Biomass.csv"))

fcomp.stab$model <- factor(fcomp.stab$model, levels=c("STEPPS", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))
biom.stab $model <- factor(biom.stab $model, levels=c("ReFAB", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))

fcomp.stab$scale <- as.factor("other")
fcomp.stab$var <- as.factor("composition")
biom.stab $var <- as.factor("biomass")

cols.bind <- c("lon", "lat", "var", "model", "type", "scale", "pft.abs", "value", "diff.abs", "diff.pdsi")
stab.env <- rbind(fcomp.stab[,cols.bind], biom.stab[,cols.bind])
stab.env$model <- factor(stab.env$model, levels=c("STEPPS", "ReFAB", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))
stab.env <- merge(stab.env, drivers.all[,c("lon", "lat", "umw")], all.x=T)
summary(stab.env)

# Doing some relativization to see if we can find a way of saying composition is MORE or LESS stable than biomass
for(VAR in unique(stab.env$var)){
  for(MOD in unique(stab.env[stab.env$var==VAR, "model"])){
    var.mean <- mean(stab.env[stab.env$var==VAR & stab.env$model==MOD, "value"], na.rm=T)
    var.max  <- max (stab.env[stab.env$var==VAR & stab.env$model==MOD, "value"], na.rm=T)
    
    stab.env[stab.env$var==VAR & stab.env$model==MOD, "diff.rel.mean"] <- stab.env[stab.env$var==VAR & stab.env$model==MOD, "diff.abs"]/var.mean
    stab.env[stab.env$var==VAR & stab.env$model==MOD, "diff.rel.max"] <- stab.env[stab.env$var==VAR & stab.env$model==MOD, "diff.abs"]/var.max
    
  }
}
summary(stab.env)

mean(stab.env[stab.env$var=="biomass" & stab.env$type=="model", "diff.rel.mean"], na.rm=T); sd(stab.env[stab.env$var=="biomass" & stab.env$type=="model", "diff.rel.mean"], na.rm=T)
mean(stab.env[stab.env$var=="composition" & stab.env$type=="model", "diff.rel.mean"], na.rm=T); sd(stab.env[stab.env$var=="composition" & stab.env$type=="model", "diff.rel.mean"], na.rm=T)

# ----------------
# Comparing relative stability of composition & biomass in and between models & data
# ----------------
hist(log(stab.env[stab.env$var=="biomass" & stab.env$type=="model", "diff.rel.mean"]))
hist(log(stab.env[stab.env$var=="composition" & stab.env$type=="model", "diff.rel.mean"]))

# Comparing biomass & composition within data types
t.test(log(stab.env[stab.env$var=="biomass" & stab.env$type=="model", "diff.rel.mean"]), log(stab.env[stab.env$var=="composition" & stab.env$type=="model", "diff.rel.mean"]), paired=T)
t.test(log(stab.env[stab.env$var=="biomass" & stab.env$type=="empirical", "diff.rel.mean"]), log(stab.env[stab.env$var=="composition" & stab.env$type=="empirical", "diff.rel.mean"]))

t.test(log(stab.env[stab.env$var=="biomass" & stab.env$type=="model" & stab.env$umw=="y", "diff.rel.mean"]), log(stab.env[stab.env$var=="biomass" & stab.env$type=="empirical", "diff.rel.mean"]))
t.test(log(stab.env[stab.env$var=="composition" & stab.env$type=="model" & stab.env$umw=="y", "diff.rel.mean"]), log(stab.env[stab.env$var=="composition" & stab.env$type=="empirical", "diff.rel.mean"]))
# ----------------

# Creating a dummy data frame since I can't do coord_cartesian with scales="free_y"
# stab.env2 <- stab.env
# stab.env2[stab.env2$var=="composition" & log(stab.env2$diff.abs) < -11 & !is.na(stab.env2$diff.abs),"diff.abs"] <- NA
# stab.env2[stab.env2$var=="biomass" & log(stab.env2$diff.abs) < -8 & !is.na(stab.env2$diff.abs),"diff.abs"] <- NA

png(file.path(path.google, "Manuscript/Figures", "4_Model_v_Data_Stability_EnvironmentalSpace_zoomed.png"), height=10, width=8, units="in", res=220)
ggplot(data=stab.env) +
  facet_grid(var~., scales="fixed", shrink=T) +
  geom_point(aes(x=log(diff.pdsi), y=log(diff.rel.mean), color=model), size=0.5, alpha=0.3) +
  stat_smooth(aes(x=log(diff.pdsi), y=log(diff.rel.mean), color=model, fill=model), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.env$model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(stab.env$model),"color"])) +
  # geom_abline(intercept=0, slope=1, color="black", linetype="dashed") +
  coord_cartesian(ylim=c(-10.0, -5.5)) +
  theme_bw() +
  theme(legend.position="top")
dev.off()
# -------------------------------------------


# -------------------------------------------
# 5. Results: Scatterplots by models (group by complexity/responses)
# -------------------------------------------
# ------------
# Most model variables (compsition excluded)
# ------------
stab.models <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))
summary(stab.models)
vars.use <- c("gpp", "npp", "nee", "lai", "bm")
type.use <- "diff."

models.long <- stack(stab.models[,paste0(type.use, vars.use)])
names(models.long) <- c("diff.abs", "var")
models.long$val.mean <- stack(stab.models[,paste0("mean.", vars.use)])[,1]
models.long[,c("lon", "lat", "Model")] <- stab.models[,c("lon", "lat", "Model")]
models.long$var <- as.factor(matrix(unlist(lapply(models.long$var, function(x) stringr::str_split(x, "[.]"))), ncol=2, byrow=T)[,2])
models.long$pft.abs <- NA # a place holder that will get filled in a minute
summary(models.long)
# ------------


# ------------
# Fcomp -- select only dominant PFT and merge into models.long
# ------------
load(file.path(path.data, "PalEON_siteInfo_all.RData"))
paleon$cell <- 1:nrow(paleon)
summary(paleon)

fcomp.diff  <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_FCOMP_DIFF.csv"))
# fcomp.deriv <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_FCOMP_DERIV.csv"))

# Fixing some column names
names(fcomp.diff)[which(names(fcomp.diff) %in% c("deriv.abs", "deriv.mean"))] <- c("diff.abs", "diff.mean")
summary(fcomp.diff)

# Adding lat/lon into the fcomp data frames
fcomp.diff  <- merge(fcomp.diff , paleon[,c("cell", "lon", "lat")], all.x=T)

# Pulling out just the dominant PFT
fcomp <- stab.models[,c("lon", "lat", "Model")]
fcomp$var <- as.factor("fcomp")
summary(fcomp)

for(mod in unique(fcomp$Model)){
  for(lon in unique(fcomp[fcomp$Model==mod, "lon"])){
    for(lat in unique(fcomp[fcomp$Model==mod & fcomp$lon==lon, "lat"])){
      df.fcomp <- fcomp.diff[fcomp.diff$model==mod & fcomp.diff$lon==lon & fcomp.diff$lat==lat,]
      if(nrow(df.fcomp)==0) next
      
      # Make sure to get rid of the "total" column
      df.fcomp <- df.fcomp[df.fcomp$pft!="Total",]
      
      ind.fcomp <- which(df.fcomp$value==max(df.fcomp$value))[1]
      ind.stab <- which(fcomp$Model==mod & fcomp$lon==lon & fcomp$lat==lat)
      
      fcomp[ind.stab, "pft.abs" ] <- paste(df.fcomp$pft[ind.fcomp])
      fcomp[ind.stab, "diff.abs"] <- df.fcomp$diff.abs[ind.fcomp]
      fcomp[ind.stab, "val.mean"] <- df.fcomp$value[ind.fcomp]
      
      # Adding the dominant PFT into the model output
      models.long[models.long$Model==mod & models.long$lon==lon & models.long$lat==lat,"pft.abs"] <- paste(df.fcomp$pft[ind.fcomp])
    }
  }
}
fcomp$pft.abs <- as.factor(fcomp$pft.abs)
models.long$pft.abs <- as.factor(models.long$pft.abs)
summary(fcomp)
summary(models.long)

# Merge the composition into models.long
models.long <- merge(models.long, fcomp, all=T)
models.long <- models.long[complete.cases(models.long),]
summary(models.long)
# ------------

# ------------
# Met -- merge into models.long
# ------------
stab.met   <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
summary(stab.met)

# Merge the met drivers into models.long
models.long <- merge(models.long, stab.met[,c("lon", "lat", "pdsi.diff", "tair.diff", "precip.diff", "tair.yr.set", "precip.yr.set")])
models.long$var <- factor(models.long$var, levels=c("gpp", "nee", "npp", "lai", "bm", "fcomp"))
summary(models.long)
# ------------


for(v in unique(models.long$var)){
  for(mod in unique(models.long$Model)){
    dat.ind <- which(models.long$var==v & models.long$Model==mod & !is.na(models.long$diff.abs))
    dat.tmp <- models.long[dat.ind, ]
    
    models.long[dat.ind, "rel.val.mean"  ] <- dat.tmp$diff.abs/mean(dat.tmp$val.mean)
    models.long[dat.ind, "rel.diff.mean" ] <- dat.tmp$diff.abs/mean(dat.tmp$diff.abs)
    models.long[dat.ind, "rel.diff.range"] <- 1 - (max(dat.tmp$diff.abs) - dat.tmp$diff.abs) / 
      (max(dat.tmp$diff.abs) - min(dat.tmp$diff.abs))
    
    models.long[dat.ind, "pdsi.rel.diff.mean" ] <- dat.tmp$pdsi.diff/mean(dat.tmp$pdsi.diff)
    models.long[dat.ind, "pdsi.rel.diff.range"] <- 1 - (max(dat.tmp$pdsi.diff) - dat.tmp$pdsi.diff) / 
      (max(dat.tmp$pdsi.diff) - min(dat.tmp$pdsi.diff))
    
  }
}
summary(models.long)

models.long$Model <- factor(models.long$Model, levels=c("ED2", "LPJ-GUESS", "LPJ-WSL", "LINKAGES", "TRIFFID"))

png(file.path(path.google, "Manuscript/Figures", "5_Model_VarComparison_EcosysStability_zoomed.png"), height=6, width=8, units="in", res=220)
ggplot(data=models.long, aes(x=log(pdsi.rel.diff.mean), y=log(rel.diff.mean), color=var, fill=var)) +
  labs(x="PDSI Stability", y="Ecosystem Stability") +
  facet_wrap(~Model) +
  geom_point(size=0.1, alpha=0.3) +
  stat_smooth(method=lm, alpha=0.5) + 
  scale_color_manual(values=c("blue4", "darkseagreen4", "turquoise4", "darkgoldenrod2", "darkorange2", "deeppink3")) +
  scale_fill_manual(values=c("blue4", "darkseagreen4", "turquoise4", "darkgoldenrod2", "darkorange2", "deeppink3")) +
  # scale_fill_brewer(palette="Dark2") +
  # scale_color_brewer(palette="Dark2") +
  coord_cartesian(ylim=c(-2, 1.5)) +
  theme_bw() +
  theme(legend.position=c(0.75, 0.25))
dev.off()
# -------------------------------------------
