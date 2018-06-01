# -------------------------------------------
# 2a. Data Availability in Space & Time; Bias of missing PDSI data
# 
# 1. Map of Data availability in upper midwest + Time series of Composition, Biomass & PDSI
# 2. Test to see if PDSI stability in UPM time period is equivalent to that of full period
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
# -------------------------------------------

# -------------------------------------------
# Is uneven PDSI representation across the region a problem?
# -------------------------------------------
pdsi.df <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_PDSI_Drivers_LBDA_time_100.csv"))
summary(pdsi.df)

# -------------------------------------------
