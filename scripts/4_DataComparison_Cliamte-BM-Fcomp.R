# -------------------------------------------
# Data-Climate Comparisons: 
# 
# H0: Climate change causes ecosystem change, as seen through shifts in 
#     biomass and composition (consistent & proportional responses; data)
#        Prediction: Locations that have experienced the most changes in 
#                    climate show the most changes in biomass and composition
#        Prediction: Places with the greatest change in composition also 
#                    have the greatest change in biomass (Zhang 2018 Nature Climate Change)
# 
# -------------------------------------------
# Steps:
# -------------------------------------------
# 0. Load libraries, set file paths etc.
# 1. Read in Datasets; Standardize Stability Metrics
# 2. Compare Biomass & Composition to Climate (climate window; do they show the same pattern?)
# 3. Compare Biomass to Composition (Full time?; are they correlated at all?)
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
# 1. Read in Datasets; Standardize Stability Metrics
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
lbda$stability <- -log(lbda$diff.abs/mean(lbda$diff.abs, na.rm=T))  # Note: Positive numbers mean MORE stable
summary(lbda)

# plot(stability ~ log(diff.abs), data=lbda)

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
stepps$stability <- -log(stepps$diff.abs/mean(stepps$diff.abs, na.rm=T))
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
refab$stability <- -log(refab$diff.abs/mean(refab$diff.abs, na.rm=T))
summary(refab)
# -----------
# -------------------------------------------

# -------------------------------------------
# 3. Compare Biomass to Composition (Full time?; are they correlated at all?)
# -------------------------------------------
# -----------
# Extracting STEPPS info for refab
# -----------

coords.stepps <- stepps[stepps$pft=="OTHER.HARDWOOD", c("lon", "lat")]

dat.sites <- data.frame(refab[,c("lat", "lon")],
                        stab.refab=refab$stability)
for(i in 1:nrow(dat.sites)){
  # Find the closest stepps site
  stepps.dist <- sqrt((stepps$lon - refab$lon[i])^2 + (stepps$lat - refab$lat[i])^2)
  stepps.ind <- which(stepps.dist==min(stepps.dist, na.rm=T) & !df.fcomp$pft %in% c("Deciduous", "Evergreen")) # Don't include generic Decid/Evergreen
  
  # Subset to just that site
  fcomp.tmp <- stepps[stepps.ind,]

  # Find the dominant PFT 
  ind.fcomp <- which(fcomp.tmp$value==max(fcomp.tmp$value))
  
  dat.sites[i, "H.prime"    ] <- - sum(fcomp.tmp$value * log(fcomp.tmp$value)) # calculate shannon-weiner index
  dat.sites[i, "dom.pft"    ] <- paste(fcomp.tmp$pft[ind.fcomp])
  dat.sites[i, "stab.stepps"] <- fcomp.tmp$stability[ind.fcomp]
}
dat.sites$dom.pft <- as.factor(dat.sites$dom.pft)
summary(dat.sites)
# -----------


# -----------
# Quantitative comparisons between Biomass & Composition Stability
# -----------
# Comparing means -- is one more stable than the other?
t.test(dat.sites$stab.refab, dat.sites$stab.stepps, paired=T)
mean(dat.sites$stab.refab - dat.sites$stab.stepps)
sd(dat.sites$stab.refab - dat.sites$stab.stepps)

# Is biomass stability correlated with composition stability?
lm.comp <- lm(stab.refab ~ stab.stepps, data=dat.sites)
summary(lm.comp)

# Is biomass stability correlated with diversity?
lm.diversity <- lm(stab.refab ~ H.prime, data=dat.sites)
summary(lm.diversity)

# Is composition stability correlated with Diversity
lm.diversity2 <- lm(stab.stepps ~ H.prime, data=dat.sites)
summary(lm.diversity2)

# Some summary figures
png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Data_Biomass_v_Composition.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites) +
  geom_point(aes(x=stab.stepps, y=stab.refab)) +
  stat_smooth(aes(x=stab.stepps, y=stab.refab), method="lm") +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Data_Biomass_v_Diversity.png"), height=6, width=6, units="in", res=320)
ggplot(data=dat.sites) +
  geom_point(aes(x=H.prime, y=stab.refab)) +
  stat_smooth(aes(x=H.prime, y=stab.refab), method="lm") +
  theme_bw()
dev.off()
# -----------

# -------------------------------------------
# -------------------------------------------
