# -------------------------------------------
# Data-Climate Comparisons: 
# 
# H0: Climate change causes ecosystem change, as seen through shifts in 
#     biomass and composition (consistent & proportional responses; data)
#        Prediction: Locations that have experienced the most changes in 
#                    climate show the most changes in biomass and composition
#        Prediction: Places with the greatest change in composition also 
#                    have the greatest change in biomass (Zhang 2018 Nature Climate Chante)
#
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

stepps[which(stepps$deriv.abs==max(stepps$deriv.abs, na.rm=T)),]

refab.stepps <- refab[refab$lon<max(stepps[stepps$lat<45, "lon"]),]
refab.stepps[which(refab.stepps$refab.mean.slope.abs==max(refab.stepps$refab.mean.slope.abs, na.rm=T)),]
