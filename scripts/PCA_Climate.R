# --------------------------------------------
# Script Description
# --------------------------------------------
# Purpose: Do an ordination analysis (PCA) on climate data to identify environmental 
#          ecotones from temperature, preciptiation. The PCs and loadings of sites 
#          from the environment ordination will be compared to ecosystem states like 
#          composition and biomass.
#
# Analysis Details:
# - Environment (Driver) Variables: 
#     - Water Holding Capactity (calculated from depth, % sand, % clay)
#        - could break down and run these 3 factors independently rather than combining
#     - Temperature: Mean annual, summer (June-July-August) 
#     - Precipitation: Mean annual, summer (JJA)
# - Spatial Extent:
#    - All sites in ecosystem model runs.
#    * maybe also run with full spatial coverage (includes sites many models didn't do)
# - Temporal Extent:
#    - Settlement Era: 1800-1850 mean
#    - Empirical Climate: 1901-1930 mean
#
#
# Workflow
# 0. Define file paths etc
# 1. Read in & format environemnt data (temp, precip, WHC)
#    - Annual/Summer Mean, 1800-1850/1901-1930
# 2. Perform Oridination (PCA)
#    - Extract loadings & PCs
# 3. Save output, 
#    - generate & save a couple figures
# --------------------------------------------

# --------------------------------------------
# 0. Define file paths etc
# --------------------------------------------
# Path to github repository/working directory
path.repo <- "~/Desktop/Research/EcosystemStability/"
setwd(path.repo)

# Path to where data are
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"
path.met  <- "~/Desktop/Research/PalEON_MIP_Region/phase2_met_regional_v2_monthly/"
path.soil <- "~/Dropbox/PalEON_CR/env_regional/phase2_env_drivers_v2/soil/"

# This is the path to my pdsi code is, which has the formula for calculating WHC
# This is currently part of my met ensemble code and code 
# and can be found here: https://github.com/PalEON-Project/modeling_met_ensemble.git
path.pdsi <- "~/Dropbox/PalEON_CR/met_ensemble/scripts/"


# Set up some time variables just to help with indexing
yrs <- 850:2010
mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)
# --------------------------------------------


# --------------------------------------------
# Reading in & extracting data
# --------------------------------------------
# loading the key to dimensions in the .RDS files; this is Christy's sheet generated
# this loads an object called "paleon"
paleon <- read.csv(file.path(path.repo, "data/paleon_site_info.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
paleon <- paleon[,c("lon", "lat", "latlon", "umw", "x", "y")]
summary(paleon)

# Extracting soil data for the appropriate sites 
# Note: We need this for the "upper" and "lower" soil layers
library(ncdf4)
tair.nc <- nc_open(file.path(path.met, "tair.nc"))
precipf.nc <- nc_open(file.path(path.met, "precipf.nc"))
sand.t <- nc_open(file.path(path.soil, "paleon_soil_t_sand.nc"))
sand.s <- nc_open(file.path(path.soil, "paleon_soil_s_sand.nc"))
clay.t <- nc_open(file.path(path.soil, "paleon_soil_t_clay.nc"))
clay.s <- nc_open(file.path(path.soil, "paleon_soil_s_clay.nc"))
depth  <- nc_open(file.path(path.soil, "paleon_soil_soil_depth.nc"))

# Getting the lat/lon index for each point to extract soil & met
# Note: for some reason the met & the soil are not in the same order
lon <- ncvar_get(sand.t, "longitude")
lat <- ncvar_get(sand.t, "latitude")
lon2 <- ncvar_get(tair.nc, "lon")
lat2 <- ncvar_get(tair.nc, "lat")


# Creating a dataframe to plug all sites into
paleon.all <- data.frame(lon=rep(lon, length=length(lon)*length(lat)), lat=rep(lat, each=length(lon)))
paleon.all$latlon <- as.factor(paste0("lat", paleon.all$lat, "lon", paleon.all$lon))
summary(paleon.all)


# Merging the subset and the full region together; we can always separate this back out later
paleon.all <- merge(paleon.all, paleon, all=T)
summary(paleon.all)

# Determining where we should be pulling to get the settlement-era & empirical climates
set.all <- which(time.mos$year>=1800 & time.mos$year<=1850)
set.jja <- which(set.all %in% which(time.mos$year>=1800 & time.mos$year<=1850 & time.mos$month>=6 & time.mos$month<=8))
cru.all <- which(time.mos$year>=1901 & time.mos$year<=1930)
cru.jja  <- which(cru.all %in% which(time.mos$year>=1901 & time.mos$year<=1930 & time.mos$month>=6 & time.mos$month<=8))

# Lets convert precip (kg/m2/s = mm/s) to total mm per yr or summer
library(lubridate)
sec_2_yr <- 60*60*24*365.25
sec_2_smr <- 60*60*24*sum(days_in_month(6:8))

# This is the ugly way of doing it; could do it prettier since we're doing the whole region, but whatever
for(i in 1:nrow(paleon.all)){
  # Extracting soil info
  x.ind  <- which(lon == paleon.all[i,"lon"])
  y.ind  <- which(lat == paleon.all[i,"lat"])
  paleon.all[i,"sand.t"] <- ncvar_get(sand.t, "t_sand", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"sand.s"] <- ncvar_get(sand.s, "s_sand", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"clay.t"] <- ncvar_get(clay.t, "t_clay", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"clay.s"] <- ncvar_get(clay.s, "s_clay", c(x.ind, y.ind), c(1,1))
  paleon.all[i,"depth"]  <- ncvar_get(depth , "soil_depth", c(x.ind, y.ind), c(1,1)) # cm

  # Extracting met data
  # Convert to 
  x.ind2 <- which(lon2 == paleon.all[i,"lon"])
  y.ind2 <- which(lat2 == paleon.all[i,"lat"])

  tair.set   <- ncvar_get(tair.nc   , "tair"   , c(x.ind2, y.ind2, min(set.all)), c(1,1,length(set.all)))
  precip.set <- ncvar_get(precipf.nc, "precipf", c(x.ind2, y.ind2, min(set.all)), c(1,1,length(set.all)))
  paleon.all[i,"tair.yr.set"   ] <- mean(tair.set)
  paleon.all[i,"tair.jja.set"  ] <- mean(tair.set[set.jja])
  paleon.all[i,"precip.yr.set" ] <- mean(precip.set)
  paleon.all[i,"precip.jja.set"] <- mean(precip.set[set.jja])

  tair.cru   <- ncvar_get(tair.nc   , "tair"   , c(x.ind2, y.ind2, min(cru.all)), c(1,1,length(cru.all)))
  precip.cru <- ncvar_get(precipf.nc, "precipf", c(x.ind2, y.ind2, min(cru.all)), c(1,1,length(cru.all)))
  paleon.all[i,"tair.yr.cru"   ] <- mean(tair.cru)
  paleon.all[i,"tair.jja.cru"  ] <- mean(tair.cru[cru.jja])
  paleon.all[i,"precip.yr.cru" ] <- mean(precip.cru)
  paleon.all[i,"precip.jja.cru"] <- mean(precip.cru[cru.jja])
}

# Get rid of the cells that are water
paleon.all <- paleon.all[!is.na(paleon.all$sand.t),]
summary(paleon.all)

# To convert precip (kg/m2/s = mm/s) to total mm per yr or summer
library(lubridate)
sec_2_yr <- 60*60*24*365.25
sec_2_smr <- 60*60*24*sum(days_in_month(6:8))



# Calculating Water Holding Capacity for each site from texture
source(file.path(path.pdsi, "calc.awc.R"))
paleon.all$awc.t <- calc.awc(paleon.all$sand.t, paleon.all$clay.t)
paleon.all$awc.s <- calc.awc(paleon.all$sand.s, paleon.all$clay.s)

# Calculating water hodling capacity **in cm**
# If our soil depth is greater than 30 cm, we assume a 30 cm topsoil; if our soil depth is less than 
# 30 cm, we assume all but 1 cm of that is topsoil (yeah, that's not great, but it's what I came up with)
paleon.all$whc.t   <- ifelse(paleon.all$depth > 30, paleon.all$awc.t * 30, paleon.all$awc.t * paleon.all$depth-1) # 30 cm depth
paleon.all$whc.s   <- ifelse(paleon.all$depth > 30, paleon.all$awc.s * (paleon.all$depth-30), paleon.all$awc.s*1) # depth - 30 cm depth 
paleon.all$whc.tot <- paleon.all$whc.t + paleon.all$whc.s
summary(paleon.all)


# Mapping our environment 
library(ggplot2); library(grid)
plot.clay.t   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=clay.t)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "coral2"), name="Clay 1")
plot.clay.s   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=clay.s)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "coral2"), name="Clay 2")
plot.sand.t   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=sand.t)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "goldenrod2"), name="Sand 1")
plot.sand.s   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=sand.s)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "goldenrod2"), name="Sand 2")
plot.wcap.t   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=whc.t)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "cyan3"), name="WHC 1")
plot.wcap.s   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=whc.s)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "cyan3"), name="WHC 2", limits=c(0,1000))
plot.depth   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=depth)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "tan3"), name="Depth")
plot.whc     <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=whc.tot)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue"), name="Total WHC", limits=c(0,1000))

tair.yr.set   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.yr.set-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red3"), name="Tair.yr")
tair.gs.set   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.jja.set-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red2"), name="Tair.jja")
precip.yr.set <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.yr.set*sec_2_yr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.yr")
precip.gs.set <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.jja.set*sec_2_smr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.jja")

tair.yr.cru   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.yr.cru-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red3"), name="Tair.yr")
tair.gs.cru   <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=tair.jja.cru-273.15)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "red2"), name="Tair.jja")
precip.yr.cru <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.yr.cru*sec_2_yr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.yr")
precip.gs.cru <- ggplot(data=paleon.all) + theme_bw() + geom_raster(aes(x=lon, y=lat, fill=precip.jja.cru*sec_2_smr)) + coord_equal(expand=0) + scale_fill_gradientn(colours=c("gray20", "blue2"), name="Precip.jja")


# ---------------
# Some exploratory graphs to see if we'd expect many differences; 
# in general, summer shows stronger spatial patterns
# ---------------
pdf("figures/Environmental_Drivers.pdf")
# The non-climate environment
grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2)))
print(plot.clay.t, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot.clay.s, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(plot.sand.t, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(plot.sand.s, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(plot.wcap.t, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(plot.wcap.s, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
print(plot.depth , vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
print(plot.whc   , vp = viewport(layout.pos.row = 4, layout.pos.col = 2))


# Compare year vs. summer: settlement era
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(tair.yr.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(tair.gs.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(precip.yr.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(precip.gs.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

# Compare year vs. summer: CRUNCEP
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(tair.yr.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(tair.gs.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(precip.yr.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(precip.gs.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

# Compare Settlement vs CRU: Year
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(tair.yr.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(tair.yr.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(precip.yr.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(precip.yr.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

# Compare Settlement vs CRU: Summer
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(tair.gs.set  , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(tair.gs.cru  , vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(precip.gs.set, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(precip.gs.cru, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()
# ---------------
# --------------------------------------------


# --------------------------------------------
# 2. Perform Oridination (PCA)
#
# We're going to do a couple versions for now:
# 2.a. All climate & env variables
# 2.b. Summer climate + top env + depth
# 2.c. Model sites only
# --------------------------------------------
library(labdsv)
set.seed(1333)

pc.all.set1 <- princomp(~ sand.t + sand.s + clay.t + clay.s + depth + awc.t + awc.s + whc.tot + 
                      tair.yr.set + tair.jja.set + precip.yr.set + precip.jja.set,
                    data=paleon.all, cor=T, scores=T)
pc.all.set2 <- princomp(~ sand.t + clay.t + depth + whc.tot + tair.jja.set + precip.jja.set,
                    data=paleon.all, cor=T, scores=T)
summary(pc.all.set1)
summary(pc.all.set2)
pc.all.set1$loadings

summary(pc.set1$scores)


pc.all.cru1 <- princomp(~ sand.t + sand.s + clay.t + clay.s + depth + awc.t + awc.s + whc.tot  + 
                      tair.yr.cru + tair.jja.cru + precip.yr.cru + precip.jja.cru,
                    data=paleon.all, cor=T, scores=T)
pc.all.cru2 <- princomp(~ sand.t + clay.t + depth + whc.tot + tair.jja.cru + precip.jja.cru,
                    data=paleon.all, cor=T, scores=T)
summary(pc.all.cru1)
summary(pc.all.cru2)
pc.all.cru1$loadings
print(pc.all.cru1$loadings, cutoff=0.00)


# Running things just at the sites the models ran at
pc.mod.set1 <- princomp(~ sand.t + sand.s + clay.t + clay.s + depth + awc.t + awc.s + whc.tot  + 
                          tair.yr.set + tair.jja.set + precip.yr.set + precip.jja.set,
                        data=paleon.all[!is.na(paleon.all$umw),], cor=T, scores=T)
pc.mod.cru1 <- princomp(~ sand.t + sand.s + clay.t + clay.s + depth + awc.t + awc.s + whc.tot  + 
                          tair.yr.cru + tair.jja.cru + precip.yr.cru + precip.jja.cru,
                        data=paleon.all[!is.na(paleon.all$umw),], cor=T, scores=T)
summary(pc.mod.set1)
summary(pc.mod.cru1)

pc.all.set1$loadings[,1:3]
pc.mod.set1$loadings[,1:3]

# Full Region
paleon.all$set1.pc1 <- pc.all.set1$scores[,1]
paleon.all$set1.pc2 <- pc.all.set1$scores[,2]
paleon.all$set1.pc3 <- pc.all.set1$scores[,3]

paleon.all$cru1.pc1 <- pc.all.cru1$scores[,1]
paleon.all$cru1.pc2 <- pc.all.cru1$scores[,2]
paleon.all$cru1.pc3 <- pc.all.cru1$scores[,3]

# Just model sites
paleon.all[!is.na(paleon.all$umw),"set2.pc1"] <- pc.mod.set1$scores[,1]
paleon.all[!is.na(paleon.all$umw),"set2.pc2"] <- pc.mod.set1$scores[,2]
paleon.all[!is.na(paleon.all$umw),"set2.pc3"] <- pc.mod.set1$scores[,3]

paleon.all[!is.na(paleon.all$umw),"cru2.pc1"] <- pc.mod.cru1$scores[,1]
paleon.all[!is.na(paleon.all$umw),"cru2.pc2"] <- pc.mod.cru1$scores[,2]
paleon.all[!is.na(paleon.all$umw),"cru2.pc3"] <- pc.mod.cru1$scores[,3]

# --------------------------------------------


# --------------------------------------------
# 3. Save output, 
#    - generate & save a couple figures
# --------------------------------------------
write.csv(paleon.all, "data/PCA_Scores_Environment_AllRegion.csv", row.names=F, eol="\r\n")

# Manually calculating the percentage of variance
set1.POV <- pc.all.set1$sdev^2/sum(pc.all.set1$sdev^2)
cru1.POV <- pc.all.cru1$sdev^2/sum(pc.all.cru1$sdev^2)
set2.POV <- pc.all.set2$sdev^2/sum(pc.all.set2$sdev^2)
cru2.POV <- pc.all.cru2$sdev^2/sum(pc.all.cru2$sdev^2)


set.pc1   <- ggplot(data=paleon.all) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(aes(x=lon, y=lat, fill=set1.pc1)) + 
  scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1", limits=quantile(paleon.all$set1.pc1, c(0.001, 0.999))) + 
  ggtitle(paste0("Region PCA, Settlement Era PC1 (", round(set1.POV["Comp.1"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
set.pc2   <- ggplot(data=paleon.all) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(aes(x=lon, y=lat, fill=set1.pc2)) + 
  scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="Set PC2", limits=quantile(paleon.all$set1.pc2, c(0.001, 0.999))) + 
  ggtitle(paste0("Region PCA, Settlement Era PC2 (", round(set1.POV["Comp.2"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
set.pc3   <- ggplot(data=paleon.all) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(aes(x=lon, y=lat, fill=set1.pc3)) + 
  scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="Set PC3", limits=quantile(paleon.all$set1.pc3, c(0.001, 0.999))) + 
  ggtitle(paste0("Region PCA, Settlement Era PC3 (", round(set1.POV["Comp.3"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))

cru.pc1   <- ggplot(data=paleon.all) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(aes(x=lon, y=lat, fill=cru1.pc1)) + 
  scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="cru PC1", limits=quantile(paleon.all$cru1.pc1, c(0.001, 0.999))) + 
  ggtitle(paste0("Region PCA, CRUNCEP PC1 (", round(cru1.POV["Comp.1"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5))
cru.pc2   <- ggplot(data=paleon.all) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(aes(x=lon, y=lat, fill=cru1.pc2)) + 
  scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="cru PC2", limits=quantile(paleon.all$cru1.pc2, c(0.001, 0.999))) + 
  ggtitle(paste0("Region PCA, CRUNCEP PC2 (", round(cru1.POV["Comp.2"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5))
cru.pc3   <- ggplot(data=paleon.all) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(aes(x=lon, y=lat, fill=cru1.pc3)) +
  scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="cru PC3", limits=quantile(paleon.all$cru1.pc3, c(0.001, 0.999))) + 
  ggtitle(paste0("Region PCA, CRUNCEP PC3 (", round(cru1.POV["Comp.3"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))


set2.pc1   <- ggplot(data=paleon.all[!is.na(paleon.all$umw),]) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(data=paleon.all, aes(x=lon, y=lat), fill="gray50") +
  geom_point(aes(x=lon, y=lat, color=set2.pc1), size=2, alpha=0.8)  + 
  scale_color_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
  ggtitle(paste0("Region PCA, Settlement Era PC1 (", round(set2.POV["Comp.1"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
set2.pc2   <- ggplot(data=paleon.all[!is.na(paleon.all$umw),]) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(data=paleon.all, aes(x=lon, y=lat), fill="gray50") +
  geom_point(aes(x=lon, y=lat, color=set2.pc2), size=2, alpha=0.8)  + 
  scale_color_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
  ggtitle(paste0("Region PCA, Settlement Era PC2 (", round(set2.POV["Comp.2"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
set2.pc3   <- ggplot(data=paleon.all[!is.na(paleon.all$umw),]) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(data=paleon.all, aes(x=lon, y=lat), fill="gray50") +
  geom_point(aes(x=lon, y=lat, color=set2.pc3), size=2, alpha=0.8)  + 
  scale_color_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC3") + 
  ggtitle(paste0("Region PCA, Settlement Era PC3 (", round(set2.POV["Comp.3"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))

cru2.pc1   <- ggplot(data=paleon.all[!is.na(paleon.all$umw),]) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(data=paleon.all, aes(x=lon, y=lat), fill="gray50") +
  geom_point(aes(x=lon, y=lat, color=cru2.pc1), size=2, alpha=0.8)  + 
  scale_color_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
  ggtitle(paste0("Region PCA, CRUNCEP Era PC1 (", round(cru2.POV["Comp.1"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
cru2.pc2   <- ggplot(data=paleon.all[!is.na(paleon.all$umw),]) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(data=paleon.all, aes(x=lon, y=lat), fill="gray50") +
  geom_point(aes(x=lon, y=lat, color=cru2.pc2), size=2, alpha=0.8)  + 
  scale_color_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
  ggtitle(paste0("Region PCA, CRUNCEP Era PC2 (", round(cru2.POV["Comp.2"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
cru2.pc3   <- ggplot(data=paleon.all[!is.na(paleon.all$umw),]) + theme_bw() + coord_equal(expand=0) + 
  geom_raster(data=paleon.all, aes(x=lon, y=lat), fill="gray50") +
  geom_point(aes(x=lon, y=lat, color=cru2.pc3), size=2, alpha=0.8)  + 
  scale_color_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC3") + 
  ggtitle(paste0("Region PCA, CRUNCEP PC3 (", round(cru2.POV["Comp.3"]*100,1),"%)")) + 
  theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))


# Graphing the site scores
pdf("figures/PC_Scores_Region.pdf", height=6, width=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
print(set.pc1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(set.pc2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(set.pc3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(set2.pc1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(set2.pc2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(set2.pc3, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))
print(cru.pc1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(cru.pc2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(cru.pc3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(cru2.pc1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(cru2.pc2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(cru2.pc3, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
dev.off()


# install.packages("ggfortify")
# library(devtools)
# install_github("ggbiplot", "vqv")

# library(ggbiplot)
# ggbiplot(pc.all.set1, var.scale=2, label.fontface=2, loadings.label.colour="blue")
# install.packages("ggbiplot")

# Making a dataframe to plot the loadigns in ggplot
# adapted from: https://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
set1.loadings <- data.frame(var=row.names(pc.all.set1$loadings), pc.all.set1$loadings[,1:4])
mult <- min( (max(paleon.all$set1.pc2)-min(paleon.all$set1.pc2))/(max(set1.loadings$Comp.2)-min(set1.loadings$Comp.2)),
             (max(paleon.all$set1.pc1)-min(paleon.all$set1.pc1))/(max(set1.loadings$Comp.1)-min(set1.loadings$Comp.1))
             )
set1.loadings$load.PC1 <- 0.7*mult*set1.loadings$Comp.1
set1.loadings$load.PC2 <- 0.7*mult*set1.loadings$Comp.2
set1.loadings$var <- c("sand (top)", "sand (low)", "clay (top)", "clay (low)", "depth", "awc (top)", "awc (low)", "Water Cap", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")

cru1.loadings <- data.frame(var=row.names(pc.all.cru1$loadings), pc.all.cru1$loadings[,1:4])
mult <- min( (max(paleon.all$cru1.pc2)-min(paleon.all$cru1.pc2))/(max(cru1.loadings$Comp.2)-min(cru1.loadings$Comp.2)),
             (max(paleon.all$cru1.pc1)-min(paleon.all$cru1.pc1))/(max(cru1.loadings$Comp.1)-min(cru1.loadings$Comp.1))
)
cru1.loadings$load.PC1 <- 0.7*mult*cru1.loadings$Comp.1
cru1.loadings$load.PC2 <- 0.7*mult*cru1.loadings$Comp.2
cru1.loadings$var <- c("sand (top)", "sand (low)", "clay (top)", "clay (low)", "depth", "awc (top)", "awc (low)", "Water Cap", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")


set2.loadings <- data.frame(var=row.names(pc.mod.set1$loadings), pc.mod.set1$loadings[,1:4])
mult <- min( (max(paleon.all$set2.pc2, na.rm=T)-min(paleon.all$set2.pc2, na.rm=T))/(max(set2.loadings$Comp.2)-min(set2.loadings$Comp.2)),
             (max(paleon.all$set2.pc1, na.rm=T)-min(paleon.all$set2.pc1, na.rm=T))/(max(set2.loadings$Comp.1)-min(set2.loadings$Comp.1))
)
set2.loadings$load.PC1 <- 0.7*mult*set2.loadings$Comp.1
set2.loadings$load.PC2 <- 0.7*mult*set2.loadings$Comp.2
set2.loadings$var <-  c("sand (top)", "sand (low)", "clay (top)", "clay (low)", "depth", "awc (top)", "awc (low)", "Water Cap", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")
cru2.loadings <- data.frame(var=row.names(pc.mod.cru1$loadings), pc.mod.cru1$loadings[,1:4])
mult <- min( (max(paleon.all$cru2.pc2,na.rm=T)-min(paleon.all$cru2.pc2,na.rm=T))/(max(cru2.loadings$Comp.2,na.rm=T)-min(cru2.loadings$Comp.2,na.rm=T)),
             (max(paleon.all$cru2.pc1,na.rm=T)-min(paleon.all$cru2.pc1,na.rm=T))/(max(cru2.loadings$Comp.1,na.rm=T)-min(cru2.loadings$Comp.1,na.rm=T))
)
cru2.loadings$load.PC1 <- 0.7*mult*cru2.loadings$Comp.1
cru2.loadings$load.PC2 <- 0.7*mult*cru2.loadings$Comp.2
cru2.loadings$var <-  c("sand (top)", "sand (low)", "clay (top)", "clay (low)", "depth", "awc (top)", "awc (low)", "Water Cap", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")

# --------------
# Manually creating a biplots
# Weights cutoff = 75th percentile for PC1 & 2 (top 25% of weights)
# --------------
set1.cutoff <- quantile(c(abs(set1.loadings$Comp.1), abs(set1.loadings$Comp.2)), 0.75)
biplot.set1 <- ggplot() +
  theme_bw() +# coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon.all, aes(x=set1.pc1, y=set1.pc2), size=0.2, color="black") +
  geom_segment(data=set1.loadings[abs(set1.loadings$Comp.1)>set1.cutoff | abs(set1.loadings$Comp.2)>set1.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=set1.loadings[abs(set1.loadings$Comp.1)>set1.cutoff | abs(set1.loadings$Comp.2)>set1.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon.all$set1.pc1, c(0.001, 0.999))) +
  scale_y_continuous(name=paste0("PC2 (", round(POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon.all$set1.pc2, c(0.001, 0.999))) +
  ggtitle("Region, Settlement-Era Climate") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))

cru1.cutoff <- quantile(c(abs(cru1.loadings$Comp.1), abs(cru1.loadings$Comp.2)), 0.75)
biplot.cru1 <- ggplot() +
  theme_bw() + #coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon.all, aes(x=cru1.pc1, y=cru1.pc2), size=0.2, color="black") +
  geom_segment(data=cru1.loadings[abs(cru1.loadings$Comp.1)>cru1.cutoff | abs(cru1.loadings$Comp.2)>cru1.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=cru1.loadings[abs(cru1.loadings$Comp.1)>cru1.cutoff | abs(cru1.loadings$Comp.2)>cru1.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=-1.5, hjust=0.5, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon.all$cru1.pc1, c(0.001, 0.999))) +
  scale_y_continuous(name=paste0("PC2 (", round(POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon.all$cru1.pc2, c(0.001, 0.999))) +
  ggtitle("Region, CRUNCEP Climate") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))

set2.cutoff <- quantile(c(abs(set2.loadings$Comp.1), abs(set2.loadings$Comp.2)), 0.75)
biplot.set2 <- ggplot() +
  theme_bw() + #coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon.all, aes(x=set2.pc1, y=set2.pc2), size=0.2, color="black") +
  geom_segment(data=set2.loadings[abs(set2.loadings$Comp.1)>set2.cutoff | abs(set2.loadings$Comp.2)>set2.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=set2.loadings[abs(set2.loadings$Comp.1)>set2.cutoff | abs(set2.loadings$Comp.2)>set2.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1)) +
  scale_y_continuous(name=paste0("PC2 (", round(POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1)) +
  ggtitle("Model Sites, Settlement-Era Climate") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))

cru2.cutoff <- quantile(c(abs(cru2.loadings$Comp.1), abs(cru2.loadings$Comp.2)), 0.75)
biplot.cru2 <- ggplot() +
  theme_bw() + #coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon.all, aes(x=cru2.pc1, y=cru2.pc2), size=0.2, color="black") +
  geom_segment(data=cru2.loadings[abs(cru2.loadings$Comp.1)>cru2.cutoff | abs(cru2.loadings$Comp.2)>cru2.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=cru2.loadings[abs(cru2.loadings$Comp.1)>cru2.cutoff | abs(cru2.loadings$Comp.2)>cru2.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=1.5, hjust=0.5, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1)) +
  scale_y_continuous(name=paste0("PC2 (", round(POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1)) +
  ggtitle("Model Sites, CRUNCEP Climate") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))


pdf("figures/PCA_Loadings.pdf", height=8, width=8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(biplot.set1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(biplot.cru1, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(biplot.set2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(biplot.cru2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()
# --------------------------------------------
