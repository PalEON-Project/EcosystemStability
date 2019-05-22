# -------------------------------------------
# Model-based analyses:
# 
# 1. Do models show the same trends as the data:
#    a. Biomass versus climate
#    b. Fcomp versus climate
#    c. Biomass versus Fcomp/Diversity
# 2. Can the models help explain the lack of correlation among 
#    climate, biomass & composition?
# -------------------------------------------
rm(list=ls())


# -------------------------------------------
# Load libraries; set file paths
# -------------------------------------------
library(ncdf4)
library(ggplot2); library(gridExtra); library(scales); library(grid)
library(mgcv)
library(plyr); library(parallel)

# Setting the path to this repository
# path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
# setwd(path.repo)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data/"

dat.colors <- data.frame(model=c("LBDA", "ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"),
                         type=c(rep("emipirical", 4), rep("model", 7)),
                         color=c(rep("#000000", 3), "#7F7F7F", rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
dat.colors

# -------------------------------------------

# -------------------------------------------
# Define some useful variables
# -------------------------------------------
# Set up some time variables just to help with indexing
yrs <- 850:2010
mos <- 1:12
time.mos <- data.frame(year=rep(yrs, each=length(mos)), month=mos)
head(time.mos)

us <- map_data("state")
# -------------------------------------------

# -------------------------------------------
# 1. Read in data; 
# -------------------------------------------
# ------------
# Most model variables (compsition excluded)
# ------------
stab.models <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_Models_100.csv"))
summary(stab.models)
vars.use <- c("gpp", "npp", "nee", "lai", "bm")
type.use <- "diff."

# Get rid of sites with no biomass because they're being weird
stab.models[!is.na(stab.models$mean.bm) & stab.models$mean.bm==0,]
stab.models[!is.na(stab.models$mean.bm) & stab.models$mean.bm==0,!names(stab.models) %in% c("lat", "lon", "latlon", "Model")] <- NA

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

fcomp.diff  <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_Models_FCOMP_100.csv"))

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
  mod <- paste(mod)
  for(lon in unique(fcomp[fcomp$Model==mod & !is.na(fcomp$lon), "lon"])){
    for(lat in unique(fcomp[fcomp$Model==mod & fcomp$lon==lon & !is.na(fcomp$lat), "lat"])){
      
      # Use the same threshold used for STEPPS
      df.fcomp <- fcomp.diff[fcomp.diff$Model==mod & fcomp.diff$lon==lon & fcomp.diff$lat==lat & fcomp.diff$pft.mean>1e-3,]
      if(nrow(df.fcomp)==0) next
      
      # Make sure to get rid of the "total" column
      df.fcomp <- df.fcomp[df.fcomp$PFT!="Total",]
      
      ind.fcomp <- which(df.fcomp$pft.mean==max(df.fcomp$pft.mean))[1]
      ind.stab <- which(fcomp$Model==mod & fcomp$lon==lon & fcomp$lat==lat)
      
      fcomp[ind.stab, "pft.abs" ] <- paste(df.fcomp$PFT[ind.fcomp])
      fcomp[ind.stab, "diff.abs"] <- df.fcomp$pft.diff.abs[ind.fcomp]
      fcomp[ind.stab, "val.mean"] <- df.fcomp$pft.mean[ind.fcomp]
      
      # Adding the dominant PFT into the model output
      models.long[models.long$Model==mod & models.long$lon==lon & models.long$lat==lat,"pft.abs"] <- paste(df.fcomp$PFT[ind.fcomp])
    }
  }
}
# fcomp$pft.abs <- as.factor(fcomp$pft.abs)
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
stab.met   <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_Drivers_100.csv"))
stab.met$var.pdsi <- stab.met$pdsi.diff/mean(stab.met$pdsi.diff)
summary(stab.met)

hist(log(stab.met$var.pdsi))

ggplot(data=stab.met) +
  coord_equal() +
  geom_tile(aes(x=lon, y=lat, fill=log(var.pdsi)))


# Merge the met drivers into models.long
models.long <- merge(models.long, stab.met[,c("lon", "lat", "pdsi.diff", "tair.diff", "precip.diff", "pdsi.mean", "tair.mean", "precip.mean")])
models.long$var <- factor(models.long$var, levels=c("gpp", "nee", "npp", "lai", "bm", "fcomp"))
summary(models.long)
# ------------

# ------------
# Calculate our final stability index value
# ------------
# Log won't work with stability=0, so we need to turn it to somethign *very* tiny
# NOTE: This will cause some very large outliers in our data, but we can trim those
models.long[models.long$diff.abs==0, "diff.abs"] <- 1e-20 
# quantile(models.long$diff.abs, 0.0001)
for(v in unique(models.long$var)){
  for(mod in unique(models.long$Model)){
    dat.ind <- which(models.long$var==v & models.long$Model==mod & !is.na(models.long$diff.abs))
    dat.tmp <- models.long[dat.ind, ]
    
    # models.long[dat.ind, "stability" ] <- -log(dat.tmp$diff.abs/abs(mean(dat.tmp$val.mean)))
    # models.long[dat.ind, "stab.rel"] <- (models.long[dat.ind, "stability" ]) / mean(models.long[dat.ind, "stability" ]) 

    models.long[dat.ind, "stab.pdsi" ] <- -log(dat.tmp$pdsi.diff/abs(mean(dat.tmp$pdsi.diff)))
    
    models.long[dat.ind, "variability" ] <- dat.tmp$diff.abs/abs(mean(dat.tmp$val.mean))
    models.long[dat.ind, "var.pdsi" ] <- dat.tmp$pdsi.diff/abs(mean(dat.tmp$pdsi.diff))
    models.long[dat.ind, "var.rel"] <- (models.long[dat.ind, "variability" ]) / mean(models.long[dat.ind, "variability" ]) 
    
  }
}
summary(models.long)
models.long$type <- "Model"

write.csv(models.long,  file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_long.csv"), row.names=F)
# ------------

# ------------
# Read in the empirical data from script #4
# ------------
dat.emp <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"))
dat.emp$var <- as.factor(ifelse(dat.emp$dataset=="ReFAB", "bm", "fcomp"))
dat.emp$type <- "Empirical"
dat.emp$dataset2 <- car::recode(dat.emp$dataset2, "'STEPPS'='STEPPS-UMW'")
summary(dat.emp)
# ------------

models.long <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_long.csv"))

# ------------
# Combining emprirical and model datasets together
# ------------
mod.v.dat <- data.frame(lat=c(models.long$lat[models.long$var %in% c("bm", "fcomp")], dat.emp$lat),
                        lon=c(models.long$lon[models.long$var %in% c("bm", "fcomp")], dat.emp$lon),
                        type=c(paste(models.long$type[models.long$var %in% c("bm", "fcomp")]), dat.emp$type),
                        var=c(paste(models.long$var[models.long$var %in% c("bm", "fcomp")]), paste(dat.emp$var)),
                        source=c(paste(models.long$Model[models.long$var %in% c("bm", "fcomp")]), paste(dat.emp$dataset2)),
                        # stab.ecosys=c(models.long$stability[models.long$var %in% c("bm", "fcomp")], dat.emp$stab.ecosys),
                        # stab.pdsi=c(models.long$stab.pdsi[models.long$var %in% c("bm", "fcomp")], dat.emp$stab.pdsi),
                        var.ecosys=c(models.long$variability[models.long$var %in% c("bm", "fcomp")], dat.emp$var.ecosys),
                        var.pdsi=c(models.long$var.pdsi[models.long$var %in% c("bm", "fcomp")], dat.emp$var.pdsi)
                        )

mod.v.dat$source <- factor(mod.v.dat$source, levels=c("ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))
summary(mod.v.dat)

# Making an additional one that will let us make a pretty scatter plot
dat.refab <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_ReFAB.csv"))
summary(dat.refab)
mod.dat2 <- data.frame(lat=c(mod.v.dat$lat, models.long$lat[models.long$var=="fcomp"], dat.refab$lat),
                       lon=c(mod.v.dat$lon, models.long$lon[models.long$var=="fcomp"], dat.refab$lon),
                       type=c(paste(mod.v.dat$type), paste(models.long$type[models.long$var=="fcomp"]), rep("Empirical", nrow(dat.refab))),
                       source=c(paste(mod.v.dat$source), paste(models.long$Model[models.long$var=="fcomp"]), rep("ReFAB", nrow(dat.refab))),
                       var1=c(rep("PDSI", nrow(mod.v.dat)), rep("Composition", nrow(models.long[models.long$var=="fcomp",])+nrow(dat.refab))),
                       var2=c(paste(mod.v.dat$var), rep("bm", nrow(models.long[models.long$var=="fcomp",])+nrow(dat.refab))),
                       # stability1 = c(mod.v.dat$stab.pdsi, models.long$stability[models.long$var=="fcomp"], dat.refab$stab.stepps.1k),
                       # stability2 = c(mod.v.dat$stab.ecosys, rep(NA, length(models.long$stability[models.long$var=="fcomp"])), dat.refab$stab.refab.1k),
                       variability1 = c(mod.v.dat$var.pdsi, models.long$variability[models.long$var=="fcomp"], dat.refab$var.stepps.1k),
                       variability2 = c(mod.v.dat$var.ecosys, rep(NA, length(models.long$variability[models.long$var=="fcomp"])), dat.refab$var.refab.1k)
                       )
summary(mod.dat2)

# Looping through to pair model BM & Fcomp stability
for(i in 1:nrow(mod.dat2)){
  if(mod.dat2$type[i]=="Empirical" | !is.na(mod.dat2$variability2[i])) next
  
  # val.fill <- models.long[models.long$var=="bm" & models.long$Model==paste(mod.dat2$source[i]) & models.long$lat==mod.dat2$lat[i] & models.long$lon==mod.dat2$lon[i], "stability"]
  val.fill2 <- models.long[models.long$var=="bm" & models.long$Model==paste(mod.dat2$source[i]) & models.long$lat==mod.dat2$lat[i] & models.long$lon==mod.dat2$lon[i], "variability"]
  
  if(length(val.fill2)==0) next 
  
  # mod.dat2[i,"stability2"] <- val.fill
  mod.dat2[i,"variability2"] <- val.fill2
}
mod.dat2$var2 <- car::recode(mod.dat2$var2, "'bm'='Biomass'; 'fcomp'='Composition'")
mod.dat2$source <- factor(mod.dat2$source, levels=c("ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))
summary(mod.dat2)

# write.csv(dat.sites, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Diversity_Data.csv"), row.names=F)

write.csv(mod.dat2, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_and_Data.csv"), row.names=F)
# ------------
# -------------------------------------------


# -------------------------------------------
# Comparing stabiltiy in models versus data
# -------------------------------------------
mod.dat2 <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_and_Data.csv"))
mod.dat2$source <- factor(mod.dat2$source, levels=c("ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "TRIFFID"))
mod.dat2$var1 <- car::recode(mod.dat2$var1, "'PDSI'='Drought'")
mod.dat2$var1 <- factor(mod.dat2$var1, levels=c("Drought", "Composition"))
mod.dat2$var2 <- factor(mod.dat2$var2, levels=c("Composition", "Biomass"))
# mod.dat2$var2
# mod.dat2[!is.na(mod.dat2$variability1) & mod.dat2$variability1>20,]
# mod.dat2[!is.na(mod.dat2$stability1) & mod.dat2$stability1>20, c("stability1", "variability1")] <- NA
# mod.dat2[!is.na(mod.dat2$variability2) & mod.dat2$stability2>20, c("stability2", c("variability2"))] <- NA

dat.colors$model <- factor(dat.colors$model, levels=c("drivers", "drivers-modified", "LBDA", "ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "TRIFFID"))
dat.colors <- dat.colors[order(dat.colors$model),]

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Model_v_Data.png"), height=6, width=6, units="in", res=320)
ggplot(dat=mod.dat2) +
  facet_grid(var2 ~ var1, scales="free", switch="both") +
  geom_point(aes(x=log(variability1), y=log(variability2), color=source), size=0.1, alpha=0.25) +
  stat_smooth(aes(x=log(variability1), y=log(variability2), color=source, fill=source), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  scale_x_continuous(name="Log Relative Variability") +
  scale_y_continuous(name="Log Relative Variability") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=rel(1)),
        # axis.title = element_blank()
        axis.title = element_text(size=rel(1.25))
        ) +
  theme(legend.position=c(0.75, 0.75),
        legend.title=element_blank())
dev.off()  

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Model_v_Data_fixedaxes.png"), height=6, width=6, units="in", res=320)
ggplot(dat=mod.dat2) +
  facet_grid(var2 ~ var1, scales="free", switch="both") +
  geom_point(aes(x=log(variability1), y=log(variability2), color=source), size=0.15, alpha=0.25) +
  stat_smooth(aes(x=log(variability1), y=log(variability2), color=source, fill=source), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  scale_x_continuous(name="Log Relative Variability") +
  scale_y_continuous(name="Log Relative Variability") +
  coord_cartesian(ylim=c(-10,-1), xlim=c(-10, 2.5)) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=rel(1)),
        # axis.title = element_blank()
        axis.title = element_text(size=rel(1.25))
  ) +
  theme(legend.position=c(0.75, 0.75),
        legend.title=element_blank())
dev.off()  
summary(models.long)

# ------------
# Biomass-PDSI Comparisons
# ------------
md.comp.bm.mean <- lm(log(var.ecosys) ~ source, data=mod.v.dat[mod.v.dat$var=="bm",])
summary(md.comp.bm.mean)

md.comp.bm.sens1 <- lm(log(var.ecosys) ~ log(var.pdsi)*source, data=mod.v.dat[mod.v.dat$var=="bm",])
summary(md.comp.bm.sens1)

md.comp.bm.sens2 <- lm(log(var.ecosys) ~ log(var.pdsi)*source-log(var.pdsi), data=mod.v.dat[mod.v.dat$var=="bm",])
summary(md.comp.bm.sens2)
# ------------

# ------------
# fcomp-PDSI Comparisons
# ------------
mod.v.dat$region <- ifelse(mod.v.dat$source=="STEPPS-UMW", "UMW", ifelse(mod.v.dat$source=="STEPPS-NEUS", "NEUS", "full"))
mod.v.dat$region <- factor(mod.v.dat$region, levels=c("UMW", "NEUS", "full"))
mod.v.dat$source2 <- as.factor(ifelse(substr(mod.v.dat$source, 1, 6)=="STEPPS", "STEPPS", paste(mod.v.dat$source)))
mod.v.dat$source2 <- factor(mod.v.dat$source2, levels=c("STEPPS", "ReFAB", "ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "TRIFFID"))
summary(mod.v.dat)

# md.comp.fcomp.mean <- lm(log(var.ecosys) ~ source, data=mod.v.dat[mod.v.dat$var=="fcomp",])
# summary(md.comp.fcomp.mean)

# This is the model that allows STEPPS regions to have different intercepts, but a common slope
md.comp.fcomp.sens1 <- lm(log(var.ecosys) ~ log(var.pdsi)*source2-source2 + source-1, data=mod.v.dat[mod.v.dat$var=="fcomp",])
summary(md.comp.fcomp.sens1)

# Effects parameterization of above
md.comp.fcomp.sens1 <- lm(log(var.ecosys) ~ log(var.pdsi)*source2-source2 + source-1 - log(var.pdsi), data=mod.v.dat[mod.v.dat$var=="fcomp",])
summary(md.comp.fcomp.sens1)

# md.comp.fcomp.sens1 <- lm(log(var.ecosys) ~ log(var.pdsi)*source2 + region-1, data=mod.v.dat[mod.v.dat$var=="fcomp",])
# summary(md.comp.fcomp.sens1)

# md.comp.fcomp.sens2 <- lm(log(var.ecosys) ~ log(var.pdsi)*source-log(var.pdsi)-1, data=mod.v.dat[mod.v.dat$var=="fcomp",])
# summary(md.comp.fcomp.sens2)
# anova(md.comp.fcomp.sens2)

# Looking at mean fcomp variability in models
mean(log(mod.v.dat$var.ecosys[mod.v.dat$source=="ED2"]), na.rm=T); sd(log(mod.v.dat$var.ecosys[mod.v.dat$source=="ED2"]), na.rm=T)
mean(log(mod.v.dat$var.ecosys[mod.v.dat$source=="LPJ-GUESS"]), na.rm=T); sd(log(mod.v.dat$var.ecosys[mod.v.dat$source=="LPJ-GUESS"]), na.rm=T)
mean(log(mod.v.dat$var.ecosys[mod.v.dat$source=="LPJ-WSL"]), na.rm=T); sd(log(mod.v.dat$var.ecosys[mod.v.dat$source=="LPJ-WSL"]), na.rm=T)
mean(log(mod.v.dat$var.ecosys[mod.v.dat$source=="LINKAGES"]), na.rm=T); sd(log(mod.v.dat$var.ecosys[mod.v.dat$source=="LINKAGES"]), na.rm=T)
mean(log(mod.v.dat$var.ecosys[mod.v.dat$source=="TRIFFID"]), na.rm=T); sd(log(mod.v.dat$var.ecosys[mod.v.dat$source=="TRIFFID"]), na.rm=T)
# ------------

# ------------
# Quick test to see if the model regions behave the same or different
# ------------
lon.umw <- range(mod.v.dat[mod.v.dat$source=="STEPPS-UMW", "lon"])
lat.umw <- range(mod.v.dat[mod.v.dat$source=="STEPPS-UMW", "lat"])
lon.neus <- range(mod.v.dat[mod.v.dat$source=="STEPPS-NEUS", "lon"])
lat.neus <- range(mod.v.dat[mod.v.dat$source=="STEPPS-NEUS", "lat"])

mod.umw <- mod.v.dat[mod.v.dat$lon>=lon.umw[1] & mod.v.dat$lon<=lon.umw[2] &
                       mod.v.dat$lat>=lat.umw[1] & mod.v.dat$lat<=lat.umw[2],]
mod.umw$region <- "UMW"

mod.neus <- mod.v.dat[mod.v.dat$lon>=lon.neus[1] & mod.v.dat$lon<=lon.neus[2] &
                       mod.v.dat$lat>=lat.neus[1] & mod.v.dat$lat<=lat.neus[2],]
mod.neus$region <- "NEUS"
mod.region <- rbind(mod.umw, mod.neus)
mod.region$source <- car::recode(mod.region$source, "'STEPPS-UMW'='STEPPS'; 'STEPPS-NEUS'='STEPPS'")
summary(mod.region)

fcomp.sens.glob <- lm(log(var.ecosys) ~ log(var.pdsi)*source-1-log(var.pdsi), data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Model",])
fcomp.sens.reg <- lm(log(var.ecosys) ~ log(var.pdsi)*source*region-1-region*log(var.pdsi) , data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Model",])
AIC(fcomp.sens.glob, fcomp.sens.reg)
summary(fcomp.sens.reg)
anova(fcomp.sens.reg)

stepps.reg <- lm(log(var.ecosys) ~ log(var.pdsi)*region, data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Empirical" & mod.region$source=="STEPPS",])
stepps.reg2 <- lm(log(var.ecosys) ~ log(var.pdsi) + region, data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Empirical" & mod.region$source=="STEPPS",])
ed.reg <- lm(log(var.ecosys) ~ log(var.pdsi)*region, data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Model" & mod.region$source=="ED2",])
link.reg <- lm(log(var.ecosys) ~ log(var.pdsi)*region, data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Model" & mod.region$source=="LINKAGES",])
lpjg.reg <- lm(log(var.ecosys) ~ log(var.pdsi)*region, data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Model" & mod.region$source=="LPJ-GUESS",])
lpjw.reg <- lm(log(var.ecosys) ~ log(var.pdsi)*region, data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Model" & mod.region$source=="LPJ-WSL",])
triff.reg <- lm(log(var.ecosys) ~ log(var.pdsi)*region, data=mod.region[mod.region$var=="fcomp" & mod.region$type=="Model" & mod.region$source=="TRIFFID",])
summary(stepps.reg)
summary(stepps.reg2)
summary(ed.reg)
summary(link.reg)
summary(lpjg.reg)
summary(lpjw.reg)
summary(triff.reg)
AIC(stepps.reg, stepps.reg2)

ggplot(data=mod.region[mod.region$var=="fcomp",]) +
  facet_wrap(~source, scales="free") +
  geom_histogram(aes(x=log(var.ecosys), fill=region)) + 
  theme_bw()

ggplot(data=mod.region[mod.region$var=="fcomp",]) +
  facet_wrap(~source, scales="free") +
  stat_smooth(aes(x=log(var.pdsi), y=log(var.ecosys)), method="lm", color="black", fill="black") +
  geom_point(aes(x=log(var.pdsi), y=log(var.ecosys), color=region), size=0.8, alpha=0.8) +
  stat_smooth(aes(x=log(var.pdsi), y=log(var.ecosys), color=region, fill=region), method="lm") +
  theme_bw()

# ------------


# ------------
# Fcomp-BM Comparisons
# ------------
summary(mod.dat2)

md.comp.v.bm1 <- lm(log(variability2) ~ log(variability1)*source, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass",])
summary(md.comp.v.bm1)

md.comp.v.bm2 <- lm(log(variability2) ~ log(variability1)*source-log(variability1)-1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass",])
summary(md.comp.v.bm2)

cvb.ed2 <- lm(variability2 ~ variability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="ED2",])
cvb.lpjg <- lm(variability2 ~ variability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="LPJ-GUESS",])
cvb.lpjw <- lm(variability2 ~ variability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="LPJ-WSL",])
cvb.link <- lm(variability2 ~ variability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="LINKAGES",])
cvb.triff <- lm(variability2 ~ variability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="TRIFFID",])
summary(cvb.ed2)
summary(cvb.lpjg)
summary(cvb.lpjw)
summary(cvb.link)
summary(cvb.triff)
# ------------
# -------------------------------------------

# -------------------------------------------
# Comparing stability sensitivity within models --> find where things fall apart
# -------------------------------------------
models.long <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_long.csv"))
summary(models.long)

models.long$Model <- factor(models.long$Model, levels=c("ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "TRIFFID"))


models.long$var <- car::recode(models.long$var, "'gpp'='GPP'; 'npp'='NPP'; 'lai'='LAI'; 'bm'='Biomass'; 'fcomp'='Composition'; 'nee'='NEE'")
models.long$var <- factor(models.long$var, levels=c("GPP", "NPP", "NEE", "LAI", "Biomass", "Composition"))


# -----
# Adding in emprical data
# -----
dat.emp <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"))
dat.emp$var <- as.factor(ifelse(dat.emp$dataset=="ReFAB", "Biomass", "Composition"))
dat.emp$type <- "Empirical"
dat.emp$dataset2 <- car::recode(dat.emp$dataset2, "'STEPPS'='STEPPS-UMW'")
dat.emp$Model <- dat.emp$dataset
summary(dat.emp)

summary(models.long)
models.long$type <- "model"
models.long[models.long$diff.abs==1e-20,c("variability")] <- NA
summary(models.long)

models.long$Model <- factor(models.long$Model, levels=c("ED2", "LPJ-GUESS", "LPJ-WSL", "LINKAGES", "TRIFFID"))


# A combined thing for graphing
var.comparison <- data.frame(lon=c(models.long$lon, dat.emp$lon),
                             lat=c(models.long$lat, dat.emp$lat),
                             type=c(models.long$type, dat.emp$type),
                             Model=c(paste(models.long$Model), rep("Pollen", nrow(dat.emp))),
                             var=c(paste(models.long$var), paste(dat.emp$var)),
                             var.ecosys=c(models.long$variability, dat.emp$var.ecosys),
                             var.pdsi=c(models.long$var.pdsi, dat.emp$var.pdsi)
)
var.comparison$Model <- factor(var.comparison$Model, levels=c("ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "TRIFFID", "Pollen"))
var.comparison$var <- factor(var.comparison$var, levels=c("GPP", "NPP", "NEE", "LAI", "Biomass", "Composition"))

summary(var.comparison)
levels(var.comparison$var)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Variability_Ecosystem_v_Climate_Models.png"), height=6, width=8, units="in", res=320)
ggplot(data=var.comparison, aes(x=log(var.pdsi), y=log(var.ecosys), color=var, fill=var)) +
  labs(x="Log Relative Drought Variability", y="Log Relative Ecosystem Variability") +
  facet_wrap(~Model) +
  geom_point(data=var.comparison[!(var.comparison$type=="Empirical" & var.comparison$var=="Composition"),], size=0.05, alpha=0.2) +
  stat_smooth(data=var.comparison[!(var.comparison$type=="Empirical" & var.comparison$var=="Composition"),], method=lm, alpha=0.5) + 
  geom_point(data=var.comparison[var.comparison$type=="Empirical" & var.comparison$var=="Composition" & var.comparison$lon< -83,], size=0.05, alpha=0.2) +
  stat_smooth(data=var.comparison[var.comparison$type=="Empirical" & var.comparison$var=="Composition" & var.comparison$lon< -83,], method=lm, alpha=0.5) + 
  geom_point(data=var.comparison[var.comparison$type=="Empirical" & var.comparison$var=="Composition" & var.comparison$lon> -83,], size=0.05, alpha=0.2) +
  stat_smooth(data=var.comparison[var.comparison$type=="Empirical" & var.comparison$var=="Composition" & var.comparison$lon> -83,], method=lm, alpha=0.5) + 
  scale_color_manual(values=c("blue4", "darkseagreen4", "turquoise4", "darkgoldenrod2", "darkorange2", "deeppink3")) +
  scale_fill_manual(values=c("blue4", "darkseagreen4", "turquoise4", "darkgoldenrod2", "darkorange2", "deeppink3")) +
  guides(fill=guide_legend(nrow=1), color=guide_legend(nrow=1)) +
  # scale_fill_brewer(palette="Dark2") +
  # scale_color_brewer(palette="Dark2") +
  coord_cartesian(ylim=c(-10.5,0), expand=0) +
  theme_bw() +
  theme(panel.background=element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size=rel(1.25)),
        axis.title = element_text(size=rel(1.25), face="bold"),
        axis.ticks = element_blank(),
        legend.position="top",
        legend.title=element_blank(),
        legend.text = element_text(size=rel(1.25)),
        strip.text = element_text(size=rel(1.5), face="bold"))
dev.off()

# ------------
# ED2
# ------------
# Comparing mean var.ecosys relative to fcomp
summary(models.long)
var.ed2 <- lm(log(var.rel) ~ relevel(var, ref="Composition"), data=models.long[models.long$Model=="ED2", ])
summary(var.ed2)

mod.ed2.comp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Composition"), data=models.long[models.long$Model=="ED2", ])
mod.ed2.bm <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Biomass"), data=models.long[models.long$Model=="ED2", ])
mod.ed2.npp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="NPP"), data=models.long[models.long$Model=="ED2", ])
mod.ed2b <- lm(log(var.rel) ~ log(var.pdsi)*var-log(var.pdsi), data=models.long[models.long$Model=="ED2", ]) # Means parameterization
summary(mod.ed2.comp)
summary(mod.ed2.bm)
summary(mod.ed2.npp)
summary(mod.ed2b)
# trend.ed2 <- emmeans::emtrends(mod.ed2, "var", var="var.pdsi")
# trend.ed2
# ------------

# ------------
# LPJ-GUESS
# ------------
# Comparing mean var.ecosys relative to fcomp
var.lpjg <- lm(log(var.rel) ~ relevel(var, ref="fcomp"), data=models.long[models.long$Model=="LPJ-GUESS", ])
summary(var.lpjg)
var.lpjg2 <- lm(log(var.rel) ~ relevel(var, ref="gpp"), data=models.long[models.long$Model=="LPJ-GUESS", ])
summary(var.lpjg2)

mod.lpjg.comp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Composition"), data=models.long[models.long$Model=="LPJ-GUESS", ])
mod.lpjg.bm <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Biomass"), data=models.long[models.long$Model=="LPJ-GUESS", ])
mod.lpjg.npp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="NPP"), data=models.long[models.long$Model=="LPJ-GUESS", ])
mod.lpjgb <- lm(log(var.rel) ~ log(var.pdsi)*var-log(var.pdsi), data=models.long[models.long$Model=="LPJ-GUESS", ]) # Means parameterization
summary(mod.lpjg.comp)
summary(mod.lpjg.bm)
summary(mod.lpjg.npp)
summary(mod.lpjgb)
# trend.lpjg <- emmeans::emtrends(mod.lpjg, "var", var="var.pdsi")
# trend.lpjg
# pairs(trend.lpjg)
# ------------

# ------------
# ------------
# Comparing mean log(var.rel) relative to fcomp
var.lpjw <- lm(log(var.rel) ~ relevel(var, ref="fcomp"), data=models.long[models.long$Model=="LPJ-WSL", ])
summary(var.lpjw)

mod.lpjw.comp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Composition"), data=models.long[models.long$Model=="LPJ-WSL", ])
mod.lpjw.bm <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Biomass"), data=models.long[models.long$Model=="LPJ-WSL", ])
mod.lpjw.npp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="NPP"), data=models.long[models.long$Model=="LPJ-WSL", ])
mod.lpjwb <- lm(log(var.rel) ~ log(var.pdsi)*var-log(var.pdsi), data=models.long[models.long$Model=="LPJ-WSL", ]) # Means parameterization
summary(mod.lpjw.comp)
summary(mod.lpjw.bm)
summary(mod.lpjw.npp)
summary(mod.lpjwb)
# trend.lpjw <- emmeans::emtrends(mod.lpjw, "var", var="var.pdsi")
# trend.lpjw
# pairs(trend.lpjw)
# ------------

# ------------
# ------------
mod.link.comp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Composition"), data=models.long[models.long$Model=="LINKAGES", ])
mod.link.bm <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Biomass"), data=models.long[models.long$Model=="LINKAGES", ])
mod.link.npp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="NPP"), data=models.long[models.long$Model=="LINKAGES", ])
mod.linkb <- lm(log(var.rel) ~ log(var.pdsi)*var-log(var.pdsi), data=models.long[models.long$Model=="LINKAGES", ]) # Means parameterization
summary(mod.link.comp)
summary(mod.link.bm)
summary(mod.link.npp)
summary(mod.linkb)
# trend.link <- emmeans::emtrends(mod.link, "var", var="var.pdsi")
# trend.link
# pairs(trend.link)
# ------------

# ------------
# ------------
mod.triff.comp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Composition"), data=models.long[models.long$Model=="TRIFFID", ])
mod.triff.bm <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="Biomass"), data=models.long[models.long$Model=="TRIFFID", ])
mod.triff.npp <- lm(log(var.rel) ~ log(var.pdsi)*relevel(var, ref="NPP"), data=models.long[models.long$Model=="TRIFFID", ])
mod.triffb <- lm(log(var.rel) ~ log(var.pdsi)*var-log(var.pdsi), data=models.long[models.long$Model=="TRIFFID", ]) # Means parameterization
summary(mod.triff.comp)
summary(mod.triff.bm)
summary(mod.triff.npp)

summary(mod.triffb)
# trend.triff <- emmeans::emtrends(mod.triff, "var", var="var.pdsi")
# trend.triff
# pairs(trend.triff)
# ------------

# -------------------------------------------
