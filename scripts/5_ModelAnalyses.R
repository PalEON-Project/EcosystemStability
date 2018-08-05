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

dat.colors <- data.frame(model=c("LBDA", "STEPPS", "ReFAB", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"),
                         type=c(rep("emipirical", 3), rep("model", 7)),
                         color=c(rep("#000000", 3), rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
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
summary(stab.met)

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
    
    models.long[dat.ind, "stability" ] <- -log(dat.tmp$diff.abs/abs(mean(dat.tmp$val.mean)))
    models.long[dat.ind, "stab.rel"] <- (models.long[dat.ind, "stability" ]) / mean(models.long[dat.ind, "stability" ]) 

    models.long[dat.ind, "stab.pdsi" ] <- -log(dat.tmp$pdsi.diff/abs(mean(dat.tmp$pdsi.diff)))
    
    models.long[dat.ind, "variability" ] <- dat.tmp$diff.abs/abs(mean(dat.tmp$val.mean))
    models.long[dat.ind, "var.pdsi" ] <- dat.tmp$pdsi.diff/abs(mean(dat.tmp$pdsi.diff))
    models.long[dat.ind, "var.rel"] <- (models.long[dat.ind, "variability" ]) / mean(models.long[dat.ind, "variability" ]) 
    
  }
}
summary(models.long)
models.long$type <- "Model"
# ------------

# ------------
# Read in the empirical data from script #4
# ------------
dat.emp <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"))
dat.emp$var <- as.factor(ifelse(dat.emp$dataset=="ReFAB", "bm", "fcomp"))
dat.emp$type <- "Empirical"
summary(dat.emp)
# ------------


# ------------
# Combining emprirical and model datasets together
# ------------
mod.v.dat <- data.frame(lat=c(models.long$lat[models.long$var %in% c("bm", "fcomp")], dat.emp$lat),
                        lon=c(models.long$lon[models.long$var %in% c("bm", "fcomp")], dat.emp$lon),
                        type=c(models.long$type[models.long$var %in% c("bm", "fcomp")], dat.emp$type),
                        var=c(paste(models.long$var[models.long$var %in% c("bm", "fcomp")]), paste(dat.emp$var)),
                        source=c(paste(models.long$Model[models.long$var %in% c("bm", "fcomp")]), paste(dat.emp$dataset)),
                        stab.ecosys=c(models.long$stability[models.long$var %in% c("bm", "fcomp")], dat.emp$stab.ecosys),
                        stab.pdsi=c(models.long$stab.pdsi[models.long$var %in% c("bm", "fcomp")], dat.emp$stab.pdsi),
                        var.ecosys=c(models.long$variability[models.long$var %in% c("bm", "fcomp")], dat.emp$var.ecosys),
                        var.pdsi=c(models.long$var.pdsi[models.long$var %in% c("bm", "fcomp")], dat.emp$var.pdsi)
                        )

mod.v.dat$source <- factor(mod.v.dat$source, levels=c("ReFAB", "STEPPS", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))
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
                       stability1 = c(mod.v.dat$stab.pdsi, models.long$stability[models.long$var=="fcomp"], dat.refab$stab.stepps.1k),
                       stability2 = c(mod.v.dat$stab.ecosys, rep(NA, length(models.long$stability[models.long$var=="fcomp"])), dat.refab$stab.refab.1k),
                       variability1 = c(mod.v.dat$var.pdsi, models.long$variability[models.long$var=="fcomp"], dat.refab$var.stepps.1k),
                       variability2 = c(mod.v.dat$var.ecosys, rep(NA, length(models.long$variability[models.long$var=="fcomp"])), dat.refab$var.refab.1k)
                       )
summary(mod.dat2)

# Looping through to pair model BM & Fcomp stability
for(i in 1:nrow(mod.dat2)){
  if(mod.dat2$type[i]=="Empirical" | !is.na(mod.dat2$stability2[i])) next
  
  val.fill <- models.long[models.long$var=="bm" & models.long$Model==paste(mod.dat2$source[i]) & models.long$lat==mod.dat2$lat[i] & models.long$lon==mod.dat2$lon[i], "stability"]
  val.fill2 <- models.long[models.long$var=="bm" & models.long$Model==paste(mod.dat2$source[i]) & models.long$lat==mod.dat2$lat[i] & models.long$lon==mod.dat2$lon[i], "variability"]
  
  if(length(val.fill)==0) next 
  
  mod.dat2[i,"stability2"] <- val.fill
  mod.dat2[i,"variability2"] <- val.fill2
}
mod.dat2$var2 <- car::recode(mod.dat2$var2, "'bm'='Biomass'; 'fcomp'='Composition'")
mod.dat2$source <- factor(mod.dat2$source, levels=c("ReFAB", "STEPPS", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "TRIFFID"))
summary(mod.dat2)

# write.csv(dat.sites, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Diversity_Data.csv"), row.names=F)

write.csv(mod.dat2, file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_and_Data.csv"), row.names=F)
# ------------
# -------------------------------------------


# -------------------------------------------
# Comparing stabiltiy in models versus data
# -------------------------------------------
mod.dat2$var1 <- factor(mod.dat2$var1, levels=c("PDSI", "Composition"))
mod.dat2$var2 <- factor(mod.dat2$var2, levels=c("Composition", "Biomass"))
mod.dat2[!is.na(mod.dat2$stability1) & mod.dat2$stability1>20, "stability1"] <- NA
mod.dat2[!is.na(mod.dat2$stability2) & mod.dat2$stability2>20, "stability2"] <- NA

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Model_v_Data.png"), height=6, width=6, units="in", res=320)
ggplot(dat=mod.dat2) +
  facet_grid(var2 ~ var1, scales="free", switch="both") +
  geom_point(aes(x=stability1, y=stability2, color=source), size=0.25, alpha=0.5) +
  stat_smooth(aes(x=stability1, y=stability2, color=source, fill=source), method="lm") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=rel(1)),
        axis.title = element_blank()) +
  theme(legend.position=c(0.75, 0.75),
        legend.title=element_blank())
dev.off()  

summary(models.long)

# ------------
# Biomass-PDSI Comparisons
# ------------
md.comp.bm.mean <- lm(stab.ecosys ~ source, data=mod.v.dat[mod.v.dat$var=="bm",])
summary(md.comp.bm.mean)

md.comp.bm.sens1 <- lm(stab.ecosys ~ stab.pdsi*source, data=mod.v.dat[mod.v.dat$var=="bm",])
summary(md.comp.bm.sens1)

md.comp.bm.sens2 <- lm(stab.ecosys ~ stab.pdsi*source-stab.pdsi, data=mod.v.dat[mod.v.dat$var=="bm",])
summary(md.comp.bm.sens2)
# ------------

# ------------
# Fcomp-PDSI Comparisons
# ------------
md.comp.fcomp.mean <- lm(stab.ecosys ~ source, data=mod.v.dat[mod.v.dat$var=="fcomp",])
summary(md.comp.fcomp.mean)

md.comp.fcomp.sens1 <- lm(stab.ecosys ~ stab.pdsi*source, data=mod.v.dat[mod.v.dat$var=="fcomp",])
summary(md.comp.fcomp.sens1)

md.comp.fcomp.sens2 <- lm(stab.ecosys ~ stab.pdsi*source-stab.pdsi, data=mod.v.dat[mod.v.dat$var=="fcomp",])
summary(md.comp.fcomp.sens2)
# ------------

# ------------
# Fcomp-BM Comparisons
# ------------
summary(mod.dat2)

md.comp.v.bm1 <- lm(stability2 ~ stability1*source, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass",])
summary(md.comp.v.bm1)

md.comp.v.bm2 <- lm(stability2 ~ stability1*source-stability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass",])
summary(md.comp.v.bm2)

cvb.link <- lm(stability2 ~ stability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="LINKAGES",])
cvb.triff <- lm(stability2 ~ stability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="TRIFFID",])
cvb.lpjw <- lm(stability2 ~ stability1, data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$source=="LPJ-WSL",])
summary(cvb.link)
summary(cvb.triff)
summary(cvb.lpjw)
# ------------
# -------------------------------------------

# -------------------------------------------
# Comparing stability sensitivity within models --> find where things fall apart
# -------------------------------------------
models.long[models.long$diff.abs==1e-20,"stability"] <- NA
summary(models.long)

models.long$Model <- factor(models.long$Model, levels=c("ED2", "LPJ-GUESS", "LPJ-WSL", "LINKAGES", "TRIFFID"))

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Models.png"), height=6, width=8, units="in", res=320)
ggplot(data=models.long, aes(x=stab.pdsi, y=stability, color=var, fill=var)) +
  labs(x="PDSI Stability", y="Ecosystem Stability") +
  facet_wrap(~Model) +
  geom_point(size=0.1, alpha=0.3) +
  stat_smooth(method=lm, alpha=0.5) + 
  scale_color_manual(values=c("blue4", "darkseagreen4", "turquoise4", "darkgoldenrod2", "darkorange2", "deeppink3")) +
  scale_fill_manual(values=c("blue4", "darkseagreen4", "turquoise4", "darkgoldenrod2", "darkorange2", "deeppink3")) +
  # scale_fill_brewer(palette="Dark2") +
  # scale_color_brewer(palette="Dark2") +
  coord_cartesian(ylim=c(3,11.5), expand=0) +
  theme_bw() +
  theme(legend.position=c(0.75, 0.25),
        legend.title=element_blank())
dev.off()

# ------------
# ED2
# ------------
# Comparing mean stability relative to FCOMP
stab.ed2 <- lm(stability ~ relevel(var, ref="fcomp"), data=models.long[models.long$Model=="ED2", ])
summary(stab.ed2)

mod.ed2 <- lm(stability ~ stab.pdsi*var, data=models.long[models.long$Model=="ED2", ])
mod.ed2b <- lm(stability ~ stab.pdsi*var-stab.pdsi, data=models.long[models.long$Model=="ED2", ]) # Means parameterization
summary(mod.ed2)
summary(mod.ed2b)
# trend.ed2 <- emmeans::emtrends(mod.ed2, "var", var="stab.pdsi")
# trend.ed2
# ------------

# ------------
# LPJ-GUESS
# ------------
# Comparing mean stability relative to FCOMP
stab.lpjg <- lm(stability ~ relevel(var, ref="fcomp"), data=models.long[models.long$Model=="LPJ-GUESS", ])
summary(stab.lpjg)
stab.lpjg2 <- lm(stability ~ relevel(var, ref="gpp"), data=models.long[models.long$Model=="LPJ-GUESS", ])
summary(stab.lpjg2)

mod.lpjg <- lm(stability ~ stab.pdsi*var, data=models.long[models.long$Model=="LPJ-GUESS", ])
mod.lpjgb <- lm(stability ~ stab.pdsi*var-stab.pdsi, data=models.long[models.long$Model=="LPJ-GUESS", ]) # Means parameterization
summary(mod.lpjg)
summary(mod.lpjgb)
# trend.lpjg <- emmeans::emtrends(mod.lpjg, "var", var="stab.pdsi")
# trend.lpjg
# pairs(trend.lpjg)
# ------------

# ------------
# ------------
# Comparing mean stability relative to FCOMP
stab.lpjw <- lm(stability ~ relevel(var, ref="fcomp"), data=models.long[models.long$Model=="LPJ-WSL", ])
summary(stab.lpjw)

mod.lpjw <- lm(stability ~ stab.pdsi*var, data=models.long[models.long$Model=="LPJ-WSL", ])
mod.lpjwb <- lm(stability ~ stab.pdsi*var-stab.pdsi, data=models.long[models.long$Model=="LPJ-WSL", ]) # Means parameterization
summary(mod.lpjw)
summary(mod.lpjwb)
# trend.lpjw <- emmeans::emtrends(mod.lpjw, "var", var="stab.pdsi")
# trend.lpjw
# pairs(trend.lpjw)
# ------------

# ------------
# ------------
mod.link <- lm(stability ~ stab.pdsi*var, data=models.long[models.long$Model=="LINKAGES", ])
mod.linkb <- lm(stability ~ stab.pdsi*var-stab.pdsi, data=models.long[models.long$Model=="LINKAGES", ]) # Means parameterization
summary(mod.link)
summary(mod.linkb)
# trend.link <- emmeans::emtrends(mod.link, "var", var="stab.pdsi")
# trend.link
# pairs(trend.link)
# ------------

# ------------
# ------------
mod.triff <- lm(stability ~ stab.pdsi*var, data=models.long[models.long$Model=="TRIFFID", ])
mod.triffb <- lm(stability ~ stab.pdsi*var-stab.pdsi, data=models.long[models.long$Model=="TRIFFID", ]) # Means parameterization
summary(mod.triff)
summary(mod.triffb)
# trend.triff <- emmeans::emtrends(mod.triff, "var", var="stab.pdsi")
# trend.triff
# pairs(trend.triff)
# ------------

# -------------------------------------------
