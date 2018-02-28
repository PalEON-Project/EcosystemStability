# -------------------------------------------
# Assessing the changes in climate sensitivity during a process cascade in differnet models
# Author: Christy Rollinson, crollinson@gmail.com
#
# 1. Figure out which climate variable correlates strongest with GPP (basic photosynthesis)
# 2. Look at how (relative) sensitivity to that climate variable as we move from upstream to downstream processes
#    GPP --> NEE --> NPP --> LAI --> BM --> Fcomp
#
# -------------------------------------------
# Workflow
# -------------------------------------------
# 1. Read in & merge data; 
# 2. relativize ecosystem variables so they can be analyzed on common scales
# 3. determine which climate predictor has the tightest correlation with GPP
# 4. compare climate sensitivity of ecosystem variables WITHIN MODELS
# 5. Determine where in the ecosystem the biggest jumps in changes in cliamte sensitivity occurr
#     - this can then be explained in the discussion by leveraging what we know about the models
# 6. compare cliamte sensitivity of ecosystem variables generally across models: are fast or slow climate variables MORE or LESS sensitive to climate
#     - relativized slopes
#     - fraction points showing significant change
# 7. Are different variables sensitive to different aspects of climate?
#     - e.g. GPP tracks temperature, but Fcomp is more correlated with precip
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
path.repo <- "~/Desktop/Research/PalEON_EcosystemStability/"
setwd(path.repo)

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


# ------------
# An exploratory graph comparing models
# ------------
png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_VarComparison_PDSI.png"), height=6, width=8, units="in", res=220)
ggplot(data=models.long, aes(x=log(pdsi.diff), y=log(diff.abs), color=Model, fill=Model)) +
  facet_wrap(~var, scales="free_y") +
  geom_point(size=0.25) +
  stat_smooth(method=lm, alpha=0.5) + 
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_VarComparison_Tair.png"), height=6, width=8, units="in", res=220)
ggplot(data=models.long, aes(x=log(tair.diff), y=log(diff.abs), color=Model, fill=Model)) +
  facet_wrap(~var, scales="free_y") +
  geom_point(size=0.25) +
  stat_smooth(method=lm, alpha=0.5) + 
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  theme_bw()
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_VarComparison_Precip.png"), height=6, width=8, units="in", res=220)
ggplot(data=models.long, aes(x=log(precip.diff), y=log(diff.abs), color=Model, fill=Model)) +
  facet_wrap(~var, scales="free_y") +
  geom_point(size=0.25) +
  stat_smooth(method=lm, alpha=0.5) + 
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  theme_bw()
dev.off()

# ------------

# -------------------------------------------


# -------------------------------------------
# 3. relativize ecosystem variables so they can be analyzed on common scales
# Notes
# - Mean centering would help with comparison among models for a given var, but if 
#   we want to compare across vars within a model, we need additional scaling
# - Concerns:
#    - Do we scale before or after log-transformation?  I'm leaning toward before, and scaling everything to 0-1 
#      to reflect Fcomp's natural scale.
# Two Options:
# 1. Mean-relativizing: obs/mean
# 2. Range-relativizing 1-(max-obs)/(max-min)
# -------------------------------------------
summary(models.long)

# relativize pdsi.diff to test
# models.long$pdsi.rel.mean  <- models.long$pdsi.diff/mean(models.long$pdsi.diff)
# models.long$pdsi.rel.range <- 1 - (max(models.long$pdsi.diff) - models.long$pdsi.diff) / 
                                  (max(models.long$pdsi.diff) - min(models.long$pdsi.diff))
summary(models.long)

for(v in unique(models.long$var)){
  for(mod in unique(models.long$Model)){
    dat.ind <- which(models.long$var==v & models.long$Model==mod & !is.na(models.long$diff.abs))
    dat.tmp <- models.long[dat.ind, ]
    
    models.long[dat.ind, "rel.mean" ] <- dat.tmp$diff.abs/mean(dat.tmp$diff.abs)
    models.long[dat.ind, "rel.range"] <- 1 - (max(dat.tmp$diff.abs) - dat.tmp$diff.abs) / 
                                             (max(dat.tmp$diff.abs) - min(dat.tmp$diff.abs))
    
    models.long[dat.ind, "pdsi.rel.mean" ] <- dat.tmp$pdsi.diff/mean(dat.tmp$pdsi.diff)
    models.long[dat.ind, "pdsi.rel.range"] <- 1 - (max(dat.tmp$pdsi.diff) - dat.tmp$pdsi.diff) / 
                                                  (max(dat.tmp$pdsi.diff) - min(dat.tmp$pdsi.diff))
    
  }
}
summary(models.long)


pdf(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_VarComposition_PDSI_RelMethods.pdf"))
print(
  ggplot(data=models.long, aes(x=log(pdsi.diff), y=log(diff.abs), color=Model, fill=Model)) +
    facet_wrap(~var, scales="free_y") +
    geom_point(size=0.25) +
    stat_smooth(method=lm, alpha=0.5) + 
    scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
    scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
    theme_bw()
)
print(
  ggplot(data=models.long, aes(x=log(pdsi.rel.mean), y=log(rel.mean), color=Model, fill=Model)) +
    facet_wrap(~var) +
    geom_point(size=0.25) +
    stat_smooth(method=lm, alpha=0.5) + 
    scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
    scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
    coord_cartesian(ylim=c(-3, 2.5)) +
    theme_bw()
)
print(
  ggplot(data=models.long, aes(x=log(pdsi.rel.range), y=log(rel.range), color=Model, fill=Model)) +
    facet_wrap(~var) +
    geom_point(size=0.25) +
    stat_smooth(method=lm, alpha=0.5) + 
    scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
    scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
    coord_cartesian(ylim=c(-5, 0)) +
    theme_bw()
)
dev.off()


png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_varComparison_PDSI_RelMean_Models.png"), height=6, width=8, units = "in", res=220)
ggplot(data=models.long, aes(x=log(pdsi.rel.mean), y=log(rel.mean), color=Model, fill=Model)) +
  labs(x="PDSI Stability", y="Ecosystem Stability") +
  facet_wrap(~var) +
  geom_point(size=0.1, alpha=0.5) +
  stat_smooth(method=lm, alpha=0.5) + 
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(models.long$Model),"color"])) +
  coord_cartesian(ylim=c(-2, 1.5)) +
  theme_bw() #+
  # theme(legend.position=c(0.75, 0.25))
dev.off()

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Model_VarComparison_PDSI_RelMean.png"), height=6, width=8, units = "in", res=220)
ggplot(data=models.long, aes(x=log(pdsi.rel.mean), y=log(rel.mean), color=var, fill=var)) +
  labs(x="PDSI Stability", y="Ecosystem Stability") +
  facet_wrap(~Model) +
  geom_point(size=0.1, alpha=0.5) +
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

# -------------------------------------------
# 4. Compare variability of PDSI sensitivity ACROSS models
#    - Is there a variable with more or less consistency?
#      Hypothesis: GPP should be more similar b/c of standard photosynthesis equations
# -------------------------------------------
summary(models.long)

mod.gpp <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*Model-log(pdsi.rel.mean), data=models.long[models.long$var=="gpp" , ])
gpp.trend <- emmeans::emtrends(mod.gpp, "Model", var="log(pdsi.rel.mean)")
summary(mod.gpp)
pairs(gpp.trend)

mod.nee <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*Model-log(pdsi.rel.mean), data=models.long[models.long$var=="nee" , ])
nee.trend <- emmeans::emtrends(mod.nee, "Model", var="log(pdsi.rel.mean)")
summary(mod.nee)
pairs(nee.trend)

mod.npp <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*Model-log(pdsi.rel.mean), data=models.long[models.long$var=="npp" , ])
npp.trend <- emmeans::emtrends(mod.npp, "Model", var="log(pdsi.rel.mean)")
summary(mod.npp)
pairs(npp.trend)

mod.lai <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*Model-log(pdsi.rel.mean), data=models.long[models.long$var=="lai" , ])
lai.trend <- emmeans::emtrends(mod.lai, "Model", var="log(pdsi.rel.mean)")
summary(mod.lai)
pairs(lai.trend)

mod.bm <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*Model-log(pdsi.rel.mean), data=models.long[models.long$var=="bm" , ])
bm.trend <- emmeans::emtrends(mod.bm, "Model", var="log(pdsi.rel.mean)")
summary(mod.bm)
pairs(bm.trend)

mod.fcomp <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*Model-log(pdsi.rel.mean), data=models.long[models.long$var=="fcomp" , ])
fcomp.trend <- emmeans::emtrends(mod.fcomp, "Model", var="log(pdsi.rel.mean)")
summary(mod.fcomp)
pairs(fcomp.trend)
# -------------------------------------------

# -------------------------------------------
# 5. compare climate sensitivity of ecosystem variables WITHIN MODELS
# -------------------------------------------
mod.ed2 <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*var, data=models.long[models.long$Model=="ED2", ])
trend.ed2 <- emmeans::emtrends(mod.ed2, "var", var="log(pdsi.rel.mean)")
summary(mod.ed2)
trend.ed2

mod.link <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*var, data=models.long[models.long$Model=="LINKAGES", ])
trend.link <- emmeans::emtrends(mod.link, "var", var="log(pdsi.rel.mean)")
summary(mod.link)
trend.link
pairs(trend.link)

mod.lpjg <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*var, data=models.long[models.long$Model=="LPJ-GUESS", ])
trend.lpjg <- emmeans::emtrends(mod.lpjg, "var", var="log(pdsi.rel.mean)")
summary(mod.lpjg)
trend.lpjg

mod.lpjw <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*var, data=models.long[models.long$Model=="LPJ-WSL", ])
trend.lpjw <- emmeans::emtrends(mod.lpjw, "var", var="log(pdsi.rel.mean)")
summary(mod.lpjw)
trend.lpjw

mod.triff <- lm(log(rel.mean) ~ log(pdsi.rel.mean)*var, data=models.long[models.long$Model=="TRIFFID", ])
trend.triff <- emmeans::emtrends(mod.triff, "var", var="log(pdsi.rel.mean)")
summary(mod.triff)
trend.triff
# -------------------------------------------

# -------------------------------------------
# 6. Determine where in the ecosystem the biggest jumps in changes in cliamte sensitivity occurr
#     - this can then be explained in the discussion by leveraging what we know about the models
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 6. compare cliamte sensitivity of ecosystem variables generally across models: are fast or slow climate variables MORE or LESS sensitive to climate
#     - relativized slopes
#     - fraction points showing significant change
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 7. Are different variables sensitive to different aspects of climate?
#     - e.g. GPP tracks temperature, but Fcomp is more correlated with precip
# -------------------------------------------
# -------------------------------------------




