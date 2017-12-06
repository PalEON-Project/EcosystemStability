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
# 1. Read in data; 
# 2. merge climate & fcomp stability with other ecosys variables
# 3. relativize ecosystem variables so they can be analyzed on common scales
# 4. determine which climate predictor has the tightest correlation with GPP
# 5. compare climate sensitivity of ecosystem variables WITHIN MODELS
# 6. Determine where in the ecosystem the biggest jumps in changes in cliamte sensitivity occurr
#     - this can then be explained in the discussion by leveraging what we know about the models
# 7. compare cliamte sensitivity of ecosystem variables generally across models: are fast or slow climate variables MORE or LESS sensitive to climate
#     - relativized slopes
#     - fraction points showing significant change
# 8. Are different variables sensitive to different aspects of climate?
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
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data/"

# path.gamm.func <- "~/Desktop/R_Functions/"  # Path to github repository of my GAMM helper functions: https://github.com/crollinson/R_Functions.git
path.gamm.func <- "~/Desktop/Research/R_Functions/"
# mip.utils <- "~/Dropbox/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git
mip.utils <- "~/Desktop/Research/PalEON_CR/MIP_Utils/" # Path to PalEON MIP Utility repository: https://github.com/PalEON-Project/MIP_Utils.git

path.co2 <- "~/Dropbox/PalEON_CR/env_regional/env_paleon/co2/paleon_annual_co2.nc"
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
stab.met   <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Drivers_100.csv"))
stab.models <- read.csv(file.path(path.google, "Current Data/Stability_GAMs", "Stability_Models_100.csv"))

summary(stab.met)
summary(stab.models)


# library(GGally)
vars.graph <- c("gpp", "npp", "nee", "lai", "bm")
type.graph <- "diff."

# npp.range <- quantile(stab.models[stab.models$Model!="LINKAGES","diff.npp"], c(0.001, 0.75), na.rm=T)
ggplot(data=stab.models[stab.models$Model!="LINKAGES",]) +
  geom_point(aes(x=diff.gpp, y=diff.npp, color=Model)) +
  stat_smooth(aes(x=diff.gpp, y=diff.npp, color=Model, fill=Model), method="lm") 

ggplot(data=stab.models[stab.models$Model!="LINKAGES",]) +
  geom_point(aes(x=diff.gpp, y=diff.nee, color=Model)) +
  stat_smooth(aes(x=diff.gpp, y=diff.nee, color=Model, fill=Model), method="lm") 


# Duping model data into something very long 
vars.use <- c("gpp", "npp", "nee", "lai", "bm")
type.use <- "diff."
dat.long <- data.frame()
pb <- txtProgressBar(min=0, max=length(vars.use)^2-1)
pb.ind=1
for(i in 1:length(vars.use)){
  while(i < length(vars.use)){
    for(j in (i+1):length(vars.use)){
      print(paste0(vars.use[i], " - ", vars.use[j]))
      
      dat.temp <- data.frame(stab.models[, c("lon", "lat", "Model")], 
                             var1=vars.use[i], value1=stab.models[,paste0(type.use, vars.use[i])],
                             var2=vars.use[j], value2=stab.models[,paste0(type.use, vars.use[j])]
                             )
      
      dat.long <- rbind(dat.long, dat.temp)
    }
  }
}

summary(dat.long)

png(file.path(path.google, "Current Figures/Stability_Synthesis", "Stability_Models_CorrelationMatrix.png"), height=9, width=10, units="in", res=220)
ggplot(data=dat.long[dat.long$Model!="LINKAGES",]) +
  facet_grid(var2 ~ var1, scales="free") +
  geom_point(aes(x=value1, y=value2, color=Model), size=0.2, alpha=0.8) +
  stat_smooth(aes(x=value1, y=value2, color=Model, fill=Model), method="lm") +
  coord_cartesian(expand=c(0,0)) +
  theme_bw()
dev.off()
# -------------------------------------------


# -------------------------------------------
# 2. merge climate & fcomp stability with other ecosys variables
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 3. relativize ecosystem variables so they can be analyzed on common scales
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 4. determine which climate predictor has the tightest correlation with GPP
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 5. compare climate sensitivity of ecosystem variables WITHIN MODELS
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 6. Determine where in the ecosystem the biggest jumps in changes in cliamte sensitivity occurr
#     - this can then be explained in the discussion by leveraging what we know about the models
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 7. compare cliamte sensitivity of ecosystem variables generally across models: are fast or slow climate variables MORE or LESS sensitive to climate
#     - relativized slopes
#     - fraction points showing significant change
# -------------------------------------------
# -------------------------------------------

# -------------------------------------------
# 8. Are different variables sensitive to different aspects of climate?
#     - e.g. GPP tracks temperature, but Fcomp is more correlated with precip
# -------------------------------------------
# -------------------------------------------




