# -------------------------------------------
# Manuscript Figures in one place to make reproduction and manipulation easy
#
# Manuscript Tables: 
# 1. [added on to trait table] variability & PDSI sensitivity of composition 
#    & biomass for models and data
#
# Manuscript Figures:
# 1. Conceptual Figure (Google Slides)
# 2. Map of Empiricial Dataset Variability
# 3. Biomass-Composition-PDSI Correllations (Models + Data)
# 4. Model Latent States
#
# Supplemental Tables:
# ST1. Sensitivity of ecosystem variability to hydroclimate (PDSI) variabiltiy
# 
# Supplemental Figures: 
# SF1. Empricial vs. Driver PDSI Map: 
#      a. LBDA 1st yr; 
#      b. LBDA Drought; 
#      c. Driver Common Period; 
#      d. Driver Full period
# SF2. Emprirical Biomass/Composition Variabilty vs. Richness
#

# -------------------------------------------
rm(list=ls())

# -------------------------------------------
# Load Libaries & Set file paths
# -------------------------------------------
library(ggplot2); library(cowplot)

# Path to where the raw output is
path.data <- "~/Dropbox/PalEON_CR/PalEON_MIP2_Region/PalEON_Regional_Extract/"

# Path to where data are; lets just pull straight from the Google Drive folder
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data"

path.figs <- file.path(path.google, "Manuscript/EcologyLetters_2020-03-30/Figures")
path.tables <- file.path(path.google, "Manuscript/EcologyLetters_2020-03-30/Tables")
if(!dir.exists(path.figs)) dir.create(path.figs)
if(!dir.exists(path.tables)) dir.create(path.tables)

# Set up a table with colors for models & data
# Using a CB palette from here: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
dat.colors <- data.frame(model=c("LBDA", "STEPPS", "ReFAB", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "JULES-TRIFFID"),
                         type=c(rep("emipirical", 3), rep("model", 7)),
                         color=c(rep("#000000", 3), rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
dat.colors

# Storing the paths for US state outlines
us <- ggplot2::map_data("state")

# -------------------------------------------


# -------------------------------------------
# Figure 2. Map of Empiricial Dataset Variability
# -------------------------------------------
climate.comparison.sp <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data_Spatial.csv"))

climate.comparison.sp$dataset <- car::recode(climate.comparison.sp$dataset, "'LBDA'='Drought'; 'STEPPS'='Composition'; 'ReFAB'='Biomass'")
climate.comparison.sp$dataset <- factor(climate.comparison.sp$dataset, levels=c("Drought", "Composition", "Biomass"))
summary(climate.comparison.sp)

plot.lbda <-  ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability) & climate.comparison.sp$dataset=="Drought",]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  # geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_tile(data=climate.comparison.sp[climate.comparison.sp$dataset=="Drought" & !is.na(climate.comparison.sp$variability),], aes(x=lon, y=lat, fill=log(variability))) +
  geom_path(data=us, aes(x=long, y=lat, group=group), color="gray15") +
  geom_text(x=-97.75, y=48, label="a)", size=5, fontface="bold") +
  coord_equal(xlim=range(climate.comparison.sp$lon[climate.comparison.sp$dataset=="Composition"], na.rm=T), ylim=range(climate.comparison.sp$lat[climate.comparison.sp$dataset=="Composition"], na.rm=T)) +
  scale_fill_gradientn(name="PDSI\nVariability", colors=c("#018571","#80cdc1", "#f5f5f5", "#dfc27d", "#a6611a"), limits=range(log(climate.comparison.sp$variability[climate.comparison.sp$dataset=="Drought"]), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

plot.stepps <- ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability) & climate.comparison.sp$dataset=="Composition",]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group), color="gray15") +
  geom_text(x=-97.75, y=48, label="b)", size=5, fontface="bold") +
  coord_equal(xlim=range(climate.comparison.sp$lon[climate.comparison.sp$dataset=="Composition"], na.rm=T), ylim=range(climate.comparison.sp$lat[climate.comparison.sp$dataset=="Composition"], na.rm=T)) +
  scale_color_gradientn(name="Comp.\nVariability", colors=c("#008837","#a6dba0", "#f7f7f7", "#c2a5cf", "#7b3294"), limits=range(log(climate.comparison.sp$variability[climate.comparison.sp$dataset=="Composition"]), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

plot.refab <- ggplot(data=climate.comparison.sp[!is.na(climate.comparison.sp$variability) & climate.comparison.sp$dataset=="Biomass",]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_path(data=us, aes(x=long, y=lat, group=group), color="gray15") +
  geom_text(x=-97.75, y=48, label="c)", size=5, fontface="bold") +
  coord_equal(xlim=range(climate.comparison.sp$lon[climate.comparison.sp$dataset=="Composition"], na.rm=T), ylim=range(climate.comparison.sp$lat[climate.comparison.sp$dataset=="Composition"], na.rm=T)) +
  scale_color_gradientn(name="Biomass\nVariability", colors=c("#4dac26","#b8e186", "#f7f7f7", "#f1b6da", "#d01c8b"), limits=range(log(climate.comparison.sp$variability[climate.comparison.sp$dataset=="Biomass"]), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

png(file.path(path.figs, "Figure2_Variability_Ecosystem_v_Climate_Data_Map_FreeColor.png"), height=6, width=6, units="in", res=220)
plot_grid(plot.lbda, plot.stepps, plot.refab, ncol=1)
dev.off()
# -------------------------------------------

# -------------------------------------------
# 3. Biomass-Composition-PDSI Correllations (Models + Data)
# -------------------------------------------
mod.dat2 <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_and_Data.csv"))
mod.dat2$source <- car::recode(mod.dat2$source, "'TRIFFID'='JULES-TRIFFID'")
mod.dat2$source <- factor(mod.dat2$source, levels=c("ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "JULES-TRIFFID"))
# mod.dat2$var1 <- car::recode(mod.dat2$var1, "'PDSI'='Drought'")
mod.dat2$var1 <- factor(mod.dat2$var1, levels=c("PDSI", "Composition"))
mod.dat2$var2 <- factor(mod.dat2$var2, levels=c("Composition", "Biomass"))

dat.colors$model <- factor(dat.colors$model, levels=c("drivers", "drivers-modified", "LBDA", "ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "JULES-TRIFFID"))
dat.colors <- dat.colors[order(dat.colors$model),]

panel.labs <- data.frame(var1=c("PDSI", "PDSI", "Composition"), 
                         var2=c("Composition", "Biomass", "Biomass"), 
                         x=c(-6.5, -6.5, -12.0),
                         y=c(-12.0, -12.5, -12.5),
                         label=c("a)", "b)", "c)"))

dat.colors <- data.frame(model=c("LBDA", "ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "drivers", "drivers-modified", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "JULES-TRIFFID"),
                         type=c(rep("emipirical", 4), rep("model", 7)),
                         color=c(rep("#000000", 3), "#7F7F7F", rep("#999999", 2),  "#009E73", "#0072B2", "#D55E00", "#D69F00", "#CC79A7"))
dat.colors


png(file.path(path.figs, "Figure3_Variability_Model_v_Data.png"), height=6, width=6, units="in", res=320)
ggplot(dat=mod.dat2) +
  facet_grid(var2 ~ var1, scales="free", switch="both") +
  geom_point(aes(x=log(variability1), y=log(variability2), color=source), size=0.1, alpha=0.25) +
  stat_smooth(aes(x=log(variability1), y=log(variability2), color=source, fill=source), method="lm") +
  geom_text(data=panel.labs, aes(x=x, y=y, label=label), fontface="bold") +
  scale_fill_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  scale_color_manual(values=paste(dat.colors[dat.colors$model %in% unique(mod.dat2$source),"color"])) +
  scale_x_continuous(name="Log Normalized Variability") +
  scale_y_continuous(name="Log Normalized Variability") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=rel(1)),
        # axis.title = element_blank()
        axis.title = element_text(size=rel(1.25))
  ) +
  theme(legend.position=c(0.75, 0.75),
        legend.title=element_blank(),
        panel.grid=element_blank())
dev.off()  
# -------------------------------------------

# -------------------------------------------
# 4. Model Latent States
# -------------------------------------------
models.long <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_long.csv"))
models.long$Model <- car::recode(models.long$Model, "'TRIFFID'='JULES-TRIFFID'")
summary(models.long)

models.long$Model <- factor(models.long$Model, levels=c("ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "JULES-TRIFFID"))
models.long$var <- car::recode(models.long$var, "'gpp'='GPP'; 'npp'='NPP'; 'lai'='LAI'; 'bm'='Biomass'; 'fcomp'='Composition'; 'nee'='NEE'")
models.long$var <- factor(models.long$var, levels=c("GPP", "NPP", "NEE", "LAI", "Biomass", "Composition"))

png(file.path(path.figs, "Figure4_Variability_Ecosystem_v_Climate_Models.png"), height=6, width=8, units="in", res=220)
ggplot(data=models.long, aes(x=log(var.pdsi), y=log(variability), color=var, fill=var)) +
  labs(x="Log Normalized PDSI Variability", y="Log Normalized Ecosystem Variability") +
  facet_wrap(~Model) +
  geom_point(size=0.05, alpha=0.2) +
  stat_smooth(method=lm, alpha=0.5) + 
  scale_color_manual(values=c("darkgreen", "darkolivegreen4", "darkslategray3", "darkseagreen3", "maroon2", "purple3")) +
  scale_fill_manual(values=c("darkgreen", "darkolivegreen4", "darkslategray3", "darkseagreen3", "maroon2", "purple3")) +
  guides(fill=guide_legend(nrow=1), color=guide_legend(nrow=1)) +
  # scale_fill_brewer(palette="Dark2") +
  # scale_color_brewer(palette="Dark2") +
  coord_cartesian(ylim=quantile(log(models.long$variability), c(0.001, 0.999)), expand=0) +
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
# -------------------------------------------

# -------------------------------------------
# Supplemental Figure 1: LBDA vs. Model Driver PDSI comparisons
# -------------------------------------------
lbda <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_LBDA_100.csv"))
lbda$var.pdsi <- lbda$diff.abs/abs(mean(lbda$lbda.mean, na.rm=T))  # Note: Positive numbers mean MORE stable
lbda <- lbda[!is.na(lbda$var.pdsi),]
summary(lbda)

stab.met   <- read.csv(file.path(path.google, "Current Data/Stability", "Stability_Drivers_100.csv"))
stab.met$var.pdsi <- stab.met$pdsi.diff/mean(stab.met$pdsi.diff)
summary(stab.met)

plot.lbda.yr <- ggplot(data=lbda[,]) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  geom_tile(aes(x=lon, y=lat, fill=1850-n.yrs+1)) +
    
  # geom_tile(aes(x=lon, y=lat, fill=cut(1850-n.yrs+1, breaks=seq(850-0.001, 1850, by=100)))) +
  geom_path(data=us, aes(x=long, y=lat, group=group), color="gray15") +
  geom_text(x=-98.5, y=48, label="a)", size=5, fontface="bold") +
  coord_equal(xlim=range(lbda$lon, na.rm=T), ylim=range(lbda$lat, na.rm=T), expand=0) +
  scale_fill_gradientn(name="LBDA\nEarliest\nYear A.D.", 
                    colors=c("#034e7b", "#0570b0", "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#ece7f2", "#fff7fb")) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))


plot.lbda.var <- ggplot(data=lbda[,]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  # geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_tile(aes(x=lon, y=lat, fill=log(var.pdsi))) +
  geom_path(data=us, aes(x=long, y=lat, group=group), color="gray15") +
  geom_text(x=-98.5, y=48, label="b)", size=5, fontface="bold") +
  coord_equal(xlim=range(lbda$lon, na.rm=T), ylim=range(lbda$lat, na.rm=T), expand=0) +
  scale_fill_gradientn(name="LBDA\nPDSI\nVariability", colors=c("#018571","#80cdc1", "#f5f5f5", "#dfc27d", "#a6611a"), limits=range(c(log(lbda$var.pdsi), log(stab.met$var.pdsi)), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

plot.drivers.var <- ggplot(data=stab.met[,]) +
  # facet_grid(dataset~.) +
  geom_polygon(data=us, aes(x=long, y=lat, group=group), fill="gray75") +
  # geom_point(aes(x=lon, y=lat, color=log(variability)), size=2) +
  geom_tile(aes(x=lon, y=lat, fill=log(var.pdsi))) +
  geom_path(data=us, aes(x=long, y=lat, group=group), color="gray15") +
  geom_text(x=-98.5, y=48, label="c)", size=5, fontface="bold") +
  coord_equal(xlim=range(lbda$lon, na.rm=T), ylim=range(lbda$lat, na.rm=T), expand=0) +
  scale_fill_gradientn(name="Driver\nPDSI\nVariability", colors=c("#018571","#80cdc1", "#f5f5f5", "#dfc27d", "#a6611a"), limits=range(c(log(lbda$var.pdsi), log(stab.met$var.pdsi)), na.rm=T)) +
  scale_y_continuous(breaks=c(40, 45, 50)) +
  theme(panel.background=element_rect(fill="gray25", color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(margin=unit(c(1,1,1,1), "lines"), color="black"),
        # axis.text.y = element_text(hjust=-1),
        axis.title=element_blank(),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.title= element_text(size=rel(0.75), face="bold"))

png(file.path(path.figs, "SupplementalFigure1_Variability_PDSI_LBDA_Models.png"), height=6, width=6, units="in", res=220)
plot_grid(plot.lbda.yr, plot.lbda.var, plot.drivers.var, ncol=1)
dev.off()
# -------------------------------------------


# -------------------------------------------
# Supplemental Figure 2: Richness vs. Variability in Empirical Datasets
# -------------------------------------------
dat.sites <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Variability_Ecosystem_v_Diversity_Data.csv"))
dat.sites <- dat.sites[!is.na(dat.sites$dom.pft),]
dat.sites$dataset2 <- car::recode(dat.sites$dataset2, "'STEPPS'='STEPPS-UMW'")
dat.sites$dom.pft <- as.character(dat.sites$dom.pft)
# Recde wasn't working, so doing this the ugly way
dat.sites[dat.sites$dom.pft=="BIRCH","dom.pft"] <- "Betula" #
dat.sites[dat.sites$dom.pft=="OTHER.HARDWOOD","dom.pft"] <- "Other Angio."
dat.sites[dat.sites$dom.pft=="OAK","dom.pft"] <- "Quercus" # 
dat.sites[dat.sites$dom.pft=="PINE","dom.pft"] <- "Pinus"
dat.sites[dat.sites$dom.pft=="HEMLOCK","dom.pft"] <- "Tsuga"
dat.sites[dat.sites$dom.pft=="BEECH","dom.pft"] <- "Fagus" #
dat.sites[dat.sites$dom.pft=="OTHER.CONIFER","dom.pft"] <- "Other Gymno."
dat.sites[dat.sites$dom.pft=="TAMARACK","dom.pft"] <- "Larix"
dat.sites[dat.sites$dom.pft=="ELM","dom.pft"] <- "Ulmus" #
dat.sites[dat.sites$dom.pft=="MAPLE","dom.pft"] <- "Acer" #
dat.sites[dat.sites$dom.pft=="ASH","dom.pft"] <- "Fraxinus" # 
dat.sites[dat.sites$dom.pft=="SPRUCE","dom.pft"] <- "Picea"
dat.sites[dat.sites$dom.pft=="CHESTNUT","dom.pft"] <- "Castanea" #
unique(dat.sites$dom.pft)
dat.sites$dom.pft <- factor(dat.sites$dom.pft, levels=c("Acer", "Fagus", "Betula", "Ulmus", "Fraxinus", "Quercus", "Castanea",  "Other Angio.", "Larix", "Picea", "Tsuga", "Pinus", "Other Gymno."))
summary(dat.sites)

colors.taxon <- data.frame(PFT=levels(dat.sites$dom.pft), color=NA)
row.names(colors.taxon) <- colors.taxon$PFT
colors.taxon["Other Gymno", "color"] <- "darkseagreen4"
colors.taxon["Larix","color"] <- "darkslategray4"
colors.taxon["Picea", "color"] <- "darkgreen"
colors.taxon["Tsuga", "color"] <- "darkolivegreen4"
colors.taxon["Pinus", "color"] <- "green3"
colors.taxon["Acer", "color"] <- "hotpink2"
colors.taxon["Fagus", "color"] <- "lightsalmon2"
colors.taxon["Betula", "color"] <- "goldenrod4"
colors.taxon["Ulmus", "color"] <- "tomato2"
colors.taxon["Fraxinus", "color"] <- "darkorange"
colors.taxon["Quercus", "color"] <- "indianred3"
colors.taxon["Castanea", "color"] <- "sienna3"
colors.taxon["Other Angio", "color"] <- "bisque2"


png(file.path(path.figs, "SupplementalFigure2_Variability_Data_Variability_v_Richness_NoTrendline.png"), height=6, width=6, units="in", res=220)
ggplot(data=dat.sites[!is.na(dat.sites$dom.pft),]) +
  facet_grid(dataset2~., scales="free_y") +
  geom_point(aes(x=richness, y=log(variability), color=dom.pft), position=position_jitter(width=0.1), size=0.75) +
  scale_x_continuous(name="Taxonomic Richness") +
  scale_y_continuous(name="Log Normalized Variability") +
  scale_color_manual(name="Dominant\nTaxa", values=colors.taxon$color) +
  theme(panel.background=element_rect(fill=NA, color="black", size=1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(-0.5, "lines"),
        # axis.ticks.margin = unit(0.5, "lines"),
        axis.text.x = element_text(size=rel(1.5), margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.text.y = element_text(size=rel(1.5), margin=unit(c(1,1,1,1), "lines"), color="black"),
        axis.title = element_text(size=rel(1.5), face="bold"),
        strip.text = element_text(size=rel(1.25), face="bold"),
        legend.position = "right",
        legend.key.height = unit(0.75, "lines"),
        legend.key = element_rect(fill=NA),
        legend.title= element_text(size=rel(0.75), face="bold"),
        legend.text=element_text(face="italic"))
dev.off()
# -------------------------------------------


# -------------------------------------------
# Table 1; results in par 1
# -------------------------------------------
# ------------
# Get needed data sets & set up empty data frame to save
# ------------
dat.emp <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Ecosystem_v_Climate_Data.csv"))
dat.emp$var <- as.factor(ifelse(dat.emp$dataset=="ReFAB", "bm", "fcomp"))
dat.emp$type <- "Empirical"
dat.emp$dataset2 <- car::recode(dat.emp$dataset2, "'STEPPS'='STEPPS-UMW'")
summary(dat.emp)

models.long <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_long.csv"))
models.long$Model <- car::recode(models.long$Model, "'TRIFFID'='JULES-TRIFFID'")

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
mod.v.dat$source <- factor(mod.v.dat$source, levels=c("ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "JULES-TRIFFID"))
mod.v.dat$source2 <- as.factor(ifelse(substr(mod.v.dat$source, 1, 6)=="STEPPS", "STEPPS", paste(mod.v.dat$source)))
mod.v.dat$source2 <- factor(mod.v.dat$source2, levels=c("STEPPS", "ReFAB", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "JULES-TRIFFID"))
summary(mod.v.dat)

ms.table1 <- data.frame(Source=c("Pollen", "ED2", "LINKAGES", "LPJ-GUESS", "LPJ-WSL", "JULES-TRIFFID"),
                        Var.Comp=NA, Sens.Comp=NA, Var.Biom=NA, Sens.Biom=NA)
# ------------

# ------------
# Composition: Mean & SD
# ------------
for(MOD in ms.table1$Source){
  if(MOD=="Pollen"){
    mod.mn <- round(mean(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & substr(mod.v.dat$source, 1, 6)=="STEPPS"]), na.rm=T), 3)
    mod.sd <- round(sd(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & substr(mod.v.dat$source, 1, 6)=="STEPPS"]), na.rm=T), 3)
    
    ms.table1[ms.table1$Source==MOD, "Var.Comp"] <- paste0(stringr::str_pad(mod.mn, ifelse(mod.mn<0, 6, 5), "right", pad=0), " (", stringr::str_pad(mod.sd, 5, "right", pad=0), ")")
    
  } else {
    mod.mn <- round(mean(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & mod.v.dat$source==MOD]), na.rm=T), 3)
    mod.sd <- round(sd(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & mod.v.dat$source==MOD]), na.rm=T), 3)
    
    ms.table1[ms.table1$Source==MOD, "Var.Comp"] <- paste0(stringr::str_pad(mod.mn, ifelse(mod.mn<0, 6, 5), "right", pad=0), " (", stringr::str_pad(mod.sd, 5, "right", pad=0), ")")
  }
  
}
# ------------

# ------------
# Composition-PDSI comparisons
# ------------
comp.sens <- lm(log(var.ecosys) ~ log(var.pdsi)*source2-source2 + source -log(var.pdsi) - 1, data=mod.v.dat[mod.v.dat$var=="fcomp",])
comp.sens.sum <- summary(comp.sens)

comp.sens.est <- round(comp.sens.sum$coefficients[paste0("log(var.pdsi):source2", levels(mod.v.dat$source2)[levels(mod.v.dat$source2) %in% unique(mod.v.dat$source2[mod.v.dat$var=="fcomp"])]),"Estimate"], 3)
comp.sens.se <- round(comp.sens.sum$coefficients[paste0("log(var.pdsi):source2", levels(mod.v.dat$source2)[levels(mod.v.dat$source2) %in% unique(mod.v.dat$source2[mod.v.dat$var=="fcomp"])]),"Std. Error"], 3)
comp.sens.p <- ifelse(comp.sens.sum$coefficients[paste0("log(var.pdsi):source2", levels(mod.v.dat$source2)[levels(mod.v.dat$source2) %in% unique(mod.v.dat$source2[mod.v.dat$var=="fcomp"])]),"Pr(>|t|)"]<0.05, "*", "")

ms.table1$Sens.Comp <- paste0(stringr::str_pad(comp.sens.est, ifelse(comp.sens.est<0, 6, 5), "right", pad=0), " (", stringr::str_pad(comp.sens.se, 5, "right", pad=0), ")", comp.sens.p)
# ms.table1$Sens.Comp.P <- paste0(comp.sens.p)

# ------------

# ------------
# Biomass: Mean & SD
# ------------
for(MOD in ms.table1$Source){
  if(MOD=="Pollen"){
    mod.mn <- round(mean(log(mod.v.dat$var.ecosys[mod.v.dat$var=="bm" & substr(mod.v.dat$source, 1, 6)=="ReFAB"]), na.rm=T), 3)
    mod.sd <- round(sd(log(mod.v.dat$var.ecosys[mod.v.dat$var=="bm" & substr(mod.v.dat$source, 1, 6)=="ReFAB"]), na.rm=T), 3)
    
    ms.table1[ms.table1$Source==MOD, "Var.Biom"] <- paste0(stringr::str_pad(mod.mn, ifelse(mod.mn<0, 6, 5), "right", pad=0), " (", stringr::str_pad(mod.sd, 5, "right", pad=0), ")")
    
  } else {
    mod.mn <- round(mean(log(mod.v.dat$var.ecosys[mod.v.dat$var=="bm" & mod.v.dat$source==MOD]), na.rm=T), 3)
    mod.sd <- round(sd(log(mod.v.dat$var.ecosys[mod.v.dat$var=="bm" & mod.v.dat$source==MOD]), na.rm=T), 3)
    
    ms.table1[ms.table1$Source==MOD, "Var.Biom"] <- paste0(stringr::str_pad(mod.mn, ifelse(mod.mn<0, 6, 5), "right", pad=0), " (", stringr::str_pad(mod.sd, 5, "right", pad=0), ")")
  }
  
}
# ------------

# ------------
# Biomass-PDSI Comparisons
# ------------
bm.sens <- lm(log(var.ecosys) ~ log(var.pdsi)*source-log(var.pdsi), data=mod.v.dat[mod.v.dat$var=="bm",])
bm.sens.sum <- summary(bm.sens)

bm.sens.est <- round(bm.sens.sum$coefficients[paste0("log(var.pdsi):source", levels(mod.v.dat$source)[levels(mod.v.dat$source) %in% unique(mod.v.dat$source[mod.v.dat$var=="bm"])]),"Estimate"], 3)
bm.sens.se <- round(bm.sens.sum$coefficients[paste0("log(var.pdsi):source", levels(mod.v.dat$source)[levels(mod.v.dat$source) %in% unique(mod.v.dat$source[mod.v.dat$var=="bm"])]),"Std. Error"], 3)
bm.sens.p <- ifelse(bm.sens.sum$coefficients[paste0("log(var.pdsi):source", levels(mod.v.dat$source)[levels(mod.v.dat$source) %in% unique(mod.v.dat$source[mod.v.dat$var=="bm"])]),"Pr(>|t|)"]<0.05, "*", "")

ms.table1$Sens.Biom <- paste0(stringr::str_pad(bm.sens.est, ifelse(bm.sens.est<0, 6, 5), "right", pad=0), " (", stringr::str_pad(bm.sens.se, 5, "right", pad=0), ")", bm.sens.p)
# ms.table1$Sens.Biom.P <- paste0(bm.sens.p)
# ------------

write.csv(ms.table1, file.path(path.tables, "Table1_Numbers.csv"), row.names=F)


# ------------
# Additional Comparisons for just the empirical data (in text)
# ------------
# Composition had the highest mean log-normalized variability, with the northeastern US (NEUS) more variable than the upper midwestern (UMW) US (NEUS: XXX SD XXX; UMW: XXX SD XXX)
round(mean(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & mod.v.dat$source=="STEPPS-NEUS"]), na.rm=T), 3); round(sd(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & mod.v.dat$source=="STEPPS-NEUS"]), na.rm=T), 3)

round(mean(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & mod.v.dat$source=="STEPPS-UMW"]), na.rm=T), 3); round(sd(log(mod.v.dat$var.ecosys[mod.v.dat$var=="fcomp" & mod.v.dat$source=="STEPPS-UMW"]), na.rm=T), 3)

summary(mod.v.dat)


# Compositional variability is weakly sensitive to drought variability in the northeastern US (slope = XXX SE XXX), but this relationship does not significantly differ from the insensitive pattern observed in the upper midwest (slope = XXX SE XXX)
comp.sens.poll <- lm(log(var.ecosys) ~ log(var.pdsi)*source -log(var.pdsi) -1, data=mod.v.dat[mod.v.dat$var=="fcomp" & mod.v.dat$type=="Empirical",])
comp.sens.poll.sum <- summary(comp.sens.poll)

comp.sens.poll.est <- round(comp.sens.poll.sum$coefficients[paste0("log(var.pdsi):source", levels(mod.v.dat$source)[levels(mod.v.dat$source) %in% unique(mod.v.dat$source[mod.v.dat$var=="fcomp" & mod.v.dat$type=="Empirical"])]),"Estimate"], 3)
comp.sens.poll.se <- round(comp.sens.poll.sum$coefficients[paste0("log(var.pdsi):source", levels(mod.v.dat$source)[levels(mod.v.dat$source) %in% unique(mod.v.dat$source[mod.v.dat$var=="fcomp" & mod.v.dat$type=="Empirical"])]),"Std. Error"], 3)
comp.sens.poll.p <- round(comp.sens.poll.sum$coefficients[paste0("log(var.pdsi):source", levels(mod.v.dat$source)[levels(mod.v.dat$source) %in% unique(mod.v.dat$source[mod.v.dat$var=="fcomp" & mod.v.dat$type=="Empirical"])]),"Pr(>|t|)"], 3)

comp.sens.poll.est
comp.sens.poll.se
comp.sens.poll.p

# Effects comparison
comp.sens.poll2 <- lm(log(var.ecosys) ~ log(var.pdsi)*source , data=mod.v.dat[mod.v.dat$var=="fcomp" & mod.v.dat$type=="Empirical",])
comp.sens.poll2.sum <- summary(comp.sens.poll2)

# In contrast to composition, drought variability was highest in the UMW (Fig 2), and mean log-normalized variability across both regions was considerably lower 
round(mean(log(mod.v.dat$var.pdsi[mod.v.dat$var=="fcomp" & mod.v.dat$type=="Empirical"]), na.rm=T), 3); round(sd(log(mod.v.dat$var.pdsi[mod.v.dat$var=="fcomp" & mod.v.dat$type=="Empirical"]), na.rm=T), 3)

round(mean(log(mod.v.dat$var.pdsi[mod.v.dat$var=="fcomp" & mod.v.dat$source=="STEPPS-UMW"]), na.rm=T), 3); round(sd(log(mod.v.dat$var.pdsi[mod.v.dat$var=="fcomp" & mod.v.dat$source=="STEPPS-UMW"]), na.rm=T), 3)

# Biomass variability correlated with compositional variability (Fig. 3c, slope = 0.48 SE 0.19;  R2=0.09;)
mod.dat2 <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_and_Data.csv"))
mod.dat2$source <- factor(mod.dat2$source, levels=c("ReFAB", "STEPPS-UMW", "STEPPS-NEUS", "ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "TRIFFID"))
# mod.dat2$var1 <- car::recode(mod.dat2$var1, "'PDSI'='Drought'")
mod.dat2$var1 <- factor(mod.dat2$var1, levels=c("PDSI", "Composition"))
mod.dat2$var2 <- factor(mod.dat2$var2, levels=c("Composition", "Biomass"))
summary(mod.dat2)

comp.v.bm1 <- lm(log(variability2) ~ log(variability1), data=mod.dat2[mod.dat2$var1=="Composition" & mod.dat2$var2=="Biomass" & mod.dat2$type=="Empirical",])
summary(comp.v.bm1)

summary(mod.dat2)
# ------------
# -------------------------------------------

# -------------------------------------------
# Supplemental Table 1: Sensitivity of model latent states to PDSI
# -------------------------------------------
models.long <- read.csv(file.path(path.google, "Current Data/Stability_Synthesis", "Stability_Models_long.csv"))
models.long$Model <- car::recode(models.long$Model, "'TRIFFID'='JULES-TRIFFID'")
models.long$Model <- factor(models.long$Model, levels=c("ED2", "LPJ-WSL", "LPJ-GUESS", "LINKAGES", "JULES-TRIFFID"))
models.long$var <- car::recode(models.long$var, "'gpp'='GPP'; 'npp'='NPP'; 'lai'='LAI'; 'bm'='Biomass'; 'fcomp'='Composition'; 'nee'='NEE'")
models.long$var <- factor(models.long$var, levels=c("GPP", "NPP", "NEE", "LAI", "Biomass", "Composition"))
summary(models.long)

supp.table1 <- data.frame(Model=levels(models.long$Model), GPP=NA, NPP=NA, NEE=NA, LAI=NA, Biomass=NA, Composition=NA)

for(VAR in unique(models.long$var)){
  for(MOD in unique(models.long$Model)){
    mod.sub <- models.long[models.long$var==VAR & models.long$Model==MOD,]
    
    if(nrow(mod.sub)==0) next
    mod.tmp <- lm(log(var.rel) ~ log(var.pdsi), data=mod.sub)
    mod.tmp.sum <- summary(mod.tmp)
    
    sens.est <- round(mod.tmp.sum$coefficients["log(var.pdsi)","Estimate"], 3)
    sens.se <- round(mod.tmp.sum$coefficients["log(var.pdsi)","Std. Error"], 3)
    sens.p <- ifelse(mod.tmp.sum$coefficients["log(var.pdsi)","Pr(>|t|)"]<0.05, "*", "")
    
    supp.table1[supp.table1$Model==MOD, VAR] <- paste0(stringr::str_pad(sens.est, ifelse(sens.est<0, 6, 5), "right", pad=0), " (", stringr::str_pad(sens.se, 5, "right", pad=0), ")", sens.p)
  }
}
supp.table1

write.csv(supp.table1, file.path(path.tables, "SupplementalTable1_Numbers.csv"), row.names=F)

# -------------------------------------------
