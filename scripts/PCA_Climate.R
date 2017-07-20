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
# 2. Correlation analysis to try and reduce # variables
# 3. Perform Oridination (PCA)
#    - Extract loadings & PCs
# 4. Save output, 
#    - generate & save a couple figures
# --------------------------------------------

# --------------------------------------------
# 0. Define file paths etc
# --------------------------------------------
# Path to github repository/working directory
path.repo <- "~/Desktop/Research/EcosystemStability/"
setwd(path.repo)

# Path to where data are
path.data <- "data/"

# This is the path to my pdsi code is, which has the formula for calculating WHC
# This is currently part of my met ensemble code and code 
# and can be found here: https://github.com/PalEON-Project/modeling_met_ensemble.git
path.pdsi <- "~/Dropbox/PalEON_CR/met_ensemble/scripts/"

# --------------------------------------------


# --------------------------------------------
# 1. Reading in & extracting data
# --------------------------------------------
# loading the key to dimensions in the .RDS files; this is Christy's sheet generated
# this loads an object called "paleon"
paleon <- read.csv(file.path(path.repo, "data/paleon_models_environment_master.csv")) 
paleon$latlon <- as.factor(paleon$latlon)
summary(paleon)
# --------------------------------------------

# --------------------------------------------
# 2. Correlation analysis to reduce dimensionality & spread out variance
# --------------------------------------------
# Making table with variables, definitions, and groups
table.vars <- data.frame(var = c("sand.t", "sand.s", "clay.t", "clay.s", "depth", "awc.t", "awc.s", "whc.t", "whc.s", "whc.tot", "tair.yr.set", "tair.yr.cru", "tair.jja.set", "tair.jja.cru", "precip.yr.set", "precip.yr.cru", "precip.jja.set", "precip.jja.cru"),
                         class = c(rep("soil",10), rep("temperature", 4), rep("precipitation", 4)))
names(paleon)[7:24]
env.cor <- cor(paleon[,paste(table.vars$var)])
env.cor[upper.tri(env.cor, diag = T)] <- NA
summary(env.cor)
dim(env.cor)

library(reshape2)
cor.melt <- melt(env.cor, na.rm = TRUE)
# cor.melt <- merge(cor.melt, table.vars)
for(v in unique(table.vars$class)){
  cor.melt[cor.melt$Var1 %in% table.vars[table.vars$class==v, "var"], "Class.Var1"] <- v
  cor.melt[cor.melt$Var2 %in% table.vars[table.vars$class==v, "var"], "Class.Var2"] <- v
}
cor.melt$Class.Var1 <- factor(cor.melt$Class.Var1, levels=c("soil", "temperature", "precipitation"))
cor.melt$Class.Var2 <- factor(cor.melt$Class.Var2, levels=c("precipitation", "temperature", "soil"))
summary(cor.melt)


library(ggplot2); library(grid)
png("figures/environment_correlation_matrix.png", height=8, width=10, units="in", res=320)
ggplot(data=cor.melt) +
  facet_grid(Class.Var2 ~ Class.Var1, scales="free") +
  geom_tile(aes(x=Var1, y=Var2, fill=value), color="white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  coord_fixed()
dev.off()

pdf("figures/environment_correlation_matrix_by_class.pdf", height=8, width=10)
for(v in unique(cor.melt$Class.Var1)){
  print(
  ggplot(data=cor.melt[cor.melt$Class.Var1==v & cor.melt$Class.Var2==v,]) +
    # facet_wrap(~class, scales="free", ncol=1)+
    geom_tile(aes(x=Var1, y=Var2, fill=value), color="white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") + 
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
    coord_fixed() + ggtitle(v)
  )
}
dev.off()

write.csv(cor.melt, "data/environment_correlation_matrix.csv", row.names=F)

# Something's being odd, so we're just going to paste to add quotes
# cor.melt$Var1 <- paste(cor.melt$Var1)
# cor.melt$Var2 <- paste(cor.melt$Var2)
# Getting the mean & max correlation of each factor
for(i in unique(cor.melt$Var1)){
  cor.mean <- mean(abs(cor.melt[cor.melt$Var1==i,"value"]),na.rm=T)
  cor.max  <- max(abs(cor.melt[cor.melt$Var1==i,"value"]),na.rm=T)
  table.vars[table.vars$var==i, "cor.mean.abs"   ] <- cor.mean
  table.vars[table.vars$var==i, "cor.max.abs"    ] <- cor.max
  table.vars[table.vars$var==i, "cor.max.id" ] <- cor.melt[cor.melt$Var1==i & abs(cor.melt$value)==cor.max,"Var2"][1] # just pick the first one
}
table.vars
# summary(table.vars)
cor.low <- table.vars[table.vars$cor.mean.abs<=quantile(abs(table.vars$cor.mean.abs), 0.33, na.rm=T),]

# Pulling the least-correlated variable on average from each class
cor.melt[cor.melt$Var1=="whc.tot" & cor.melt$Class.Var2=="soil",]
table.vars[!table.vars$var %in% c("whc.tot", "whc.s", "awc.s") & table.vars$class=="soil",]

vars.low <- c("whc.tot", "tair.yr.set", "precip.yr.set")
vars.low2 <- c("awc.t", "tair.jja.set", "precip.jja.set")

summary(cor.melt[cor.melt$Var1 %in% c(vars.low, vars.low2) & cor.melt$Var2 %in% c(vars.low, vars.low2),])
# --------------------------------------------


# --------------------------------------------
# 3. Perform Oridination (PCA)
#
# We're going to do a couple versions for now:
# 2.a. All climate & env variables
# 2.b. Summer climate + top env + depth
# 2.c. Model sites only
# --------------------------------------------
library(labdsv)
set.seed(1333)

pc.all <- princomp(~ whc.tot + awc.t + tair.yr.set + tair.jja.set + precip.yr.set + precip.jja.set,
                   data=paleon, cor=T, scores=T)
pc.setveg <- princomp(~ whc.tot + awc.t + tair.yr.set + tair.jja.set + precip.yr.set + precip.jja.set,
                      data=paleon[!is.na(paleon$domain.paleon),], cor=T, scores=T)
pc.setveg.mod <- princomp(~ whc.tot + awc.t + tair.yr.set + tair.jja.set + precip.yr.set + precip.jja.set,
                          data=paleon[!is.na(paleon$domain.paleon) & !is.na(paleon$umw),], cor=T, scores=T)
pc.umw <- princomp(~ whc.tot + awc.t + tair.yr.set + tair.jja.set + precip.yr.set + precip.jja.set,
                   data=paleon[paleon$domain.paleon %in% c("MN", "WI", "MI - upper"),], cor=T, scores=T)
pc.umw.mod <- princomp(~ whc.tot + awc.t + tair.yr.set + tair.jja.set + precip.yr.set + precip.jja.set,
                          data=paleon[paleon$umw=="y",], cor=T, scores=T)
summary(pc.all)
summary(pc.setveg)
summary(pc.setveg.mod)
summary(pc.umw)
summary(pc.umw.mod)

# Looking at & loadings
pc.all$loadings
pc.setveg$loadings
pc.setveg.mod$loadings
pc.umw$loadings
pc.umw.mod$loadings


# Extracting site scores
paleon[,"pc.all1"] <- pc.all$scores[,1]
paleon[,"pc.all2"] <- pc.all$scores[,2]
paleon[!is.na(paleon$domain.paleon),"pc.setveg1"] <- pc.setveg$scores[,1]
paleon[!is.na(paleon$domain.paleon),"pc.setveg2"] <- pc.setveg$scores[,2]
paleon[!is.na(paleon$domain.paleon) & !is.na(paleon$umw),"pc.setveg.mod1"] <- pc.setveg.mod$scores[,1]
paleon[!is.na(paleon$domain.paleon) & !is.na(paleon$umw),"pc.setveg.mod2"] <- pc.setveg.mod$scores[,2]
paleon[paleon$domain.paleon %in% c("MN", "WI", "MI - upper"),"pc.umw1"] <- pc.umw$scores[,1]
paleon[paleon$domain.paleon %in% c("MN", "WI", "MI - upper"),"pc.umw2"] <- pc.umw$scores[,2]
paleon[paleon$umw=="y" & !is.na(paleon$umw),"pc.umw.mod1"] <- pc.umw.mod$scores[,1]
paleon[paleon$umw=="y" & !is.na(paleon$umw),"pc.umw.mod2"] <- pc.umw.mod$scores[,2]
# --------------------------------------------


# --------------------------------------------
# 3. Save output, 
#    - generate & save a couple figures
# --------------------------------------------
write.csv(paleon[,names(paleon)[!names(paleon) %in% table.vars$var]], "data/PCA_Scores_Environment_AllRegion.csv", row.names=F, eol="\r\n")

  # Manually calculating the percentage of variance
  pc.all.POV <- pc.all$sdev^2/sum(pc.all$sdev^2)
  pc.setveg.POV <- pc.setveg$sdev^2/sum(pc.setveg$sdev^2)
  pc.setveg.mod.POV <- pc.setveg.mod$sdev^2/sum(pc.setveg.mod$sdev^2)
  pc.umw.POV <- pc.umw$sdev^2/sum(pc.umw$sdev^2)
  pc.umw.mod.POV <- pc.umw.mod$sdev^2/sum(pc.umw.mod$sdev^2)
  
# Graphing spatial scores
{
  pc.setveg.mod1   <- ggplot(data=paleon) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(aes(x=lon, y=lat, fill=pc.all1)) + 
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
    ggtitle(paste0("Region PCA, PC1 (", round(pc.all.POV["Comp.1"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  pc.all2   <- ggplot(data=paleon) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(aes(x=lon, y=lat, fill=pc.all2)) + 
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC2") + 
    ggtitle(paste0("Region PCA, PC2 (", round(pc.all.POV["Comp.2"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  
  pc.setveg1   <- ggplot(data=paleon[!is.na(paleon$pc.setveg1),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon, aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.setveg1)) + 
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
    ggtitle(paste0("Region PCA, PC1 (", round(pc.setveg.POV["Comp.1"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  pc.setveg2   <- ggplot(data=paleon[!is.na(paleon$pc.setveg2),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon, aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.setveg2)) + 
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC2") + 
    ggtitle(paste0("Region PCA, PC2 (", round(pc.setveg.POV["Comp.2"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  
  pc.setveg.mod1   <- ggplot(data=paleon[!is.na(paleon$pc.setveg.mod1),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon[!is.na(paleon$pc.setveg1),], aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.setveg.mod1)) + 
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
    ggtitle(paste0("Region PCA, PC1 (", round(pc.setveg.mod.POV["Comp.1"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  pc.setveg.mod2   <- ggplot(data=paleon[!is.na(paleon$pc.setveg.mod2),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon[!is.na(paleon$pc.setveg1),], aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.setveg.mod2)) + 
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC2") + 
    ggtitle(paste0("Region PCA, PC2 (", round(pc.setveg.mod.POV["Comp.2"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  
  pc.umw1   <- ggplot(data=paleon[!is.na(paleon$pc.umw1),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon[,], aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.umw1)) + 
    scale_x_continuous(limits=range(paleon[!is.na(paleon$pc.umw1),"lon"], na.rm=T)) +
    scale_y_continuous(limits=range(paleon[!is.na(paleon$pc.umw1),"lat"], na.rm=T)) +
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
    ggtitle(paste0("Region PCA, PC1 (", round(pc.umw.POV["Comp.1"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  
  pc.umw2   <- ggplot(data=paleon[!is.na(paleon$pc.umw2),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon[,], aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.umw2)) + 
    scale_x_continuous(limits=range(paleon[!is.na(paleon$pc.umw2),"lon"], na.rm=T)) +
    scale_y_continuous(limits=range(paleon[!is.na(paleon$pc.umw2),"lat"], na.rm=T)) +
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC2") + 
    ggtitle(paste0("Region PCA, PC2 (", round(pc.umw.POV["Comp.2"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  
  pc.umw.mod1   <- ggplot(data=paleon[!is.na(paleon$pc.umw.mod1),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon[!is.na(paleon$pc.umw1),], aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.umw.mod1)) + 
    scale_x_continuous(limits=range(paleon[!is.na(paleon$pc.umw1),"lon"], na.rm=T)) +
    scale_y_continuous(limits=range(paleon[!is.na(paleon$pc.umw1),"lat"], na.rm=T)) +
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
    ggtitle(paste0("Region PCA, PC1 (", round(pc.umw.mod.POV["Comp.1"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  
  pc.umw.mod2   <- ggplot(data=paleon[!is.na(paleon$pc.umw.mod2),]) + theme_bw() + coord_equal(expand=0) + 
    geom_raster(data=paleon[!is.na(paleon$pc.umw2),], aes(x=lon, y=lat), fill="gray50") +
    geom_raster(aes(x=lon, y=lat, fill=pc.umw.mod2)) + 
    scale_x_continuous(limits=range(paleon[!is.na(paleon$pc.umw2),"lon"], na.rm=T)) +
    scale_y_continuous(limits=range(paleon[!is.na(paleon$pc.umw2),"lat"], na.rm=T)) +
    scale_fill_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC2") + 
    ggtitle(paste0("Region PCA, PC2 (", round(pc.umw.mod.POV["Comp.2"]*100,1),"%)")) + 
    theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
  
  
  # Graphing the site scores
  png("figures/PC_Scores_Env_Region_All.png", height=5, width=5, units="in", res=320)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(pc.all1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc.all2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  dev.off()
  
  png("figures/PC_Scores_Env_SetVeg_All.png", height=5, width=5, units="in", res=320)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(pc.umw.mod1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc.setveg2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  dev.off()
  
  png("figures/PC_Scores_Env_SetVeg_Mod.png", height=5, width=5, units="in", res=320)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(pc.setveg.mod1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc.setveg.mod2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  dev.off()
  
  png("figures/PC_Scores_Env_UMW_All.png", height=5, width=5, units="in", res=320)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(pc.umw1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc.umw2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  dev.off()
  
  png("figures/PC_Scores_Env_UMW_Mod.png", height=5, width=5, units="in", res=320)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1)))
  print(pc.umw.mod1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(pc.umw.mod2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  dev.off()

}

# Making a dataframe to plot the loadings in ggplot
# adapted from: https://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
# (I'm doing it the long way because I want more control over customization thatn I had previously)
pc.all.loadings <- data.frame(var=row.names(pc.all$loadings), pc.all$loadings[,1:4])
mult <- min( (max(paleon$pc.all2)-min(paleon$pc.all2))/(max(pc.all.loadings$Comp.2)-min(pc.all.loadings$Comp.2)),
             (max(paleon$pc.all1)-min(paleon$pc.all1))/(max(pc.all.loadings$Comp.1)-min(pc.all.loadings$Comp.1))
             )
pc.all.loadings$load.PC1 <- 0.7*mult*pc.all.loadings$Comp.1
pc.all.loadings$load.PC2 <- 0.7*mult*pc.all.loadings$Comp.2
pc.all.loadings$var <- c("Water Cap", "AWC (top)", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")

pc.setveg.loadings <- data.frame(var=row.names(pc.setveg$loadings), pc.setveg$loadings[,1:4])
mult <- min( (max(paleon$pc.setveg2, na.rm=T)-min(paleon$pc.setveg2, na.rm=T))/(max(pc.setveg.loadings$Comp.2, na.rm=T)-min(pc.setveg.loadings$Comp.2, na.rm=T)),
             (max(paleon$pc.setveg1, na.rm=T)-min(paleon$pc.setveg1, na.rm=T))/(max(pc.setveg.loadings$Comp.1, na.rm=T)-min(pc.setveg.loadings$Comp.1, na.rm=T))
             )
pc.setveg.loadings$load.PC1 <- 0.7*mult*pc.setveg.loadings$Comp.1
pc.setveg.loadings$load.PC2 <- 0.7*mult*pc.setveg.loadings$Comp.2
pc.setveg.loadings$var <- c("Water Cap", "AWC (top)", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")

pc.setveg.mod.loadings <- data.frame(var=row.names(pc.setveg.mod$loadings), pc.setveg.mod$loadings[,1:4])
mult <- min( (max(paleon$pc.setveg.mod2, na.rm=T)-min(paleon$pc.setveg.mod2, na.rm=T))/(max(pc.setveg.mod.loadings$Comp.2, na.rm=T)-min(pc.setveg.mod.loadings$Comp.2, na.rm=T)),
             (max(paleon$pc.setveg.mod1, na.rm=T)-min(paleon$pc.setveg.mod1, na.rm=T))/(max(pc.setveg.mod.loadings$Comp.1, na.rm=T)-min(pc.setveg.mod.loadings$Comp.1, na.rm=T))
)
pc.setveg.mod.loadings$load.PC1 <- 0.7*mult*pc.setveg.mod.loadings$Comp.1
pc.setveg.mod.loadings$load.PC2 <- 0.7*mult*pc.setveg.mod.loadings$Comp.2
pc.setveg.mod.loadings$var <- c("Water Cap", "AWC (top)", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")

pc.umw.loadings <- data.frame(var=row.names(pc.umw$loadings), pc.umw$loadings[,1:4])
mult <- min( (max(paleon$pc.umw2, na.rm=T)-min(paleon$pc.umw2, na.rm=T))/(max(pc.umw.loadings$Comp.2, na.rm=T)-min(pc.umw.loadings$Comp.2, na.rm=T)),
             (max(paleon$pc.umw1, na.rm=T)-min(paleon$pc.umw1, na.rm=T))/(max(pc.umw.loadings$Comp.1, na.rm=T)-min(pc.umw.loadings$Comp.1, na.rm=T))
)
pc.umw.loadings$load.PC1 <- 0.7*mult*pc.umw.loadings$Comp.1
pc.umw.loadings$load.PC2 <- 0.7*mult*pc.umw.loadings$Comp.2
pc.umw.loadings$var <- c("Water Cap", "AWC (top)", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")

pc.umw.mod.loadings <- data.frame(var=row.names(pc.umw.mod$loadings), pc.umw.mod$loadings[,1:4])
mult <- min( (max(paleon$pc.umw.mod2, na.rm=T)-min(paleon$pc.umw.mod2, na.rm=T))/(max(pc.umw.mod.loadings$Comp.2, na.rm=T)-min(pc.umw.mod.loadings$Comp.2, na.rm=T)),
             (max(paleon$pc.umw.mod1, na.rm=T)-min(paleon$pc.umw.mod1, na.rm=T))/(max(pc.umw.mod.loadings$Comp.1, na.rm=T)-min(pc.umw.mod.loadings$Comp.1, na.rm=T))
)
pc.umw.mod.loadings$load.PC1 <- 0.7*mult*pc.umw.mod.loadings$Comp.1
pc.umw.mod.loadings$load.PC2 <- 0.7*mult*pc.umw.mod.loadings$Comp.2
pc.umw.mod.loadings$var <- c("Water Cap", "AWC (top)", "Tair (yr)", "Tair (JJA)", "Precip (yr)", "Precip (JJA)")

# --------------
# Manually creating a biplots
# Weights cutoff = 75th percentile for PC1 & 2 (top 25% of weights)
# --------------
# pc.all.cutoff <- quantile(c(abs(pc.all.loadings$Comp.1), abs(pc.all.loadings$Comp.2)), 0.75, na.rm=T)
pc.all.cutoff <- 0
biplot.pc.all <- ggplot() +
  theme_bw() +# coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon, aes(x=pc.all1, y=pc.all2), size=0.2, color="black") +
  geom_segment(data=pc.all.loadings[abs(pc.all.loadings$Comp.1)>pc.all.cutoff | abs(pc.all.loadings$Comp.2)>pc.all.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=pc.all.loadings[abs(pc.all.loadings$Comp.1)>pc.all.cutoff | abs(pc.all.loadings$Comp.2)>pc.all.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(pc.all.POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.all.pc1, c(0.001, 0.999))) +
  scale_y_continuous(name=paste0("PC2 (", round(pc.all.POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.all.pc2, c(0.001, 0.999))) +
  ggtitle("Region, All") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))

# pc.setveg.cutoff <- quantile(c(abs(pc.setveg.loadings$Comp.1), abs(pc.setveg.loadings$Comp.2)), 0.75, na.rm=T)
pc.setveg.cutoff <- 0
biplot.pc.setveg <- ggplot() +
  theme_bw() +# coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon[,], aes(x=pc.setveg1, y=pc.setveg2), size=0.2, color="black") +
  geom_segment(data=pc.setveg.loadings[abs(pc.setveg.loadings$Comp.1)>pc.setveg.cutoff | abs(pc.setveg.loadings$Comp.2)>pc.setveg.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=pc.setveg.loadings[abs(pc.setveg.loadings$Comp.1)>pc.setveg.cutoff | abs(pc.setveg.loadings$Comp.2)>pc.setveg.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(pc.setveg.POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.setveg.pc1, c(0.001, 0.999))) +
  scale_y_continuous(name=paste0("PC2 (", round(pc.setveg.POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.setveg.pc2, c(0.001, 0.999))) +
  ggtitle("Settlement Vegetation, All") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))

# pc.setveg.mod.cutoff <- quantile(c(abs(pc.setveg.mod.loadings$Comp.1), abs(pc.setveg.mod.loadings$Comp.2)), 0.75, na.rm=T)
pc.setveg.mod.cutoff <- 0
biplot.pc.setveg.mod <- ggplot() +
  theme_bw() +# coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon[,], aes(x=pc.setveg.mod1, y=pc.setveg.mod2), size=0.2, color="black") +
  geom_segment(data=pc.setveg.mod.loadings[abs(pc.setveg.mod.loadings$Comp.1)>pc.setveg.mod.cutoff | abs(pc.setveg.mod.loadings$Comp.2)>pc.setveg.mod.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=pc.setveg.mod.loadings[abs(pc.setveg.mod.loadings$Comp.1)>pc.setveg.mod.cutoff | abs(pc.setveg.mod.loadings$Comp.2)>pc.setveg.mod.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(pc.setveg.mod.POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.setveg.mod.pc1, c(0.001, 0.999))) +
  scale_y_continuous(name=paste0("PC2 (", round(pc.setveg.mod.POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.setveg.mod.pc2, c(0.001, 0.999))) +
  ggtitle("Settlement Vegetation, Models") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))


# pc.umw.cutoff <- quantile(c(abs(pc.umw.loadings$Comp.1), abs(pc.umw.loadings$Comp.2)), 0.75, na.rm=T)
pc.umw.cutoff <- 0
biplot.pc.umw <- ggplot() +
  theme_bw() +# coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon[,], aes(x=pc.umw1, y=pc.umw2), size=0.2, color="black") +
  geom_segment(data=pc.umw.loadings[abs(pc.umw.loadings$Comp.1)>pc.umw.cutoff | abs(pc.umw.loadings$Comp.2)>pc.umw.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=pc.umw.loadings[abs(pc.umw.loadings$Comp.1)>pc.umw.cutoff | abs(pc.umw.loadings$Comp.2)>pc.umw.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(pc.umw.POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.umw.pc1, c(0.001, 0.999))) +
  scale_y_continuous(name=paste0("PC2 (", round(pc.umw.POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.umw.pc2, c(0.001, 0.999))) +
  ggtitle("Upper Midwest, All") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))

# pc.umw.mod.cutoff <- quantile(c(abs(pc.umw.mod.loadings$Comp.1), abs(pc.umw.mod.loadings$Comp.2)), 0.75, na.rm=T)
pc.umw.mod.cutoff <- 0
biplot.pc.umw.mod <- ggplot() +
  theme_bw() +# coord_equal() +
  geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
  geom_point(data=paleon[,], aes(x=pc.umw.mod1, y=pc.umw.mod2), size=0.2, color="black") +
  geom_segment(data=pc.umw.mod.loadings[abs(pc.umw.mod.loadings$Comp.1)>pc.umw.mod.cutoff | abs(pc.umw.mod.loadings$Comp.2)>pc.umw.mod.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
               arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
  geom_text(data=pc.umw.mod.loadings[abs(pc.umw.mod.loadings$Comp.1)>pc.umw.mod.cutoff | abs(pc.umw.mod.loadings$Comp.2)>pc.umw.mod.cutoff,], 
            aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
  scale_x_continuous(name=paste0("PC1 (", round(pc.umw.mod.POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.umw.mod.pc1, c(0.001, 0.999))) +
  scale_y_continuous(name=paste0("PC2 (", round(pc.umw.mod.POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$pc.umw.mod.pc2, c(0.001, 0.999))) +
  ggtitle("Upper Midwest, Models") +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))



png("figures/PCA_Loadings.png", height=10, width=8, units="in", res=320)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
print(biplot.pc.all       , vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(biplot.pc.setveg    , vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(biplot.pc.setveg.mod, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(biplot.pc.umw       , vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(biplot.pc.umw.mod   , vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()
# --------------------------------------------
