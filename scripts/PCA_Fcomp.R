library(reshape2)
library(labdsv)
library(raster)
library(geosphere)

set.seed(1333)

#######################################################################################################################
## Load meta data and fcomp values
#######################################################################################################################

# loading the key to dimensions in the .RDS files; this is Christy's sheet generated
# this loads an object called "paleon"
paleon <- read.csv("data/paleon_site_info.csv") 
paleon$latlon <- as.factor(paleon$latlon)
paleon <- paleon[,c("X", "num", "lon", "lat", "latlon", "umw", "x", "y")]

# store all pca scores
fcomp_pca = data.frame(lon=numeric(0), lat=numeric(0), model=character(0))

# get coords for full paleon domain
paleon.all <- read.csv('data/PCA_Scores_Environment_AllRegion.csv')

# determine which grid cells from 0.5 degree grid are in the UMW
# base = raster("data/paleon.unit.ll_01.tif")
base = raster("data/paleon_full_ll_v0.1.tif")
base_matrix = rasterToPoints(base)
colnames(base_matrix) = c('lon', 'lat', 'state')
base_matrix = base_matrix[which(base_matrix[,'state'] %in% c(5, 6, 11, 12)),] # state codes for UMW; see paleon wiki for code list
dist_mat = distm(as.matrix(cbind(paleon.all$lon, paleon.all$lat)), as.matrix(cbind(base_matrix[,1], base_matrix[,2])), fun = distHaversine)/1000
in_umw   = apply(dist_mat, 1, function(x) any(x<56))
paleon.all$umwR = in_umw

# check UMW domain
plot(paleon.all$lon[which(paleon.all$umwR)], paleon.all$lat[which(paleon.all$umwR)])

# read in Fcomp for models and STEPPS
fcomp  = readRDS('data/fcomp-all.RDS')
load('data/input.rdata')

centers_veg_alb = centers_veg*1e6
coordinates(centers_veg_alb) <- ~ x + y
proj4string(centers_veg_alb) <- CRS('+init=epsg:3175')#CRS('+proj=longlat +ellps=WGS84')

centers_veg_ll <- spTransform(centers_veg_alb, CRS('+proj=longlat +ellps=WGS84'))
centers_veg_ll <- data.frame(centers_veg_ll)
colnames(centers_veg_ll) = c('lon', 'lat')

#######################################################################################################################
## PCA for models for full domain
#######################################################################################################################

models = unique(fcomp$model)
models = models[which(models != 'STEPPS')]
pc_out = vector("list", length=length(models))
names(pc_out) = models

for (i in 1:length(models)){
  model=models[i]
  
  fcomp_mod = fcomp[which((fcomp$model == model) & (fcomp$year == 1750)),]
  
  pfts = unique(fcomp_mod$pft)
  
  fcomp_cast = dcast(fcomp_mod, model + x + y + cell_id ~ pft)
  
  pfts_drop = names(which(colSums(fcomp_cast[,pfts]) == 0))
  fcomp_cast = fcomp_cast[,!(colnames(fcomp_cast) %in% pfts_drop)]
  fcomp_cast = fcomp_cast[with(fcomp_cast, order(cell_id)),]
  
  pfts = pfts[!(pfts %in% pfts_drop)]
  
  pc <- princomp(fcomp_cast[,pfts], cor=T, scores=T)
  pc[['cell_id']] = fcomp_cast$cell_id
  
  pc[['coords']] = data.frame(lat = paleon$lat[pc$cell_id], lon=paleon$lon[pc$cell_id])
 
  pc_out[[i]] = pc
}

# --------------------------------------------
# Plot the scores
# --------------------------------------------
plot_scores <- function(pc, paleon, paleon.all, umw){
  POV <- pc$sdev^2/sum(pc$sdev^2)
  
  if (ncol(pc$scores) < 3){nscores=ncol(pc$scores)} else {nscores=3}
  
  if (umw) {paleon.all = paleon.all[which(paleon.all$umwR),]}
  
  # store first three principal components
  scores_plot = vector("list", length=3)
  for (i in 1:nscores) {
    pca_scores = data.frame(lat=pc$coords$lat, lon=pc$coords$lon, pc = pc$scores[,i])  
    p   <- ggplot(data=pca_scores) + theme_bw() + coord_equal(expand=0) + 
      geom_raster(data=paleon.all, aes(x=lon, y=lat), fill="gray50") +
      geom_point(aes(x=lon, y=lat, color=pc), size=2, alpha=0.8)  + 
      scale_color_gradient2(low="blue2", high="red3", mid="gray80", midpoint=0, name="PC1") + 
      ggtitle(paste0(model, " , Settlement Era PC", i, " (", round(POV[paste0("Comp.",i)]*100,1),"%)")) + 
      theme(plot.title=element_text(face="bold", hjust=0.5, size=rel(1)))
    
    scores_plot[[i]] = p
  }
  
  return(scores_plot)
}

scores_plot = vector("list", length=length(models))
for (i in 1:length(models)){
  model = as.vector(models[i])

  scores_plot[[i]] = plot_scores(pc_out[[model]], paleon, paleon.all, umw=FALSE)
}

# Graphing the site scores
pdf("figures/PC_Scores_Fcomp_Full.pdf", height=6, width=12)
for (page in 1:ceiling(length(models)/2)){
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3)))
  idx_model = page*2 - 1
  for (i in idx_model:(idx_model+1)){
    if (i %% 2 == 0){layout = 2} else {layout = 1}
    for (j in 1:3){
      print(scores_plot[[i]][[j]], vp = viewport(layout.pos.row = layout, layout.pos.col = j))
    }
  }
}
dev.off()

#######################################################################################################################
## PCA for models and STEPPS for UMW
#######################################################################################################################

models = unique(fcomp$model)
# models = models[which(models != 'STEPPS')]
pc_out = vector("list", length=length(models))
names(pc_out) = models

for (i in 1:length(models)){
  model=models[i]
  fcomp_mod = fcomp[which((fcomp$model == model) & (fcomp$year == 1750)),]
  
  if (model != 'STEPPS'){
    fcomp_mod = fcomp_mod[which(paleon$umw[fcomp_mod$cell_id] == "y"),]
  }
  
  pfts = unique(fcomp_mod$pft)
  
  fcomp_cast = dcast(fcomp_mod, model + x + y + cell_id ~ pft)
  
  pfts_drop = names(which(colSums(fcomp_cast[,pfts]) == 0))
  fcomp_cast = fcomp_cast[,!(colnames(fcomp_cast) %in% pfts_drop)]
  fcomp_cast = fcomp_cast[with(fcomp_cast, order(cell_id)),]
  
  pfts = pfts[!(pfts %in% pfts_drop)]
  
  pc <- princomp(fcomp_cast[,pfts], cor=T, scores=T)
  
  pc[['cell_id']] = fcomp_cast$cell_id
  
  if (model == 'STEPPS'){
    pc[['coords']] = data.frame(lat = centers_veg_ll$lat[pc$cell_id], lon=centers_veg_ll$lon[pc$cell_id])
  } else {
    pc[['coords']] = data.frame(lat = paleon$lat[pc$cell_id], lon=paleon$lon[pc$cell_id])
  }
  
  pc_out[[i]] = pc
}

# --------------------------------------------
# Plot the scores
# --------------------------------------------

scores_plot = vector("list", length=length(models))
for (i in 1:length(models)){
  model = as.vector(models[i])
  scores_plot[[i]] = plot_scores(pc_out[[model]], paleon, paleon.all, umw=TRUE)
}

# Graphing the site scores
pdf("figures/PC_Scores_Fcomp_UMW.pdf", height=6, width=12)
for (page in 1:ceiling(length(models)/2)){
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3)))
  idx_model = page*2 - 1
  for (i in idx_model:(idx_model+1)){
    if (i > length(models)) {break}
    if (i %% 2 == 0){layout = 2} else {layout = 1}
    for (j in 1:3){
      print(scores_plot[[i]][[j]], vp = viewport(layout.pos.row = layout, layout.pos.col = j))
    }
  }
}
dev.off()


#######################################################################################################################
## loadings
#######################################################################################################################
# 
# # Making a dataframe to plot the loadigns in ggplot
# # adapted from: https://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
# set1.loadings <- data.frame(var=row.names(pc$loadings), pc$loadings[,1:5])
# mult <- min( (max(fcomp_pca$set1.pc2)-min(fcomp_pca$set1.pc2))/(max(set1.loadings$Comp.2)-min(set1.loadings$Comp.2)),
#              (max(fcomp_pca$set1.pc1)-min(fcomp_pca$set1.pc1))/(max(set1.loadings$Comp.1)-min(set1.loadings$Comp.1))
# )
# set1.loadings$load.PC1 <- 0.7*mult*set1.loadings$Comp.1
# set1.loadings$load.PC2 <- 0.7*mult*set1.loadings$Comp.2
# set1.loadings$var <- pfts
# 
# # --------------
# # Manually creating a biplots
# # Weights cutoff = 75th percentile for PC1 & 2 (top 25% of weights)
# # --------------
# set1.cutoff <- quantile(c(abs(set1.loadings$Comp.1), abs(set1.loadings$Comp.2)), 0.75)
# biplot.set1 <- ggplot() +
#   theme_bw() +# coord_equal() +
#   geom_hline(yintercept=0, size=0.3, color="gray50") + geom_vline(xintercept = 0, size=0.3, color="gray50") +
#   geom_point(data=fcomp_pca, aes(x=set1.pc1, y=set1.pc2), size=0.2, color="black") +
#   geom_segment(data=set1.loadings[abs(set1.loadings$Comp.1)>set1.cutoff | abs(set1.loadings$Comp.2)>set1.cutoff,], aes(x=0, y=0, xend=load.PC1, yend=load.PC2), 
#                arrow=arrow(length=unit(0.3,"cm")), size=1, alpha=0.9, color="red") +
#   geom_text(data=set1.loadings[abs(set1.loadings$Comp.1)>set1.cutoff | abs(set1.loadings$Comp.2)>set1.cutoff,], 
#             aes(x=load.PC1, y=load.PC2, label=var), size = 4, vjust=2, color="red", fontface="bold") +
#   scale_x_continuous(name=paste0("PC1 (", round(set1.POV["Comp.1"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$set1.pc1, c(0.001, 0.999))) +
#   scale_y_continuous(name=paste0("PC2 (", round(set1.POV["Comp.2"]*100,1),"%)"), expand=c(0.1,0.1), limits=quantile(paleon$set1.pc2, c(0.001, 0.999))) +
#   ggtitle(paste0("Region, Settlement-Era Fcomp, ", model)) +
#   theme(plot.title=element_text(hjust=0.5, face="bold", size=rel(1.25)))
# 
# print(biplot.set1)
