library(sp)
library(ggplot2)
library(grid)

# read in stability index table
stepps_si = read.csv('data/Stability_STEPPS.csv')
path.google <- "/Volumes/GoogleDrive/My Drive/PalEON_ecosystem-change_models-vs-data"

# read in stability model annual
fcomp_annual = read.csv('data/Stability_FCOMP_annual.csv')
fcomp_annual$lat = paleon[fcomp_annual$cell, 'lat'] 
fcomp_annual$lon = paleon[fcomp_annual$cell, 'lon'] 
write.csv(fcomp_annual, 'data/Stability_FCOMP_annual.csv')

# make some maps
us <- map_data("state")
models = as.vector(unique(fcomp_annual$model))

# annual
# png(file.path("figures", "Stability_FCOMP.png"), height=4, width=6, units="in", res=320)
pdf(file.path("figures", "Stability_FCOMP.pdf"), height=4, width=6)
for (model in models){
  pfts = unique(fcomp_annual[which(fcomp_annual$model == model),'pft'])
  for (pft in pfts){
    fcomp_annual_sub = fcomp_annual[which((fcomp_annual$model==model)&(fcomp_annual$pft== pft)),]
    
    fcomp.deriv <- ggplot(data=fcomp_annual_sub) +
      geom_tile(aes(x=lon, y=lat, fill=deriv.abs)) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
      coord_equal(xlim=range(fcomp_annual$lon), ylim=range(fcomp_annual$lat), expand=0) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(fcomp_annual$deriv.abs, na.rm=T)) +
      theme_bw() +
      ggtitle(paste0("Mean Absolute Rate of Change; ", model,": ", pft))
    
    fcomp.sig <- ggplot(data=fcomp_annual_sub) +
      geom_tile(aes(x=lon, y=lat, fill=n.yrs.sig)) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
      coord_equal(xlim=range(fcomp_annual$lon), ylim=range(fcomp_annual$lat), expand=0) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(fcomp_annual$n.yrs.sig, na.rm=T)) +
      theme_bw() +
      ggtitle(paste0("Number of Years Showing Change; ", model,": ", pft))
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    print(fcomp.deriv, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(fcomp.sig, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  }
}
dev.off()


# centennial
fcomp_100 = read.csv('data/Stability_FCOMP_centennial.csv')
pdf(file.path("figures", "Stability_FCOMP_100.pdf"), height=4, width=6)
for (model in models){
  pfts = unique(fcomp_100[which(fcomp_100$model == model),'pft'])
  for (pft in pfts){
    fcomp_sub = fcomp_100[which((fcomp_100$model==model)&(fcomp_100$pft== pft)),]
    
    fcomp.deriv <- ggplot(data=fcomp_sub) +
      geom_tile(aes(x=lon, y=lat, fill=deriv.abs)) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
      coord_equal(xlim=range(fcomp_100$lon), ylim=range(fcomp_100$lat), expand=0) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(fcomp_100$deriv.abs, na.rm=T)) +
      theme_bw() +
      ggtitle(paste0("Mean Absolute Rate of Change; ", model,": ", pft))
    
    fcomp.sig <- ggplot(data=fcomp_sub) +
      geom_tile(aes(x=lon, y=lat, fill=n.yrs.sig)) +
      geom_path(data=us,aes(x=long, y=lat, group=group), color="gray50") + 
      coord_equal(xlim=range(fcomp_100$lon), ylim=range(fcomp_100$lat), expand=0) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mean(fcomp_100$n.yrs.sig, na.rm=T)) +
      theme_bw() +
      ggtitle(paste0("Number of Years Showing Change; ", model,": ", pft))
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    print(fcomp.deriv, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(fcomp.sig, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  }
}
dev.off()