library(sp)
library(ggplot2)
library(grid)

# read in stability index table
stepps_si = read.csv('data/Stability_STEPPS.csv')
path.google <- "~/Google Drive/PalEON_ecosystem-change_models-vs-data"

# convert coordinates to lat/long
coordinates(stepps_si) <- ~ x + y
proj4string(stepps_si) <- CRS('+init=epsg:3175')

stepps_ll <- spTransform(stepps_si, CRS('+proj=longlat +ellps=WGS84'))
stepps_ll <- data.frame(stepps_ll)


# read in stability model annual
fcomp_annual = read.csv('data/Stability_FCOMP_annual.csv')
fcomp_annual$lat = paleon[fcomp_annual$cell, 'lat'] 
fcomp_annual$lon = paleon[fcomp_annual$cell, 'lon'] 

# make some maps
us <- map_data("state")
models = as.vector(unique(fcomp_annual$model))

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