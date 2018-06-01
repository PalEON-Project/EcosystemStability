mip_fcomp = readRDS('data/mip_fcomp.RDS')
mip2coarse = read.csv('data/mip2coarse.csv')

prob_sig <- function(x, prob){
  n_side = max(sum(x>0), sum(x<0))
  
  if (n_side/length(x) > prob){
    return(mean(x))
  } else {
    return(NA)
  }
}

load('data/input.rdata')
rIts = readRDS('data/rIts_small.RDS')
niter = dim(rIts)[4]

rIts_decid = apply(rIts[,which(mip2coarse$pft.coarse[match(taxa, mip2coarse$pft.short)]=='deciduous'),,], c(1,3,4), sum, na.rm=TRUE)
rIts_ever = apply(rIts[,which(mip2coarse$pft.coarse[match(taxa, mip2coarse$pft.short)]=='evergreen'),,], c(1,3,4), sum, na.rm=TRUE)

decid_deriv_all = apply(rIts_decid, c(1,3), function(x) diff(x, na.rm=TRUE)/100)
decid_deriv = apply(decid_deriv_all, c(1,2), prob_sig, prob=0.8)
decid_sig  = apply(decid_deriv, c(2), function(x) sum(!is.na(x)))
decid_deriv_abs = apply(decid_deriv, c(2), function(x) mean(abs(x), na.rm=TRUE))
decid_deriv_mean = apply(decid_deriv, c(2), function(x) mean(x, na.rm=TRUE))


ever_deriv_all = apply(rIts_ever, c(1,3), function(x) diff(x, na.rm=TRUE)/100)
ever_deriv = apply(ever_deriv_all, c(1,2), prob_sig, prob=0.8)
ever_sig  = apply(ever_deriv, c(2), function(x) sum(!is.na(x)))
ever_deriv_abs = apply(ever_deriv, c(2), function(x) mean(abs(x), na.rm=TRUE))
ever_deriv_mean = apply(ever_deriv, c(2), function(x) mean(x, na.rm=TRUE))


taxon_deriv_all = apply(rIts, c(1,2,4), function(x) diff(x, na.rm=TRUE)/100)
taxon_deriv = apply(taxon_deriv_all, c(1,2,3), prob_sig, prob=0.8)
taxon_sig  = apply(taxon_deriv, c(2,3), function(x) sum(!is.na(x)))
taxon_deriv_abs = apply(taxon_deriv, c(2,3), function(x) mean(abs(x), na.rm=TRUE))
taxon_deriv_mean = apply(taxon_deriv, c(2,3), function(x) mean(x, na.rm=TRUE))

library(reshape2)
taxon_si = melt(taxon_deriv_mean)
colnames(taxon_si) = c('cell', 'pft', 'deriv.mean')
taxon_si$pft = taxa[taxon_si$pft]
taxon_si$deriv.abs= melt(taxon_deriv_abs)[,3]
taxon_si$sig = melt(taxon_sig)[,3]

stepps_si = data.frame(x=rep(centers_veg[,1], times=2), 
                       y=rep(centers_veg[,2], times=2), 
                       cell=rep(seq(1,nrow(centers_veg)), 2),
                       deriv.mean=c(decid_deriv_mean, ever_deriv_mean), 
                       deriv.abs=c(decid_deriv_abs, ever_deriv_abs), 
                       pft=c(rep('Deciduous', nrow(centers_veg)), rep('Evergreen', nrow(centers_veg))),
                       sig = c(decid_sig, ever_sig))
stepps_si = rbind(stepps_si, data.frame(x=rep(centers_veg[,1], times=12),
                                        y=rep(centers_veg[,2], times=12),
                                        cell=seq(1,nrow(centers_veg)),
                                        deriv.mean = taxon_si$deriv.mean,
                                        deriv.abs = taxon_si$deriv.abs,
                                        pft = taxon_si$pft,
                                        sig = taxon_si$sig))

stepps_si$x = stepps_si$x*1e6 
stepps_si$y = stepps_si$y*1e6 

# convert coordinates to lat/long
coordinates(stepps_si) <- ~ x + y
proj4string(stepps_si) <- CRS('+init=epsg:3175')

stepps_ll <- spTransform(stepps_si, CRS('+proj=longlat +ellps=WGS84'))
stepps_ll <- data.frame(stepps_ll)

stepps_save = stepps_ll[,c('cell', 'y', 'x', 'pft', 'deriv.mean', 'deriv.abs', 'sig')]
colnames(stepps_save) = c('cell', 'lat', 'lon', 'pft', 'deriv.mean', 'deriv.abs', 'sig')

rIts_mean = apply(rIts, c(1,2,3), mean, na.rm=TRUE)
rIts_mean2 = apply(rIts_mean, c(1,2), mean, na.rm=TRUE)

rIts_decid_mean = apply(rIts_decid, c(1,2), mean, na.rm=TRUE)
rIts_decid_mean2 = apply(rIts_decid_mean, c(1), mean, na.rm=TRUE)

rIts_ever_mean = apply(rIts_ever, c(1,2), mean, na.rm=TRUE)
rIts_ever_mean2 = apply(rIts_ever_mean, c(1), mean, na.rm=TRUE)

rIts_mean_melt = melt(rIts_mean2)
colnames(rIts_mean_melt) = c('cell', 'pft', 'value')
rIts_mean_melt$pft = taxa[rIts_mean_melt$pft]
rIts_mean_melt = rbind(rIts_mean_melt,
                       data.frame(cell=seq(1,nrow(centers_veg)), pft=rep('Deciduous', nrow(centers_veg)), value=rIts_decid_mean2),
                       data.frame(cell=seq(1,nrow(centers_veg)), pft=rep('Evergreen', nrow(centers_veg)), value=rIts_decid_mean2))

stepps_si = merge(stepps_save, rIts_mean_melt, by=c("cell", "pft"))

# now save a csv file with columns: lon, lat, deriv.mean, deriv.abs, n.yrs, n.yrs.sig, fract.sig
write.csv(stepps_si, 'data/Stability_STEPPS.csv')


lims = list(xlims=c(min(paleon$x), max(paleon$x)), ylims=c(min(paleon$y), max(paleon$y)))

p <- ggplot(data=stepps_si[which(stepps_si$pft=='Deciduous'),]) + geom_point(aes(x=x, y=y, colour=abs(deriv.mean)), size=3) 
# p <- p + scale_colour_gradientn(colours=tim.colors(10))
p <- p + scale_colour_gradientn(colours=brewer.pal(8,"YlOrRd"), name="Stability index") + coord_fixed()
p <- add_map_albers(p, us.shp, lims)
p <- p + theme_bw() + theme(axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank())
print(p)
ggsave('figures/maps_deriv_si_coarsePFT.pdf')


p <- ggplot(data=stepps_si[which(stepps_si$pft != 'Evergreen'),]) + geom_point(aes(x=x, y=y, colour=abs(deriv.mean)), size=0.1) 
# p <- p + scale_colour_gradientn(colours=tim.colors(10))
p <- p + scale_colour_gradientn(colours=brewer.pal(8,"YlOrRd"), name="Stability index") + coord_fixed()
p <- add_map_albers(p, us.shp, lims)
p <- p + facet_wrap(~pft) + theme_bw() + theme(axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank())
print(p)
ggsave('figures/maps_deriv_si_taxon.pdf')


# now do this for all of the models
# for models, use gams
# then do this at annual and centennial
# calc.stability <- function(x, width){
#   dat.tmp <- data.frame(Y=x, Year=1:length(x))
#   k.use=round(length(x)/width, 0)
#   mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
#   mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp, vars="Year")
#   return(mod.deriv)
# }

calc.stability <- function(x, width, res){
  dat.tmp <- data.frame(Y=x, Year=seq(1, length(x)) )
  if (res=="cent"){
    idx = seq(1, length(x), by=100) 
    dat.tmp.pred <- data.frame(Y=x[idx], Year=seq(1, length(x))[idx] )
  } else {
    dat.tmp.pred = dat.tmp
  }
  k.use=round(length(x)/width, 0)
  mod.gam <- gam(Y ~ s(Year, k=k.use), data=dat.tmp)
  mod.deriv <- calc.derivs(mod.gam, newdata=dat.tmp.pred, vars="Year")
  return(mod.deriv)
}


library(mgcv)

paleon = read.csv('data/paleon_site_info.csv')

mip_fcomp$x = paleon[match(mip_fcomp$cell_id, paleon$num), 'x']
mip_fcomp$y = paleon[match(mip_fcomp$cell_id, paleon$num), 'y']

# only want 850-1850
mip_fcomp = mip_fcomp[which((mip_fcomp$year>=850)&(mip_fcomp$year<=1850)),]

models = unique(mip_fcomp$model)

# annual-scale stability indices
fcomp_si = data.frame(model=character(0), 
                      cell=numeric(0), 
                      x=numeric(0), 
                      y=numeric(0), 
                      pft=character(0), 
                      year = numeric(0),
                      deriv=numeric(0),
                      sig = character(0))
for (model in models){
  print(model)
  model_idx = which(mip_fcomp$model == model)
  pfts = unique(mip_fcomp[model_idx, 'pft'])
  cells = unique(mip_fcomp[model_idx, 'cell_id'])
  N_cells = length(cells)
  for (pft in pfts){
    print(pft)
    
    x_pft=mip_fcomp[which((mip_fcomp$model==model)&(mip_fcomp$pft==pft)), ]
    if (all(x_pft$value==0)){
      print("PFT does not occur in domain! Skipping.") 
      next
    }
    
    for (n in 1:N_cells){
      if ((n %% 50) == 0){
        print(paste0('Cell ', n, ' of ', N_cells))
      }
        
      cell=cells[n]
      x=x_pft[which((x_pft$model==model)&(x_pft$cell_id==cell)&(x_pft$pft==pft)), ]
      
      if (all(x$value==0)){
        print("PFT does not occur in grid cell! Skipping.") 
        next
      }
      
      stability = calc.stability(x$value, width=25, res="annual")
      nyears = nrow(stability)
      
      fcomp_si = rbind(fcomp_si, data.frame(model=rep(model, nyears), 
                                            cell=rep(cell, nyears), 
                                            x=rep(x$x[1], nyears), 
                                            y=rep(x$y[1], nyears), 
                                            pft=rep(pft, nyears), 
                                            year = stability$Year,
                                            deriv = stability$mean,
                                            sig = stability$sig))
    }
  }
}

fcomp_si$sig.bool = !is.na(fcomp_si$sig)
write.csv(fcomp_si, 'data/Stability_FCOMP_annual_by_year.csv')

deriv.abs  = aggregate(deriv ~ cell + pft + model, fcomp_si, function(x) mean(abs(x), na.rm=TRUE))
deriv.mean = aggregate(deriv ~ cell + pft + model, fcomp_si, function(x) mean(x, na.rm=TRUE))
n.yrs.sig  = aggregate(sig.bool ~ cell + pft + model, fcomp_si, function(x) sum(x))
n.yrs      = aggregate(year ~ cell + pft + model, fcomp_si, function(x) length(x))
fract.sig  = n.yrs.sig$sig.bool/n.yrs$year

fcomp_si_summary = data.frame(deriv.abs[,1:3], 
                              deriv.abs = deriv.abs[,4], 
                              deriv.mean = deriv.mean[,4], 
                              n.yrs = n.yrs$year, 
                              n.yrs.sig = n.yrs.sig$sig.bool, 
                              fract.sig = fract.sig)

mip_fcomp_mean = aggregate(value~model+cell_id+pft, mip_fcomp, mean, na.rm=TRUE)
colnames(mip_fcomp_mean) = c('model', 'cell', 'pft', 'value')

# library(reshape2)
# mip_fcomp_mean_cast = dcast(mip_fcomp_mean, model+cell_id~pft)

fcomp_si_summary = merge(fcomp_si_summary, mip_fcomp_mean)

write.csv(fcomp_si_summary, 'data/Stability_FCOMP_annual.csv')


# centennial-scale stability indices
fcomp_si_100 = data.frame(model=character(0), 
                      cell=numeric(0), 
                      x=numeric(0), 
                      y=numeric(0), 
                      pft=character(0), 
                      year = numeric(0),
                      deriv=numeric(0),
                      sig = character(0))
for (model in models){
  print(model)
  model_idx = which(mip_fcomp$model == model)
  pfts = unique(mip_fcomp[model_idx, 'pft'])
  cells = unique(mip_fcomp[model_idx, 'cell_id'])
  N_cells = length(cells)
  for (pft in pfts){
    print(pft)
    
    x_pft=mip_fcomp[which((mip_fcomp$model==model)&(mip_fcomp$pft==pft)), ]
    if (all(x_pft$value==0)){
      print("PFT does not occur in domain! Skipping.") 
      next
    }
    
    for (n in 1:N_cells){
      print(paste0('Cell ', n, ' of ', N_cells))
      cell=cells[n]
      
      x=x_pft[which((x_pft$model==model)&(x_pft$cell_id==cell)&(x_pft$pft==pft)), ]
      
      if (all(x$value==0)){
        print("PFT does not occur in grid cell! Skipping.") 
        next
      }
      
      stability = calc.stability(x$value, width=100, res="cent")
      nyears = nrow(stability)
      
      fcomp_si_100 = rbind(fcomp_si_100, data.frame(model=rep(model, nyears), 
                                            cell=rep(cell, nyears), 
                                            x=rep(x$x[1], nyears), 
                                            y=rep(x$y[1], nyears), 
                                            pft=rep(pft, nyears), 
                                            year = stability$Year,
                                            deriv = stability$mean,
                                            sig = stability$sig))
    }
  }
}

fcomp_si_100$sig.bool = !is.na(fcomp_si_100$sig)
write.csv(fcomp_si_100, 'data/Stability_FCOMP_centennial_by_year.csv')

deriv.abs  = aggregate(deriv ~ cell + pft + model, fcomp_si_100, function(x) mean(abs(x), na.rm=TRUE))
deriv.mean = aggregate(deriv ~ cell + pft + model, fcomp_si_100, function(x) mean(x, na.rm=TRUE))
n.yrs.sig  = aggregate(sig.bool ~ cell + pft + model, fcomp_si_100, function(x) sum(x))
n.yrs      = aggregate(year ~ cell + pft + model, fcomp_si_100, function(x) length(x))
fract.sig  = n.yrs.sig$sig.bool/n.yrs$year

fcomp_si_summary = data.frame(deriv.abs[,1:3], 
                              deriv.abs = deriv.abs[,4], 
                              deriv.mean = deriv.mean[,4], 
                              n.yrs = n.yrs$year, 
                              n.yrs.sig = n.yrs.sig$sig.bool, 
                              fract.sig = fract.sig)

fcomp_si_summary$lat = paleon[fcomp_si_summary$cell, 'lat'] 
fcomp_si_summary$lon = paleon[fcomp_si_summary$cell, 'lon'] 

mip_fcomp_mean = aggregate(value~model+cell_id+pft, mip_fcomp, mean, na.rm=TRUE)
colnames(mip_fcomp_mean) = c('model', 'cell', 'pft', 'value')

# library(reshape2)
# mip_fcomp_mean_cast = dcast(mip_fcomp_mean, model+cell_id~pft)

fcomp_si_summary = merge(fcomp_si_summary, mip_fcomp_mean)

write.csv(fcomp_si_summary, 'data/Stability_FCOMP_centennial.csv')

## to do:
## add lat and long to fcomp_si
## add centennial fcomp analysis
## convert stepps output to lat long

## add deciduous and conifer to fcomp analysis
## add LINKAGES
## stepps only for overlapping interval as well as full