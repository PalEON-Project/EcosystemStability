
load('data/refab.biomass.CI13.Rdata')
load('data/refab.biomass.x.meta.w.settle.Rdata')
names(biomassCI) <- unique(x.meta$site.id)[1:182]

refab.diff <- lapply(biomassCI,function(x) diff(rev(x[2,1:11])/100))

refab.signif <- lapply(refab.diff,function(x) length(which(x > .01)))
refab.signif1 <- lapply(refab.diff,function(x) length(which(x <  -.01)))

n.signif = as.vector(unlist(refab.signif) + unlist(refab.signif1))


refab.diff.mean <- lapply(biomassCI,function(x) mean(diff(rev(x[2,1:11])/100)))
x <- unlist(refab.diff.mean)

lat <- long <- list()
for(i in 1:182){
  lat[[i]] <- x.meta[x.meta$site.id==names(biomassCI)[i],'lat'][1]
  long[[i]] <- x.meta[x.meta$site.id==names(biomassCI)[i],'long'][1]
}

refab.mean.slope=data.frame(lat = unlist(lat,use.names = F), long = unlist(long,use.names = F),
                            refab.mean.slope = unlist(refab.diff.mean,use.names=F),
                            n.signif = n.signif)

write.csv(refab.mean.slope,file = 'refab.mean.slope.csv')
