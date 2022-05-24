####ILandscape genetic analysis of PMI samples
#load env and metadata
pkg = c("tidyverse", "gridExtra", "RColorBrewer", "moments", "viridis" , "data.table", "sf", "geosphere", 'UpSetR', 'raster', 'igraph', 'adegenet', 'vegan')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#load metadata and clean

path='~/Projects/MOVE/anopheles_pmi_cnrfp/'
#read metadata and subset to useful stuff, read relatedness data and join metadata to each individual 'side'
metadata = fread('~/Projects/MOVE/anopheles_pmi_cnrfp/metadata/pmi_metadata.csv')
metadata = metadata[, c('sample_id','lat', 'long', 'Pop', 'angsdno', 'collection.year')]
res = fread('~/Projects/MOVE/anopheles_pmi_cnrfp/data/relatedness/AgamP4_3L.res')
res = left_join(res, metadata, by=c('a' = 'angsdno')) %>% left_join(.,metadata, by=c('b' = 'angsdno'))

#create comparison columns and new df with median relatedness
res$loc_comp = paste0(res$Pop.x, res$Pop.y)
median_relatedness_by_site_year = aggregate(rab ~ Pop.x + Pop.y + lat.x + long.x + lat.y + long.y+ collection.year.x + collection.year.y, data=res,FUN=median)
res$yearcomp = paste0(res$collection.year.x, '_', res$collection.year.y)
res$sameyearyesno = ifelse(res$collection.year.x == res$collection.year.y, 1, 0)

##-----------------------------------------------------------------------------------------------------------##
##                                            Test for isolation-by-distance                                 ##
##-----------------------------------------------------------------------------------------------------------##


#calculate distance between each point/sample
res$dist = distGeo(res[,c('lat.x','long.x')], res[,c('lat.y','long.y')])


res_2008 = res[res$collection.year.x == 2008 & res$collection.year.y ==2008]
res_2014 = res[res$collection.year.x == 2014 & res$collection.year.y ==2014]
res_2008 = res[res$collection.year.x == 2016 & res$collection.year.y ==2016]
res_2018 = res[res$collection.year.x == 2018 & res$collection.year.y ==2018]
res_2020 = res[res$collection.year.x == 2020 & res$collection.year.y ==2020]

newres = res[res$lat.x < 9.8 & lat.y < 9.8] #excluding outgroup sample

geogendist = ggplot(res[res$sameyearyesno == 1], aes(x=dist, y=rab))+
  geom_point()+
  theme_classic()+
  facet_wrap(~yearcomp)

#repeat this for each year (TURN INTO LAPPLY OR LOOP)

#test for isobd for each year of sampling
calc_isobd <- function(yr) {
  resdf = res
  resdf = resdf[resdf$collection.year.x == yr & collection.year.y == yr]
  rabmat = pivot_wider(resdf[,c('a', 'b', 'rab')], names_from = a, values_from = rab)
  rabmat = rabmat[,-1] #remove row ids
  rabmat = as.dist(rabmat, upper = FALSE) #take lower triangle to dist matrix object
  
  distmat = pivot_wider(resdf[,c('a', 'b', 'dist')], names_from = a, values_from = dist)
  distmat = distmat[,-1] #remove row ids
  distmat = as.dist(distmat, upper = FALSE) #take lower triangle to dist matrix object
  ibd = mantel.randtest(rabmat, distmat)
  return(ibd)
}

#run on all years
years = c(2008, 2014, 2016, 2018, 2020)
isos = lapply(years, calc_isobd)

#plot ibd plots
par(mfrow = c(2,3))
x = lapply(seq_along(1:length(isos)), function(i){plot(isos[[i]], main=years[[i]])})





#negative all around for isobd
#what happens if we exclude the 'outlier' samples?





#get Burkina Faso map
bf<-getData('GADM', country="GHA", level=1)
bf2<-getData('GADM', country="GHA", level=2)
bf3<-getData('GADM', country="GHA", level=3)

