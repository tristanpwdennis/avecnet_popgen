pkg = c("tidyverse", "rnaturalearth", "rnaturalearthdata", "sf", "cowplot", "scatterpie")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#read all our data in
res = read.delim('~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/wg_filts_glf3/testres.res', header = F)
bamlist = read.csv('~/Projects/MOVE/data/anopheles_03_21_DW/angsd-output/bamlist.txt', header =F)
bamlist$no = seq(0, (length(bamlist$V1) - 1))
colnames(bamlist) = c('sample_id', 'bamname', 'samplenum')
metadata = read.csv('~/Projects//MOVE/data/anopheles_03_21_DW/metadata/location_ecozone_metadata.csv')
oldmeta = read.csv('~/Projects//MOVE/td_je_angam_2022/metadata/metadata_with_insecticide_info.csv')
newmeta = left_join(metadata, bamlist, by=c('sample' = 'sample_id'))

metastrip = metadata %>% select(Lat, Long,Site)
fullstrip = oldmeta %>% select(-Lat, -long)



s = left_join(fullstrip, metastrip, by = 'Site') %>% unique()

#download earth data from rnaturalrearth
world <- ne_countries(scale = 'large', returnclass = "sf", continent = 'africa')
lakes <- ne_download(scale = 'large', type = 'lakes', category = 'physical', returnclass = "sf")
rivers <- ne_download(scale = 'large', type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
ocean <- ne_download(scale = 'large', type = 'ocean', category = 'physical', returnclass = "sf")
metadata
sites = metadata %>% select(Site, Lat, Long, Ecozone) %>% unique() #%>% full_join(colesites)
habitatocolours = c('#e41a1c','#377eb8','#4daf4a','#984ea3')
darkerhabitatcolours3 = colorspace::darken(habitatocolours, amount=0.3)
iseccolours = c('#a6cee3','#1f78b4','#b2df8a','#33a02c')
darkeriseccolours = colorspace::darken(iseccolours, amount=0.3)
#isecusage = metadata %>% group_by(Form, Site, Lat, Long, INSECTICIDE_USE, habitat) %>% unique() 
#isecusage
#plot eas africa
maps = ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  geom_point(data=s, aes(y=Lat,x=Long, color=habitat, fill=habitat, shape=INSECTICIDE_USE), size = 6)+
  scale_color_manual(values=darkerhabitatcolours3,guide='none')+
  scale_fill_manual(values=habitatocolours)+
  scale_shape_manual(values=c(21, 22, 23, 24))+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-3.6, 1.5), ylim = c(4.5, 11.5), expand = FALSE)+
  labs(x='Longitude', y='Latitude', color='Ecozone',fill = 'Ecozone', shape='Insecticide Usage Pattern')+
  guides(fill=guide_legend(override.aes=list(colour=habitatocolours)))+
  theme_classic()+
  facet_wrap(~Form)
maps

ggsave('~/OneDrive - University of Glasgow/MOVE/data/dw_anopheles_03_21/maps.pdf', width = 11, height = 8, units = 'in', device = 'pdf')

#load preseq metadata and join to gps data for ghana samples
ghana_library_metadata = fread('~/OneDrive - University of Glasgow/MOVE/samples/pmi/complete_preseq_metadata.csv')
colesites = read.csv('~/OneDrive - University of Glasgow/MOVE/samples/pmi/ghana_site_gps.csv')
ghana_library_metadata = ghana_library_metadata %>% left_join(colesites, by = c('Pop'= 'site'))
gsites = unique(ghana_library_metadata[`sine PCR result` == 'gambiae'][,c('Pop', 'lat', 'long')])
gsites$schedule = 'PMI - IRS'
anet= fread('~/OneDrive - University of Glasgow/MOVE/samples/cnrfp/avecnet_db_updated.csv')
anet_locs = fread('~/OneDrive - University of Glasgow/MOVE/samples/cnrfp/cluster_coords.csv')

anet_clusters = c(1,2,3,4,9,10,11,12,18,20,26,29,36,37,39)
bsites = anet_locs[anet_locs$cluster %in% anet_clusters]
colnames(bsites) = c('Pop', 'lat', 'long')
bsites$schedule = 'AvecNet'
dsites = sites[c('Site', 'Lat', 'Long')]
colnames(dsites) = c('Pop', 'lat', 'long')
dsites$schedule = 'LSTM'
allsites = rbind(bsites, dsites, gsites)
dsites
irs_colours = c('#e41a1c', '#377eb8', '#4daf4a')
darker_irs = colorspace::darken(irs_colours, amount=0.4)
  
  
ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  geom_point(data=allsites, aes(y=lat,x=long, colour = schedule, fill=schedule), size = 3, shape = 23)+
  scale_color_manual(values=darker_irs)+
  scale_fill_manual(values=irs_colours)+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-5, 1.5), ylim = c(4.5, 11.5), expand = FALSE)+
  labs(x='Latitude', y='Longitude', fill = 'Programme')+guides(colour = "none")+
  theme(text = element_text(size = 15))

ggplot()+
  geom_sf(data=world, color = '#dfc27d', fill='#f6e8c3')+
  geom_sf(data = rivers, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = lakes, colour = '#4575b4', fill = '#a6bddb')+
  geom_sf(data = ocean, colour = '#4575b4', fill = '#a6bddb')+
  geom_point(data=metadata, aes(y=Lat,x=long, colour = INSECTICIDE_USE, fill=INSECTICIDE_USE), size = 3, shape = 23)+
  scale_color_manual(values=darker_irs)+
  scale_fill_manual(values=irs_colours)+
  ggspatial::annotation_scale(location = "br", width_hint = 0.4) +
  coord_sf(xlim = c(-5, 1.5), ylim = c(4.5, 11.5), expand = FALSE)+
  labs(x='Latitude', y='Longitude', fill = 'Programme')+guides(colour = "none")+
  theme(text = element_text(size = 15))


pietwo = metadata %>% group_by(Form, Site) %>% count()
pieone


ggplot()+
  geom_scatterpie(aes(x=Lat, y=long, cols = INSECTICIDE_USE), data = metadata)

