library(tidyverse)
library(data.table)

metadata <- fread('Projects/avecnet_popgen/metadata/sequenced_metadata_avecnet.csv')

#read sample names
files <- list.files('~/Projects/avecnet_popgen/data/roh/', pattern = 'summary')
files <- gsub('_EKD.*','',files)

#read rohlem
av_rohlen <- fread('Projects/avecnet_popgen/data/roh/avglenfix.txt', col.names = c('MeanRohLen','CIL_RohLen','CIU_RohLen'))
#read & in ROH
fracroh <- fread('Projects/avecnet_popgen/data/roh/fracrohfix.txt', col.names = c('FracRoh','CIL_FracRoh','CIU_FracRoh'))

#bind
rohstats <- cbind(files, av_rohlen, fracroh)

ggplot(rohstats, aes(x=MeanRohLen,y=FracRoh))+
  geom_point()

hist(rohstats$FracRoh)

rohstats <- left_join(rohstats, metadata,by =c('files' = 'sample_name'))
ggplot(rohstats, aes(x=as.factor(year),y=FracRoh, colour=as.factor(treated)))+
  geom_boxplot()
