######################################################################
#identifying wolbachia and plasmodium reads in anopheles samples#
######################################################################
#load our shit
pkg = c('tidyverse', 'data.table')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

setwd('~/Projects/anopheles_pmi/analysis/read_mapping_data_14_04_22/')

idxfiles = list.files('.', pattern = 'idxstats')
idxdata = lapply(seq_along(idxfiles), function(i){fread(idxfiles[[i]]) %>% mutate(fn=idxfiles[[i]])}) #fread all, adding filename as col
idxdata = do.call(rbind, idxdata)
idxdata$sample = gsub('_220322_L001.srt.dp.bam.idxstats', '', idxdata$fn)
unique(idxdata$V1)

rbind(idxdata[grepl('Pf3', idxdata$V1)],idxdata[grepl('NZ_CP050531_Wolbachia_pipientis', idxdata$V1)])

hist(idxdata[grepl('NZ_CP050531_Wolbachia_pipientis', idxdata$V1)]$V3)
