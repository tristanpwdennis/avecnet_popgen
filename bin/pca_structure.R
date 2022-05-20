library(data)
#ls pkg
pkg = c("ggtree", "ape", "BiocManager", 'phangorn', 'phytools', 'data.table')
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#
#BiocManager::install('ggtree')

#devtools::install_version("dplyr", version = "1.0.5")
#remove.packages('ggtree')
#BiocManager::install("YuLab-SMU/treedataverse")
#import metadata, join to bamlist order
setwd('~/Projects/anopheles_pmi/analysis/')

metadata = read.csv('metadata/allsamples_qcpass_metadata.csv')

cov = fread('data/pca/AgamP4_3L.pcangsd.cov')
pca = eigen(cov)

s = data.frame(pca$vectors)
metadata$

s = cbind(s, metadata$sample_id)
names(s)[names(s) == 'metadata$seq_id'] <- 'seq_id'
metadata$collection.year
#plot by treatment
plot(pca$vectors[,2:3],col=as.factor(metadata$collection.year))
legend("topleft",fill=1:length(levels(as.factor(metadata$collection.year))),levels(as.factor(metadata$collection.year)))

#plot by pop
plot(pca$vectors[,3:2],col=as.factor(metadata$Pop))
legend("topright",fill=1:9,levels(as.factor(metadata$Pop)))  

#plot by pop
plot(pca$vectors[,1:2],col=as.factor(metadata$Pop))
legend("topright",fill=1:9,levels(as.factor(metadata$Pop)))  
