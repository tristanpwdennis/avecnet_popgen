#cnrfp pca, admixture, selection, covariance matrix analysis

library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
#metadata
metadata <- fread("~/Projects/anopheles_cnrfp/sequenced_metadata_avecnet.csv")
metadata$yearmo = paste0(metadata$year, '_', metadata$month)
metadata$yearmo <- factor(metadata$yearmo, levels = c("2014_7","2014_8","2014_9","2014_11","2015_7","2015_9","2015_11"))
#all samples analysis
allsample_cov <- fread('~/Projects/anopheles_cnrfp/pca/pcangsd.cov')
allsample_admixtue <- fread('~/Projects/anopheles_cnrfp/pca/pcangsd.admix.2.Q')
allsample_inbreeding <- fread('~/Projects/anopheles_cnrfp/pca/pcangsd.inbreed.samples')
allsample_selection <- fread('~/Projects/anopheles_cnrfp/pca/pcangsd.inbreed.samples')


allsample_pca <- eigen(allsample_cov)
allsample_pcavecs <- data.frame(allsample_pca$vectors)
allsample_pcavecs <- cbind(metadata, allsample_pcavecs)


#by treatment
pca = ggplot(allsample_pcavecs, aes(x=X1,y=X2, colour=as.factor(treated)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC1",y="PC2", colour="Treated")

pcb = ggplot(allsample_pcavecs, aes(x=X3,y=X4, colour = as.factor(treated)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC3",y="PC4", colour="Treated")

pcc = ggplot(allsample_pcavecs, aes(x=X5,y=X6, colour=as.factor(treated)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC5",y="PC6", colour="Treated")

plot_grid(pca, pcb, pcc)

#by cluster
pcd = ggplot(allsample_pcavecs, aes(x=X1,y=X2, colour=as.factor(Cluster)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC1",y="PC2", colour="Cluster")

pce = ggplot(allsample_pcavecs, aes(x=X3,y=X4, colour = as.factor(Cluster)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC3",y="PC4", colour="Cluster")

pcf = ggplot(allsample_pcavecs, aes(x=X5,y=X6, colour=as.factor(Cluster)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC5",y="PC6", colour="Cluster")

plot_grid(pcd, pce, pcf)

#by yearmonth
pcg = ggplot(allsample_pcavecs, aes(x=X1,y=X2, colour=as.factor(yearmo)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC1",y="PC2", colour="Year-Month")

pch = ggplot(allsample_pcavecs, aes(x=X3,y=X4, colour = as.factor(yearmo)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC3",y="PC4", colour="Year-Month")

pci = ggplot(allsample_pcavecs, aes(x=X5,y=X6, colour=as.factor(yearmo)))+
  geom_point()+
  theme_tufte()+
  labs(x="PC5",y="PC6", colour="Year-Month")

plot_grid(pcg, pch, pci)


###inbreeding
colnames(allsample_inbreeding) <- "Fis"
inbreeddat <- cbind(metadata, allsample_inbreeding)

ggplot(inbreeddat, aes(y=Fis, x=as.factor(treated)))+
  geom_violin()+
  geom_boxplot(width=0.5)+
  theme_tufte()+
  labs(x="Treated")

ggplot(inbreeddat, aes(y=Fis, x=yearmo, colour=as.factor(treated)))+
  geom_boxplot()+
  theme_tufte()+
  labs(x="YearMo", colour = "Treated")+
  scale_color_brewer(palette = "Dark2")



