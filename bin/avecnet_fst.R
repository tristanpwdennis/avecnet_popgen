library(tidyverse)
library(data.table)

#read fst
fst_avecnet <- fread('/Users/dennistpw/Projects/anopheles_pmi_cnrfp/fst/treated_untreated.win10step1')
colnames(fst_avecnet) <- c('region','chrom','midpos','nsites','fst')
annotation <- fread('')
ann <- fread('~/Projects/td_je_angam_2022/data/annotation/VectorBase-55_AgambiaePEST.gff', col.names = c('chrom','db','type','ann_start','ann_end','pt','strand','pt2','name'))

#ir genes amirite
name <- c('kdr','rdl','coeae','ace1','cyp6','gst','cyp6m2','cyp9k1')
chrom <- c('AgamP4_2L','AgamP4_2L','AgamP4_2L','AgamP4_2R','AgamP4_2R','AgamP4_3R','AgamP4_3R','AgamP4_X')
start <- c(2358158,25363652,28548433,3484107,28463000,28580000,6900000,15222000)
end <- c(2431617,25434556,28550748,3495790,28568000,28605000,7030000,15257000)
irdf <- data.frame(name, chrom, start, end)
irdf$chrom <- factor(irdf$chrom)
irdf$meanpos <- (irdf$start + irdf$end) / 2

ggplot(fst_avecnet, aes(x=midpos, y=fst,colour=chrom))+
  geom_point()+
  facet_wrap(~chrom, nrow = 1, scales = "free_x")+
  geom_vline(data=irdf, aes(xintercept=meanpos))+
  theme_classic()



fst_avecnet$start = fst_avecnet$midpos - 5000
fst_avecnet$end = fst_avecnet$midpos + 5000

x = fst_avecnet %>% filter(fst > 0.001)

intersected_fst = tidygenomics::genome_intersect(x, ann, by = c('chrom', 'start' = 'ann_start', 'end' = 'ann_end'))

protein_coding_gene_fst = intersected_fst %>% filter(type == 'protein_coding_gene')


thetas_treated_10 <- fread('/Users/dennistpw/Projects/avecnet_popgen/data/thetas/treated_avecnet_samples.list.wg.nomaffilter.10kb.thetasWindow.gz.pestPG')
thetas_treated_20 <- fread('/Users/dennistpw/Projects/avecnet_popgen/data/thetas/treated_avecnet_samples.list.wg.nomaffilter.20kb.thetasWindow.gz.pestPG')

thetas_untreated_10 <- fread('/Users/dennistpw/Projects/avecnet_popgen/data/thetas/untreated_avecnet_samples.list.wg.nomaffilter.10kb.thetasWindow.gz.pestPG')
thetas_untreated_20 <- fread('/Users/dennistpw/Projects/avecnet_popgen/data/thetas/untreated_avecnet_samples.list.wg.nomaffilter.20kb.thetasWindow.gz.pestPG')


piuntx10 = ggplot(thetas_untreated_10, aes(x=WinCenter,y=tP/nSites))+
  geom_line(colour="dodgerblue2", alpha=0.5)+
  ylim(0,0.4)+
  facet_wrap(~Chr, nrow = 1, scales = "free_x")+
  theme_tufte()

pitx10 = ggplot(thetas_treated_10, aes(x=WinCenter,y=tP/nSites))+
  geom_line(colour="dodgerblue2", alpha=0.5)+
  ylim(0,0.4)+
  facet_wrap(~Chr, nrow = 1, scales = "free_x")+
  theme_tufte()

plot_grid(piuntx10, pitx10, ncol=1)



tduntx10 = ggplot(thetas_untreated_10, aes(x=WinCenter,y=Tajima))+
  geom_line(colour="dodgerblue2", alpha=0.5)+
  facet_wrap(~Chr, nrow = 1, scales = "free_x")+
  theme_tufte()

tdtx10 = ggplot(thetas_treated_10, aes(x=WinCenter,y=Tajima))+
  geom_line(colour="dodgerblue2", alpha=0.5)+
  facet_wrap(~Chr, nrow = 1, scales = "free_x")+
  theme_tufte()

plot_grid(tduntx10, tdtx10, ncol=1)
