#-----------------------------------------------------------#
#avecnet genomics analysis
#-----------------------------------------------------------#
#setup env
#load and.or install packages we need
pkg = c("tidyverse", "sjPlot", "cowplot", "data.table", "DHARMa", "lme4", "MASS", "ggthemes")
#install.packages(pkg) #install packages if you need them and load
new.packages <- pkg[!(pkg %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(pkg, require, character.only = TRUE)#

#set cowplot theme
theme_set(theme_cowplot(font_size=17, font_family = "serif") + 
            theme(text = element_text(colour = "black")))

#-----------------------------------------------------------#
#load and wrangle data
#-----------------------------------------------------------#
#metadata
metadata <- fread('~/Projects/avecnet_popgen/metadata/sequenced_metadata_avecnet.csv')
#load fsts
afst_df =  fread("~/Projects/avecnet_popgen/data/fst/fst_xyearmo_xcluster.txt")
fst_df <- fread('/Users/dennistpw/Library/Mobile\ Documents/com~apple~CloudDocs/avecnet\ paper/fst_bysiteby.txt')
morefst = cbind(fst_df, afst_df[,4:10])
hist(morefst$V1)
colnames(morefst) <- c('unweighted_fst','weighted_fst',colnames(afst_df[,4:10]))

hist(fst_df$fst_weighted, breaks=50)


#between each group and every other group
fst =morefst %>% group_by(month_a, cluster_a, year_a, treated) %>% summarise(mean_fst = mean(weighted_fst))

#make timepoints and cap negative fst at zero
fst$fst = pmax(fst$mean_fst,0)
fst$timep = paste0(fst$year_a,'_',fst$month_a)
fst$timep <- factor(fst$timep, levels = c("2014_7","2014_8","2014_9","2014_11","2015_7","2015_9","2015_11"))

#read thetas
thetas=fread('~/Projects/avecnet_popgen/data/thetas/avecnet_thetas.csv')

#inbreeding coefficients
inbreedin_coeffs = fread('~/Projects/avecnet_popgen/data/fis/cnrfp_3l.34856_seed.final', col.names = c('fis'))

#bind iunbreeding coeffs to data
metadata_seq_fis = cbind(metadata, inbreedin_coeffs)

#-----------------------------------------------------------#
#models of difference in fst, pi and fis by year and tx
#-----------------------------------------------------------#
#model (insignificant)

#uncapped fsts make gaussian distribution (though for allele freq)
m1 = glm(formula = mean_fst ~ treated * as.factor(year_a), data = fst, family="gaussian")
m1simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
residuals(m1simulationOutput, quantileFunction = qnorm, outlierValues = c(0,1))
plot(m1simulationOutput)
drop1(m1)
summary(m1)
#model for pi
#pi also follows gaussian
m2 = glm(formula = log(pi) ~ treated * as.factor(year), data = thetas,family="gaussian")
summary(m2)
m2simulationOutput <- simulateResiduals(fittedModel = m2, plot = F)
residuals(m2simulationOutput, quantileFunction = qnorm, outlierValues = c(0,1))
plot(m2simulationOutput)
plot(m2)
#fis - gaussian
m3 = glm(formula = fis ~ treated * as.factor(year), data = metadata_seq_fis, family="gaussian")
summary(m3)
m3simulationOutput <- simulateResiduals(fittedModel = m3, plot = F)
residuals(m3simulationOutput, quantileFunction = qnorm, outlierValues = c(0,1))
plot(m3simulationOutput)
drop1(m3)
#hist(metadata_seq_fis$fis_ngsf)

#-----------------------------------------------------------#
#plots of fst, fis, pi and pca by time and treatment
#-----------------------------------------------------------#

#dual colour palette yellow and blue
pal = c('#f2be54','#153e5c')

#pca
x =  fread('~/Projects/avecnet_popgen/data/pca/pcangsd.cov')
pca = eigen(x)
pcaframe = cbind(data.frame(pca$vectors), metadata_seq_fis)

#add pc variance explained
varex = pca$values/sum(pca$values)
barplot(varex[1:20])

pc12plot = ggplot(pcaframe, aes(x=X1, y=X2, colour=as.factor(treated)))+
  scale_color_manual(values = pal)+
  labs(x="PC1 (0.7%)", y="PC2 (0.1%)", colour = 'Treated')+
  geom_point(alpha=0.7)+
  theme_tufte()+
  theme(text = element_text(size = 17)) 

#fsts
fstplt = ggplot(fst, aes(x=as.factor(year_a), y=mean_fst, colour=as.factor(treated)))+
  geom_jitter(position = position_jitterdodge(),alpha=0.7)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  labs(x='', y='Fst')+
  theme_tufte()+
  scale_color_manual(values = pal)+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

#theta treated vs untreated all
piplt = ggplot(thetas, aes(x=as.factor(year), y=pi,color=as.factor(treated)))+
  geom_jitter(position = position_jitterdodge(),alpha=0.7)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  labs(x='', y='π')+
  theme_tufte()+
  scale_color_manual(values = pal)+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

#inbreeding coeffs
fisplt = ggplot(metadata_seq_fis, aes(x=as.factor(year), y=fis,color=as.factor(treated)))+
  geom_jitter(position = position_jitterdodge(),alpha=0.7)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  labs(x='', y='Fis')+
  theme_tufte()+
  scale_color_manual(values = pal)+
  theme(text = element_text(size = 17)) +
  theme(legend.position = "none")

pcalegend = get_legend(pc12plot) #extract legend

#plot all together
divplots = cowplot::plot_grid(fstplt, 
                              piplt, 
                              fisplt, 
                              pcalegend,
                              labels=c('B','C','D'))


countplots = cowplot::plot_grid(pc12plot + theme(legend.position = "none"),
                                divplots,
                                labels = 'A')

#-----------------------------------------------------------#
#do population connectivity and size show a relationship? include this or not? wait till mv sends data, do analysis, then report to her
#-----------------------------------------------------------#
#get monthly abu,ndance data
monthlyabundance <- fread('~/Projects/anopheles_pmi_cnrfp/MonthyData_TD.csv', col.names = c('v1','cluster','year','month','counts'))
#wrangle everything into df where  mean fis, pi and fst are in df by treatment, cluster, site

#subset thetas and fst
thetasub = thetas[,c('month','year','timep','treated','cluster','pi')]
fstsub = fst[,c('month_a','year_a','timep','treated','cluster_a','mean_fst')]
colnames(fstsub) <- c('month','year','timep','treated','cluster','mean_fst')

fissum = aggregate(fis ~cluster+ month + year + treated, data = metadata_seq_fis, FUN = mean)
fissum$timep = as.character(paste0(fissum$year, '_', fissum$month))
colnames(fissum) <- c('cluster','month','year','treated','mean_fis','timep')

f = full_join(thetasub, fstsub)
g = left_join(f, monthlyabundance)
mean_divs = full_join(g, fissum)

#plot connectivity vs counts
fstcounts = ggplot(data=mean_divs, aes(x=mean_fst,y=counts, colour=as.factor(treated)))+
  scale_color_manual(values = pal)+
  geom_point(size=2, alpha=0.7)+
  theme_tufte()+
  labs(x="Fst",y='Counts', colour='Treated')+
  theme(text = element_text(size = 15),
        legend.position = "none")

#plot inbreeding vs counts
fiscounts = ggplot(data=mean_divs, aes(x=mean_fis,y=counts, colour=as.factor(treated)))+
  scale_color_manual(values = pal)+
  geom_point(size=2, alpha=0.7)+
  theme_tufte()+
  labs(x="Fis",y='', colour='Treated')+
  theme(text = element_text(size = 15),
        legend.position = "none")

picounts = ggplot(data=mean_divs, aes(x=pi,y=counts, colour=as.factor(treated)))+
  scale_color_manual(values = pal)+
  geom_point(size=2, alpha=0.7)+
  theme_tufte()+
  labs(x="π",y='Counts', colour='Treated')+
  theme(text = element_text(size = 15)) 

txlegend = get_legend(picounts) #extract legend

#plot final nice plot
countplots = cowplot::plot_grid(
  fstcounts,
  fiscounts,
  picounts + theme(legend.position = "none"),
  txlegend,
  labels = c('A','B','C')
)
countplots
#-----------------------------------------------------------#
#models of relationship between diffs and divs and counts
#-----------------------------------------------------------#

####NBINOM
#look at residuals for each mode,  in drop1, the year doesn't add anything, so will remove
#ggplot predict GLM


hist(mean_divs$mean_fis)
c1 = glm.nb(formula = counts ~ mean_fis* as.factor(treated) + year * as.factor(treated), data = mean_divs)


mean_divs$cluster <- as.factor(mean_divs$cluster)

c1 = lme4::glmer(formula = mean_fis ~ counts+ as.factor(treated) + year * as.factor(treated) + (1|cluster), data = mean_divs, family=binomial(link='logit'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



c1simulationOutput <- simulateResiduals(fittedModel = c1, plot = F)
residuals(c1simulationOutput, quantileFunction = qnorm, outlierValues = c(0,1))
plot(c1simulationOutput)
drop1(c1)
MuMIn::r.squaredGLMM(c1)
modc1 = glm.nb(formula = counts ~ mean_fis* as.factor(treated), data = mean_divs)
summary(modc1)

#look at residuals for each mode,  in drop1, the year doesn't add anything, so will remove
#ggplot predict GLM
c2 = glm.nb(formula = counts ~ mean_fst* as.factor(treated) + year * as.factor(treated), data = mean_divs)
summary(c2)
c2simulationOutput <- simulateResiduals(fittedModel = c2, plot = F)
residuals(c2simulationOutput, quantileFunction = qnorm, outlierValues = c(0,1))
plot(c2simulationOutput)
drop1(c2)
MuMIn::r.squaredGLMM(c2)
modc2 = glm.nb(formula = counts ~ mean_fst* as.factor(treated), data = mean_divs)
summary(modc2)

#look at residuals for each model, in drop1, the year doesn't add anything, so will remove
#ggplot predict GLM
c3 = glm.nb(formula = counts ~ pi* as.factor(treated) + year * as.factor(treated), data = mean_divs)
summary(c3)
c3simulationOutput <- simulateResiduals(fittedModel = c3, plot = F)
residuals(c3simulationOutput, quantileFunction = qnorm, outlierValues = c(0,1))
plot(c3simulationOutput)
drop1(c3) 
MuMIn::r.squaredGLMM(c3)

#-----------------------------------------------------------#
#genome scans
#-----------------------------------------------------------#

#fst





#-----------------------------------------------------------#
#isolation by distance and relatedness
#-----------------------------------------------------------#

res = fread('~/Projects/avecnet_popgen/data/relatedness/AgamP4_3L.res')
res = left_join(res, metadata_seq_fis, by=c('a'='bamorder')) %>% left_join(., metadata, by=c('b'='bamorder'))
res$pdist = geosphere::distVincentyEllipsoid(res[,c('lat.x','long.x')], res[,c('lat.y','long.y')])
ggplot(res, aes(x=pdist, y=rab))+
  geom_point()


res$txpair <- paste0(res$treated.x, res$treated.y)

ggplot(res,aes(x=pdist, y=rab))+
  geom_point()


####genome scans
#ann <- fread('/Users/dennistpw/Library/Mobile Documents/com~apple~CloudDocs/Projects/anopheles_cnrfp/VectorBase-62_AgambiaePEST.gff', col.names = c("chrom",'source','feature','start','end','pt','sense','pt2','note'))
#fsttxuntx <- fread('/Users/dennistpw/Projects/anopheles_pmi_cnrfp/treated_untreated.win5step1')
#colnames(fsttxuntx) <- c('gubbins','chrom','pos','nsites','fst')
#ggplot(fsttxuntx, aes(x=pos,y=fst,colour=chrom))+
#  theme_classic()+
#  geom_point()+
#  facet_grid(~chrom, scales = "free_x")+
#  scale_color_brewer(palette="Set2")
#
#newfst = fsttxuntx[fsttxuntx$fst > 0.01]
#
#
#
#fst_bed = data.frame(cbind(newfst$chrom,
#                           newfst$pos - 2500,
#                           newfst$pos + 2500))
#fst_bed$start <- as.numeric(fst_bed$start)
#fst_bed$end <- as.numeric(fst_bed$end)
#
#
#colnames(fst_bed) <- c('chrom', 'start', 'end')
#x = valr::bed_intersect(fst_bed, ann)
#t = x %>% select(feature.y, note.y) %>% filter(feature.y == 'protein_coding_gene')
#
#
#