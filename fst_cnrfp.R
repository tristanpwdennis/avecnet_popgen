library(tidyverse)
library(data.table)
fsttab = fread('~/Projects/anopheles_cnrfp/fst_bysiteby.txt', header=F)
fstnam = fread('~/Projects/anopheles_cnrfp/filenames.txt', header=F)
colnames(fsttab) <- c('fst_unweighted', 'fst_weighted')
fstnam = separate(fstnam, col = V1, into = c('cluster_a', 'month_a', 'year_a', 'cluster_b', 'month_b', 'year_b'))
fstframe = cbind(fstnam, fsttab)


#mean fst per site per tp with respect to all other sites and tps

fst_sub = select(fstframe, cluster_a, month_a, year_a, fst_weighted)

#split dataframe into. comparisons only across timepoints (so we are only looking at comparisons between 
#contemporary popuylations