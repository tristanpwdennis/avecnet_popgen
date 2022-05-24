######
#Metadata wrangling for anopheles_pmi project
######
library(data.table)
library(tidyverse)
setwd('~/Projects/MOVE/anopheles_pmi_cnrfp/')
lib_metadata = setDT(read.csv('metadata/complete_library_metadata.csv'))

#read irs history and clean, melt to long
sentinel_site_history = setDT(fread('metadata/pmi_sentinel_site_history.csv', header = T))
sentinel_site_history = data.table::melt.data.table(sentinel_site_history, id.vars= c('District', 'Site'))
sentinel_site_history$variable = gsub('X', '', sentinel_site_history$variable) #date needs a clean
colnames(sentinel_site_history) = c('District', 'Pop', 'collection.year', 'Treatment')
sentinel_site_history$collection.year = as.integer(sentinel_site_history$collection.year)

unique(lib_metadata$Pop)
unique(sentinel_site_history$Pop)

#setkey(sentinel_site_history, Site, Date)
#setkey(lib_metadata, Pop, collection.year, format = 'y')
lib_metadata = left_join(lib_metadata, sentinel_site_history, by = c('Pop', 'collection.year'))

qualimap_qc = read.delim('data/sequencing_qc_pmi_samples/multiqc_data/multiqc_qualimap_bamqc_genome_results.txt')
flagstat_qc = read.delim('data/sequencing_qc_pmi_samples//multiqc_data/multiqc_samtools_flagstat.txt')

flagstat_qc$cgr_id = as.numeric(gsub("\\_.*","",flagstat_qc$Sample))
qualimap_qc$cgr_id = as.numeric(gsub("\\_.*","",qualimap_qc$Sample))


fullmeta = left_join(lib_metadata, qualimap_qc, by=c('cgr_num' = 'cgr_id'))


sequenced_metadata = fullmeta[!is.na(fullmeta$bam_file)]
hist(sequenced_metadata$mean_coverage, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
qc_init_pass = sequenced_metadata[sequenced_metadata$mean_coverage >= 1]
qc_pass_minus_daves_tests =  qc_init_pass[!qc_init_pass$Pop == 'DW_test']

dw_points = filter(qc_init_pass, Pop == 'DW_test')

qc_init_pass$testornot = ifelse(qc_init_pass$Pop == 'DW_test', 'Test', 'Not_test')

a = qc_init_pass %>% select(percentage_aligned, mean_coverage, median_insert_size, Pop, testornot) %>% pivot_longer(cols=1:3) %>% 
  filter(name == 'percentage_aligned') %>% 
  ggplot(aes(x=name, y=value))+
  geom_jitter(aes(color = testornot))+
  geom_boxplot(width = 0.5, alpha=0.5)+
  theme_classic()+
  labs(x='', y='Percentage')

b = qc_init_pass %>% select(percentage_aligned, mean_coverage, median_insert_size, Pop, testornot) %>% pivot_longer(cols=1:3) %>% 
  filter(name == 'mean_coverage') %>% 
  ggplot(aes(x=name, y=value))+
  geom_jitter(aes(color = testornot))+
  geom_boxplot(width = 0.5, alpha=0.5)+
  theme_classic()+
  labs(x='', y='Mean DOC')

c = qc_init_pass %>% select(percentage_aligned, mean_coverage, median_insert_size, Pop, testornot) %>% pivot_longer(cols=1:3) %>% 
  filter(name == 'median_insert_size') %>% 
  ggplot(aes(x=name, y=value))+
  geom_jitter(aes(color = testornot))+
  geom_boxplot(width = 0.5, alpha=0.5)+
  theme_classic()+
  labs(x='', y='Median insert Size')


cowplot::plot_grid(a,b,c, align = 'h')


hist(qc_init_pass$mean_coverage, breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

#write_csv(sequenced_metadata, 'metadata/cov_over_0.csv')

bampath = '/export/projects/III-data/lamberton/tdd3v/anopheles_pmi/bam/'
cgrsuf = '_220322_L001.srt.dp.bam'

bamloc = paste0(bampath, qc_pass_minus_daves_tests$cgr_num, cgrsuf)
qc_pass_minus_daves_tests = cbind(qc_pass_minus_daves_tests, bamloc)

#add angsd bam ordering (zero indexed rownum)
qc_pass_minus_daves_tests$angsdno = seq(0, nrow(qc_pass_minus_daves_tests)-1)
#add sample loc
loc = fread('metadata/pmi/ghana_site_gps.csv')
qc_pass_minus_daves_tests = left_join(qc_pass_minus_daves_tests, loc, by=c('Pop' = 'site'))
#write pmi metadata
write_csv(x = qc_pass_minus_daves_tests, 'metadata/pmi_metadata.csv')
#all samples with cov >= 1
#write(paste0(bampath, qc_pass_minus_daves_tests$cgr_num, cgrsuf), file = 'poplists/allsamples')

true_untreated = qc_pass_minus_daves_tests[Pop == 'Kulaa' | Pop == 'Tugu']
treated_all = qc_pass_minus_daves_tests[Treatment == 'ACy' | Treatment == 'PM' | Treatment == 'CLD' | Treatment == 'CLD+DM']
acy = qc_pass_minus_daves_tests[Treatment == 'ACy']
pm = qc_pass_minus_daves_tests[Treatment == 'PM']
cld = qc_pass_minus_daves_tests[Treatment == 'CLD' | Treatment == 'CLD+DM']

write(paste0(bampath, treated_all$cgr_num, cgrsuf), file = 'poplists/treated_all')
write(paste0(bampath, acy$cgr_num, cgrsuf), file = 'poplists/treated_acy')
write(paste0(bampath, pm$cgr_num, cgrsuf), file = 'poplists/treated_pm')
write(paste0(bampath, cld$cgr_num, cgrsuf), file = 'poplists/treated_cld')

for (y in unique(qc_pass_minus_daves_tests$Pop)) {
  s = qc_pass_minus_daves_tests[Pop == y]
  write(paste0(bampath, s$cgr_num, cgrsuf), file = paste0('poplists/by_village_', y))
}

hist(qc_pass_minus_daves_tests$mean_mapping_quality)

write_csv(qc_pass_minus_daves_tests, 'metadata/allsamples_qcpass_metadata.csv')


#all vs and tps
for (s in unique(qc_pass_minus_daves_tests$collection.year)) {
  for (x in unique(qc_pass_minus_daves_tests$Pop)) {
    t = qc_pass_minus_daves_tests[Pop == x & collection.year == s]
    write(paste0(bampath, t$cgr_num, cgrsuf), file = paste0('poplists/by_v_tp_', s,'_',x))
  }
}

#all untreated
for (s in unique(qc_pass_minus_daves_tests$collection.year)) {
    t = qc_pass_minus_daves_tests[Pop == 'Kulaa' | Pop == 'Tugu']
    v = t[collection.year == s]
    write(paste0(bampath, v$cgr_num, cgrsuf), file = paste0('poplists/untreated_',s))
}



#post tx 2014 
x = qc_pass_minus_daves_tests[Pop == 'Gbullung' | Pop == 'Dimabi' | Pop == 'Woribugu']
post_tx_2014 = x[collection.year == 2014]
write(paste0(bampath, post_tx_2014$cgr_num, cgrsuf), file = 'poplists/condition5_post_tx_2014')
#write(paste0(bampath, all_treated$cgr_num, cgrsuf), file = 'poplists/all_treated')

#post tx 2016
x = qc_pass_minus_daves_tests[Pop == 'Tarikpaa' | Pop == 'Nanton' | Pop == 'Woribugu' | Pop == 'Dimabi']
post_tx_2016 = x[collection.year == 2016]
write(paste0(bampath, post_tx_2016$cgr_num, cgrsuf), file = 'poplists/condition5_post_tx_2016')

#post tx 2018
x = qc_pass_minus_daves_tests[Pop == 'Tarikpaa' | Pop == 'Nanton' | Pop == 'Woribugu' | Pop == 'Dimabi']
post_tx_2018 = x[collection.year == 2018]
write(paste0(bampath, post_tx_2016$cgr_num, cgrsuf), file = 'poplists/condition5_post_tx_2018')

#newtx_2014
x = qc_pass_minus_daves_tests[Pop == 'Nanton' | Pop == 'Bunbuna']
tx_2014 = x[collection.year == 2014]
write(paste0(bampath, tx_2014$cgr_num, cgrsuf), file = 'poplists/condition6_tx_2014')

#post tx 2016
x = qc_pass_minus_daves_tests[Pop == 'Tarikpaa' | Pop == 'Nanton' | Pop == 'Woribugu' | Pop == 'Dimabi']
tx_2014 = x[collection.year == 2016]
write(paste0(bampath, tx_2014$cgr_num, cgrsuf), file = 'poplists/condition6_posttx_2016')

x = qc_pass_minus_daves_tests[Pop == 'Kulaa' | Pop == 'Tugu']
write(paste0(bampath, x$cgr_num, cgrsuf), file = 'poplists/untreated_all')

#tx 2018
x = qc_pass_minus_daves_tests[Pop == 'Gbullung' | Pop == 'Bunbuna' | Pop == 'Binduli']
tx_2014 = x[collection.year == 2018]
write(paste0(bampath, tx_2014$cgr_num, cgrsuf), file = 'poplists/condition6_tx_2018')

#post tx 2018
x = qc_pass_minus_daves_tests[Pop == 'Tarikpaa' | Pop == 'Nanton' | Pop == 'Woribugu' | Pop == 'Dimabi']
tx_2014 = x[collection.year == 2018]
write(paste0(bampath, tx_2014$cgr_num, cgrsuf), file = 'poplists/condition6_posttx_2018')


#tx 2020
x = qc_pass_minus_daves_tests[Pop == 'Gbullung' | Pop == 'Bunbuna' | Pop == 'Binduli']
tx_2014 = x[collection.year == 2020]
write(paste0(bampath, tx_2014$cgr_num, cgrsuf), file = 'poplists/condition6_tx_2020')


write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Binduli' & collection.year == 2008]$cgr_num, cgrsuf), file = 'poplists/condition7_binduli_2008')
write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Nanton' & collection.year == 2008]$cgr_num, cgrsuf), file = 'poplists/condition7_nanton_2008')
write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Nanton' & collection.year == 2014]$cgr_num, cgrsuf), file = 'poplists/condition7_nanton_2014')
write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Bunbuna' & collection.year == 2014]$cgr_num, cgrsuf), file = 'poplists/condition7_bunbuna_2014')
write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Bunbuna' & collection.year == 2018]$cgr_num, cgrsuf), file = 'poplists/condition7_Bunbuna_2018')
write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Bunbuna' & collection.year == 2020]$cgr_num, cgrsuf), file = 'poplists/condition7_Bunbuna_2020')
write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Nanton' & collection.year == 2018]$cgr_num, cgrsuf), file = 'poplists/condition7_Nanton_2018')
write(paste0(bampath, qc_pass_minus_daves_tests[Pop == 'Kulaa' & collection.year == 2020]$cgr_num, cgrsuf), file = 'poplists/condition7_Kulaa_2020')







