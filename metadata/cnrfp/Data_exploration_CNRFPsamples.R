
rm(list=ls())

#libraries
library(raster)
library(lattice)
library(R2jags)
library(coda)
library(runjags)
library(boot)
#library(ggplot2)
library(grid)
library(lubridate)
library(stringi)
library(jagsUI)
library(sp)   
library(geoR)
library(fields)
library(maptools)
library(RColorBrewer)
library(tidyverse)
library(chron)

#par(mar=c(4.5,4.5,1,1))

setwd("/Users/tristanpwdennis/OneDrive - University of Glasgow/MOVE/samples/")

#-------------#
# Import data #
#-------------#

density_data14<- read.csv('cnrfp/MosquitoesDensity_2014_MV.csv', header=T)
density_data15<- read.csv('cnrfp/MosquitoesDensity_2015_MV.csv', header=T)
#charact_data<- read.csv('HouseCharacteristics_2014_MV.csv', header=T)
rollout_data<- read.csv('cnrfp/NetsRollout_MV.csv', header=T)
nhouses_data<- read.csv('cnrfp/NHouses_perCluster_MV.csv', header=T)
tubedata <- read.csv('cnrfp/avecnet_wp6_tubes_mosquitoes.csv')


#-------------#
# Format data #
#-------------#

#Make a column with julian week time


density_data14$jweek<- lubridate::week(dmy((density_data14$Date)))
density_data14$jmonth<- lubridate::month(dmy((density_data14$Date)))

density_data15$jweek<- lubridate::week(dmy((density_data15$Date)))+52
density_data15$jmonth<- lubridate::month(dmy((density_data15$Date)))+12

#combine density data
density_data<- rbind(density_data14,density_data15)

#add jweek and jmonth to tube data as above
tubedata <- tubedata %>% mutate(jweek = case_when(
  Year == 14 ~ lubridate::week(dmy((Date.collecte))),
  Year == 15 ~ 52 + lubridate::week(dmy((Date.collecte)))
))
tubedata <- tubedata %>% mutate(jmonth = case_when(
  Year == 14 ~ lubridate::month(dmy((Date.collecte))),
  Year == 15 ~ 12 + lubridate::month(dmy((Date.collecte)))
))





#---------------------------------------#
# MOSQUITO DATA MATRICES for JAGS model #
#---------------------------------------#

Nclusters<-40
jul.time<- c(sort(unique(density_data$jweek))[1] : 104) # from week 19 to 104.
tmax<- length(jul.time)





## TOTAL DENSITY MATRIX ##

#Anopheles gambiae TOTAL density per cluster and time point
agTotal_table<- with(density_data, aggregate(density_data$TotalAG, by=list(cluster=cluster, jweek=jweek ), FUN='sum') )

agTotal_matrix<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )

for(i in 1:nrow(agTotal_table)){
  agTotal_matrix[rownames(agTotal_matrix)==agTotal_table[i,1],colnames(agTotal_matrix)==agTotal_table[i,2]] <- agTotal_table[i,3] }

agTotal_matrix
#agTotal_matrix[is.na(agTotal_matrix)] <- 0

## PAROUS & SAMPLED MATRIX ##

#Number of PAROUS individuals per cluster and time point
agParous_table<- with(density_data, aggregate(density_data$TotalParous, by=list(cluster=cluster, jweek=jweek ), FUN='sum') )

agParous_matrix<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )

for(i in 1:nrow(agParous_table)){
  agParous_matrix[rownames(agParous_matrix)==agParous_table[i,1],colnames(agParous_matrix)==agParous_table[i,2]] <- agParous_table[i,3] }

agParous_matrix
#agParous_matrix[is.na(agParous_matrix)] <- 0

#Number of DISSECTED individuals per cluster and time point
agDissec_table<- with(density_data, aggregate(density_data$TotalDissecFemales, by=list(cluster=cluster, jweek=jweek ), FUN='sum') )

agDissec_matrix<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )

for(i in 1:nrow(agDissec_table)){
  agDissec_matrix[rownames(agDissec_matrix)==agDissec_table[i,1],colnames(agDissec_matrix)==agDissec_table[i,2]] <- agDissec_table[i,3] }

agDissec_matrix  # !!!! there are 2 Parous indiv not counted in the Dissected! @ Cluster 33 month 9
agDissec_matrix[is.na(agDissec_matrix)] <- 0

#Number of NON_PAROUS
agNonParous_matrix<- agDissec_matrix - agParous_matrix
agNonParous_matrix<- ifelse(agNonParous_matrix<0, 0, agNonParous_matrix)
#agNonParous_matrix[is.na(agNonParous_matrix)] <- 0

## GRAVID MATRIX ##

agGravid_table<- with(density_data, aggregate(density_data$TotalGravAG, by=list(cluster=cluster, jweek=jweek ), FUN='sum') )

agGravid_matrix<-matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )

for(i in 1:nrow(agGravid_table)){
  agGravid_matrix[rownames(agGravid_matrix)==agGravid_table[i,1],colnames(agGravid_matrix)==agGravid_table[i,2]] <- agGravid_table[i,3] }

agGravid_matrix
#agGravid_matrix[is.na(agGravid_matrix)] <- 0

## FED/UNFED MATRIX ##

#Fed matrix
agFed_table<- with(density_data, aggregate(density_data$TotalFedAG, by=list(cluster=cluster, jweek=jweek ), FUN='sum') )

agFed_matrix<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )

for(i in 1:nrow(agFed_table)){
  agFed_matrix[rownames(agFed_matrix)==agFed_table[i,1],colnames(agFed_matrix)==agFed_table[i,2]] <- agFed_table[i,3] }

agFed_matrix
#agFed_matrix[is.na(agFed_matrix)] <- 0

#Unfed matrix
agUnfed_table<- with(density_data, aggregate(density_data$TotalUnfedAG, by=list(cluster=cluster, jweek=jweek ), FUN='sum') )

agUnfed_matrix<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )

for(i in 1:nrow(agUnfed_table)){
  agUnfed_matrix[rownames(agUnfed_matrix)==agUnfed_table[i,1],colnames(agUnfed_matrix)==agUnfed_table[i,2]] <- agUnfed_table[i,3] }

agUnfed_matrix
#agUnfed_matrix[is.na(agUnfed_matrix)] <- 0


#---------------------------------------#
# DUO COVARIATE MATRICES for JAGS model #
#---------------------------------------#

#Extract month and year of DUO rollout in each cluster
for(i in 1:nrow(rollout_data)){ rollout_data$Day[i]<- as.numeric(strsplit(as.character(rollout_data$StartDuoNets), "/")[[i]][1]) }
for(i in 1:nrow(rollout_data)){ rollout_data$Month[i]<- as.numeric(strsplit(as.character(rollout_data$StartDuoNets), "/")[[i]][2]) }
for(i in 1:nrow(rollout_data)){ rollout_data$Year[i]<- as.numeric(strsplit(as.character(rollout_data$StartDuoNets), "/")[[i]][3]) }

rollout_data$Year[89]<-2014

#add julian week
rollout_data$jweek<- ifelse(rollout_data$Year==2014, week(dmy((rollout_data$StartDuoNets))), week(dmy((rollout_data$StartDuoNets)))+52 ) 
    
#Switch for DUO nets in each cluster. Once it's delivered, it stays on
rollout_table<- with(rollout_data, aggregate(rollout_data$Year, by=list(Cluster=Cluster, jweek=jweek ), FUN='unique') )
rollout_matrix<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )

for(i in 1:nrow(rollout_table)){
  rollout_matrix[rownames(rollout_matrix)==rollout_table[i,1],colnames(rollout_matrix)==rollout_table[i,2]] <- 1 }

rollout_table[order(rollout_table$Cluster),]

rollout_matrix[1,(34-18):ncol(rollout_matrix)]<-rollout_matrix[2,(30-18):ncol(rollout_matrix)]<-rollout_matrix[3,(78-18):ncol(rollout_matrix)]<-
  rollout_matrix[4,(30-18):ncol(rollout_matrix)]<- rollout_matrix[5,(38-18):ncol(rollout_matrix)]<-rollout_matrix[6,(78-18):ncol(rollout_matrix)]<-
  rollout_matrix[7,(30-18):ncol(rollout_matrix)]<-rollout_matrix[8,(86-18):ncol(rollout_matrix)]<-rollout_matrix[9,(82-18):ncol(rollout_matrix)]<-
  rollout_matrix[10,(38-18):ncol(rollout_matrix)]<-rollout_matrix[11,(86-18):ncol(rollout_matrix)]<-rollout_matrix[12,(29-18):ncol(rollout_matrix)]<-
  rollout_matrix[13,(77-18):ncol(rollout_matrix)]<-rollout_matrix[14,(91-18):ncol(rollout_matrix)]<-rollout_matrix[15,(25-18):ncol(rollout_matrix)]<-
  rollout_matrix[16,(86-18):ncol(rollout_matrix)]<-rollout_matrix[17,(86-18):ncol(rollout_matrix)]<-rollout_matrix[18,(38-18):ncol(rollout_matrix)]<-
  rollout_matrix[19,(26-18):ncol(rollout_matrix)]<-rollout_matrix[20,(34-18):ncol(rollout_matrix)]<-rollout_matrix[21,(34-18):ncol(rollout_matrix)]<-
  rollout_matrix[22,(82-18):ncol(rollout_matrix)]<-rollout_matrix[23,(27-18):ncol(rollout_matrix)]<-rollout_matrix[24,(86-18):ncol(rollout_matrix)]<-
  rollout_matrix[25,(92-18):ncol(rollout_matrix)]<-rollout_matrix[26,(77-18):ncol(rollout_matrix)]<-rollout_matrix[27,(82-18):ncol(rollout_matrix)]<-
  rollout_matrix[28,(38-18):ncol(rollout_matrix)]<-rollout_matrix[29,(82-18):ncol(rollout_matrix)]<-rollout_matrix[30,(82-18):ncol(rollout_matrix)]<-
  rollout_matrix[31,(38-18):ncol(rollout_matrix)]<-rollout_matrix[32,(26-18):ncol(rollout_matrix)]<- rollout_matrix[33,(24-18):ncol(rollout_matrix)]<-
  rollout_matrix[34,(29-18):ncol(rollout_matrix)]<-rollout_matrix[35,(91-18):ncol(rollout_matrix)]<-rollout_matrix[36,(77-18):ncol(rollout_matrix)]<-
  rollout_matrix[37,(34-18):ncol(rollout_matrix)]<-rollout_matrix[38,(92-18):ncol(rollout_matrix)]<-rollout_matrix[39,(92-18):ncol(rollout_matrix)]<-
  rollout_matrix[40,(34-18):ncol(rollout_matrix)]<-1

rollout_matrix[is.na(rollout_matrix)] <- 0

#-------------------#
# Time decay matrix #
#-------------------#

decay_tmax<- rowSums(rollout_matrix)
decay_matrix<- rollout_matrix

for(i in 1:Nclusters){
decay_matrix[i,((ncol(decay_matrix)-decay_tmax[i])+1):ncol(decay_matrix)]<- c(0:(decay_tmax[i]-1)) }
  

#--------------------------------#
# OTHER VARIABLES for JAGS model #
#--------------------------------#



##NUMBER OF TRAPS PER CLUSTER PER TIMEPOINT
Ntraps_withNA<- as.data.frame(table(density_data$cluster,density_data$jweek))

Ntraps<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) )
for(i in 1:nrow(Ntraps_withNA)){
  Ntraps[rownames(Ntraps)==Ntraps_withNA[i,1],colnames(Ntraps)==Ntraps_withNA[i,2]] <- Ntraps_withNA[i,3] }

Ntraps[is.na(Ntraps)] <- 0

cluster<- 1:40





julweek<- as.numeric(colnames(agFed_matrix))
julmonth<-julweek

for(i in 1:length(julmonth)){
  julmonth[i]<- density_data$jmonth[which(julweek[i]==density_data$jweek)][1]
}




julmonth[3]<- 5
julmonth[5:7]<- 6
julmonth[35:38]<- 13
julmonth[39:42]<-14
julmonth[43:47]<- 15
julmonth[72:73]<- 21
julmonth[84:86]<- 24

month.matrix<- matrix(NA, nrow=Nclusters, ncol=tmax, dimnames=list(c(1:Nclusters),c(jul.time)) ) 
for(i in 1:Nclusters) {month.matrix[i,] <- julmonth }

year.matrix<- ifelse(month.matrix<=12,1,2)

##-----------------------------------------------------------------------------------------------------------##
##                                             RUN JAGS model                                                ##
##-----------------------------------------------------------------------------------------------------------##

for(i in 1:nrow(density_data)){
density_data$rollout[i]<- as.character(rollout_data$StartDuoNets)[which(
                        density_data$cluster[i]==rollout_data$Cluster & density_data$jweek[i]>=rollout_data$jweek)][1]
}

density_data$rollout2<- ifelse(is.na(density_data$rollout)==T,0,1)


density_data$Season<- ifelse(density_data$jmonth<=6, 'Y1.S1', 
                      ifelse(density_data$jmonth==7 | density_data$jmonth==8, 'Y1.S2',
                      ifelse(density_data$jmonth==9 | density_data$jmonth==10, 'Y1.S3',
                      ifelse(density_data$jmonth==11 | density_data$jmonth==12, 'Y1.S4',
                         ifelse(density_data$jmonth==17 | density_data$jmonth==18, 'Y2.S1',
                         ifelse(density_data$jmonth==19 | density_data$jmonth==20, 'Y2.S2',
                         ifelse(density_data$jmonth==21 | density_data$jmonth==22, 'Y2.S3',
                         ifelse(density_data$jmonth==23 | density_data$jmonth==24, 'Y2.S4', NA) )))))))


tubedata$Season<- ifelse(tubedata$jmonth<=6, 'Y1.S1', 
                      ifelse(tubedata$jmonth==7 |  tubedata$jmonth==8, 'Y1.S2',
                      ifelse(tubedata$jmonth==9 |  tubedata$jmonth==10, 'Y1.S3',
                      ifelse(tubedata$jmonth==11 | tubedata$jmonth==12, 'Y1.S4',
                      ifelse(tubedata$jmonth==17 | tubedata$jmonth==18, 'Y2.S1',
                      ifelse(tubedata$jmonth==19 | tubedata$jmonth==20, 'Y2.S2',
                      ifelse(tubedata$jmonth==21 | tubedata$jmonth==22, 'Y2.S3',
                      ifelse(tubedata$jmonth==23 | tubedata$jmonth==24, 'Y2.S4', NA) )))))))



tubedata %>% filter(jmonth == 24)

t = tubedata %>% dplyr::select(jmonth, jweek,Season) %>% distinct()

s = density_data %>% dplyr::select(jmonth, jweek, Season) %>% distinct()
density_data$Year = lubridate::year(dmy((density_data$Date)))


season.table1 <- aggregate(density_data$TotalUnfedAG, by=list(density_data$Season, density_data$rollout2, density_data$cluster), FUN='sum')
colnames(season.table1)<- c('season','rollout','cluster','density')

season.table2 <- aggregate(density_data$TotalUnfedAG, by=list(density_data$Season, density_data$cluster), FUN='sum')
colnames(season.table2)<- c('season','cluster','density')

#create table that we base our picks on
x = season.table2 %>% pivot_wider(names_from =season, values_from = density)










counttable
t = counttable %>% dplyr::select(cluster, season, enough_samples, density) %>% group_by(cluster) %>% summarise(totalgood = sum(enough_samples), totalAg = sum(density)) %>% filter(totalgood > 3) 
clusterpicks  = c(t$cluster, 18)
sum(t$totalAg)

clusterpicks


counttable %>% dplyr::select(season, cluster, density) %>% pivot_wider(names_from =season, values_from = density)

counttable


counttable
dplyr::select(season, cluster, density)


eligible_tubes = density_data %>% filter(Season == 'Y1.S2' | Season == 'Y1.S4' | Season == 'Y2.S2' | Season == 'Y2.S4') %>% filter(cluster %in% clusterpicks) 
sum(eligible_tubes$TotalUnfedAG)

season.table


picks %>% dplyr::select(cluster, Season, TotalUnfedAG) %>% group_by(Season, cluster) %>% summarise(sum(TotalUnfedAG))
sum(by_cluster_and_timepoint[,3])


density_data %>% mutate(yearmonth = paste0(Year,"-", Month)) %>% dplyr::select(yearmonth, cluster, rollout2) %>% pivot_wider(names_from = yearmonth, values_from = rollout2)


tube.table = tubedata %>% dplyr::select(Cluster, Season) %>% group_by(Cluster, Season) %>% count()
tube.counttable <- rbind(
  tube.table %>% dplyr::filter(Season == 'Y1.S2' & n > 20) %>% mutate(enough_samples = 1),
  tube.table %>% dplyr::filter(Season == 'Y1.S4' & n > 20) %>% mutate(enough_samples = 1),
  tube.table %>% dplyr::filter(Season == 'Y2.S2' & n > 20) %>% mutate(enough_samples = 1),
  tube.table %>% dplyr::filter(Season == 'Y2.S4' & n > 20) %>% mutate(enough_samples = 1)
)


tubebreakdown = tubedata %>% group_by(Cluster, Month, Year) %>% count()

density_breakdown = density_data %>% group_by(cluster, Month, Year) %>% summarise(AgTotal = sum(TotalUnfedAG))
density_breakdown




for(i in 1:nrow(season.table)){
season.table$lat[i]<- cluster_coords$lat[which(season.table$cluster[i]==cluster_coords$cluster)][1]
season.table$lon[i]<- cluster_coords$lon[which(season.table$cluster[i]==cluster_coords$cluster)][1]
}

llin.stable<- season.table[which(season.table$rollout==0),]
duo.stable<- season.table[which(season.table$rollout==1),]

llin.y1s1<- llin.stable[which(llin.stable$season=='Y1.S1'),]
llin.y1s2<- llin.stable[which(llin.stable$season=='Y1.S2'),]
llin.y1s3<- llin.stable[which(llin.stable$season=='Y1.S3'),]
llin.y1s4<- llin.stable[which(llin.stable$season=='Y1.S4'),]
llin.y2s1<- llin.stable[which(llin.stable$season=='Y2.S1'),]
llin.y2s2<- llin.stable[which(llin.stable$season=='Y2.S2'),]
llin.y2s3<- llin.stable[which(llin.stable$season=='Y2.S3'),]
llin.y2s4<- llin.stable[which(llin.stable$season=='Y2.S4'),]

duo.y1s1<- duo.stable[which(duo.stable$season=='Y1.S1'),]
duo.y1s2<- duo.stable[which(duo.stable$season=='Y1.S2'),]
duo.y1s3<- duo.stable[which(duo.stable$season=='Y1.S3'),]
duo.y1s4<- duo.stable[which(duo.stable$season=='Y1.S4'),]
duo.y2s1<- duo.stable[which(duo.stable$season=='Y2.S1'),]
duo.y2s2<- duo.stable[which(duo.stable$season=='Y2.S2'),]
duo.y2s3<- duo.stable[which(duo.stable$season=='Y2.S3'),]
duo.y2s4<- duo.stable[which(duo.stable$season=='Y2.S4'),]


write.csv(season.table, '~/OneDrive - University of Glasgow/MOVE/samples/rollout_mv.csv')

boxplot(llin.stable$density~llin.stable$season, xlim=c(0.5,8))
boxplot(duo.stable$density~duo.stable$season, add=T, offset=0.2)

#get Burkina Faso map
bf<-getData('GADM', country="BFA", level=1)
bf2<-getData('GADM', country="BFA", level=2)
bf3<-getData('GADM', country="BFA", level=3)


par(mfrow=c(1,1), mar=c(5,5,2,1), mgp=c(2.2,0.8,0), cex.lab=1.4,cex.axis=1.2)
plot(bf3, col='grey90', bg='grey', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2, col=NULL, bg='grey', border='blue' , add=T)
points(hh_coords$LAT~hh_coords$LONG, col=rgb(red=160/255, green=205/255, blue=90/255, alpha=0.1), pch=20, cex=0.8)
points(cluster_coords$lat~cluster_coords$lon, col='darkred',cex=2, lwd=1, pch=20)
#text(cluster_coords$lon[1:5]+0.03, cluster_coords$lat[1:5]+0.01, cluster_coords$cluster[1:5], col='darkblue',cex=1)
axis(1, labels=T)
axis(2, labels=T)


#######
#LLIN
#######

par(mfrow=c(2,4), mar=c(5,5,2,1), mgp=c(2.2,0.8,0), cex.lab=1.4,cex.axis=1.2)
plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y1s1$lat~llin.y1s1$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y1s2$lat~llin.y1s2$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y1s3$lat~llin.y1s3$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y1s4$lat~llin.y1s4$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     
#-
plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y2s1$lat~llin.y2s1$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y2s2$lat~llin.y2s2$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y2s3$lat~llin.y2s3$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(llin.y2s4$lat~llin.y2s4$lon, col='darkgoldenrod2', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     



#######
#DUO
#######


par(mfrow=c(2,4), mar=c(5,5,2,1), mgp=c(2.2,0.8,0), cex.lab=1.4,cex.axis=1.2)
plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y1s1$lat~duo.y1s1$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y1s2$lat~duo.y1s2$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y1s3$lat~duo.y1s3$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y1s4$lat~duo.y1s4$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     
#-
plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y2s1$lat~duo.y2s1$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y2s2$lat~duo.y2s2$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y2s3$lat~duo.y2s3$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     

plot(bf3, col='grey98', bg='white', border='grey40', ylim=c(9.5,11), xlim=c(-5,-3.8), mar=c(0,0,0,0), xlab='Longitude', ylab='Latitude')
plot(bf2[which(bf2$NAME_1=='Cascades'),], col=NULL, bg='grey', border='blue', add=T)
points(duo.y2s4$lat~duo.y2s4$lon, col='darkred', pch=20, cex=2)
scalebar(50, divs=2, type='bar', lonlat=TRUE, below='km', xy=c(-4.7,9.5), cex=1, xpd='n')  #click()     




#points(hh_coords$LAT~hh_coords$LONG, col=rgb(red=160/255, green=205/255, blue=90/255, alpha=0.1), pch=20, cex=0.8)
#points(cluster_coords$lat~cluster_coords$lon, col='darkred',cex=2, lwd=1, pch=20)
axis(1, labels=T)
axis(2, labels=T)



mydat<- list(Nclusters=40, tmax=tmax, NObs_G=agGravid_matrix, NObs_Fed=agFed_matrix,  NObs_N=agNonParous_matrix, 
              NObs=agTotal_matrix, NSampled=agDissec_matrix, Ntraps=Ntraps, DUO=rollout_matrix, 
              
              cluster=cluster, year=year.matrix, Nid=sp_matrix, NspS=spS_matrix )



