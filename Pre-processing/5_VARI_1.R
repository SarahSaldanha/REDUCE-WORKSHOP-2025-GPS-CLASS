#-------------------------------------------------------------------------------
#
# Title: Environmental data extraction
# Course: Analysis of GPS/PTT Tracking Data 
#
# Author: Sarah Saldanha
# Email: sarahsaldanha9@gmail.com
# Last revision: 2025-04-08
#
#-------------------------------------------------------------------------------




#rm(list = ls())

# list.of.packages <- c("tidyverse", "tidylog", "plyr", "viridis", "sp", "sf", "adehabitatHR", "janitor",
#                       "lubridate", "fossil", "emmeans")
# 
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if (length(new.packages)) install.packages(new.packages)
# lapply(list.of.packages, require, character.only = T)
# remove(list.of.packages, new.packages)

library(rgeos)
library(lubridate)
library(dplyr)
library(sp)
library(rgdal)
library(geosphere)
library(maps)
library(tibble)

###GLMM packages
library(nlme)
library(lme4)

###GAMM packages
library(gamm4)
library(mgcv)
library(raster)

Sys.setenv(TZ='UTC') #

setwd("C:/Users/benat/OneDrive - Universitat de Barcelona (1)/Escriptori/GPS-ANALYSIS")

birdTrack <- read.csv("./output/Trips_cut.csv")

birdTrack$date_time<-lubridate::ymd_hms(as.character(birdTrack$time))
birdTrack$date <- as.Date(birdTrack$date_time)

############## Static variables

prof <- raster("C:/Users/Fernando/Dropbox/Doctorat/CURSOS/Curs_R_2023/04 GPS data filtering; how to calculate metrics and kernels from tracking data; how to extract corresponding environmental variables/01_Environmental_data/bathymetry.asc")
plot(prof)
hist(prof)

prof2 <- calc(prof, fun=function(x){ x[x > 0] <- 0; return(x)} )

hist(prof2)

#### Slope calculation

slope <- terrain(prof2, opt=c("slope"), unit='degrees')
prof2$layer

crs(prof2) <- CRS("+proj=longlat +datum=WGS84")
slope <- terrain(prof2, opt=c("slope"), unit='degrees')

plot(prof2)


#### Second step. Data extraction
masterBirds_rep <- subset(birdTrack, (birdTrack$Latitude == birdTrack$nestLoc_lat) == FALSE)
birdTrack$prof <- rep(NA, nrow(birdTrack))
birdTrack$slope <- rep(NA, nrow(birdTrack))
birdTrack$prof <- raster::extract(prof2, data.frame(x = birdTrack$longitude, y = birdTrack$latitude), method = 'simple')
birdTrack$slope <- raster::extract(slope, data.frame(x = birdTrack$longitude, y = birdTrack$latitude), method = 'simple')

hist(birdTrack$prof)
hist(birdTrack$slope)


###Analisis prof
ggplot(birdTrack, aes(x=sex, y= (prof))) +   
  geom_violin(aes(fill = factor(sex)), size = 1, alpha = 0.3, border = 0.01, colour = NA) +
  geom_boxplot(aes(factor(sex)), outlier.alpha = 0.01, coef = 0.5, color = "gray40", width = .15) +
  coord_flip()+
  theme(axis.text = element_text(size = 15)) +
  scale_y_continuous(name= "Depth in trips (meters)")


ggplot(birdTrack, aes(x=prof, y= slope)) +   
  geom_point(color ="gray40")+
  geom_smooth()
  




#########Variables dinamicas
library(raster)

sst <- raster::brick("C:/Users/Fernando/Dropbox/Doctorat/CURSOS/Curs_R_2023/04 GPS data filtering; how to calculate metrics and kernels from tracking data; how to extract corresponding environmental variables/01_Environmental_data/SST.nc")

plot(sst$X2021.06.01)

birds_matchDate300 <- as.Date(birdTrack$date)
dateTimes <- getZ(sst)

birdTrack$sst <- rep(NA, nrow(birdTrack))


for(II in 1:length(dateTimes)){
  
  #II<-133
  
  indREF_300 <- which(birds_matchDate300 == dateTimes[II])
  
  if(length(indREF_300) > 0){
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## open netCDF file and prep data
    sst_current <- subset(sst, paste0("X", year(dateTimes[II]), ".", ifelse(month(dateTimes[II])>9, month(dateTimes[II]), paste0("0",month(dateTimes[II]))), ".", ifelse(day(dateTimes[II])>9, day(dateTimes[II]), paste0("0",day(dateTimes[II])))))
    #CHLA <- rotate(CHLA)
    
    #ncREF[ncREF > 10000] = NA
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## extract data
    if(length(indREF_300) > 0){
      birdTrack$sst[indREF_300] <- raster::extract(sst_current, data.frame(x = birdTrack$longitude[indREF_300], y = birdTrack$latitude[indREF_300]), method = 'simple')
      
    }
    rm(sst_current)
  }
  rm(indREF_300)
  print(paste0(" - loop ", II))
  
}

birdTrack$sst_celsius <- birdTrack$sst - 273


###Analisis SST
ggplot(birdTrack, aes(x=sex, y= (sst_celsius))) +   
  geom_violin(aes(fill = factor(sex)), size = 1, alpha = 0.3, border = 0.01, colour = NA) +
  geom_boxplot(aes(factor(sex)), outlier.alpha = 0.01, coef = 0.5, color = "gray40", width = .15) +
  #geom_jitter(color="gray40", size=0.8, alpha=0.5)+
  coord_flip()+
  theme(axis.text = element_text(size = 15)) +
  scale_y_continuous(name= "SST in trips (celsius)")


