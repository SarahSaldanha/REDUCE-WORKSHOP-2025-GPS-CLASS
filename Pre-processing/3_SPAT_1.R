
#-------------------------------------------------------------------------------
#
# Title: Manual Metrics
# Course: Analysis of GPS/PTT Tracking Data 
#
# Author: Sarah Saldanha
# Email: sarahsaldanha9@gmail.com
# Last revision: 2025-04-08
#
#-------------------------------------------------------------------------------
rm(list = ls())

# list.of.packages <- c("tidyverse", "tidylog", "plyr", "viridis", "sp", "sf", "adehabitatHR", "janitor", "lubridate", "fossil", "emmeans")
# 
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if (length(new.packages)) install.packages(new.packages)
# lapply(list.of.packages, require, character.only = T)
# remove(list.of.packages, new.packages)


library(lubridate)
library(dplyr)
library(sp)
library(geosphere)
library(maps)
library(tibble)
library(fossil)

###GLMM packages
library(nlme)
library(lme4)

###GAMM packages
library(gamm4)
library(mgcv)


setwd("C:/Users/benat/OneDrive - Universitat de Barcelona (1)/Escriptori/GPS-ANALYSIS")

GPS <- read.csv("output/Trips_cut.csv")

table(GPS$tripID)

######## 1. Spatial analysis

glimpse(GPS)

GPS$lat <- as.numeric(as.character(GPS$Latitude))
GPS$lon <- as.numeric(as.character(GPS$Longitude))
GPS$date_time <- ymd_hms(as.character(GPS$DateTime))
GPS$longitude_colony <- as.numeric(as.character(GPS$nestLoc_lon))
GPS$latitude_colony <- as.numeric(as.character(GPS$nestLoc_lat))





###calcualte distance between consecutive positions and difference in time 
GPS_list <- split(GPS, GPS$tripID)

for (i in seq_along(GPS_list)) {
  x <- GPS_list[[i]]
  p2p_dist <- deg.dist(x$lon, x$lat, x$lon[-1], x$lat[-1])
  p2p_dist[length(p2p_dist)] <- deg.dist(x$lon, x$lat, x$longitude_colony, x$latitude_colony)
  x$p2p_dist <- p2p_dist
  time_interval_h <- as.numeric(difftime(x$date_time[-1],
                                         x$date_time,
                                         units = "hours"))
  time_interval_h[length(time_interval_h)] <- NA
  x$time_interval_h <- time_interval_h
  GPS_list[i] <- list(x)
}

##bring everything together
GPS_data <- do.call(rbind.data.frame, GPS_list)

#clean up workshpace
rm(x, GPS_list, i, p2p_dist, time_interval_h, GPS)


###distance from the colony
GPS_data$colony_dist = deg.dist(GPS_data$lon, GPS_data$lat, GPS_data$longitude_colony, GPS_data$latitude_colony)

##For trip metrics, we need to select only the complete trips, that returned to the colony
GPS_data <- subset(GPS_data, GPS_data$Returns== "Yes")

GPS_data$tripID <- factor(as.character(GPS_data$tripID))

###add in bird Sex (make them up for now!)

IDs<-as.character(unique(GPS_data$tripID))
sex<-c("Male", "Female", "Female", "Female", "Male")
dat<-as.data.frame(cbind(IDs, sex))


GPS_data<-merge(GPS_data, dat, by.x = "tripID", by.y = "IDs")

GPS_data %>%
  dplyr::group_by(tripID, sex)%>%
  dplyr::summarize(max_dist      = max(colony_dist, na.rm=TRUE),
            duration      = sum(time_interval_h, na.rm = T),
            duration_days = sum(time_interval_h, na.rm = T)/24,
            trip_length   = sum(p2p_dist, na.rm = T),
            n_pos = length(lat)) %>%
  ungroup() -> sum_trips # this contains all trips longer than 3 hours







ggplot(sum_trips) + 
  geom_violin(aes(x = sex, y = duration_days, fill = sex), size = 1, alpha = 0.3, border = 0.01, color = NA) +
  geom_boxplot(aes(x = sex, y = duration_days), outlier.alpha = 0.01, coef = 0.5, color = "gray40", width = .15) +
  geom_jitter(aes(x = sex, y = duration_days), color="black", size=1, alpha=0.3)+
  theme(axis.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_y_continuous(name= "Duration of GPS trips (days)", breaks = c(0, 2, 4, 6, 8, 10))+theme_bw()


###practice with the other parameters










