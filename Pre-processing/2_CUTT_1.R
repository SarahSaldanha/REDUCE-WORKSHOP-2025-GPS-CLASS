#-------------------------------------------------------------------------------
#
# Title: Cutting GPS trips with track2kba
# Course: Analysis of GPS/PTT Tracking Data 
#
# Author: Sarah Saldanha
# Email: sarahsaldanha9@gmail.com
# Last revision: 2025-04-08
#
#-------------------------------------------------------------------------------


#Set working directory
setwd("C:/Users/benat/OneDrive - Universitat de Barcelona (1)/Escriptori/GPS-ANALYSIS")


##clean environment
rm(list = ls())

library(lubridate)
library(track2KBA)
library(dplyr)
library(rnaturalearth)
library(ggplot2)
library(stringr)



gps <- readRDS("output/GPS_data_joined_pathtrack_int.rds")

###Now lets get ready to split trips with Track2KBA
#for this we will need to split the dateTime into two columsn for the functions of the package to work 

gps$Date<-ymd(str_sub(gps$Date_Time, 1, 10))
gps$Time<-str_sub(gps$Date_Time, 12, 20)

datagroup <- formatFields(
  dataGroup = gps, 
  fieldID   = "tracking_event", 
  fieldDate = "Date", 
  fieldTime = "Time",
  fieldLon  = "Longitude", 
  fieldLat  = "Latitude",
  formatDT = "ymd_HMS"
)



wrld <- ne_countries(scale = "Medium", returnclass = "sf")


Coord_colony <- data.frame(Latitude = 14.970429, Longitude = -24.638775)

datagroup$nestLoc_lat <- c(Coord_colony$Latitude)
datagroup$nestLoc_lon <- c(Coord_colony$Longitude)

# 2 set colony location when trips out from a centrally-located place 
colony <- datagroup %>% 
  summarise(
    Latitude  = first(nestLoc_lat),
    Longitude = first(nestLoc_lon)
  )

# 3 set some parameters to decide what constitutes a trip (from the colony)

 Trips <- tripSplit(
   dataGroup  = datagroup,
   colony     = colony,
   innerBuff  = 35,      # the minimum distance from the colony, in km
   returnBuff = 40,  # can be set further out in order to catch incomplete trips
   duration   = 15,      # hours
   rmNonTrip  = TRUE # remove the periods when the animals were not on trips.
 )
 
 mapTrips(trips = Trips, colony = colony)

cut_data <- as.data.frame(Trips)

ggplot()+
  geom_point(data = cut_data, aes(x=Longitude, y = Latitude, color = tripID)) +
  geom_point(data = Coord_colony, aes(x=Longitude, y = Latitude), fill = "black", size = 2)+
  geom_sf(data = wrld, fill = "gray75", color = NA, size = 0.8)+
  scale_x_continuous(breaks = seq(-35, -20, 5), name= "Longitude")+
  scale_y_continuous(breaks = seq(10, 20, 2), name= "Latitude")+
  coord_sf(xlim=c(-35, -20), ylim=c(10, 20), expand=F) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_light()
mapTrips(trips = trips, colony = colony)

##Select only complete trips
trips_compelte <- subset(cut_data,Returns == "Yes" )
trips<- cut_data

trips$depID<- trips$ID
trips$ID<- trips$tripID

sumTrips <- tripSummary(trips = trips, colony = colony, nests=T )

sumTrips



##project the trips
tracks <- projectTracks( dataGroup = trips, projType = 'azim', custom=T )
class(tracks)


###find h-value for kernels
hVals <- findScale(
  tracks   = tracks,
  scaleARS = TRUE,
  sumTrips = sumTrips)

hVals



KDE <- estSpaceUse(
  tracks = tracks, 
  scale = hVals$mag, 
  levelUD = 50, 
  polyOut = TRUE,
)



DF_50_KDE<-KDE$UDPolygons

KDE_95 <- estSpaceUse(
  tracks = tracks, 
  scale = hVals$step_length, 
  levelUD = 95, 
  polyOut = TRUE,
)

DF_95_KDE<-KDE_95$UDPolygons


mapKDE(KDE = KDE$UDPolygons)
mapKDE(KDE = KDE_95$UDPolygons)


###too small a sample size to run these for now and they also take a long time on a larger sample size, but here is the code for calculating the representativeness and the KBA's
repr <- repAssess(
  tracks    = tracks, 
  KDE       = KDE$KDE.Surface,
  levelUD   = 50,
  iteration = 100, 
  bootTable = FALSE)



Site <- findSite(
  KDE = KDE$KDE.Surface,
  represent = repr$out,
  levelUD = 50,
  polyOut = TRUE
)
Sitemap <- mapSite(Site, colony = colony)

write.csv(cut_data, "output/Trips_cut.csv")

###E 



