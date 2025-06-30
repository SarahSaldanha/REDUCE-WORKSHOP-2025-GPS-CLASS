rm(list = ls())

setwd("E:/Dropbox/PhD/RClasses_2022/EmBC and HMMs")
Sys.setenv(TZ='UTC')

library('ggplot2')
library("dplyr")
library('marmap') # mapping tools
library('lubridate') # datetime functions
library("viridis")# viridis colour map
library("trip")
library("oce")

options(expressions = 20000)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## folders

if (file.exists(paste("./", "1_identifiedTrips_PLOTS", sep = "")) == FALSE){
  dir.create(paste("./", "1_identifiedTrips_PLOTS", sep = ""))
}

if (file.exists(paste("./", "1_assignedBlocks", sep = "")) == FALSE){
  dir.create(paste("./", "1_assignedBlocks", sep = ""))
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Read in data

masterBirds <- readRDS(paste("./extractedTrips_allBirds.rds", sep = ""))

head(masterBirds)




#####Take a look at the trips
##########The trips have already been cut into seperate foraging trips (the refID)

theme_base <- theme(axis.title = element_text(size=12), 
                    axis.text = element_text(size=10),
                    legend.title = element_text(size=12),
                    legend.text= element_text(size=10),
                    legend.key.size = unit(0.75, "cm"),
                    panel.grid = element_blank(), panel.border = element_blank(),
                    panel.background = element_rect(fill = "white", colour = "black"))
  

ggplot(data = masterBirds, aes(lat, lon)) + 
     geom_point()+ 
    geom_path()+
facet_wrap(~refID) +
    geom_point(data = masterBirds, aes(nestLoc_lat, nestLoc_lon, col = "red", size= 2)) +
  theme(legend.position="none") +
    theme_base 




######order by trip, burst and time

masterBirds <- masterBirds[order(masterBirds$refID, masterBirds$tripID, masterBirds$dateTime), ]

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  mark points of new bird, new trip or time gap of more than 20 minutes...

cutPointsBird <- which(c(NA, diff(masterBirds$refID)) > 0)
cutPointsBird ###this should have the length of the number of trips
cutPointsTrip <- which(c(NA, diff(masterBirds$tripID)) > 0)
cutPointsTrip 
cutPointsDuration <- which(c(NA, diff(masterBirds$fixInterval)) > 1200) ###
cutPointsDuration

cutPoints <- sort(unique(c(cutPointsBird, cutPointsTrip, cutPointsDuration)))
rm(cutPointsBird, cutPointsTrip, cutPointsDuration)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## assign continuous blocks

JJ <- 1
contBlock <- rep(NA, nrow(masterBirds)) #####make a contBlock one for the first segments of the trips
contBlock[1:(cutPoints[1]-1)] <- JJ
for(II in 2:length(cutPoints)){
##II<-2  
  JJ <- JJ+1
  contBlock[cutPoints[II-1]:(cutPoints[II]-1)] <- JJ
  #contBlock ####Continue to fill in this dataframe
}
contBlock[cutPoints[II] : length(contBlock)] <- JJ +1
#contBlock ####Continue to fill in this dataframe
rm(JJ, II)
contBlock

identical(cutPoints, which(c(NA, diff(contBlock)) > 0)) ####want this to be true

masterBirds$contBlock <- contBlock


masterBirds <- masterBirds[order(masterBirds$refID, masterBirds$tripID, masterBirds$dateTime), ]


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## write out data


write.table(masterBirds, paste("./1_assignedBlocks/extractedTrips_allBirds_withContBlocks.csv", sep = ""), col.names = T, row.names = F, sep = ",")
saveRDS(masterBirds, paste("./1_assignedBlocks/extractedTrips_allBirds_withContBlocks.rds", sep = ""))



