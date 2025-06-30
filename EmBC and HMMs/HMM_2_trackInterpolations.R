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

if (file.exists(paste("./", "2_trackInterpolations", sep = "")) == FALSE){
  dir.create(paste("./", "2_trackInterpolations", sep = ""))
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Read in data to test interpolations on

masterBirds <- readRDS(paste("./1_assignedBlocks/extractedTrips_allBirds_withContBlocks.rds", sep = ""))
length(unique(masterBirds$refID))


######select columns that I want in the interpolation
track <- select(masterBirds, lon, lat, dateTime, refID, tripID, contBlock, island)
  

head(track)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  Interpolate function


#interpStep = 300


interpFunc <- function(track, interpStep, method = c( "linear", "cubic")){

  if(difftime(max(track$dateTime), min(track$dateTime), units = "secs") > 5*(interpStep) & nrow(track) > 5){
    
  newTime <- seq(from = ceiling_date(min(track$dateTime, na.rm = T), unit = "minute"),
                 to = floor_date(max(track$dateTime, na.rm = T), unit = "minute"),
                 by = interpStep)
  
  if(method == "cubic"){
    trackNEW <- data.frame(lon = spline(track$dateTime, track$lon, xout = newTime, method = "natural")$y,
                                lat = spline(track$dateTime, track$lat, xout = newTime, method = "natural")$y,
                                dateTime = newTime)  
    }
  
  if(method == "linear"){
    trackNEW <- data.frame(lon = approx(track$dateTime, track$lon, xout = newTime)$y,
                           lat = approx(track$dateTime, track$lat, xout = newTime)$y,
                           dateTime = newTime)  
  }
  
  if(length(track$refID) > 0){
    trackNEW$refID <- rep(track$refID[1], nrow(trackNEW))
  }
  if(length(track$tripID) > 0){
    trackNEW$tripID <- rep(track$tripID[1], nrow(trackNEW))
  }
  if(length(track$contBlock) > 0){
    trackNEW$contBlock <- rep(track$contBlock[1], nrow(trackNEW))
  }
  if(length(track$region) > 0){
    trackNEW$region <- rep(as.character(track$region[1]), nrow(trackNEW))
  }
  
  return(trackNEW)
  rm(newTime, trackNEW)
  
  }
  
}
    
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  Interpolate data at 30 and 60 seconds using linear and cubic spline



linearTracks300 <- lapply(split(track, as.factor(masterBirds$contBlock)), function(x) interpFunc(x, interpStep = 300, method = "linear"))

cubicTracks300 <- lapply(split(track, as.factor(masterBirds$contBlock)), function(x) interpFunc(x, interpStep = 300, method = "cubic"))

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## unlist to giant data frame

linearTracks300_ALL <- do.call("rbind", linearTracks300)
cubicTracks300_ALL <- do.call("rbind", cubicTracks300)



length(unique(linearTracks300_ALL$refID))
length(unique(cubicTracks300_ALL$refID))#####check if you loose any short tracks here
unique(track$refID)-unique(cubicTracks300_ALL$refID)
unique(track$refID)-unique(cubicTracks300_ALL$refID) ####make sure they still line up


linearTracks300_ALL <- linearTracks300_ALL[order(linearTracks300_ALL$refID,linearTracks300_ALL$tripID, linearTracks300_ALL$dateTime), ]

cubicTracks300_ALL <- cubicTracks300_ALL[order(cubicTracks300_ALL$refID,cubicTracks300_ALL$tripID, cubicTracks300_ALL$dateTime), ]
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## save out

##%%%%%%%%%%%%%%%%%%%
##  linear
write.table(linearTracks300_ALL, paste("./2_trackInterpolations/interpolatedTRACKS_linear300sec.csv", sep = ""), col.names = T, row.names = F, sep = ",")
saveRDS(linearTracks300_ALL, paste("./2_trackInterpolations/interpolatedTRACKS_linear300sec.rds", sep = ""))


##%%%%%%%%%%%%%%%%%%%
##  cubic
write.table(cubicTracks300_ALL, paste("./2_trackInterpolations/interpolatedTRACKS_cubic300sec.csv", sep = ""), col.names = T, row.names = F, sep = ",")
saveRDS(cubicTracks300_ALL, paste("./2_trackInterpolations/interpolatedTRACKS_cubic300sec.rds", sep = ""))
