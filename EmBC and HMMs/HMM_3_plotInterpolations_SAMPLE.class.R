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

if (file.exists(paste("./", "3_trackInterpolations_samplePlots", sep = "")) == FALSE){
  dir.create(paste("./", "3_trackInterpolations_samplePlots", sep = ""))
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Read in data for plotting

birdInfo <- read.table("./birdInfo.csv", sep = ",", header = T, as.is = T)

masterBirds <- readRDS(paste("./1_assignedBlocks/extractedTrips_allBirds_withContBlocks.rds", sep = ""))
length(unique(masterBirds$refID))#####

linearTrack300 <- readRDS(paste("./2_trackInterpolations/interpolatedTRACKS_linear300sec.rds", sep = ""))
length(unique(linearTrack300$refID))

cubicTrack300 <- readRDS(paste("./2_trackInterpolations/interpolatedTRACKS_cubic300sec.rds", sep = ""))
length(unique(cubicTrack300$refID))


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  plot interpolations

plotInterpolation <- function(newTRACK, oldTRACK, folderOut, fileOut){
  
  if (file.exists(paste("./", "3_trackInterpolations_samplePlots/", folderOut, sep = "")) == FALSE){
    dir.create(paste("./", "3_trackInterpolations_samplePlots/", folderOut, sep = ""))
  }
  
  png(paste("./3_trackInterpolations_samplePlots/", folderOut, fileOut, ".png", sep=""), width = 35, height = 25, res = 150, units = "cm")
  par(mfrow = c(2, 1), mar = c(5, 4, 2, 4), oma = c(2, 2, 1, 1))
  
  plot(oldTRACK$lon ~ oldTRACK$dateTime, xlab= "time", col ="blue", ylab = "lon", main = "blue = raw | red = interp", pch = 20, cex = 0.6)
  points(oldTRACK$lon ~ oldTRACK$dateTime, type = "l", cex = 0.6, col = "blue")
  points(newTRACK$lon ~ newTRACK$dateTime, pch = 20, cex = 1, col = "red")
  points(newTRACK$lon ~ newTRACK$dateTime, type = "l", cex = 1, col = "red")
  
  plot(oldTRACK$lat ~ oldTRACK$dateTime, xlab= "time", col ="blue", ylab = "lat", main = " ", pch = 20, cex = 0.6)
  points(oldTRACK$lat ~ oldTRACK$dateTime, type = "l", cex = 0.6, col = "blue")
  points(newTRACK$lat ~ newTRACK$dateTime, pch = 20, cex = 1, col = "red")
  points(newTRACK$lat ~ newTRACK$dateTime, type = "l", cex = 1, col = "red")
  dev.off()
  
  birdPlot <- ggplot(oldTRACK, aes(lon, lat)) +
    
    geom_point(col = "blue", size = 1, alpha = 1) +   
    geom_path(col = "blue", size = 1, alpha = 0.6) +  
    
    geom_point(data = newTRACK, col = "red", size = 1, alpha = 0.5) +   
    geom_path(data = newTRACK, col = "red", size = 1, alpha = 0.5) +   
    
    xlab("lon") +
    ylab("lat") +
    ggtitle("blue = raw | red = interp") +
    coord_fixed(ratio = 1) +
    scale_x_continuous(limits = c(min(oldTRACK$lon) -abs(diff(range(oldTRACK$lon)))/5, max(oldTRACK$lon) +abs(diff(range(oldTRACK$lon)))/5), expand = c(-0.1, 0)) +
    scale_y_continuous(limits = c(min(oldTRACK$lat) -abs(diff(range(oldTRACK$lat)))/5, max(oldTRACK$lat) +abs(diff(range(oldTRACK$lat)))/5), expand = c(-0.1, 0)) +
    theme(axis.title = element_text(size=12), 
          axis.text = element_text(size=10),
          legend.title = element_text(size=12),
          legend.text= element_text(size=10),
          legend.key.size = unit(0.75, "cm"),
          panel.grid = element_blank(), panel.border = element_blank())
  
  ggsave(paste("./3_trackInterpolations_samplePlots/", folderOut, "/spatialComp_", fileOut, ".png", sep=""), birdPlot, width = 20, height = 20, units = "cm", dpi = 150)

}
    
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  select samples to plot

 set.seed(345)
 
#####take sample of birds to plot... not really usefull here but I kept in just in case you wanted to use it later. 
 
r<-unique(linearTrack300$refID)
birdInfo<-birdInfo[birdInfo$refID %in% r,]
birdsToPlot.all <- sample(birdInfo$refID, 8)

# 
birdsToPlot <- c( birdsToPlot.all)
# 
blockID <- rep(NA, length(birdsToPlot))


for(II in 1:length(birdsToPlot)){
  
   blockID[II] <- sample(unique(linearTrack300$contBlock[linearTrack300$refID == birdsToPlot[II]]), 1)
 }



table(linearTrack300$refID)
 saveRDS(birdsToPlot, paste("./3_trackInterpolations_samplePlots/sampled_birdIDs.rds", sep = ""))
 saveRDS(blockID, paste("./3_trackInterpolations_samplePlots/sampled_blockIDs.rds", sep = ""))

birdsToPlot <- readRDS(paste("./3_trackInterpolations_samplePlots/sampled_birdIDs.rds", sep = ""))
blockID <- readRDS(paste("./3_trackInterpolations_samplePlots/sampled_blockIDs.rds", sep = ""))

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  apply function to samples - linear

##%%%%%%%%%%%%%%%%%%%%%%%%
## 300 seconds

for(II in 1:length(birdsToPlot)){
  
  fileOut <- paste(linearTrack300$region[linearTrack300$contBlock == blockID[II]][1],
                   "_bird", linearTrack300$refID[linearTrack300$contBlock == blockID[II]][1], 
                  "_trip",  linearTrack300$tripID[linearTrack300$contBlock == blockID[II]][1], 
                  "_block",  linearTrack300$contBlock[linearTrack300$contBlock == blockID[II]][1], sep = "")
                   
  plotInterpolation(newTRACK = linearTrack300[linearTrack300$contBlock == blockID[II], ], oldTRACK = masterBirds[masterBirds$contBlock == blockID[II], ], folderOut = "interpLINEAR_300sec", fileOut)
    
  rm(fileOut)
  
}
rm(II)




##%%%%%%%%%%%%%%%%%%%%%%%%
## 300 seconds

for(II in 1:length(birdsToPlot)){
  
  fileOut <- paste(cubicTrack300$region[cubicTrack300$contBlock == blockID[II]][1],
                   "_bird", cubicTrack300$refID[cubicTrack300$contBlock == blockID[II]][1], 
                   "_trip",  cubicTrack300$tripID[cubicTrack300$contBlock == blockID[II]][1], 
                   "_block",  cubicTrack300$contBlock[cubicTrack300$contBlock == blockID[II]][1], sep = "")
  
  plotInterpolation(newTRACK = cubicTrack300[cubicTrack300$contBlock == blockID[II], ], oldTRACK = masterBirds[masterBirds$contBlock == blockID[II], ], folderOut = "interpCUBIC_300sec", fileOut)
  
  rm(fileOut)
  
}
rm(II)






