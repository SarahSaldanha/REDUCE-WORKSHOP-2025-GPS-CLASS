

#-------------------------------------------------------------------------------
#
# Title: Population Kernels
# Course: Analysis of GPS/PTT Tracking Data 
#
# Author: Sarah Saldanha
# Email: sarahsaldanha9@gmail.com
# Last revision: 2025-04-08
#
#-------------------------------------------------------------------------------

######## Remove all the previous data
rm(list = ls())

######## Loading packages

# list.of.packages <- c("tidyverse", "tidylog", "plyr", "viridis", "sp", "sf", "adehabitatHR", "janitor",
#                       "lubridate", "fossil", "emmeans", "rgeos", "dplyr", "rgdal", "geosphere", "maps")
# 
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if (length(new.packages)) install.packages(new.packages)
# lapply(list.of.packages, require, character.only = T)
# remove(list.of.packages, new.packages)

setwd("C:/Users/benat/OneDrive - Universitat de Barcelona (1)/Escriptori/GPS-ANALYSIS")


load("./projections.Rdata")

selection <- read.csv("./output/Trips_cut.csv")


###Add in sex again

IDs<-as.character(unique(selection$tripID))
sex<-c("Male", "Female", "Female", "Female", "Female","Male","Male")
dat<-as.data.frame(cbind(IDs, sex))

selection<-merge(selection, dat, by.x = "tripID", by.y = "IDs")


####### 2. Analisis de kernel

selection$lon <- as.numeric(as.character(selection$X))
selection$lat <- as.numeric(as.character(selection$Y))


selection_df<-selection ##keep this for later


coordinates(selection) <- ~lon + lat
plot(selection, col = factor(selection$sex), pch = ".")
selection@proj4string <- sp::CRS(projections$WGS84)
selection_WGS <- spTransform(selection, CRS(projections$WGS84))


###Here we want to make population kernels for each Sex
selection_WGS@data <- data.frame(ID = selection_WGS$sex) #variable to overlap

kuds_sum <- kernelUD(selection_WGS, h = "href", extent = 0.5, grid = 300, same4all = TRUE)

plot(kuds_sum[[1]])

homerange  <-  kerneloverlaphr(kuds_sum, method = "BA", percent = 95, conditional = T) #overlap among home ranges
core <-  kerneloverlaphr(kuds_sum, method = "BA", percent = 50, conditional = T) #overlap among core areas

test1 <- as.data.frame(homerange)
test2 <- as.data.frame(core)

test3 <- rbind(test1, test2)

write.table(test3, "output/overlap_sex.csv", col.names = TRUE, row.names = TRUE)


### Now, let's plot these overlaps!
##take original dataset
selection_df


selection_df %>%
  filter(sex == "Male") -> selection_male

selection_df %>%
  filter(sex == "Female") -> selection_female

coordinates(selection_male) <- ~lon + lat
selection_male@proj4string <- sp::CRS(projections$WGS84)

coordinates(selection_female) <- ~lon + lat
selection_female@proj4string <- sp::CRS(projections$WGS84)

selection_WGS_male <- spTransform(selection_male, CRS(projections$WGS84))
selection_WGS_female <- spTransform(selection_female, CRS(projections$WGS84))

kuds_male <- kernelUD(selection_WGS_male, h = "href", extent = 0.5, grid = 300, same4all = TRUE)
kuds_female <- kernelUD(selection_WGS_female, h = "href", extent = 0.5, grid = 300, same4all = TRUE)

hranges_male <- getverticeshr(kuds_male, percent = 95) # home ranges
hrange_WGS_male <- spTransform(hranges_male, CRS(projections$WGS84))

careas_male <- getverticeshr(kuds_male, percent = 50) # core areas
careas_WGS_male <- spTransform(careas_male, CRS(projections$WGS84))

hranges_female <- getverticeshr(kuds_female, percent = 95) # home ranges
hrange_WGS_female <- spTransform(hranges_female, CRS(projections$WGS84))

careas_female <- getverticeshr(kuds_female, percent = 50) # core areas
careas_WGS_female <- spTransform(careas_female, CRS(projections$WGS84))

library("rnaturalearth")
wrld <- ne_countries(scale = "Medium", returnclass = "sf")

ggplot()+
  geom_polygon(data = hrange_WGS_male, aes(x=long, y = lat, group = group), fill = NA, color = "orangered1") +
  geom_polygon(data = hrange_WGS_female, aes(x=long, y = lat, group = group), fill = NA, color = "deepskyblue1") +
  geom_polygon(data = careas_WGS_male, aes(x=long, y = lat, group = group), fill = NA, color = "orangered3") +
  geom_polygon(data = careas_WGS_female, aes(x=long, y = lat, group = group), fill = NA, color = "deepskyblue4") +
  geom_sf(data = wrld, fill = "gray75", color = NA, size = 0.8)+
  scale_x_continuous(name= "Longitude")+
  scale_y_continuous(name= "Latitude")+
  coord_sf(xlim=c(-31, -21), ylim=c(11, 19), expand=F) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_light()


###Now lets bring in an env dataset to match up with these kernels


