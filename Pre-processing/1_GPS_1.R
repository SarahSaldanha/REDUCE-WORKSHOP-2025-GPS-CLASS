

#-------------------------------------------------------------------------------
#
# Title: Importing and Processing
# Course: Analysis of GPS/PTT Tracking Data 
#
# Author: Sarah Saldanha
# Email: sarahsaldanha9@gmail.com
# Last revision: 2025-04-08
#
#-------------------------------------------------------------------------------

###First lets load in and combine some GPS tracks from seabirds & turtles and combine them, this will facilitate data processing and analysis 


#Set working directory
setwd("C:/Users/benat/OneDrive - Universitat de Barcelona (1)/Escriptori/GPS-ANALYSIS")

###packages
library(stringr)
library(tidyverse)
library(rworldmap)
library(tidylog)
library(track2KBA)
library(tidyr)
library(lubridate)
library(assertthat)
library(stringr)
library(Rmisc)
library(readxl)
library(terra)
library(tidyterra)
library(aniMotum)
library(rnaturalearth)



Sys.setenv(TZ = "GMT") #####SET tz of session

### This code is developed to work with GPS of Pathtrack. You need to make changes for importing the GPS of Catlog (those are in a csv file, and could be imported changing the pattern in the list.files, but also check the separator and how many rows to skip). A second option would be to read the files separately and then to bind them with rbind.


txt_files <- list.files("./Data/pathtrack/", pattern = '.pos')

#Function to import a GPS text file and add a new column with the file name for a unique id.
gps_import <- function(file_name){
  # First import the data from the file, separator set to \t for tab deliminated data.
  gps_data <- read.table( paste0("./Data/pathtrack/",file_name) , sep = ",",  skip = 5)
  # Remove file extension - so we just have the main file name, not the .txt extension.
  file_name_only <- strsplit(file_name,".pos")
  #Make this into a unit length vector
  tracking_event <- unlist(file_name_only)
  #Add the file names as a new column to the GPS data
  gps_data <- cbind(tracking_event, gps_data)
  #Output this new table.
  return(gps_data)
}

# Run the above function for all text files, produces a list of dataframes, with a dataframe produced for each file.
data_list <- lapply(txt_files, gps_import)
# Merge list 'data_list' into a single dataframe called 'GPS_all'
GPS_all <- NULL
for(i in seq(along = txt_files)){
  GPS_all <- rbind(GPS_all, as.data.frame((data_list[i])))
}

##add in column names (since the .pos files from pathtrack don't have them!)
GPS_all <- dplyr::rename(GPS_all, dia = V1, mes = V2, year = V3, hora = V4, minutos = V5, segundos = V6, Satellites = V8, Latitude = V9, Longitude = V10)
GPS_all <- dplyr::select(GPS_all, tracking_event, dia, mes, year, hora, minutos, segundos, Satellites, Latitude, Longitude)

glimpse(GPS_all)

GPS_all$Latitude <- as.numeric(as.character(GPS_all$Latitude))
GPS_all$Longitude <- as.numeric(as.character(GPS_all$Longitude))
GPS_all$dia <- as.numeric(as.character(GPS_all$dia))
GPS_all$mes <- as.numeric(as.character(GPS_all$mes))
GPS_all$year <- as.numeric(as.character(GPS_all$year))
GPS_all$hora <- as.numeric(as.character(GPS_all$hora))
GPS_all$minutos <- as.numeric(as.character(GPS_all$minutos))
GPS_all$segundos <- as.numeric(as.character(GPS_all$segundos))

###bring all ymd_hms into one column

GPS_all <- unite(GPS_all, Date, c(2:4),  sep = "/", remove = TRUE) #many rows in one
GPS_all <- unite(GPS_all, Time, c(3:5),  sep = ":", remove = TRUE)
GPS_all <- unite(GPS_all, Date_Time, c(2,3),  sep = "_", remove = TRUE)

str(GPS_all)

##correct date format
GPS_all$Date_Time <- lubridate::dmy_hms(as.character(GPS_all$Date_Time))

###trun trip and bird ID into 
GPS_all$BirdId <- as.factor(substr(GPS_all$tracking_event, 1, 6))
GPS_all$trip_id <- as.factor(substr(GPS_all$tracking_event, 1, 20))

rm(data_list, gps_import, i, txt_files)

levels(as.factor(GPS_all$tracking_event)) #checking we imported the data properly 


################################################################################################
####bring in and join catlog tracks
###############################################################################################

Catlog_files <- list.files("./Data/catlog/",pattern = ".csv")
gps_Catlog_files <- function(file_name){
  # First import the data from the file, sepperator set to \t for tab deliminated data.
  file.error<- try(read.csv(paste0("./Data/catlog/",file_name), header=T, sep=","))
  if(is.error(file.error)==FALSE){gps_data<-read.csv(paste0("./Data/catlog/",file_name), header=T, sep=",")}
  if(is.error(file.error)==TRUE){gps_data<-read.csv(paste0("./Data/catlog/",file_name), header=T, sep=",",skip = 7)}

  data.error<- try(dplyr::select(gps_data, Date, Time,Latitude,Longitude, Altitude, Satellites))
  if(is.error(data.error)==TRUE){
    gps_data<-dplyr::select(gps_data, Date, Time,Latitude,Longitude, Altitude)
    gps_data$Satellites<-"NA"}  else {gps_data<-dplyr::select(gps_data, Date, Time,Latitude,Longitude, Altitude, Satellites)}
  
  # Remove file extension - so we just have the main file name, not the .txt extension.
  file_name_only <- strsplit(file_name,".csv")
  print(file_name_only)
  #   #Make this into a unit length vector
  TrackId <- unlist(file_name_only)
  #   #Add the file names as a new column to the GPS data
  gps_data <- cbind(TrackId, gps_data)
  #   #Output this new table.
  return(gps_data)
}

data_list <- lapply(Catlog_files, gps_Catlog_files)
# Merge list 'data_list' into a single dataframe called 'GPS_all'
GPS_catlog <- NULL
for(i in seq(along = Catlog_files)){
  GPS_catlog<- rbind(GPS_catlog, as.data.frame((data_list[i])))
}

GPS_catlog


write.csv(GPS_catlog, "output/GPS_data_joined_catlog.csv", row.names = TRUE)
saveRDS(GPS_catlog, "output/GPS_data_joined_catlog.rds")

###########################################################################################
###Explore the trips
#####################################################################################

glimpse(GPS_all)


wrld <- ne_countries(scale = "Medium", returnclass = "sf")
Coord_colony <- data.frame(Latitude = 14.970429, Longitude = -24.638775)

ggplot()+
  geom_point(data = GPS_all, aes(x=Longitude, y = Latitude, color = tracking_event)) +
  geom_point(data = Coord_colony, aes(x=Longitude, y = Latitude), fill = "black", size = 2)+
  theme_bw()

GPS_all <- subset(GPS_all,  Longitude < -10)

ggplot()+
  geom_point(data = GPS_all, aes(x=Longitude, y = Latitude, color = tracking_event)) +
  geom_point(data = Coord_colony, aes(x=Longitude, y = Latitude), fill = "black", size = 2)+
  theme_bw()

ggplot()+
  geom_point(data = GPS_all, aes(x=Longitude, y = Latitude, color = tracking_event)) +
  geom_point(data = Coord_colony, aes(x=Longitude, y = Latitude), fill = "black", size = 2)+
  geom_sf(data = wrld, fill = "gray75", color = NA, size = 0.8)+
  scale_x_continuous(breaks = seq(-35, -20, 5), name= "Longitude")+
  scale_y_continuous(breaks = seq(10, 20, 2), name= "Latitude")+
  coord_sf(xlim=c(-35, -20), ylim=c(10, 20), expand=F) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_light()


write.csv(GPS_all, "output/GPS_data_joined_pathtrack.csv", row.names = TRUE)
saveRDS(GPS_all, "output/GPS_data_joined_pathtrack.rds")


###########################################################################################
###Processing with aniMotum
#####################################################################################
GPS<-GPS_all


length(unique(GPS$tracking_event))

BirdaniMotum <- GPS %>%
  mutate(lc = "G",
         median_time_diff = 7200) %>% #set this for all since they are GPS
  dplyr::select(id = tracking_event, date = Date_Time , lc, lon =Longitude , lat =Latitude , median_time_diff)



noInter <- fit_ssm(BirdaniMotum,
                   model = "crw",
                   vmax = 20, #10 m/s
                   min.dt = 2/60, #to define the temporal window in hours use to consider near duplicates
                   time.step = NA, #No interpoation
                   map = list(rho_o = factor(NA)))
table(noInter$converged)


##run the model with interpolation takes much longer 


#I just want to regularize the positions if
#the difference in time between positions is maximum up to
#3 times the median time between position of each trip
#for that I will create a segment variable to establish
#in which segments should I regularize the positions and
#in which not

# Calculate the segment variable
data <- BirdaniMotum %>%
  mutate(median_time_diff = as.numeric(median_time_diff)) %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(time_diff = difftime(date, lag(date), units = "hours"),
         segment = cumsum(ifelse(is.na(time_diff) | time_diff > median_time_diff*3/3600, 1, 0))) %>%
  ungroup()%>%
  dplyr::select(-time_diff)

# Print the result
print(data)


# Get the beginning and ending of each segment
segments <- data %>%
  group_by(id, segment) %>%
  summarise(start = min(date), end = max(date)) %>%
  ungroup()

print(segments)

# Merge segments with the median time intervals between positions


segments$median_time_diff<-7200

# Function to create 2-hour sequences
create_sequences <- function(start, end, interval) {
  seq(from = start, to = end, by = interval)
}

# Apply the function to each row of the segments dataframe
sequences <- segments %>%
  rowwise() %>%
  mutate(datetime_seq = list(create_sequences(start, end, median_time_diff))) %>%
  ungroup() %>%
  dplyr::select(id, segment, datetime_seq)

# Expand the sequences into a long format dataframe
sequences_long <- sequences %>%
  unnest(datetime_seq)

# Print the result
print(sequences_long)

time_steps <- sequences_long %>%
  dplyr::select(-segment)%>%
  dplyr::rename(date = datetime_seq)


### 2. fit a correlated random walk model for all individuals ####
### Interpolating positions every 6 hours depending but only when the gap between positions
### is < 3*median time diff between positions

inter<- fit_ssm(BirdaniMotum,
                model = "crw",
                vmax = 20, #20 m/s
                min.dt = 2/60, #to define the temporal window in hours use to consider near duplicates
                #time.step = 6, #6h
                time.step = time_steps, # interpolation every 2 hours except when data gaps > 3*median time diff between positions
                map = list(rho_o = factor(NA)))
table(inter$converged)


###check model: 
map(inter, what = "fitted")
plot(inter,what = "predicted")

## grab predicted locations
ploc <- grab(inter, what = "predicted")
ploc[1:5,]


## add GPS id to predicted locations
ploc <- ploc %>%
  dplyr::select(-id)
id <- time_steps$id
id <- factor(id)

ploc <- cbind(id, ploc)
str(ploc)

## grab observed (prefiltered) data
fp.data <- grab(inter, what = "data")
fp.data[1:5,]
str(fp.data$id)

unique(fp.data$id)
str(id)

GPSID <- unique(id)

dev.off()


### format and save 
names(GPS_all)


ploc_to_save<- ploc %>%
  dplyr::select(tracking_event = id, Date_Time =  date, Latitude =lat, Longitude = lon)
  


###add in trip and bird ID into 
ploc_to_save$BirdId <- as.factor(substr(ploc_to_save$tracking_event, 1, 6))
ploc_to_save$trip_id <- as.factor(substr(ploc_to_save$tracking_event, 1, 20))

write.csv(ploc_to_save, "output/GPS_data_joined_pathtrack_int.csv", row.names = TRUE)
saveRDS(ploc_to_save, "output/GPS_data_joined_pathtrack_int.rds")


###run the same for the catlog GPS!
GPS_catlog$Date_Time<- ymd_hms(paste0(GPS_catlog$Date, " ", GPS_catlog$Time))



BirdaniMotum <- GPS_catlog%>%
  mutate(lc = "G",
         median_time_diff = 300) %>% #set this for all since they are GPS
  dplyr::select(id = TrackId, date = Date_Time , lc, lon =Longitude , lat =Latitude , median_time_diff)




noInter <- fit_ssm(BirdaniMotum,
                   model = "crw",
                   vmax = 20, #10 m/s
                   min.dt = 1/60, #to define the temporal window in hours use to consider near duplicates
                   time.step = NA, #No interpoation
                   map = list(rho_o = factor(NA)))
table(noInter$converged)


##run the model with interpolation takes much longer 


#I just want to regularize the positions if
#the difference in time between positions is maximum up to
#3 times the median time between position of each trip
#for that I will create a segment variable to establish
#in which segments should I regularize the positions and
#in which not

# Calculate the segment variable
data <- BirdaniMotum %>%
  mutate(median_time_diff = as.numeric(median_time_diff)) %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(time_diff = difftime(date, lag(date), units = "hours"),
         segment = cumsum(ifelse(is.na(time_diff) | time_diff > median_time_diff*3/3600, 1, 0))) %>%
  ungroup()%>%
  dplyr::select(-time_diff)

# Print the result
print(data)


# Get the beginning and ending of each segment
segments <- data %>%
  group_by(id, segment) %>%
  summarise(start = min(date), end = max(date)) %>%
  ungroup()

print(segments)

# Merge segments with the median time intervals between positions


segments$median_time_diff<-300

# Function to create 2-hour sequences
create_sequences <- function(start, end, interval) {
  seq(from = start, to = end, by = interval)
}

# Apply the function to each row of the segments dataframe
sequences <- segments %>%
  rowwise() %>%
  mutate(datetime_seq = list(create_sequences(start, end, median_time_diff))) %>%
  ungroup() %>%
  dplyr::select(id, segment, datetime_seq)

# Expand the sequences into a long format dataframe
sequences_long <- sequences %>%
  unnest(datetime_seq)

# Print the result
print(sequences_long)

time_steps <- sequences_long %>%
  dplyr::select(-segment)%>%
  dplyr::rename(date = datetime_seq)


### 2. fit a correlated random walk model for all individuals ####
### Interpolating positions every 6 hours depending but only when the gap between positions
### is < 3*median time diff between positions

inter<- fit_ssm(BirdaniMotum,
                model = "crw",
                vmax = 20, #20 m/s
                min.dt = 1/60, #to define the temporal window in hours use to consider near duplicates
                #time.step = 6, #6h
                time.step = time_steps, # interpolation every  5 mins except when data gaps > 3*median time diff between positions
                map = list(rho_o = factor(NA)))
table(inter$converged)


###check model: 
map(inter, what = "fitted")
plot(inter,what = "predicted")






