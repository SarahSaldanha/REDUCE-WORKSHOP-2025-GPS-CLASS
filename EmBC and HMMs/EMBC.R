
##### ~~~~~~~~~~~~~~~~~~~   EMBC
# install.packages("EMbC")
library(EMbC)
# install.packages("reshape2")
library(reshape2)
# install.packages("plyr")
library(plyr)
# install.packages("ggplot2")
library(ggplot2)

setwd("E:/Dropbox/PhD/RClasses_2022/EmBC and HMMs")
load("./gps.RData")
names(gps)

str(gps)
gps$Date_Time <- as.POSIXct(gps$Date_Time, format = "%d/%m/%Y %H:%M", tz = "GMT")

# plot
unique(gps$tracking_event)
example <- subset(gps, gps$tracking_event== "8000015_19012016_22012016_c")
plot(example$Longitude, example$Latitude, asp= 1, type = "o", ylab= "Latitude", xlab= "Longitude")


#### 1.   Run EMbC algorithm  ####

# split by tracking event, so the program will know when a trip starts and ends
BM_POB_list <- split (gps, as.factor(as.character(gps$tracking_event))) 
length(unique(gps$tracking_event))

BBMM_POB <- EMbC::stbc(BM_POB_list) # run the EMBC algorithm
gps$label_EMBC_POB <- BBMM_POB@bC@A  ## 1:LL, 2:LH, 3:HL, and 4:HH, 5:NC 
sctr(BBMM_POB) # we get the graph.

head(gps)
table(gps$label_EMBC_POB)

# lblp(BBMM_POB, lims=c(100, 300))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#        interpretation of last column: behavioral mode                   #
# yellow     --> 1 (low speed & low turning angles) - Resting             #
# red        --> 2 (low speed & high turning angles) - Intensive search   #
# light blue --> 3 (high speed & low turning angles) - Relocating         #
# blue       --> 4 (high speed & high turning angles) - Extensive search  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# grey       --> 5 Unknown                                                #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

save(gps, file = "./gps_BM_POB.RData")

# Remove unknown (label=5)
gps2 <- gps[gps$label_EMBC_POB != 5,]

#### 2. Visualizing output ####
# R
temp <- slct(BBMM_POB,2) # the number corresponds to the track you want to see
view(temp)

# Google earth
setwd("E:/Dropbox/PhD/RClasses_2022/EmBC and HMMs/")
pkml(slct(BBMM_POB, 1), display=TRUE) # to plot a trajectory in kml (you can open it in google earth to see the path)


#### 3. Smoothing ####

BBMM_POB_smth <- EMbC::smth(BBMM_POB)

# visualization differences
temp <- slct(BBMM_POB,10)
view(temp)

temp2 <- slct(BBMM_POB_smth,10)
view(temp2)

#### 3.3 some possible calculations #### 

data_wide <- dcast(gps2, tracking_event + Ring ~ label_EMBC_POB)

# calculate proportions for each behaviour
names(data_wide) <- c("Tracking_event", "Ring", "Resting", "Foraging", "Relocating", "ExtensiveSearch")
props <- ddply(data_wide, ~ Ring+Tracking_event, summarize, 
               prop.foraging = Foraging/sum(Foraging, Resting, Relocating, ExtensiveSearch), 
               prop.resting = Resting/sum(Foraging, Resting, Relocating, ExtensiveSearch), 
               prop.relocating = Relocating/sum(Foraging, Resting, Relocating, ExtensiveSearch))

# boxplots by Ring

plot <- ggplot(data = props)  +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.y  = element_text(size=14),
        panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank())


prop.foraging.plot <- plot + geom_boxplot(aes(x = factor(Ring), y = prop.foraging, fill = Ring)) + 
  labs(list(title = "Proportion Foraging", x = "Ring", y = "%")) + 
  guides(fill=FALSE)

plot(prop.foraging.plot)

