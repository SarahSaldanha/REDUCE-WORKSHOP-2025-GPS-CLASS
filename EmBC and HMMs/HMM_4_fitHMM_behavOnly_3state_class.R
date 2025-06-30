rm(list = ls())

setwd("E:/Dropbox/PhD/RClasses_2022/EmBC and HMMs")
#setwd("C:/Users/Sarah Saldanha/Dropbox/PhD/RClasses_2022/EmBC and HMMs")
Sys.setenv(TZ='UTC')

library("momentuHMM")
library("ggplot2")
library("circular") ####to set the angles 
library("parallel")
options(expressions = 20000)

# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#folders out
if (file.exists(paste("./", "dataOut", sep = "")) == FALSE){
  dir.create(paste("./", "dataOut", sep = ""))
 }
 
 if (file.exists(paste("./", "trackStates", sep = "")) == FALSE){
   dir.create(paste("./", "trackStates", sep = ""))
 }
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## read in data
birdTrack <- readRDS("./2_trackInterpolations/interpolatedTRACKS_linear300sec.rds")
class(birdTrack$dateTime)
head(birdTrack)
# 
length(unique(birdTrack$refID))
length(unique(birdTrack$tripID))

###take a subset of 3 trips to make things run fast

birdTrack<-subset(birdTrack, tripID %in% c(2:3))



# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ## calculate step lengths and turning angles
birdPrep <- data.frame(ID = birdTrack$contBlock, time = birdTrack$dateTime, lon = birdTrack$lon, lat = birdTrack$lat)
birdPrep <- prepData(birdPrep, type = "LL", coordNames =c("lon", "lat"))
class(birdPrep)
head(birdPrep)

# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ## replace block ID to have the trip IDs rather than the continuous block IDs.
 birdPrep$ID <- birdTrack$refID
 class(birdPrep)
 head(birdPrep)
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## decide on number of states and set initial starting values

# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ## try k-means clustering to get starting values for steps.
 set.seed(234)
 clusterBird_step <- kmeans(na.omit(data.frame(birdPrep$step)), 3)
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## assign initial values for STEP - look at shape and scale parametres alongside distribution types (gamma versus weibull)
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## intial values for step length - from kmeans
clusterBird_step$centers
# Set starting parameters of the mean

muS_1 <- 0.1#sort(clusterBird_step$centers)[1]
muS_2 <- sort(clusterBird_step$centers)[2]
muS_3 <- sort(clusterBird_step$centers)[3]

###set starting values for the standard deviation

sdS_1 <- sd(na.omit(birdPrep$step)[clusterBird_step[[1]] == which(clusterBird_step$centers == sort(clusterBird_step$centers)[1])])
sdS_2 <- sd(na.omit(birdPrep$step)[clusterBird_step[[1]] == which(clusterBird_step$centers == sort(clusterBird_step$centers)[2])])
sdS_3 <- sd(na.omit(birdPrep$step)[clusterBird_step[[1]] == which(clusterBird_step$centers == sort(clusterBird_step$centers)[3])])




# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## distribution - gamma
xSeq_step <- hist(birdPrep$step, breaks = 200)$breaks # create indices to plot against (the range of data we have)
 xSeq_step[1] <- 0.001 # zero returns infinty so we replace this.
# 
hist(birdPrep$step, freq = FALSE, ylab = "density", breaks = 60, xlab = "step length", main = "GAMMA")
 lines(xSeq_step,
       dgamma(xSeq_step,
              shape = muS_1^2/sdS_1^2,
              rate = muS_1/sdS_1^2), col = "gold", lwd = 3)
 lines(xSeq_step,
       dgamma(xSeq_step,
              shape = muS_2^2/sdS_2^2,
              rate = muS_2/sdS_2^2), col = "red", lwd = 3)
 lines(xSeq_step,
       dgamma(xSeq_step,
              shape = muS_3^2/sdS_3^2,
              rate = muS_3/sdS_3^2), col = "cyan", lwd = 3)

 
 rm(xSeq_step)
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## assign initial values for angle
# ## look at concentrations (high is small angels and low is larger variable angles)
# ## assume mean to be zero (see estAngleMean within HMM info to estimate this as well)
# ## descide distributions - here go with von mises
kappaA_1 <- 25
kappaA_2 <- 1
kappaA_3 <- 50
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## von mises
 xSeq_angle <- hist(birdPrep$angle, breaks = 200)$breaks # create indices to plot against (the range of data we have)
# 

hist(birdPrep$angle, freq = FALSE, ylab = "density", breaks = 60, xlab = "angle", main = "VON MISES")
lines(xSeq_angle,
      dvonmises(xSeq_angle,
                 mu = circular(0),
                 kappa = kappaA_1), col = "gold", lwd = 3)
 lines(xSeq_angle,
       dvonmises(xSeq_angle,
                 mu = circular(0),
                 kappa = kappaA_2), col = "red", lwd = 3)
 lines(xSeq_angle,
       dvonmises(xSeq_angle,
                 mu = circular(0),
                 kappa = kappaA_3), col = "cyan", lwd = 3)

# 
rm(xSeq_angle)
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## check for step lengths of zero - we need to deal with this in the model if we have step lengths of zero
sum(birdPrep$step == 0, na.rm = TRUE) #if more than zero need to define zero-mass parameter for zero inflation

# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## zero-mass step parameters - required when there are step lengths of zero - take this as proportion of zeros in dataset
# length(which(birdPrep$step == 0))/nrow(birdPrep)
# length(which(birdPrep$step[clusterBird_step[[1]] == which(clusterBird_step$centers == max(clusterBird_step$centers))] == 0))/nrow(birdPrep[clusterBird_step[[1]] == which(clusterBird_step$centers == max(clusterBird_step$centers)), ])
# length(which(birdPrep$step[clusterBird_step[[1]] == whic
# 
 #zeroMass_1 <- length(which(birdPrep$step[clusterBird_step[[1]] == which(clusterBird_step$centers == max(clusterBird_step$centers))] == 0))/nrow(birdPrep[clusterBird_step[[1]] == which(clusterBird_step$centers == max(clusterBird_step$centers)), ])
 #zeroMass_2 <- length(which(birdPrep$step[clusterBird_step[[1]] == which(clusterBird_step$centers == min(clusterBird_step$centers))] == 0))/nrow(birdPrep[clusterBird_step[[1]] == which(clusterBird_step$centers == max(clusterBird_step$centers)), ])
# 
 rm(clusterBird_step)
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## compile initial starting values - put all starting values togather to input into the model
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## expected values
 stepPar0 <- c(muS_1, muS_2, muS_3, sdS_1, sdS_2, sdS_3)
 stepPar0
 anglePar0 <- c(kappaA_1, kappaA_2, kappaA_3) # for von mises
anglePar0
# 
rm(muS_1, muS_2, sdS_1, sdS_2,sdS_3, kappaA_1, kappaA_2, kappaA_3)
#rm(zeroMass_1, zeroMass_2)

####
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## fit model - for expected starting values
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ## expected model
 fitTrack_behavOnlyall_3states <- fitHMM(data = birdPrep, nbStates = 3,dist = list(step ="gamma", angle = "vm"), #"weibull" for weibull, "vm" for von mises, "wrpcauchy" for wrapped cauchy
                             Par0 = list(step = stepPar0, angle = anglePar0),
                             formula = ~ 1,  optMethod="Nelder-Mead")
 
saveRDS(  fitTrack_behavOnlyall_3states, paste0("./", "dataOut/3state_300ec_fittedHMM_forageTrips_3states.rds"))
save.image(paste0("./", "dataOut/backUp_modelFitting.Rdata"))
# 
load(paste0("./", "dataOut/backUp_modelFitting.Rdata"))
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## have a quick look at model outputs
 plot(fitTrack_behavOnlyall_3states, plotCI = TRUE)
# dev.off()
# 
plotStates(fitTrack_behavOnlyall_3states)

# 
 fitTrack_behavOnlyall_3states$mle$gamma
# 
 fitTrack_behavOnlyall_3states$mle$step
 fitTrack_behavOnlyall_3states$mle$angle
#
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## refit model using varying starting values that are based on those from the first candidate model (to check convergence of this model to global rather than local maxima)
 muS_1 <-  fitTrack_behavOnlyall_3states$mle$step[1, 1]
 muS_2 <-  fitTrack_behavOnlyall_3states$mle$step[1, 2]
 muS_3 <-  fitTrack_behavOnlyall_3states$mle$step[1, 3]
 sdS_1 <-  fitTrack_behavOnlyall_3states$mle$step[2, 1]
 sdS_2 <-  fitTrack_behavOnlyall_3states$mle$step[2, 2]
 sdS_3 <-  fitTrack_behavOnlyall_3states$mle$step[2, 3]
 
 
 #zeroMass_1 <-  fitTrack_behavOnlyall_3states$mle$step[3, 1]
# zeroMass_2 <-  fitTrack_behavOnlyall_3states$mle$step[3, 2]
# 
 kappaA_1 <-  fitTrack_behavOnlyall_3states$mle$angle[2, 1]
 kappaA_2 <-  fitTrack_behavOnlyall_3states$mle$angle[2, 2]
 kappaA_3 <-  fitTrack_behavOnlyall_3states$mle$angle[2, 3]
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## sequence of values to try - step
 set.seed(123)
 muS_1seq <- abs(rnorm(10, muS_1, muS_1))
 muS_2seq <- abs(rnorm(10, muS_2, muS_2))
 muS_3seq <- abs(rnorm(10, muS_3, muS_3))
# 
 set.seed(123)
 sdS_1seq <- abs(rnorm(10, sdS_1, sdS_1))
 sdS_2seq <- abs(rnorm(10, sdS_2, sdS_2))
 sdS_3seq <- abs(rnorm(10, sdS_3, sdS_3))
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## plot out sequences - gamma
 xSeq_step <- hist(birdPrep$step, breaks = 120)$breaks
 xSeq_step[1] <- 0.001 # zero returns infinty so we replace this
# 

hist(birdPrep$step, freq = FALSE, ylab = "density", breaks = 60, xlab = "step length", main = "GAMMA")
# 
 for(II in 1:length(muS_1seq)){
# 
   lines(xSeq_step,
         dgamma(xSeq_step,
                shape = muS_1seq[II]^2/sdS_1seq[II]^2,
                rate = muS_1seq[II]/sdS_1seq[II]^2), col = "gold", lwd = 1)
   lines(xSeq_step,
         dgamma(xSeq_step,
                shape = muS_2seq[II]^2/sdS_2seq[II]^2,
                rate = muS_2seq[II]/sdS_2seq[II]^2), col = "red", lwd = 1)
   
   
   lines(xSeq_step,
         dgamma(xSeq_step,
                shape = muS_3seq[II]^2/sdS_3seq[II]^2,
                rate = muS_3seq[II]/sdS_3seq[II]^2), col = "cyan", lwd = 1)
 
 }
 rm(II)
# 
 lines(xSeq_step,
       dgamma(xSeq_step,
              shape = muS_1^2/sdS_1^2,
              rate = muS_1/sdS_1^2), col = "gold", lwd = 3)
 lines(xSeq_step,
       dgamma(xSeq_step,
              shape = muS_2^2/sdS_2^2,
              rate = muS_2/sdS_2^2), col = "red", lwd = 3)
 
 lines(xSeq_step,
       dgamma(xSeq_step,
              shape = muS_3^2/sdS_3^2,
              rate = muS_3/sdS_3^2), col = "cyan", lwd = 3)
 
 

 rm(xSeq_step)
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## sequence of values to try - angle
 set.seed(471355)
 kappaA_1seq <- abs(rnorm(10, kappaA_1, kappaA_1))
 kappaA_2seq <- abs(rnorm(10, kappaA_2, kappaA_2))
 kappaA_3seq <- abs(rnorm(10, kappaA_3, kappaA_3))
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## plot out sequences - von mises
 xSeq_angle <- hist(birdPrep$angle, breaks = 120)$breaks
# 
 
 hist(birdPrep$angle, freq = FALSE, ylab = "density", breaks = 100, xlab = "angle", main = "VON MISES")
# 
 for(II in 1:length(kappaA_1seq)){
# 
   lines(xSeq_angle,
         dvonmises(xSeq_angle,
                   mu = circular(0),
                   kappa = kappaA_1seq[II]), col = "gold", lwd = 1)
   lines(xSeq_angle,
         dvonmises(xSeq_angle,
                   mu=circular(0),
                   kappa = kappaA_2seq[II]), col = "red", lwd = 1)
   lines(xSeq_angle,
         dvonmises(xSeq_angle,
                   mu=circular(0),
                   kappa = kappaA_3seq[II]), col = "cyan", lwd = 1)
# 
}
 rm(II)
# 
 lines(xSeq_angle,
       dvonmises(xSeq_angle,
                 mu = circular(0),
                 kappa = kappaA_1), col = "gold", lwd = 3)
 lines(xSeq_angle,
       dvonmises(xSeq_angle,
                 mu = circular(0),
                 kappa = kappaA_2), col = "red", lwd = 3)
 lines(xSeq_angle,
       dvonmises(xSeq_angle,
                 mu = circular(0),
                 kappa = kappaA_3), col = "cyan", lwd = 3)
 
 dev.off()
# 
 rm(xSeq_angle)
# 
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ## compile sequence values
 stepPar0_seq <- Map(function(w, x, y, z, a, b) c(w, x, y, z, a, b), w = muS_1seq, x = muS_2seq, y = muS_3seq, z = sdS_1seq, a = sdS_2seq, b = sdS_3seq)
 anglePar0_seq <- Map(function(x, y, z) c(x, y, z), x = kappaA_1seq, y = kappaA_2seq, z = kappaA_3seq)
# 
 rm(muS_1seq, muS_2seq, sdS_1seq, sdS_2seq, kappaA_1seq, kappaA_2seq)
 rm(zeroMass_1, zeroMass_2)
# 
 rm(muS_1, muS_2, sdS_1, sdS_2, kappaA_1, kappaA_2)
# 

# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ## fit sequence models to test for local versus global maxima
.f <- function(x, y, z, iterID){
   fittedHMM <- fitHMM(data = z, nbStates = 3,
                       dist = list(step ="gamma", angle = "vm"),
                       Par0 = list(step = x, angle = y),
                       formula = ~1)
   saveRDS(fittedHMM, paste0("./", "dataOut/modelIter_3state_300sec_fittedHMM_forageTrips_iter", iterID, ".rds"))
   return(fittedHMM)
 }
 
modelIter_startPar0 <- mcMap(.f, x = stepPar0_seq[c(1:10)], y = anglePar0_seq[c(1:10)], z = list(birdPrep), iterID = c(1:10))


saveRDS(modelIter_startPar0, paste0("./", "dataOut/modelIter_3state_300sec_fittedHMM_forageTrips_allBirdsNoColony.rds"))
 
save.image(paste0("./", "dataOut/backUp_modelFitting_iteration.Rdata"))

rm(list = ls())
load(paste0("./", "dataOut/backUp_modelFitting_iteration.Rdata"))

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# look at negative log likelihoods and AIC - we wont these to be small
 fitTrack_behavOnlyall_3states$mod$minimum
AIC( fitTrack_behavOnlyall_3states)

unlist(unname(lapply(modelIter_startPar0, function(x) x$mod$minimum)))
min(unlist(unname(lapply(modelIter_startPar0, function(x) x$mod$minimum)))) # best model is model 1 to 8

unlist(unname(lapply(modelIter_startPar0, function(x) AIC(x))))
min(unlist(unname(lapply(modelIter_startPar0, function(x) AIC(x))))) # best model is model 1 to 8

 fitTrack_behavOnlyall_3states_OLD <-  fitTrack_behavOnlyall_3states
rm( fitTrack_behavOnlyall_3states)

 fitTrack_behavOnlyall_3states <- modelIter_startPar0[[1]]

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## investigate validated model - add info to dataframe, save and plot out
identical(birdPrep$x, birdTrack$lon)

plot( fitTrack_behavOnlyall_3states, plotCI = TRUE) 
dev.off()

plotStates( fitTrack_behavOnlyall_3states)
dev.off()

birdTrack$state <- viterbi( fitTrack_behavOnlyall_3states)
head(birdTrack$state)
table(birdTrack$state)

birdProbs <- stateProbs(fitTrack_behavOnlyall_3states)
head(birdProbs)
birdTrack$stateProb_state1 <- birdProbs[, 1]
birdTrack$stateProb_state2 <- birdProbs[, 2]
birdTrack$stateProb_state3 <- birdProbs[, 3]
rm(birdProbs)

## look at CI around states
birdCI <- CIreal(fitTrack_behavOnlyall_3states)

fitTrack_behavOnlyall_3states$mle$step
birdCI$step$lower
birdCI$step$upper

fitTrack_behavOnlyall_3states$mle$angle
birdCI$angle$lower
birdCI$angle$upper

fitTrack_behavOnlyall_3states$mle$gamma
birdCI$gamma$lower
birdCI$gamma$upper

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## if all looking good - save out model
saveRDS(birdTrack, paste0("./", "dataOut/trackInterp_300sec_forageTrips_withBehav_allBirdsNoColony.rds"))
write.table(birdTrack, paste0("./", "dataOut/trackInterp_¡300sec_forageTrips_withBehav_allBirdsNoColony.csv"), sep = ',', col.names = T, row.names = F)


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## make some plots of the model outputs
theme_base <- theme(axis.title = element_text(size=12), 
                    axis.text = element_text(size=10),
                    legend.title = element_text(size=12),
                    legend.text= element_text(size=10),
                    legend.key.size = unit(0.75, "cm"),
                    panel.grid = element_blank(), panel.border = element_blank(),
                    panel.background = element_rect(fill = "white", colour = "black"))
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## state distributions
breaksPlot_step = hist(birdPrep$step, breaks = 60)$breaks
breaksPlot_angle = hist(birdPrep$angle, breaks = 60)$breaks

##%% step
muS_1 <- fitTrack_behavOnlyall_3states$mle$step[1, 1]
sdS_1 <- fitTrack_behavOnlyall_3states$mle$step[2, 1]
muS_2 <- fitTrack_behavOnlyall_3states$mle$step[1, 2]
sdS_2 <- fitTrack_behavOnlyall_3states$mle$step[2, 2]
muS_3 <- fitTrack_behavOnlyall_3states$mle$step[1, 3]
sdS_3 <- fitTrack_behavOnlyall_3states$mle$step[2, 3]
kappaA_1 <- fitTrack_behavOnlyall_3states$mle$angle[2, 1]
kappaA_2 <- fitTrack_behavOnlyall_3states$mle$angle[2, 2]
kappaA_3 <- fitTrack_behavOnlyall_3states$mle$angle[2, 3]

dev.off()
png(filename = paste("./dataPlots/", "fittedStepLengthsAngles_histogramStates3_states.png", sep=""), width = 25, height = 20, res = 150, units = "cm")
par(mfrow = c(2, 1), mar = c(5, 4, 2, 4), oma = c(2, 2, 1, 1))
hist(birdPrep$step, ylab = "density", breaks = 100, xlab = "step length", main = "", freq = FALSE)
lines(breaksPlot_step,
      dgamma(breaksPlot_step,
             shape = muS_1^2/sdS_1^2,
             rate = muS_1/sdS_1^2), col = "orange", lwd = 1)
lines(breaksPlot_step,
      dgamma(breaksPlot_step,
             shape = muS_2^2/sdS_2^2,
             rate = muS_2/sdS_2^2), col = "royalblue", lwd = 1)

lines(breaksPlot_step,
      dgamma(breaksPlot_step,
             shape = muS_3^2/sdS_3^2,
             rate = muS_3/sdS_3^2), col = "cyan", lwd = 1)

hist(birdPrep$angle, ylab = "density", breaks = 100, xlab = "turning angle", main = "", freq = FALSE)
lines(breaksPlot_angle,
      dvonmises(breaksPlot_angle,
                mu = circular(0),
                kappa = kappaA_1), col = "orange", lwd = 1)
lines(breaksPlot_angle,
      dvonmises(breaksPlot_angle,
                mu=circular(0),
                kappa = kappaA_2), col = "royalblue", lwd = 1)
lines(breaksPlot_angle,
      dvonmises(breaksPlot_angle,
                mu=circular(0),
                kappa = kappaA_3), col = "cyan", lwd = 1)
dev.off()

rm(breaksPlot_step, breaksPlot_angle, muS_1, sdS_1, muS_2, sdS_2, kappaA_1, kappaA_2)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## tracks with states
set.seed(93873912)
birdsToSample <- sample(unique(birdTrack$tripBlockID), 200)


theme_base <- theme(axis.title = element_text(size=12), 
                    axis.text = element_text(size=10),
                    legend.title = element_text(size=12),
                    legend.text= element_text(size=10),
                    legend.key.size = unit(0.75, "cm"),
                    panel.grid = element_blank(), panel.border = element_blank(),
                    panel.background = element_rect(fill = "white", colour = "black"))


for(II in 1:length(birdsToSample)){
  
  birdTemp <- birdTrack[birdTrack$tripBlockID == birdsToSample[II], ]
  mapExtendLon <- abs(diff(range(birdTemp$lon)))/5
  mapExtendLat <- abs(diff(range(birdTemp$lat)))/5
  
  birdPlot <- ggplot(data = birdTemp, aes(lon, lat)) +
    geom_path(data = birdTemp, aes(lon, lat), size = 0.5, alpha = 0.2) +  
    geom_point(data = birdTemp[birdTemp$state == 1, ], 
               aes(lon, lat), size = 1, alpha = 1, color = "orange") +   
    geom_point(data = birdTemp[birdTemp$state == 2, ], 
               aes(lon, lat), size = 1, alpha = 1, color = "red") +
    geom_point(data = birdTemp[birdTemp$state == 3, ], 
               aes(lon, lat), size = 1, alpha = 1, color = "cyan") +
    geom_point(data = birdTemp[1, ], aes(nestLoc_lon, nestLoc_lat), pch = "\u2605", size = 6, color = "black") +
    xlab("lon") +
    ylab("lat") +
    coord_fixed(ratio = 1) +
    scale_x_continuous(limits = c(min(birdTemp$lon) -mapExtendLon, max(birdTemp$lon) +mapExtendLon), expand = c(-0.1, 0)) +
    scale_y_continuous(limits = c(min(birdTemp$lat) -mapExtendLat, max(birdTemp$lat) +mapExtendLat), expand = c(-0.1, 0)) +
    theme_base
  ggsave(paste("./", "trackStates/bird", birdTemp$refID[1], "_trip", birdTemp$tripID[1], ".png", sep=""), birdPlot, width = 20, height = 20, units = "cm", dpi = 150)

  rm(birdTemp, birdPlot, mapExtendLat, mapExtendLon)
  
}



for(II in  unique(birdTrack$tripBlockID)){
  
  birdTemp <- birdTrack[birdTrack$tripBlockID == II,]
  mapExtendLon <- abs(diff(range(birdTemp$lon)))/5
  mapExtendLat <- abs(diff(range(birdTemp$lat)))/5
  
  birdPlot <- ggplot(data = birdTemp, aes(dateTime,nestDistance)) +
    geom_path(data = birdTemp, aes(dateTime,nestDistance), size = 0.5, alpha = 0.2) +  
    geom_point(data = birdTemp[birdTemp$state == 1, ], 
               aes(dateTime,nestDistance), size = 1, alpha = 1, color = "orange") +   
    geom_point(data = birdTemp[birdTemp$state == 2, ], 
               aes(dateTime,nestDistance), size = 1, alpha = 1, color = "red") +
    
    geom_point(data = birdTemp[birdTemp$state == 3, ], 
               aes(dateTime,nestDistance), size = 1, alpha = 1, color = "cyan") +
    xlab("DateTime") +
    ylab("Distance from the Colony") +
    coord_fixed(ratio = 1) +
    theme_base
  ggsave(paste("./dataplots/", birdTemp$refID[1], "_timeseries", birdTemp$tripBlockID[1], ".png", sep=""), birdPlot, width = 20, height = 20, units = "cm", dpi = 150)
  
  rm(birdTemp, birdPlot, mapExtendLat, mapExtendLon)
  












rm(II)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## residuals
modRes <- pseudoRes(fitTrack_behavOnlyall_3states)

png(filename = paste("./", "modelResiduals_plotsACF3_states.png", sep=""), width = 25, height = 20, res = 150, units = "cm")
plotPR(fitTrack_behavOnlyall_3states)
dev.off()




for(II in  unique(birdTrack$tripBlockID)){
  
  birdTemp <- birdTrack[birdTrack$tripBlockID == II,]
  mapExtendLon <- abs(diff(range(birdTemp$lon)))/5
  mapExtendLat <- abs(diff(range(birdTemp$lat)))/5
  
  birdPlot <- ggplot(data = birdTemp, aes(lon, lat)) +
    geom_path(data = birdTemp, aes(lon, lat), size = 0.5, alpha = 0.2) +  
    geom_point(data = birdTemp[birdTemp$state == 1, ], 
               aes(lon, lat), size = 1, alpha = 1, color = "orange") +   
    geom_point(data = birdTemp[birdTemp$state == 2, ], 
               aes(lon, lat), size = 1, alpha = 1, color = "red") +
    
    geom_point(data = birdTemp[birdTemp$state == 3, ], 
               aes(lon, lat), size = 1, alpha = 1, color = "blue") +   
    geom_point(data = birdTemp[birdTemp$state== 4, ], 
               aes(lon, lat), size = 1, alpha = 1, color = "cyan") +
    
    
    geom_point(data = birdTemp[1, ], aes(nestLoc_lon, nestLoc_lat), pch = "\u2605", size = 6, color = "black") +
    xlab("lon") +
    ylab("lat") +
    coord_fixed(ratio = 1) +
    scale_x_continuous(limits = c(min(birdTemp$lon) -mapExtendLon, max(birdTemp$lon) +mapExtendLon), expand = c(-0.1, 0)) +
    scale_y_continuous(limits = c(min(birdTemp$lat) -mapExtendLat, max(birdTemp$lat) +mapExtendLat), expand = c(-0.1, 0)) +
    theme_base
  ggsave(paste("./dataplots/", birdTemp$refID[1], "_triplabels", birdTemp$tripBlockID[1], ".png", sep=""), birdPlot, width = 20, height = 20, units = "cm", dpi = 150)
  
  rm(birdTemp, birdPlot, mapExtendLat, mapExtendLon)
  
}




















































































