
## ----packages, include=FALSE, message = FALSE, waring = FALSE------------
.libPaths("Packages")

library(raster)
library(sp)
library(TwGeos)
library(MASS)
library(mth, lib.loc = "C:/Users/hallworthm/R_Library")
library(SGAT)
#library(maptools)
#library(FLightR)

## ----SpatialData, warning = FALSE, echo = FALSE--------------------------
# Read in spatial layers 
Americas <- raster::shapefile("Spatial_Layers/Americas.shp")

OSFLDist <- raster::shapefile("Spatial_Layers/OliveSidedFlycatcher.shp")

States <- raster::shapefile("Spatial_Layers/st99_d00.shp")

Canada <- raster::shapefile("Spatial_Layers/Canada.shp")

AK <- subset(States,NAME=="Alaska")

# Define projections used in the analysis
WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
NAEA <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
NAEA <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform Canada into WGS84

Canada <- spTransform(Canada,CRS = WGS84)

# Set the projection for shapefiles
crs(States) <- WGS84
crs(Americas) <- WGS84
crs(OSFLDist) <- WGS84

# Create a land mask without Galapagos Islands so birds are not able to migrate to and from the galapagos
Land <- subset(Americas, SUBREGION == 5 | SUBREGION == 21 | SUBREGION == 13)

E<-disaggregate(subset(Land, NAME == "Ecuador"))

Ecuador <- E[14,]

Land <- (subset(Land, NAME != "Ecuador"))
Land <-(rgeos::gUnion(Land,Ecuador))

LandMask <- shapefile("Spatial_Layers/LandMask.shp")

Land <- crop(LandMask,Land)

## ----birdId--------------------------------------------------------------
AncFiles <- list.files(path = "Data/Anchorage_recoveries",
                       pattern = ".lux",
                       full.names = TRUE)

FairFiles <- list.files(path = "Data/Fairbanks_recoveries",
                        pattern = "driftadj.lux",
                        full.names=TRUE)

TetlinFiles <- list.files(path = "Data/Tetlin_recoveries",
                          pattern = ".lux", 
                          full.names = TRUE)

# Read just the file names for Bird ID

AncNames <- list.files(path = "Data/Anchorage_recoveries",
                       pattern=".lux")

FairNames <- list.files(path = "Data/Fairbanks_recoveries",
                        pattern="driftadj.lux")

TetlinNames <- list.files(path = "Data/Tetlin_recoveries",
                          pattern = ".lux")

# Combine the files into a single object

OSFLFiles<-c(AncFiles,FairFiles,TetlinFiles)

# Combine the bird ID
BirdId <- c(AncNames,FairNames,TetlinNames)

# More readable BirdId 
BirdId <- substr(x = BirdId,
                 start = 1,
                 stop = 4)

# Determine the number of birds
nBirds<-length(BirdId)

BirdId <- c(BirdId,paste0(BirdId[4],"_Yr2"))

## ----readLux-------------------------------------------------------------
OSFLdata <- lapply(X = OSFLFiles, 
                   FUN = readMTlux) 

## ----loglight------------------------------------------------------------
for(i in 1:nBirds){
OSFLdata[[i]]$Light<-log(OSFLdata[[i]]$Light)
}

## ----caplocs-------------------------------------------------------------
# Set the capture coordinates for each bird #
CapLocs<-array(NA,c(nBirds,2))

CapLocs[1,] <- c(-149.60425,61.38084)    # F202
CapLocs[2,] <- c(-149.74919,61.28071)    # F203
CapLocs[3,] <- c(-149.80309,61.29737)    # F311
CapLocs[4,] <- c(-146.832960,65.335090)  # K288 
CapLocs[5,] <- c(-147.1002,64.71122)     # K648
CapLocs[6,] <- c(-146.800030,65.33139)   # K654
CapLocs[7,] <- c(-146.943180,65.34078)   # K657
CapLocs[8,] <- c(-142.30559,63.13592)    # K659
CapLocs[9,] <- c(-146.832960,65.335090)  # K660
CapLocs[10,] <- c(-147.1002,64.71122)    # K665
CapLocs[11,] <- c(-147.758300,64.936900) # K670
CapLocs[12,] <- c(-146.832960,65.335090) # S050 
CapLocs[13,] <- c(-146.832960,65.335090) # S052 
CapLocs[14,] <- c(-141.833333,62.666667) # Q330
CapLocs[15,] <- c(-141.833333,62.666667) # Q334

## ----seed----------------------------------------------------------------
seed <- as.POSIXct(c(rep("2014-01-01 04:00:00",3), 
                     rep("2015-01-01 04:00:00",8),
                     rep("2016-01-01 04:00:00",4)),format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

# ## ----echo = FALSE--------------------------------------------------------
# lightImage(OSFLdata[[1]], offset = 19, zlim = c(0,1))
# points(as.POSIXct("2014-01-01",format = "%Y-%m-%d",tz = "GMT"), 24+4, col = "red",pch = 19)

## ------------------------------------------------------------------------
# Create an empty list to store the results
twl <- vector('list',nBirds+1)

# Loop through the files and assign twilights using a threshold value of 1. 

for(i in 1:nBirds){
  twl[[i]] <- findTwilights(tagdata = OSFLdata[[i]],
                            threshold = 1,
                            include = seed[i],
                            dark.min = 0) # 0 hours minimum dark period
}

# K288 has two years worth of data, therefore we need to set two seeds for that bird
twl[[4]] <- findTwilights(tagdata = OSFLdata[[4]],
                          threshold = 1,
                          include = c(seed[4],seed[nBirds]),
                          dark.min = 0) # 0 hours minimum dark period

twl[[16]] <- twl[[4]][595:nrow(twl[[4]]),]
twl[[4]] <- twl[[4]][1:594,]


for(i in 1:(nBirds+1)){
twl[[i]] <- twilightEdit(twilights = twl[[i]], 
                           window = 4,           # two days before and two days after
                           outlier.mins = 45,    # difference in mins
                           stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                           plot = FALSE)
}

twl <- Map(cbind, twl, threshold = 1)


## ------------------------------------------------------------------------
# Create an empty list to store the results
twl_Ext <-twl_Edit2<- vector('list',(nBirds+1))

for(i in 1:nBirds){
  twl_Ext[[i]] <- findTwilights(tagdata = OSFLdata[[i]],
                                threshold = 4,
                                include = seed[i],
                                dark.min = 0) # 0 hours minimum dark period
}

# K288 has two years worth of data, therefore we need to set two seeds for that bird
twl_Ext4 <- findTwilights(tagdata = OSFLdata[[4]],
                              threshold = 4,
                              include = c(seed[4],seed[nBirds]),
                              dark.min = 0) # 0 hours minimum dark period

twl_Ext[[16]] <- twl_Ext4[701:nrow(twl_Ext[[4]]),]
twl_Ext[[4]] <- twl_Ext4[1:700,]

for(i in 1:(nBirds+1)){
twl_Ext[[i]] <- twilightEdit(twilights = twl_Ext[[i]], 
             		    window = 4,           # two days before and two days after
                            outlier.mins = 45,    # difference in mins
                            stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                            plot = FALSE)
}

twl_Ext <- Map(cbind, twl_Ext, threshold = 4)

## ------------------------------------------------------------------------
twlNew <- twls <- vector('list',(nBirds+1))

for(i in 1:(nBirds+1)){
  
if( i %in% c(1:3,5:15)){
newRows <- subset(twl_Ext[[i]],twl_Ext[[i]][,1] < twl[[i]][1,1])
 
twlNew[[i]] <- rbind(newRows,twl[[i]])

newRows <- subset(twl_Ext[[i]],twl_Ext[[i]][,1] > twl[[i]][nrow(twl[[i]]),1])

twlNew[[i]] <- rbind(twlNew[[i]], newRows) 
}

if(i == 4){
newRows <- subset(twl_Ext[[4]],twl_Ext[[4]][,1] < twl[[4]][1,1])

twlNew[[4]] <- rbind(newRows, twl[[4]])

newRows <- subset(twl_Ext[[4]], (twl_Ext[[4]][,1] > twl[[4]][nrow(twl[[4]]),1] & 
                                 twl_Ext[[4]][,1] < as.POSIXct("2015-06-21",tz = "GMT")))

twlNew[[4]] <- rbind(twlNew[[4]],newRows)
}

if(i == 16){
newRows <- subset(twl_Ext[[16]],(twl_Ext[[16]][,1] > as.POSIXct("2015-06-20",tz = "GMT") &
                                twl_Ext[[16]][,1] < twl[[16]][1,1]))

twlNew[[16]] <- rbind(newRows, twl[[16]])
}
}

## ------------------------------------------------------------------------
#twls <- vector('list',(nBirds+1))

#for(i in 1:(nBirds+1)){
#twls[[i]] <- twilightEdit(twilights = twlNew[[i]], 
#             window = 4,           # two days before and two days after
#             outlier.mins = 45,    # difference in mins
#             stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
#             plot = FALSE)

# Make sure to keep track of threshold
#twls[[i]]$Threshold <- twlNew[[i]]$threshold
#}

# ## ----exampleTWlEdit, echo = FALSE, fig.cap= "Red dots = Sunset, blue = Sunrise, arrows indicate the transtions that were moved and X represents transitions that were removed."----
# twlEdit2 <- twilightEdit(twilights = twl_ex, 
                    # window = 4,           # two days before and two days after
                    # outlier.mins = 45,    # difference in mins
                    # stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    # plot = TRUE)

# ## ----showdata, echo = FALSE----------------------------------------------
# head(twls[[1]])

## ----adjusttime----------------------------------------------------------
twlEdit <- lapply(X = twlNew, FUN = twilightAdjust, interval = 300)

twlEdit <- lapply(twlEdit,subset,Deleted == FALSE)

for(i in 1:(nBirds+1)){
lightImage(OSFLdata[[i]], zlim = c(0,5),offset = 19)
tsimagePoints(twlEdit[[i]]$Twilight, 
              offset = 12, 
              pch = 16, 
              cex = 0.5,
              col = ifelse(twlEdit[[i]]$Rise, "dodgerblue", "firebrick"))
}

## ----saveresults, eval = FALSE-------------------------------------------
 birds <- c(BirdId,paste0(BirdId[4],"_Yr2"))

## 
 for(i in 1:(nBirds+1)){
 write.csv(twlEdit[[i]], paste0("Twilights/",birds[[i]],"_Threshold_Combined.csv"), row.names = FALSE)
 }

## ----calibrartionDates---------------------------------------------------
# Create a vector with the dates known to be at deployment #
calibration.dates <- vector('list',3) # Using only the anchorage birds

for(i in 1:3){
calibration.dates[[i]] <- c(OSFLdata[[i]][1,1],as.POSIXct("2013-07-30",tz="GMT"))
}

## ----CalibrationData-----------------------------------------------------
calibration.data<-vector('list',3)

for(i in 1:3){
  calibration.data[[i]]<- subset(twl[[i]],
                                 subset = (twl[[i]]$Twilight>=calibration.dates[[i]][1] & 
                                           twl[[i]]$Twilight<=calibration.dates[[i]][2]))
}

## ----Zeniths-------------------------------------------------------------
# create empty vectors to store data #
sun<-z<-zenith0<-zenith1<-twl_t<-twl_deviation<-alpha<-fitml<-vector("list",3)

zenith0<-zenith1<-rep(NA,3)

# loop through each of the three individuals #
for(i in 1:3){
  
  # Calculate solar time from calibration data 
  sun[[i]]  <- solar(calibration.data[[i]][,1])
  
  # Adjust the solar zenith angle for atmospheric refraction
  z[[i]]<- refracted(zenith(sun = sun[[i]],
                            lon = CapLocs[i,1],
                            lat = CapLocs[i,2]))
  
  twl_t[[i]]   <- twilight(tm = calibration.data[[i]][,1], 
                           lon = CapLocs[i,1],
                           lat = CapLocs[i,2],
                           rise = calibration.data[[i]][,2],
                           zenith = quantile(z[[i]],probs=0.5))
  
  # Determine the difference in minutes from true twilight and the geolocator defined twilight
  twl_deviation[[i]] <- ifelse(calibration.data[[i]]$Rise, # If Rise
                    as.numeric(difftime(calibration.data[[i]][,1], twl_t[[i]], units = "mins")), # Do this
                    as.numeric(difftime(twl_t[[i]], calibration.data[[i]][,1], units = "mins"))) # If set, do this
  
  twl_deviation[[i]]<-subset(twl_deviation[[i]],twl_deviation[[i]]>=0)
  
  # Describe the distribution of the error 
  fitml[[i]] <- fitdistr(twl_deviation[[i]], "log-Normal")
  # save the Twilight model parameters
  alpha[[i]] <- c(fitml[[i]]$estimate[1], fitml[[i]]$estimate[2]) 
  
  zenith0[i] <-quantile(z[[i]],prob=0.5)
  zenith1[i]<-quantile(z[[i]],prob=0.95)
}

## ----plot1, echo = FALSE, fig.cap = "*Left* The deviation in twilights from the true twilight at the capture site. *Right* The mean (point estimate) and 95% CI for the zenith angle at the capture location."----
meanAlpha1<-mean(c(alpha[[1]][1],alpha[[2]][1],mean(alpha[[3]][1])))
meanAlpha2<-mean(c(alpha[[1]][2],alpha[[2]][2],mean(alpha[[3]][2])))
ALPHA<-alpha[[1]]
ALPHA[1]<-meanAlpha1
ALPHA[2]<-meanAlpha2


# b<-unlist(twl_deviation)
# cols<-c("red","blue","green")
# seq <- seq(0,60, length = 100)
# par(mfrow=c(1,2),mar=c(4,4,0,0))
# hist(b, freq = F,
     # yaxt="n",
     # ylim = c(0, 0.15),
     # xlim = c(0, 60),
     # breaks=15,
     # col="gray",
     # main = "",
     # xlab = "Twilight error (mins)")
# axis(2,las=2)
# for(i in 1:3){
# lines(seq, dlnorm(seq, alpha[[i]][1], alpha[[i]][2]), col = cols[i], lwd = 3, lty = 2)
# }

# #Zenith angle plot
# par(bty="l")
# plot(median(z[[1]],na.rm=TRUE),xlim=c(1,3),ylim=c(85,100),pch=19,ylab="Zenith Angle",xlab="OSFL",col=cols[1])
# segments(1,quantile(z[[1]],probs=0.025),1,quantile(z[[1]],probs=0.975),col=cols[1])
# for(i in 2:3){
  # par(new = TRUE)
  # plot(median(z[[i]],na.rm=TRUE)~i,xlim=c(1,3),ylim=c(85,100),pch=19,yaxt="n",ylab="",xlab="",col=cols[i])
  # segments(i,quantile(z[[i]],probs=0.025),i,quantile(z[[i]],probs=0.975),col=cols[i])
# }

## ----AncZ----------------------------------------------------------------
# use mean Anchorage sun angle for Fairbanks birds 
AncZ <- quantile( c(z[[1]], z[[2]], z[[3]]), probs=0.95)

## ----addZenith-----------------------------------------------------------
# Add the Zenith value to the twilight files
for(i in 1:(nBirds+1)){
  
# Create empty variable to hold the zenith angle
twlEdit[[i]]$Zenith <- NA

# assign the appropriate zenith angle 
twlEdit[[i]]$Zenith[which(twlEdit[[i]]$threshold == 1)] <- 95.76416

# calculation of zenith angle for threshold of 4 is not shown 
# but was conducted in the same manner as above.

twlEdit[[i]]$Zenith[which(twlEdit[[i]]$threshold == 4)] <- 93.05437
}

## ----scrubtwl, echo = FALSE----------------------------------------------
twlEdit[[1]] <- twlEdit[[1]][3:nrow(twlEdit[[1]]),]
twlEdit[[2]] <- twlEdit[[2]][3:nrow(twlEdit[[2]]),]
twlEdit[[3]] <- twlEdit[[3]][3:nrow(twlEdit[[3]]),]
twlEdit[[4]] <- twlEdit[[4]][1:(nrow(twlEdit[[4]])-8),]
twlEdit[[8]] <- twlEdit[[8]][1:724,]
twlEdit[[15]] <- twlEdit[[15]][3:(nrow(twlEdit[[15]])-4),]
twlEdit[[16]] <- twlEdit[[16]][7:nrow(twlEdit[[16]]),]

## ----setTols-------------------------------------------------------------
tolvalues <- array(NA,c((nBirds+1),2))
tolvalues[,1]<-0
tolvalues[,2]<-0.08

# Manual adjustments for Fall Equinox period 
tolvalues[1,1]<-0.13
tolvalues[2,1]<-0.101
tolvalues[3,1]<-0.2
tolvalues[4,1]<-0.17
tolvalues[5,1]<-0.22
tolvalues[6,1]<-0.2
tolvalues[7,1]<-0.168
tolvalues[8,1]<-0.195
tolvalues[9,1]<-0.105
tolvalues[10,1]<-0.2
tolvalues[11,1]<-0.115
tolvalues[12,1]<-0.13
tolvalues[13,1]<-0.165
tolvalues[14,1]<-0.185
tolvalues[15,1]<-0.215
tolvalues[16,1]<-0.23

# Manual adjustments for Spring Equinox period 
tolvalues[1,2]<-0.275
tolvalues[2,2]<-0.239
tolvalues[3,2]<-0.2
tolvalues[4,2]<-0.23
tolvalues[5,2]<-0.22
tolvalues[6,2]<-0.27
tolvalues[7,2]<-0.2285
tolvalues[8,2]<-0.2
tolvalues[9,2]<-0.2
tolvalues[10,2]<-0.162
tolvalues[11,2]<-0.188
tolvalues[12,2]<-0.17
tolvalues[13,2]<-0.22
tolvalues[14,2]<-0.24
tolvalues[15,2]<-0.2
tolvalues[16,2]<-0.22

## ----estPath-------------------------------------------------------------
path <- vector('list',nBirds+1)

for(i in 1:(nBirds+1)){  
   path[[i]] <- SGAT::thresholdPath(twilight = twlEdit[[i]]$Twilight,
                             rise = twlEdit[[i]]$Rise,
                             zenith = twlEdit[[i]]$Zenith,
                             tol = tolvalues[i,])
}

## ----pathsPlot, eval = FALSE, echo = FALSE-------------------------------
## ## ----echo = FALSE, fig.cap="**Figure 5** The initial annual cycle path of Olive-sided Flycatchers captured breeding in Alaska - *blue* = Fall, *green* = Spring, *red vertical lines* spring and fall equniox"----
for(i in 1:(nBirds+1)){
    layout(matrix(c(1,3,
                    2,3), 2, 2, byrow = TRUE))
 par(mar=c(2,4,2,0))
 plot(path[[i]]$time, path[[i]]$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '',xaxt="n")
 abline(h = ifelse(i != 16, CapLocs[i,2],CapLocs[4,2]))
 abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
 abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
 par(mar=c(2,4,2,0))
 plot(path[[i]]$time, path[[i]]$x[, 1], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '')
 abline(h = ifelse(i != 16, CapLocs[i,1],CapLocs[4,1]))
 abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
 abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
## 
## 
 plot(Americas, col = "grey95",xlim = c(-170,-60),ylim=c(0,65))
 box()
 lines(path[[i]]$x, col = "blue")
 points(path[[i]]$x, pch = 16, cex = 0.5, col = "blue")
Sys.sleep(2)
}

## ----intitalPaths--------------------------------------------------------
x0 <- z0 <- fixedx <- vector('list',nBirds+1)

for(i in 1:(nBirds+1)){
  # Take the location estimates created above
x0[[i]]<- path[[i]]$x

  # the model also needs the mid-points - generate those here
z0[[i]]<- trackMidpts(x0[[i]])

fixedx[[i]] <- rep(FALSE,nrow(x0[[i]]))
}

for(i in 1:3){
fixedx[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2013-07-30",format = "%Y-%m-%d",tz = "GMT"))] <- TRUE
x0[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2013-07-30",format = "%Y-%m-%d",tz = "GMT")),1] <- CapLocs[i,1]
x0[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2013-07-30",format = "%Y-%m-%d",tz = "GMT")),2] <- CapLocs[i,2]
}
for(i in 4:11){
fixedx[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2014-07-10",format = "%Y-%m-%d",tz = "GMT"))] <- TRUE
x0[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2014-07-10",format = "%Y-%m-%d",tz = "GMT")),1] <- CapLocs[i,1]
x0[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2014-07-10",format = "%Y-%m-%d",tz = "GMT")),2] <- CapLocs[i,2]
}
for(i in 12:nBirds){
fixedx[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2015-07-10",format = "%Y-%m-%d",tz = "GMT"))] <- TRUE
x0[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2015-07-10",format = "%Y-%m-%d",tz = "GMT")),1] <- CapLocs[i,1]
x0[[i]][which(twlEdit[[i]][,1] < as.POSIXct("2015-07-10",format = "%Y-%m-%d",tz = "GMT")),2] <- CapLocs[i,2]
}


## ----beta----------------------------------------------------------------
beta <- c(0.7, 0.08)

## ------------------------------------------------------------------------
xlim <- c(-170, -60)
ylim <- c(-90, 90)

## ------------------------------------------------------------------------
distribution.mask <- function(xlim, ylim, land = TRUE, shape) {
  
  r <- raster(res = c(0.25,0.25),
              xmn = xlim[1], 
              xmx = xlim[2], 
              ymn = ylim[1], 
              ymx = ylim[2], 
              crs = proj4string(shape))
  
  r <- cover(rasterize(shape, shift = c(-360, 0), r,1, silent = TRUE), 
             rasterize(shape, r, 1, silent = TRUE), rasterize(shape, r, 1, silent = TRUE))
  
  r <- as.matrix(is.na(r))[nrow(r):1, ]
  
  if (land) 
    r <- !r
  xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
  ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)
  
  function(p) {
    r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
  }
}

## ----distmask------------------------------------------------------------
is.dist <- distribution.mask(shape=Land,
                             xlim = xlim,
                             ylim = ylim,
                             land = TRUE)

## ----logprior------------------------------------------------------------
log.prior <- function(p) {
    f <- is.dist(p)
    ifelse(f | is.na(f), 0, -10)
}

## ----model---------------------------------------------------------------
model <-  vector('list', nBirds+1)

for(i in 1:(nBirds+1)){
model[[i]]<- thresholdModel(twilight = twlEdit[[i]]$Twilight,
                            rise = twlEdit[[i]]$Rise,
                            twilight.model = "ModifiedLogNormal",
                            alpha = ALPHA,
                            beta = beta,
                            # Here is where we set the constraints for land
                            logp.x = log.prior, 
                            logp.z = log.prior, 
                            x0 = x0[[i]],
                            z0 = z0[[i]],
                            fixedx = fixedx[[i]],
                            zenith = twlEdit[[i]]$Zenith)
}

# saveRDS(model,"OSFL_model.rds")
## ----error.def-----------------------------------------------------------
proposal.x <- proposal.z <- vector('list',nBirds+1)

for(i in 1:(nBirds+1)){
proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0[[i]]))
proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0[[i]]))
}

## ----eval = FALSE--------------------------------------------------------
  fit <- xsum <- zsum <- vector('list', nBirds+1)
 
  for(i in 1:(nBirds+1)){
  cat("\n", BirdId[i],"\n")
  fit[[i]] <- estelleMetropolis(model = model[[i]],
                                proposal.x = proposal.x[[i]],
                                proposal.z = proposal.z[[i]],
                                iters = 5000, # This value sets the number of iterations to run
                                thin = 5,
                                chains = 3)
## 
  xsum[[i]] <- locationSummary(fit[[i]]$x,collapse = TRUE)
  zsum[[i]] <- locationSummary(fit[[i]]$z,collapse = TRUE)
## 
 proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%')))
 proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%')))
## 
  fit[[i]] <- estelleMetropolis(model = model[[i]],
                               proposal.x = proposal.x[[i]],
                                proposal.z = proposal.z[[i]],
                                x0 = cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                                z0 = cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                                iters= 5000, # This value sets the number of iterations to run
                                thin= 5,
                                chains=3)
## 
## # Final Run
   xsum[[i]] <- locationSummary(fit[[i]]$x,collapse = TRUE)
   zsum[[i]] <- locationSummary(fit[[i]]$z,collapse = TRUE)
## 
  proposal.x[[i]] <- mvnorm(chainCov(fit[[i]]$x),s=0.1)
  proposal.z[[i]] <- mvnorm(chainCov(fit[[i]]$z),s=0.1)
## 
##  # Note the increase in number of interations - this takes a bit longer to run
 fit[[i]] <- estelleMetropolis(model = model[[i]],
                                proposal.x = proposal.x[[i]],
                                proposal.z = proposal.z[[i]],
                                x0=cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                                z0=cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                                iters=10000,  # This value sets the number of iterations to run
                                thin = 10,
                                chains=3)
 }
## 
## # Save the fit object #
## 
#saveRDS(fit,paste0("OSFL_fit_",format(Sys.Date(),"%b_%d_%Y"),".rds"))

#fit <- readRDS("OSFL_fit_Mar_07_2017.rds")


# This step makes an empty raster #
r <- raster(xmn=xlim[1],
            xmx=xlim[2],
            ymn=ylim[1],
            ymx=ylim[2],
            res = c(0.25,0.25))

## ------------------------------------------------------------------------
S <- Sp <- vector('list',(nBirds+1))

for(i in 1:(nBirds+1)){
S[[i]] <- slices(type="intermediate",
                 weights = rep(0.5,length(fit[[i]][[1]]$time)),
                 breaks="day",
                 mcmc=fit[[i]],
                 grid=r,
                 include.lowest = FALSE, 
                 right = FALSE)

Sp[[i]] <- slices(type="primary",
                 weights = rep(0.5,length(fit[[i]][[1]]$time)),
                 breaks="day",
                 mcmc=fit[[i]],
                 grid=r,
                 include.lowest = FALSE, 
                 right = FALSE)
}


## Migration Schedules 

May1 <- July <- Nov1<-Feb1<-rep(NA,(nBirds + 1))
for(i in 1:3){
May1[i] <- "2013-05-01"
July[i] <- "2013-07-30"
Nov1[i]<-"2013-11-01"
Feb1[i]<-"2014-02-01"
}
for(i in c(4:11)){
May1[i] <- "2014-05-01"
July[i] <- "2014-07-10"
Nov1[i]<-"2014-11-01"
Feb1[i]<-"2015-02-01"
}
for(i in c(12:(nBirds+1))){
May1[i] <- "2015-05-01"
July[i] <- "2015-07-10"
Nov1[i]<-"2015-11-01"
Feb1[i]<-"2016-02-01"
}
May1[16] <- "2015-06-20"


schedules <- vector('list',(nBirds+1))

detach("package:mth", unload = TRUE)
library(mth, lib.loc = "C:/Users/hallworthm/R_Library")

for(i in 1:(nBirds+1)){
schedules[[i]]<- mth::MigSchedule(MCMC = S[[i]], 
                             prob = 0.95, 
                             known.breed = c(May1[i],July[i]),
                             known.winter = c(Nov1[i],Feb1[i]),
                             plot = TRUE,
                             plot.legend = FALSE,
                             rm.lat.equinox = TRUE, 
                             days.omit = 5,
                             latAllow = c(0,5))
}

# saveRDS(schedules,paste0("Schedules_",format(Sys.Date(),"%b_%d_%Y"),".rds"))

# -- readSchedules --------------------------------------------------------------- 
schedules <- readRDS("Schedules_Mar_07_2017.rds")

# Write KML files #
#for(b in 1:(nBirds+1)){
#KML(schedules[[b]]$movements, 
#    filename = paste0(substr(BirdId[b],start = 1, stop = 16),".kml"),
#    col = rev(bpy.colors(100)),
#    alpha = 0.5,
#    overwrite = TRUE)
#}

## -- Further analyses ----------------------------------------------------------
# Read in protected areas in the Americas
ProtectedAm <- shapefile("Spatial_Layers/ProtectedAm.shp")

# Define ISO3 as Canada to subset only the protected areas in Canada
ISO3 <- unique(subset(Americas,SUBREGION == 13)$ISO3)

ProtectedCA <- ProtectedAm[which(ProtectedAm$ISO3 %in% ISO3),]

# Define ISO3 as North America to subset only the protected areas in North America
ISO3 <- unique(subset(Americas,SUBREGION == 21)$ISO3)

ProtectedNA <- ProtectedAm[which(ProtectedAm$ISO3 %in% ISO3),]

# Define ISO3 as South America to subset only the protected areas in south America
ISO3 <- unique(subset(Americas,SUBREGION == 5)$ISO3)

ProtectedSA <- ProtectedAm[which(ProtectedAm$ISO3 %in% ISO3),]


BIRDSyr <- BIRDSspring <- BIRDSfall <- vector('list',(nBirds+1))

for(i in 1:(nBirds+1)){
# move the movements into a new vector
BIRDSyr[[i]] <- schedules[[i]]$movements
BIRDSfall[[i]]<-schedules[[i]]$movements[[c(1:which(schedules[[i]]$Schedule$duration == (max(schedules[[i]]$Schedule$duration))))]]
BIRDSspring[[i]]<-schedules[[i]]$movements[[c(which(schedules[[i]]$Schedule$duration == (max(schedules[[i]]$Schedule$duration))):
                                                    nlayers(schedules[[i]]$movements))]]


# if value is non-zero set to 1. These are already 95% credible intervals 
BIRDSyr[[i]][BIRDSyr[[i]]>0]<- 1
BIRDSfall[[i]][BIRDSfall[[i]]>0]<-1
BIRDSspring[[i]][BIRDSspring[[i]]>0]<-1

# If value is na - set to 0
BIRDSyr[[i]][is.na(BIRDSyr[[i]])] <- 0
BIRDSyr[[i]] <- sum(BIRDSyr[[i]])
BIRDSyr[[i]][BIRDSyr[[i]]>=1]<-1

BIRDSfall[[i]][is.na(BIRDSfall[[i]])]<-0
BIRDSfall[[i]] <- sum(BIRDSfall[[i]])
BIRDSfall[[i]][BIRDSfall[[i]]>=1]<-1

BIRDSspring[[i]][is.na(BIRDSspring[[i]])]<-0
BIRDSspring[[i]] <- sum(BIRDSspring[[i]])
BIRDSspring[[i]][BIRDSspring[[i]]>=1]<-1
}

BirdsUse<-sum(stack(BIRDSyr[c(1:11,13,15)]))
BirdsUse[BirdsUse==0]<-NA

BirdsFall <- sum(stack(BIRDSfall[c(1:11,13,15)]))
BirdsFall[BirdsFall==0]<-NA

BirdsSpring<-sum(stack(BIRDSspring[c(1:11,13,15)]))
BirdsSpring[BirdsSpring==0]<-NA

# Protected Areas 

# get area in km2 of raster area #
a<-raster::area(BirdsUse)

protected <- mask(a,ProtectedAm)

birdusearea <- birdprotectedArea <- rep(NA,(nBirds+1))

### Calculate the area used by each bird ###
for(i in c(1:11,13,15)){
birdusearea[i] <- cellStats(BIRDSyr[[i]]*a,sum,na.rm = TRUE)
birdprotectedArea[i] <- cellStats(BIRDSyr[[i]]*protected,sum,na.rm = TRUE)
}


The mean area protected for each bird is `r round(mean((birdprotectedArea/birdusearea)*100,na.rm = TRUE),3)` &plusmn; `r round(sd((birdprotectedArea/birdusearea)*100,na.rm = TRUE)/sqrt(13),3)` (mean &plusmn; 1 SE)

### Protected areas by region and by season 


# mask BirdsUse area by Protected in the Americas 

#  mask BirdUse area by Region
NorthAm <- mask(a,subset(Americas,SUBREGION == 21))
CenAm <- mask(a,subset(Americas,SUBREGION == 13))
SouthAm <- mask(a,subset(Americas,SUBREGION == 5))

# protected areas by Region
protectedNA <- mask(protected,NorthAm)
protectedCenAm <- mask(protected,CenAm)
protectedSA <- mask(protected,SouthAm)

areaBirds <- protectedArea <- areaCenAm <- protectedCenAmArea <- rep(NA,13)
SpringNA <- SpringCenAm <- SpringSA <- rep(NA,13)
protectedSpringNA <- protectedSpringCA <- protectedSpringSA <- rep(NA,13)

FallNA <- FallCenAm <- FallSA <- rep(NA, 13)
protectedFallNA <- protectedFallCA <- protectedFallSA <- rep(NA,13)

areaSAbirds <- protectedAreaSA <- rep(NA,13)
areaNAbirds <- protectedAreaNA <- rep(NA, 13)

### Calculate the area used by all individuals birds ###

for(i in 1:13){
  
BirdsUse<-sum(stack(BIRDSyr[c(1:11,13,15)]))

BirdsUse[BirdsUse!=i]<-NA

areaBirds[i] <- cellStats(BirdsUse*a,sum)

protectedArea[i] <- cellStats(BirdsUse*protected,sum)

# Entire Year 

# North America
areaNAbirds[i] <- cellStats(BirdsUse*NorthAm,sum)
protectedAreaNA[i] <- cellStats(BirdsUse*protectedNA,sum)


# Central America
areaCenAm[i] <- cellStats(BirdsUse*CenAm,sum)
protectedCenAmArea[i] <- cellStats(BirdsUse*protectedCenAm,sum)

# South America 
areaSAbirds[i] <- cellStats(BirdsUse*SouthAm,sum)
protectedAreaSA[i] <- cellStats(BirdsUse*protectedSA,sum)

# Spring 

BirdsSpring<-sum(stack(BIRDSspring[c(1:11,13,15)]))

BirdsSpring[BirdsSpring!=i]<-NA

BirdsSpring[!is.na(BirdsSpring)]<-1

# North America
SpringNA[i]<-cellStats(BirdsSpring*NorthAm,sum)
protectedSpringNA[i] <- cellStats(BirdsSpring*protectedNA,sum)

# Central America
SpringCenAm[i] <- cellStats(BirdsSpring*CenAm,sum)
protectedSpringCA[i] <- cellStats(BirdsSpring*protectedCenAm, sum)

# South America 
SpringSA[i] <- cellStats(BirdsSpring*SouthAm,sum)
protectedSpringSA[i] <- cellStats(BirdsSpring*protectedSA,sum)

# Fall

BirdsFall<-sum(stack(BIRDSfall[c(1:11,13,15)]))

BirdsFall[BirdsFall!=i]<-NA

BirdsFall[!is.na(BirdsFall)]<-1

# North America
FallNA[i]<-cellStats(BirdsFall*NorthAm,sum)
protectedFallNA[i] <- cellStats(BirdsFall*protectedNA,sum)

# Central America
FallCenAm[i] <- cellStats(BirdsFall*CenAm,sum)
protectedFallCA[i] <- cellStats(BirdsFall*protectedCenAm, sum)

# South America 
FallSA[i] <- cellStats(BirdsFall*SouthAm,sum)
protectedFallSA[i] <- cellStats(BirdsFall*protectedSA,sum)

}

```
    
### North America - Spring    

```{r echo = FALSE}
BirdsSpring<-sum(stack(BIRDSspring[c(1:11,13,15)]))
BirdsSpring[BirdsSpring == 0]<-NA

par(mar = c(2,0,0,0))
plot(subset(Americas,SUBREGION == 21),col = "gray88",border = "gray50",ylim = c(30.32167,65.73015),xlim = c(-176.1937,-75.86459))
plot(Americas,add = TRUE,border = "gray50")
plot(ProtectedNA,col = "gray60",border = "gray60",add = TRUE)
plot(States,add = TRUE, border = "gray50")
plot(Canada,add = TRUE, border = "gray50")
plot(BirdsSpring,add = TRUE,alpha = 0.3,col = rev(bpy.colors(13)),legend = FALSE)
par(new = TRUE,fig = c(0,0.8,0,0.7),mar = c(4,4.5,4,4),bty = "l")
M<-barplot(SpringNA,col = rev(bpy.colors(13)),ylab = expression("Area 1000 km"^2),
     xlab = "Number of Birds",ylim = c(0,max(SpringNA)+60000),yaxt ="n",xlim = c(0,15),names.arg = c(1:13),axis.lty=1)
axis(2,las =2,at = seq(0,1400000,100000),labels = (seq(0,1400000,100000)/1000))
barplot(protectedSpringNA,add = TRUE,col = "gray60",ylab = "",xlab = "",ylim = c(0,max(SpringNA)+60000),yaxt ="n",xlim = c(0,15),
names.arg = c(rep("",10),11,"",13))
precent<-(protectedSpringNA/SpringNA)*100
precent[precent=="NaN"]<-0
text(x = M[1:13,1],y = SpringNA+35000, labels = paste0(round(precent,1),"%"),cex = 0.55)
```
    
### North America - Fall     


```{r echo = FALSE}
BirdsFall<-sum(stack(BIRDSfall[c(1:11,13,15)]))

BirdsFall[BirdsFall == 0]<-NA

par(mar = c(2,0,0,0))
plot(subset(Americas,SUBREGION == 21),col = "gray88",border = "gray50",ylim = c(30.32167,65.73015),xlim = c(-176.1937,-75.86459))
plot(Americas,add = TRUE,border = "gray50")
plot(ProtectedNA,col = "gray60",border = "gray60",add = TRUE)
plot(States,add = TRUE, border = "gray50")
plot(Canada,add = TRUE, border = "gray50")
plot(BirdsFall,add = TRUE,alpha = 0.3,col = rev(bpy.colors(13)),legend = FALSE)
par(new = TRUE,fig = c(0,0.8,0,0.7),mar = c(4,4.5,4,4),bty = "l")
M<-barplot(FallNA,col = rev(bpy.colors(13)),ylab = expression("Area 1000 km"^2),
     xlab = "Number of Birds",ylim = c(0,max(FallNA)+60000),yaxt ="n",xlim = c(0,15),names.arg = c(1:13),axis.lty=1)
axis(2,las =2,at = seq(0,1400000,100000),labels = (seq(0,1400000,100000)/1000))
barplot(protectedFallNA,add = TRUE,col = "gray60",ylab = "",xlab = "",ylim = c(0,max(FallNA)+60000),yaxt ="n",xlim = c(0,15),
names.arg = c(rep("",10),11,"",13))
precent<-(protectedFallNA/FallNA)*100
precent[precent=="NaN"]<-0
text(x = M[1:13,1],y = FallNA+35000, labels = paste0(round(precent,1),"%"),cex = 0.55)
```


### Central America
```{r echo = FALSE, fig.width = 6, fig.height = 6}
BirdsUse<-sum(stack(BIRDSyr[c(1:11,13,15)]))
BirdsUse[BirdsUse == 0]<-NA
 par(mar = c(0,0,0,0))
 plot(subset(Americas,SUBREGION == 13),col = "gray88",border = "gray50")
plot(Americas,add = TRUE,border = "gray50")
plot(ProtectedCA,col = "gray60",border = "gray60",add = TRUE)
plot(BirdsUse,add = TRUE,alpha = 0.3,col = rev(bpy.colors(13)),legend = FALSE)
par(new = TRUE,fig = c(0,1,0,0.7),mar = c(4,4.5,4,4),bty = "l")
M<-barplot(areaCenAm,col = rev(bpy.colors(13)),ylab = expression("Area 1000 km"^2),
     xlab = "Number of Birds",ylim = c(0,max(areaCenAm)+20000),yaxt ="n",xlim = c(0,15),names.arg = c(1:13),axis.lty=1)
 axis(2,las =2,at = seq(0,400000,40000),labels = (seq(0,400000,40000)/1000))
 barplot(protectedCenAmArea,add = TRUE,col = "gray60",ylab = "",xlab = "",ylim = c(0,max(areaCenAm)+20000),yaxt ="n",xlim = c(0,15),
 names.arg = c(rep("",10),11,"",13))
 text(x = M[1:13,1],y = areaCenAm+10000, labels = paste0(round((protectedCenAmArea/areaCenAm)*100,1),"%"),cex = 0.55)
```


### South America
```{r echo = FALSE, fig.width = 6, fig.height = 6}
par(mar = c(0,0,0,0))
plot(Americas, xlim = c(-95.06555,-65),ylim = c(-24.0403,9.32701),col = "gray88")
plot(protectedSA,add = TRUE,col = "gray60",legend = FALSE)
plot(BirdsUse,col =rev(bpy.colors(13)), alpha = 0.6,add = TRUE,legend = FALSE)
par(new = TRUE,fig = c(0,0.85,0,0.7),mar = c(4,4.5,4,4),bty = "l")
M<-barplot(areaSAbirds,col = rev(bpy.colors(13)),ylab = expression("Area 1000 km"^2),
     xlab = "Number of Birds",ylim = c(0,max(areaSAbirds)+40000),yaxt ="n",xlim = c(0,15),names.arg = c(1:13),axis.lty=1)
axis(2,las =2,at = seq(0,10800000,100000),labels = (seq(0,10800000,100000)/1000))
barplot(protectedAreaSA,add = TRUE,col = "gray60",ylab = "",xlab = "",ylim = c(0,max(protectedAreaSA)+20000),yaxt ="n",xlim = c(0,15),
names.arg = c(rep("",10),11,"",13))
text(x = M[1:13,1],y = areaSAbirds+30000, labels = paste0(round((protectedAreaSA/areaSAbirds)*100,1),"%"),cex = 0.55)
```

Second year breeding location of K288 
```{r eval = FALSE}
source("Functions/Lisovski_functions.R")
K288yr1<-subset(OSFLdata[[4]],Date < as.POSIXct("2014-07-15",format = "%Y-%m-%d",tz = "GMT"))

K288yr2<-subset(OSFLdata[[4]],subset(Date > as.POSIXct("2015-05-20",format = "%Y-%m-%d",tz = "GMT") & 
                                     Date < as.POSIXct("2015-07-01",format = "%Y-%m-%d",tz = "GMT")))

calib <- calibration(datetime = K288yr1$Date,
            light = log(K288yr1$Light),
            crds = CapLocs[4,], 
            max.light = NULL,
            aggr = 3, plot = TRUE) 

lon <- estim.lon(datetime = K288yr2$Date,
           light = log(K288yr2$Light))

lat <- estim.lat (datetime = K288yr2$Date,
           light =  log(K288yr2$Light), 
           calib = calib, 
           lonlim = c(-165,-135),
           latlim = c(55,75),
            plot = T,res=c(0.25,0.25)) 

breed2 <- projectRaster(lat[[1]],crs = NAEA)
breed2[breed2<0.7]<-NA
```


### Overlap during the winter
```{r}
# Overlapping 95% CI of all the birds during the winter 

birdwinter<-vector('list',(nBirds+1))

#for(i in c(1:11,13,15)){
for(i in 1:(nBirds+1)){
# Pull out the location where birds spent the longest time 
birdwinter[[i]]<- schedules[[i]]$movements[[which(schedules[[i]]$Schedule$duration == max(schedules[[i]]$Schedule$duration))]]

# Set that location equal to 1 #
birdwinter[[i]][!is.na(birdwinter[[i]])]<-1
#birdwinter[[i]]<-birdwinter[[i]]/cellStats(birdwinter[[i]],max,na.rm=TRUE)
}

# Sum up the birds to get a map for the number of birds
WinSum <- sum(stack(birdwinter[c(1:11,13,15)]),na.rm = TRUE)
```
      
      
Figure     
     
     
```{r , echo = FALSE, fig.width = 6, fig.height = 6, collapse = TRUE}
plot(Americas, xlim = c(-88.06555,-56.1876),ylim = c(-20.0403,11.32701),col = "gray88")
plot(WinSum,col = c("transparent",rev(bpy.colors(5))),add = TRUE,legend = FALSE)
plot(Americas,add = TRUE,border = "gray44")
plot(WinSum,useRaster = TRUE,legend.only = TRUE, horizontal = TRUE,col = c("transparent",rev(bpy.colors(5))),
       legend.width=0.5, legend.shrink=0.35,
       axis.args=list(at=pretty(0:5),
                    labels=pretty(0:5), 
                    cex.axis=1),
     legend.args=list(text='Number of birds', side=3, font=2, line=0, cex=1.25))
raster::scalebar(1000,xy = c(-90,-15),divs = 4, type = "bar", label = c(0,500,1000),below = "km") 

```
       

```{r echo = FALSE, fig.width = 6, fig.height = 6, eval = FALSE}
###### Map to show a "less pixelated" map of where they are in the winter - NOTE this is a nice visual but shouldn't be used for stats ######
wins<-vector('list',(nBirds+1))
for(i in 1:(nBirds+1)){

Dates <- as.character(strptime(S[[i]]$mcmc[[1]]$time[which(S[[i]]$mcmc[[1]]$rise == TRUE)],format = "%Y-%m-%d",tz = "GMT"))
wins[[i]] <- slice(S[[i]],
                 k = c(which(Dates==as.character(schedules[[i]]$Schedule[which(schedules[[i]]$Schedule$duration == max(schedules[[i]]$Schedule$duration)),1])):
                       which(Dates==as.character(schedules[[i]]$Schedule[which(schedules[[i]]$Schedule$duration == max(schedules[[i]]$Schedule$duration)),2]))))

wins[[i]]<-wins[[i]]/cellStats(wins[[i]],max, na.rm = TRUE)
}

WinSum <- sum(stack(wins[c(1:11,13,15)]),na.rm = TRUE)
plot(Americas, xlim = c(-88.06555,-56.1876),ylim = c(-20.0403,11.32701),col = "gray88")
plot(WinSum,col = c("transparent",rev(bpy.colors(50))),add = TRUE,legend = FALSE)
plot(Americas,add = TRUE,border = "gray44")
blankraster<-raster(nrow = 1,ncol = 5)
values(blankraster)<-0:4
plot(blankraster,useRaster = TRUE,legend.only = TRUE, horizontal = TRUE,col = c("transparent",rev(bpy.colors(50))),
       legend.width=0.5, legend.shrink=0.35,
       axis.args=list(at=c(0,1,2,3,4),
                    labels=c(0,1,2,3,4), 
                    cex.axis=1),
     legend.args=list(text='Number of birds', side=3, font=2, line=0, cex=1.25))
raster::scalebar(1000,xy = c(-90,-15),divs = 4, type = "bar", label = c(0,500,1000),below = "km") 
```


## Identifying important stop-over locations 
```{r echo = FALSE}
CapLocs <- rbind(CapLocs,c(CapLocs[4,]))
```

```{r stopovers}
# create empty vectors to store results
stops <- maxdays <- springstops <- fallstops <- vector("list",(nBirds+1))

for(i in 1:(nBirds+1)){
# store the movements in a new vector
stops[[i]]<-schedules[[i]]$movements

# Change the value of the raster to the stop duration
for(n in 1:nlayers(stops[[i]])){
stops[[i]][[n]][stops[[i]][[n]] > 0] <- as.numeric(as.character(schedules[[i]]$Schedule$duration[n]))

if(!is.na(extract(stops[[i]][[n]], SpatialPoints(cbind(CapLocs[i,1],CapLocs[i,2]))))){
stops[[i]][[n]][stops[[i]][[n]]] <- 0
 }
}

# create a single raster with the stop durations - excluding winter! 
# Fall durations
fallstops[[i]] <- sum(stack(stops[[i]][[c(1:(which(schedules[[i]]$Schedule$duration == (max(schedules[[i]]$Schedule$duration)))-1))]]),na.rm = TRUE)


# Spring durations
springstops[[i]] <- sum(stack(stops[[i]][[c((which(schedules[[i]]$Schedule$duration == (max(schedules[[i]]$Schedule$duration)))+1):
                                              nlayers(schedules[[i]]$movements))]]),na.rm = TRUE)

maxdays[[i]] <- max(stack(stops[[i]]),na.rm = TRUE)
}

# Determine the number of birds that use each stop
# Spring
birdsspring <- springstops
birdsspring <- stack(birdsspring[c(1:11,13,15)])
birdsspring[birdsspring>0]<-1

# Fall 
birdsfall<- fallstops
birdsfall <- stack(fallstops[c(1:11,13,15)])

# birdsfall[birdsfall>0]<-1
```

### Calculate Importance of each stop-over area
$$Importance = \frac{\left(\sum(Stop Duration) * proportion(BirdsUsed) \right )}{max(\sum(Stop Duration) *  proportion(BirdsUsed))}$$

SpringImportance <- (sum(stack(springstops))*(sum(birdsspring)/max(birdsspring))) / 
                    cellStats(sum(stack(springstops))*(sum(birdsspring)/max(birdsspring)),max)			
					
FallImportance <- (sum(stack(fallstops))*(sum(birdsfall)/max(birdsfall))) / 
                  cellStats(sum(stack(fallstops))*(sum(birdsfall)/max(birdsfall)),max)

     



# subset the fall importance raster into fall regions #
# named from North to South #
FallRegion1 <- crop(FallImportance,extent(c(-133,-125,57,61.5)))   # Northern Canada
FallRegion2 <- crop(FallImportance,extent(c(-123,-111,47,53)))     # BC/Washington
FallRegion3 <- crop(FallImportance,extent(c(-104,-96,25.5,32.5)))  # Texas
FallRegion4 <- crop(FallImportance,extent(c(-100,-88,15.3,20.75))) # Mexico 
FallRegion5 <- crop(FallImportance,extent(c(-89,-83.5,11.8,15.2))) # Nicaruaga/Hondoras
FallRegion6 <- crop(FallImportance,extent(c(-86,-74.5,5.3,11.55))) # Panama 
FallRegion7 <- crop(FallImportance,extent(c(-78,-74.5,-3.5,1)))    # Ecuador

FallRegion1[!is.na(FallRegion1)] <- 1
FallRegion2[!is.na(FallRegion2)] <- 2
FallRegion3[!is.na(FallRegion3)] <- 3
FallRegion4[!is.na(FallRegion4)] <- 4
FallRegion5[!is.na(FallRegion5)] <- 5
FallRegion6[!is.na(FallRegion6)] <- 6
FallRegion7[!is.na(FallRegion7)] <- 7

FallRegions <- sum(stack(extend(FallRegion1,extent(FallImportance),value = NA),
      			extend(FallRegion2,extent(FallImportance),value = NA),
      			extend(FallRegion3,extent(FallImportance),value = NA),
      			extend(FallRegion4,extent(FallImportance),value = NA),
      			extend(FallRegion5,extent(FallImportance),value = NA),
      			extend(FallRegion6,extent(FallImportance),value = NA),
      			extend(FallRegion7,extent(FallImportance),value = NA)),na.rm = TRUE)

FallRegions[FallRegions == 0]<-NA

### When moving through important areas 

FallRegion1pass <- FallRegion2pass <- FallRegion3pass <- FallRegion4pass <- FallRegion5pass <- FallRegion6pass <- FallRegion7pass <- array(NA,c(16,54))
FallRegion1Dates <- FallRegion2Dates <- FallRegion3Dates <- FallRegion4Dates <- FallRegion5Dates <- FallRegion6Dates <- FallRegion7Dates <- array(NA,c((nBirds+1),4))

winterNum <- rep(NA,(nBirds+1))

for(i in 1:(nBirds+1)){
# store the movements in a new vector
stops[[i]]<-schedules[[i]]$movements

# Change the value of the raster to the stop duration
for(n in 1:nlayers(stops[[i]])){
stops[[i]][[n]][stops[[i]][[n]] > 0] <- as.numeric(as.character(schedules[[i]]$Schedule$duration[n]))
FallRegion1pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(FallRegion1,extent(FallImportance),value = NA)),sum)
FallRegion2pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(FallRegion2,extent(FallImportance),value = NA)),sum)
FallRegion3pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(FallRegion3,extent(FallImportance),value = NA)),sum)
FallRegion4pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(FallRegion4,extent(FallImportance),value = NA)),sum)
FallRegion5pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(FallRegion5,extent(FallImportance),value = NA)),sum)
FallRegion6pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(FallRegion6,extent(FallImportance),value = NA)),sum)
FallRegion7pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(FallRegion7,extent(FallImportance),value = NA)),sum)
}


# Identify when non-breeding winter starts
winterNum[i] <- which(schedules[[i]]$Schedule$duration == (max(schedules[[i]]$Schedule$duration)))

FallRegion1Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(FallRegion1pass[i,1:winterNum[i]] > 0)),1]),"%j"))
FallRegion1Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(FallRegion1pass[i,1:winterNum[i]] > 0)),2]),"%j"))

FallRegion2Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(FallRegion2pass[i,1:winterNum[i]] > 0)),1]),"%j"))
FallRegion2Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(FallRegion2pass[i,1:winterNum[i]] > 0)),2]),"%j"))

FallRegion3Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(FallRegion3pass[i,1:winterNum[i]] > 0)),1]),"%j"))
FallRegion3Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(FallRegion3pass[i,1:winterNum[i]] > 0)),2]),"%j"))

FallRegion4Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(FallRegion4pass[i,1:winterNum[i]] > 0)),1]),"%j"))
FallRegion4Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(FallRegion4pass[i,1:winterNum[i]] > 0)),2]),"%j"))

FallRegion5Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(FallRegion5pass[i,1:winterNum[i]] > 0)),1]),"%j"))
FallRegion5Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(FallRegion5pass[i,1:winterNum[i]] > 0)),2]),"%j"))

FallRegion6Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(FallRegion6pass[i,1:winterNum[i]] > 0)),1]),"%j"))
FallRegion6Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(FallRegion6pass[i,1:winterNum[i]] > 0)),2]),"%j"))

FallRegion7Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(FallRegion7pass[i,1:winterNum[i]] > 0)),1]),"%j"))
FallRegion7Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(FallRegion7pass[i,1:winterNum[i]] > 0)),2]),"%j"))
}

# Important Areas during the fall
    
FallImportance <- (sum(stack(fallstops))*(sum(birdsfall)/max(birdsfall))) / 
                  cellStats(sum(stack(fallstops))*(sum(birdsfall)/max(birdsfall)),max)

tiff("Figures/FallImportance_stops.tiff",res = 600, width = 3600, height = 3600)
par(mar = c(0,0,0,0))
plot(Americas, xlim = c(-178.6541,-78), ylim= c(0,70))
plot(FallImportance, add = TRUE, col = rev(gray.colors(100, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL)), legend = FALSE)
plot(Canada, add = TRUE, border = "gray50")
plot(States, add = TRUE, border = "gray50")
plot(Americas,add = TRUE,border = "gray50")
FallImportance[FallImportance<0.1]<-NA
plot(FallImportance, add = TRUE, col = rev(bpy.colors(100)), legend = FALSE)
rect(-133,57,-125,61.5, lty = 1,lwd = 2) # Region1
rect(-123,47,-111,53, lty = 1,lwd = 2) # Region2
rect(-104,27.5,-96,32.5, lty = 1,lwd = 2) # Region3
rect(-100,15.3,-88,20.75, lty = 1,lwd = 2) # Region4
rect(-89,11.8,-83.5,15.2, lty = 1,lwd = 2) # Region5
rect(-86,5.3,-74.5,11.55, lty = 1,lwd = 2) # Region6
rect(-78,-1.5,-74.5,1, lty = 1,lwd = 2) # Region7
blankraster<-raster(nrow = 1,ncol = 10)
values(blankraster)<-seq(0.1,1,0.1)
plot(blankraster,useRaster = TRUE,legend.only = TRUE, horizontal = TRUE,col = rev(bpy.colors(100)),
       legend.width=0.25, legend.shrink=0.35,
       smallplot=c(0.2,0.625,0.09,0.1),
       axis.args=list(at=seq(0.1,1,0.1),
                    labels=seq(0.1,1,0.1), 
                    cex.axis=0.9),
     legend.args=list(text='Importance', side=3, font=2, line=0, cex=1))

par(bty = "l", fig = c(0,0.7,0.125,0.65),new = TRUE,mar = c(4,4,0,0))
plot(NA, xlim = c(200,300),ylim = c(0.5,7.5),xlab = "",xaxt = "n", yaxt = "n",ylab = "")
axis(1,labels = seq(200,300,5), at = seq(200,300,5))
axis(2,las = 2, at = 1:7, labels = 7:1)
mtext("Ordinal Day", 1, at = 250, line = 1.75)
mtext("Fall Stop-over Region",2, at = 4, line = 1.75)
segments(apply(FallRegion1Dates,2,min,na.rm = TRUE)[1],7,apply(FallRegion1Dates,2,max,na.rm = TRUE)[2],7, col = "black", lwd = 1, lty = 2)
rect(apply(FallRegion1Dates,2,mean,na.rm = TRUE)[1],6.9,apply(FallRegion1Dates,2,mean,na.rm = TRUE)[2],7.1, col = "gray", lwd = 1)

segments(apply(FallRegion2Dates,2,min,na.rm = TRUE)[1],6,apply(FallRegion2Dates,2,max,na.rm = TRUE)[2],6, col = "black", lwd = 1, lty = 2)
rect(apply(FallRegion2Dates,2,mean,na.rm = TRUE)[1],5.9,apply(FallRegion2Dates,2,mean,na.rm = TRUE)[2],6.1, col = "gray", lwd = 1)

segments(apply(FallRegion3Dates,2,min,na.rm = TRUE)[1],5,apply(FallRegion3Dates,2,max,na.rm = TRUE)[2],5, col = "black", lwd = 1, lty = 2)
rect(apply(FallRegion3Dates,2,mean,na.rm = TRUE)[1],4.9,apply(FallRegion3Dates,2,mean,na.rm = TRUE)[2],5.1, col = "gray", lwd = 1)

segments(apply(FallRegion4Dates,2,min,na.rm = TRUE)[1],4,apply(FallRegion4Dates,2,max,na.rm = TRUE)[2],4, col = "black", lwd = 1, lty = 2)
rect(apply(FallRegion4Dates,2,mean,na.rm = TRUE)[1],3.9,apply(FallRegion4Dates,2,mean,na.rm = TRUE)[2],4.1, col = "gray", lwd = 1)

segments(apply(FallRegion5Dates,2,min,na.rm = TRUE)[1],3,apply(FallRegion5Dates,2,max,na.rm = TRUE)[2],3, col = "black", lwd = 1, lty = 2)
rect(apply(FallRegion5Dates,2,mean,na.rm = TRUE)[1],2.9,apply(FallRegion5Dates,2,mean,na.rm = TRUE)[2],3.1, col = "gray", lwd = 1)

# These start to cross the year mark - most likely birds wintering in those areas #
FallRegion6Dates[which(FallRegion6Dates[,2]<200),] <- NA
segments(apply(FallRegion6Dates,2,min,na.rm = TRUE)[1],2,apply(FallRegion6Dates,2,max,na.rm = TRUE)[2],2, col = "black", lwd = 1, lty = 2)
rect(apply(FallRegion6Dates,2,mean,na.rm = TRUE)[1],1.9,apply(FallRegion6Dates,2,mean,na.rm = TRUE)[2],2.1, col = "gray", lwd = 1)

# These start to cross the year mark - most likely birds wintering in those areas #
FallRegion7Dates[which(FallRegion7Dates[,2]<200),] <- NA
segments(apply(FallRegion7Dates,2,min,na.rm = TRUE)[1],1,apply(FallRegion7Dates,2,max,na.rm = TRUE)[2],1, col = "black", lwd = 1, lty = 2)
rect(apply(FallRegion7Dates,2,mean,na.rm = TRUE)[1],0.9,apply(FallRegion6Dates,2,mean,na.rm = TRUE)[2],1.1, col = "gray", lwd = 1)
dev.off()




# subset the spring importance raster into spring regions #
# named from South to North #
SpringRegion1 <- crop(SpringImportance,extent(c(-80,-74,-5.5,2.5)))   # Ecuador
SpringRegion2 <- crop(SpringImportance,extent(c(-83.5,-72,2.5,10.1))) # Panama/Colombia
SpringRegion3 <- crop(SpringImportance,extent(c(-87,-81,10,14.7)))    # Nicaragua
SpringRegion4 <- crop(SpringImportance,extent(c(-100,-88,14,23)))     # Mexico 
SpringRegion5 <- crop(SpringImportance,extent(c(-124,-118,39.5,48.31)))# P. NW 
SpringRegion6 <- crop(SpringImportance,extent(c(-126,-118,48.5,53))) # BC

SpringRegion1[!is.na(SpringRegion1)] <- 8
SpringRegion2[!is.na(SpringRegion2)] <- 9
SpringRegion3[!is.na(SpringRegion3)] <- 10
SpringRegion4[!is.na(SpringRegion4)] <- 11
SpringRegion5[!is.na(SpringRegion5)] <- 12
SpringRegion6[!is.na(SpringRegion6)] <- 13


SpringRegions <- sum(stack(extend(SpringRegion1,extent(SpringImportance),value = NA),
      			extend(SpringRegion2,extent(SpringImportance),value = NA),
      			extend(SpringRegion3,extent(SpringImportance),value = NA),
      			extend(SpringRegion4,extent(SpringImportance),value = NA),
      			extend(SpringRegion5,extent(SpringImportance),value = NA),
      			extend(SpringRegion6,extent(SpringImportance),value = NA)),na.rm = TRUE)

SpringRegions[SpringRegions == 0]<-NA

SpringRegions[SpringRegions > 1] <- 1

SpringAreaImpArea <- raster::area(SpringRegions)

PercentSpringProtect <- cellStats(protected*SpringRegions,sum)/cellStats(SpringAreaImpArea*SpringRegions,sum)

par(mar = c(0,0,0,0),bty = "n")
plot(SpringRegion6, col = "gray", legend = FALSE, axes = FALSE)
plot(Americas, add = TRUE)
plot(ProtectedAm, col = "forestgreen", border = "forestgreen",add = TRUE)
plot(SpringRegion6, col = rgb(190/255,190/255,190/255,100/255),legend = FALSE, axes = FALSE,add = TRUE)

### When moving through important areas 

SpringRegion1pass <- SpringRegion2pass <- SpringRegion3pass <- SpringRegion4pass <- SpringRegion5pass <- SpringRegion6pass  <- array(NA,c(16,54))
SpringRegion1Dates <- SpringRegion2Dates <- SpringRegion3Dates <- SpringRegion4Dates <- SpringRegion5Dates <- SpringRegion6Dates <- array(NA,c((nBirds+1),2))

winterNum <- rep(NA,(nBirds+1))

for(i in 1:(nBirds+1)){
# store the movements in a new vector
stops[[i]]<-schedules[[i]]$movements

# Change the value of the raster to the stop duration
for(n in 1:nlayers(stops[[i]])){
stops[[i]][[n]][stops[[i]][[n]] > 0] <- as.numeric(as.character(schedules[[i]]$Schedule$duration[n]))
SpringRegion1pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(SpringRegion1,extent(SpringImportance),value = NA)),sum)
SpringRegion2pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(SpringRegion2,extent(SpringImportance),value = NA)),sum)
SpringRegion3pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(SpringRegion3,extent(SpringImportance),value = NA)),sum)
SpringRegion4pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(SpringRegion4,extent(SpringImportance),value = NA)),sum)
SpringRegion5pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(SpringRegion5,extent(SpringImportance),value = NA)),sum)
SpringRegion6pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(SpringRegion6,extent(SpringImportance),value = NA)),sum)
}


# Identify when non-breeding winter starts
winterNum[i] <- which(schedules[[i]]$Schedule$duration == (max(schedules[[i]]$Schedule$duration)))

SpringRegion1Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(SpringRegion1pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
SpringRegion1Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(SpringRegion1pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

SpringRegion2Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(SpringRegion2pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
SpringRegion2Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(SpringRegion2pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

SpringRegion3Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(SpringRegion3pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
SpringRegion3Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(SpringRegion3pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

SpringRegion4Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(SpringRegion4pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
SpringRegion4Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(SpringRegion4pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

SpringRegion5Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(SpringRegion5pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
SpringRegion5Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(SpringRegion5pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

SpringRegion6Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(SpringRegion6pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
SpringRegion6Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(SpringRegion6pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

}

# Important areas in the Spring    
SpringImportance <- (sum(stack(springstops))*(sum(birdsspring)/max(birdsspring))) / 
                    cellStats(sum(stack(springstops))*(sum(birdsspring)/max(birdsspring)),max)	
		
tiff("Figures/SpringImportance_stops.tiff",res = 600, width = 3600, height = 3600)	
par(mar = c(0,0,0,0))
plot(Americas, xlim = c(-178.6541,-76), ylim= c(-10,70))
plot(SpringImportance, add = TRUE, legend = FALSE, col = rev(gray.colors(100, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL)))
plot(Canada, add = TRUE, border = "gray50")
plot(States, add = TRUE, border = "gray50")
plot(Americas,add = TRUE,border = "gray50")
SpringImportance[SpringImportance < 0.1] <- NA
rect(-80,-5.5,-74,2.5, lty = 1,lwd = 2) # Region1
rect(-83.5,2.5,-72,10.1, lty = 1,lwd = 2) # Region2
rect(-87,10,-81,14.7, lty = 1,lwd = 2) # Region3
rect(-100,14,-88,23, lty = 1,lwd = 2) # Region4
rect(-124,39.5,-118,48.31, lty = 1,lwd = 2) # Region5
rect(-126,48.5,-118,53, lty = 1,lwd = 2) # Region6
plot(SpringImportance, add = TRUE, legend = FALSE, col = rev(bpy.colors(100)))
plot(blankraster,useRaster = TRUE,legend.only = TRUE, horizontal = TRUE,col = rev(bpy.colors(100)),
       legend.width=0.25, legend.shrink=0.35,
       smallplot=c(0.125,0.57,0.09,0.1),
       axis.args=list(at=seq(0.1,1,0.1),
                    labels=seq(0.1,1,0.1), 
                    cex.axis=0.9),
     legend.args=list(text='Importance', side=3, font=2, line=0, cex=1))

par(bty = "l", fig = c(0,0.6,0.1,0.65),new = TRUE,mar = c(4,3,0,0))
plot(NA, xlim = c(35,165),ylim = c(0.5,6.5),xlab = "",xaxt = "n", yaxt = "n",ylab = "")
axis(1,labels = seq(35,165,5), at = seq(35,165,5))
axis(2,las = 2, at = 1:6, labels = 8:13)
mtext("Ordinal Day", 1, at = 100, line = 1.75)
mtext("Spring Stop-over Region",2, at = 3.5, line = 2)

# These start to cross the year mark - most likely birds wintering in those areas #
SpringRegion1Dates[which(SpringRegion1Dates[,1] > 200),] <- NA

segments(apply(SpringRegion1Dates,2,min,na.rm = TRUE)[1],1,apply(SpringRegion1Dates,2,max,na.rm = TRUE)[2],1, col = "black", lwd = 1, lty = 2)
rect(apply(SpringRegion1Dates,2,mean,na.rm = TRUE)[1],0.9,apply(SpringRegion1Dates,2,mean,na.rm = TRUE)[2],1.1, col = "gray", lwd = 1)

# These start to cross the year mark - most likely birds wintering in those areas #
SpringRegion2Dates[which(SpringRegion2Dates[,1] > 200),] <- NA
segments(apply(SpringRegion2Dates,2,min,na.rm = TRUE)[1],2,apply(SpringRegion2Dates,2,max,na.rm = TRUE)[2],2, col = "black", lwd = 1, lty = 2)
rect(apply(SpringRegion2Dates,2,mean,na.rm = TRUE)[1],1.9,apply(SpringRegion2Dates,2,mean,na.rm = TRUE)[2],2.1, col = "gray", lwd = 1)

segments(apply(SpringRegion3Dates,2,min,na.rm = TRUE)[1],3,apply(SpringRegion3Dates,2,max,na.rm = TRUE)[2],3, col = "black", lwd = 1, lty = 2)
rect(apply(SpringRegion3Dates,2,mean,na.rm = TRUE)[1],2.9,apply(SpringRegion3Dates,2,mean,na.rm = TRUE)[2],3.1, col = "gray", lwd = 1)

segments(apply(SpringRegion4Dates,2,min,na.rm = TRUE)[1],4,apply(SpringRegion4Dates,2,max,na.rm = TRUE)[2],4, col = "black", lwd = 1, lty = 2)
rect(apply(SpringRegion4Dates,2,mean,na.rm = TRUE)[1],3.9,apply(SpringRegion4Dates,2,mean,na.rm = TRUE)[2],4.1, col = "gray", lwd = 1)

segments(apply(SpringRegion5Dates,2,min,na.rm = TRUE)[1],5,apply(SpringRegion5Dates,2,max,na.rm = TRUE)[2],5, col = "black", lwd = 1, lty = 2)
rect(apply(SpringRegion5Dates,2,mean,na.rm = TRUE)[1],4.9,apply(SpringRegion5Dates,2,mean,na.rm = TRUE)[2],5.1, col = "gray", lwd = 1)

segments(apply(SpringRegion6Dates,2,min,na.rm = TRUE)[1],6,apply(SpringRegion6Dates,2,max,na.rm = TRUE)[2],6, col = "black", lwd = 1, lty = 2)
rect(apply(SpringRegion6Dates,2,mean,na.rm = TRUE)[1],5.9,apply(SpringRegion6Dates,2,mean,na.rm = TRUE)[2],6.1, col = "gray", lwd = 1)
dev.off()


# Get maximum number of 'stops'
nlayersMoves <- rep(NA, nBirds+1)
for(i in 1:16){
nlayersMoves[i] <- nlayers(schedules[[i]]$movements)
}

StopRegion <- BirdsUsedRegion <- array(0,c(max(nlayersMoves),16))

for(i in 1:16){
for(n in 1:winterNum[i]){
schedules[[i]]$movements[[n]][!is.na(schedules[[i]]$movements[[n]])] <- 1
schedules[[i]]$movements[[n]][schedules[[i]]$movements[[n]] == 0] <- NA
StopRegion[n,i] <- zonal(FallRegions,schedules[[i]]$movements[[n]],na.rm = TRUE, fun = max)[2]
BirdsUsedRegion[n,i] <- zonal(BirdsFall,schedules[[i]]$movements[[n]],na.rm = TRUE, fun = max)[2]
}
for(n in (winterNum[i]+1):nlayersMoves[i]){
schedules[[i]]$movements[[n]][!is.na(schedules[[i]]$movements[[n]])] <- 1
schedules[[i]]$movements[[n]][schedules[[i]]$movements[[n]] == 0] <- NA
StopRegion[n,i] <- zonal(SpringRegions,schedules[[i]]$movements[[n]],na.rm = TRUE, fun = max)[2]
BirdsUsedRegion[n,i] <- zonal(BirdsSpring,schedules[[i]]$movements[[n]],na.rm = TRUE, fun = max)[2]
}
}


StopRegion[StopRegion == "-Inf"] <- 0
BirdsUsedRegion[BirdsUsedRegion == "-Inf"] <- 1

#Durations <- vector('list',nBirds+1)
#for(i in 1:(nBirds+1)){
#Durations[[i]] <- cbind(schedules[[i]]$Schedule,StopRegion = StopRegion[1:nlayersMoves[i],i], BirdsUsed = BirdsUsedRegion[1:nlayersMoves[i],i])
#write.csv(Durations[[i]],paste0("Schedule_",BirdId[i],".csv"),row.names = FALSE)
#}


### Important spring and fall - Importance values > 0.1 

SpringImportance[SpringImportance < 0.1] <- NA
FallImportance[FallImportance < 0.1] <- NA

SpringFallImport <- mask(FallImportance,SpringImportance)

par(mar = c(0,0,0,0))
plot(Americas, xlim = c(-178.6541,-66.02323), ylim= c(-10,70))
plot(SpringFallImport, add = TRUE, legend = FALSE, col = rev(gray.colors(100, start = 0.1, end = 0.9, gamma = 2.2, alpha = NULL)))
plot(Canada, add = TRUE, border = "gray50")
plot(States, add = TRUE, border = "gray50")
plot(Americas,add = TRUE,border = "gray50")
SpringFallImport[SpringFallImport > 0 ] <- 1
plot(SpringFallImport, add = TRUE, legend = FALSE, col = "firebrick")

Region1 <- crop(SpringFallImport,extent(c(-100,-87,15,21)))
Region2 <- crop(SpringFallImport,extent(c(-87,-82,10,15)))
Region3 <- crop(SpringFallImport,extent(c(-84,-73,4,10)))
Region4 <- crop(SpringFallImport,extent(c(-78,-74,-2,1.5)))

par(mar = c(0,0,0,0))
plot(subset(Americas, SUBREGION == 13), xlim = c(-130,-70.02323), ylim= c(-3,20), col = "gray88",border = "gray50")
plot(Americas, add = TRUE, col = "gray88")
plot(SpringFallImport, add = TRUE, legend = FALSE, col = "firebrick")
rect(-100,15,-87,21, lty = 1) # Region1
rect(-87,10,-82,15, lty = 1) # Region2
rect(-84,4,-73,10, lty = 1) # Region3
rect(-78,-2,-74,1.5, lty = 1) # Region4
plot(Canada, add = TRUE, border = "gray50")
plot(States, add = TRUE, border = "gray50")
plot(Americas,add = TRUE,border = "gray50")
 

### When moving through important areas 

Region1pass <- Region2pass <- Region3pass <- Region4pass <- array(NA,c(16,54))
Region1Dates <- Region2Dates <- Region3Dates <- Region4Dates <- array(NA,c((nBirds+1),4))
winterNum <- rep(NA,(nBirds+1))

for(i in 1:(nBirds+1)){
# store the movements in a new vector
stops[[i]]<-schedules[[i]]$movements

# Change the value of the raster to the stop duration
for(n in 1:nlayers(stops[[i]])){
stops[[i]][[n]][stops[[i]][[n]] > 0] <- as.numeric(as.character(schedules[[i]]$Schedule$duration[n]))
Region1pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(Region1,extent(SpringFallImport),value = NA)),sum)
Region2pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(Region2,extent(SpringFallImport),value = NA)),sum)
Region3pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(Region3,extent(SpringFallImport),value = NA)),sum)
Region4pass[i,n] <- cellStats(mask(stops[[i]][[n]],extend(Region4,extent(SpringFallImport),value = NA)),sum)
}

# Identify when non-breeding winter starts
winterNum[i] <- which(schedules[[i]]$Schedule$duration == (max(schedules[[i]]$Schedule$duration)))

Region1Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region1pass[i,1:winterNum[i]] > 0)),1]),"%j"))
Region1Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region1pass[i,1:winterNum[i]] > 0)),2]),"%j"))

Region1Dates[i,3]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region1pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
Region1Dates[i,4]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region1pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

Region2Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region2pass[i,1:winterNum[i]] > 0)),1]),"%j"))
Region2Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region2pass[i,1:winterNum[i]] > 0)),2]),"%j"))

Region2Dates[i,3]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region2pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
Region2Dates[i,4]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region2pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

Region3Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region3pass[i,1:winterNum[i]] > 0)),1]),"%j"))
Region3Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region3pass[i,1:winterNum[i]] > 0)),2]),"%j"))

Region3Dates[i,3]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region3pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
Region3Dates[i,4]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region3pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))

Region4Dates[i,1]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region4pass[i,1:winterNum[i]] > 0)),1]),"%j"))
Region4Dates[i,2]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region4pass[i,1:winterNum[i]] > 0)),2]),"%j"))

Region4Dates[i,3]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[min(which(Region4pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),1]),"%j"))
Region4Dates[i,4]<-as.numeric(format(as.Date(schedules[[i]]$Schedule[max(which(Region4pass[i,winterNum[i]:54] > 0))+(winterNum[i]-1),2]),"%j"))
}


par(bty = "l")
plot(NA, xlim = c(225,500),ylim = c(0.5,4.5),xlab = "Ordinal Day",xaxt = "n", yaxt = "n",ylab = "Stop-over Region")
axis(1,labels = c(seq(225,365,10),seq(10,140,10)), at = seq(225,505,10))
axis(2,las = 2, at = 1:4, labels = 1:4)
segments(apply(Region1Dates,2,min,na.rm = TRUE)[1],1,apply(Region1Dates,2,max,na.rm = TRUE)[2],1, col = "black", lwd = 1, lty = 2)
rect(apply(Region1Dates,2,mean,na.rm = TRUE)[1],0.9,apply(Region1Dates,2,mean,na.rm = TRUE)[2],1.1, col = "gray", lwd = 1)
segments(apply(Region1Dates,2,min,na.rm = TRUE)[3]+365,1,apply(Region1Dates,2,max,na.rm = TRUE)[4]+365,1, col = "black", lwd = 1, lty = 2)
rect(apply(Region1Dates,2,mean,na.rm = TRUE)[3]+365,0.9,apply(Region1Dates,2,mean,na.rm = TRUE)[4]+365,1.1, col = "gray", lwd = 1)

segments(apply(Region2Dates,2,min,na.rm = TRUE)[1],2,apply(Region2Dates,2,max,na.rm = TRUE)[2],2, col = "black", lwd = 1, lty = 2)
rect(apply(Region2Dates,2,mean,na.rm = TRUE)[1],1.9,apply(Region2Dates,2,mean,na.rm = TRUE)[2],2.1, col = "gray", lwd = 1)
segments(apply(Region2Dates,2,min,na.rm = TRUE)[3]+365,2,apply(Region2Dates,2,max,na.rm = TRUE)[4]+365,2, col = "black", lwd = 1, lty = 2)
rect(apply(Region2Dates,2,mean,na.rm = TRUE)[3]+365,1.9,apply(Region2Dates,2,mean,na.rm = TRUE)[4]+365,2.1, col = "gray", lwd = 1)

segments(apply(Region3Dates[c(2,4:6,8,11:16),],2,min,na.rm = TRUE)[1],3,
         apply(Region3Dates[c(2,4:6,8,11:16),],2,max,na.rm = TRUE)[2],3, col = "black", lwd = 1, lty = 2)
rect(apply(Region3Dates[c(2,4:6,8,11:16),],2,mean,na.rm = TRUE)[1],2.9,apply(Region3Dates[c(2,4:6,8,11:16),],2,mean,na.rm = TRUE)[2],3.1, col = "gray", lwd = 1)
segments(apply(Region3Dates[c(2,4:6,8,11:16),],2,min,na.rm = TRUE)[3]+365,3,
         apply(Region3Dates[c(2,4:6,8,11:16),],2,max,na.rm = TRUE)[4]+365,3, col = "black", lwd = 1, lty = 2)
rect(apply(Region3Dates[c(2,4:6,8,11:16),],2,mean,na.rm = TRUE)[3]+365,2.9,
     apply(Region3Dates[c(2,4:6,8,11:16),],2,mean,na.rm = TRUE)[4]+365,3.1, col = "gray", lwd = 1)


segments(apply(Region4Dates[c(2,4:6,8,11,13:16),],2,min,na.rm = TRUE)[1],4,
         apply(Region4Dates[c(2,4:6,8,11,13:16),],2,max,na.rm = TRUE)[2],4, col = "black", lwd = 1, lty = 2)
rect(apply(Region4Dates[c(2,4:6,8,11,13:16),],2,mean,na.rm = TRUE)[1],3.9,
     apply(Region3Dates[c(2,4:6,8,11:16),],2,mean,na.rm = TRUE)[2],4.1, col = "gray", lwd = 1)
segments(apply(Region4Dates[c(2,4:6,8,11,13:16),],2,min,na.rm = TRUE)[3]+365,4,
         apply(Region4Dates[c(2,4:6,8,11,13:16),],2,max,na.rm = TRUE)[4]+365,4, col = "black", lwd = 1, lty = 2)
rect(apply(Region4Dates[c(2,4:6,8,11,13:16),],2,mean,na.rm = TRUE)[3]+365,3.9,
     apply(Region4Dates[c(2,4:6,8,11,13:16),],2,mean,na.rm = TRUE)[4]+365,4.1, col = "gray", lwd = 1)
```

############################### COLOMBIA ######################################################

Colombia <- subset(Americas, NAME == "Colombia")
COcover <- raster("Spatial_Layers/COL_msk_cov.grd")
waterlines <- shapefile("Spatial_Layers/COL_water_lines_dcw.shp")
waterarea <- shapefile("Spatial_Layers/COL_water_areas_dcw.shp")
roads <- shapefile("Spatial_Layers/COL_roads.shp")

tiff("Figures/Colombia_BirdsUse.tiff", res = 600, width = 3600, height = 3600)
par(mar = c(0,0,0,0))
plot(Colombia)
#plot(waterlines, add = TRUE, col = "blue")
plot(waterarea, add = TRUE, col = "blue",border = "blue")
plot(COcover, add = TRUE,alpha = 0.5,legend = FALSE)
plot(Colombia, add = TRUE)
plot(Americas, add = TRUE, border = "gray88")
BirdsUseCrop <- crop(BirdsUse,extent(Colombia))
plot(BirdsUseCrop, add = TRUE, col = rev(bpy.colors(13)), breaks = seq(1,13,1),alpha = 0.3,legend = FALSE)
plot(BirdsUseCrop,useRaster = TRUE,legend.only = TRUE, horizontal = TRUE,col = rev(bpy.colors(13)),
       legend.width=0.25, legend.shrink=0.35,
       smallplot=c(0.6,0.9,0.8,0.8125),
       axis.args=list(at=seq(1,13,1),
                    labels=seq(1,13,1), 
                    cex.axis=0.9),
     legend.args=list(text='Birds', side=3, font=2, line=0, cex=1))
dev.off()




#Code for creating GIFs (not run here)
```{r eval = FALSE}
DATES <- vector('list',(nBirds+1))

library(animation)

# this needs to be set for your own machine. Note - ImageMagick needs to be installed
ani.options(convert="C:/Program Files/ImageMagick-6.8.9-Q8/convert.exe",interval=0.3)

for(i in 1:(nBirds+1)){ 
# Store the dates for each bird 
DATES[[i]] <- S[[i]]$mcmc[[1]]$time[S[[i]]$mcmc[[1]]$rise == TRUE]

# Here is the wrapper function to save the gifs
saveGIF({

for(d in 1:length(DATES[[i]])){
  
tm<-sliceInterval(S[[i]],k=d)

sk<-slice(S[[i]],k=d)

cols <- colorRampPalette(c("gray50","red"))

# plot only the 95% confidence interval
sk[sk<quantile(sk,probs = 0.95)]<-NA

#sk[!is.na(sk)]<-1

  plot(Americas,col="gray74",border = "gray88",
       ylim=c(-20,70),xlim=c(-160,-60))
  plot(sk,
       ylim=c(-20,70),xlim=c(-160,-60),
       axes=FALSE, add=TRUE,
       legend=FALSE,
       col=cols(20),
       cex.axis=0.7)
  plot(Americas,border="gray88",add=TRUE)
  plot(Canada, border = "gray88",add = TRUE)
  plot(States, border = "gray88",add = TRUE)
  if(is.null(tm[1])){
  legend(-120,70,legend="",bty = "n")}
  else{
  legend(-120,70,legend=paste(strptime(tm[1],format = "%Y-%m-%d")),cex=1.25,bty="n")
  if(tm[1]> as.POSIXct("2013-09-04") & tm[1]<as.POSIXct("2013-10-05")){
  legend(-120,60,legend="Equinox",bty="n",cex=1.25)}
  if(tm[1]> as.POSIXct("2014-03-04") & tm[1]<as.POSIXct("2014-04-10")){
  legend(-120,60,legend="Equinox",bty="n",cex=1.25)}
  if(tm[1]> as.POSIXct("2014-09-04") & tm[1]<as.POSIXct("2014-10-05")){
  legend(-120,60,legend="Equinox",bty="n",cex=1.25)}
  if(tm[1]> as.POSIXct("2015-03-04") & tm[1]<as.POSIXct("2015-04-10")){
  legend(-120,60,legend="Equinox",bty="n",cex=1.25)}
  if(tm[1]> as.POSIXct("2015-09-04") & tm[1]<as.POSIXct("2015-10-05")){
  legend(-120,60,legend="Equinox",bty="n",cex=1.25)}
  if(tm[1]> as.POSIXct("2016-03-04") & tm[1]<as.POSIXct("2016-04-10")){
  legend(-120,60,legend="Equinox",bty="n",cex=1.25)}
}
}
}, movie.name=paste0(birdnames[i],".gif"))
}
```
