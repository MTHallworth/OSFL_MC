## ----eval = FALSE--------------------------------------------------------
## # Check to make sure the required packages are installed on your machine
## # If not, they will be installed
## 
## reqPackages <- c("devtools","raster","sp","maptools","rgeos","MASS")
## get.packages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]
## if(length(get.packages)>0) install.packages(get.packages)
## 
## # Install necessary packages from Github using the devtools library #
## library(devtools)
## install_github("SWotherspoon/SGAT")
## install_github("SLisovski/TwGeos")
## install_github("eldarrak/FLightR")

## ---- warning = FALSE, message = FALSE-----------------------------------
library(raster)
library(sp)
library(rgeos)
library(geosphere)
library(SGAT)
library(TwGeos)
library(MASS)
library(maptools)
library(FLightR)

## ------------------------------------------------------------------------
# read in a simple world map from the maptools package #
Americas<-raster::shapefile("Spatial_Layers/Americas.shp")
OSFLDist<-raster::shapefile("Spatial_Layers/OliveSidedFlycatcher.shp")
States <- raster::shapefile("Spatial_Layers/st99_d00.shp")
Canada <- raster::shapefile("Spatial_Layers/Canada.shp")

WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Canada <- spTransform(Canada,CRS = WGS84)


WinterFiles <- list.files("Data/WinterCalibration", pattern = "driftadj.lux", full.names = TRUE)
WinterData <- lapply(WinterFiles, readMTlux)

WinterCalibDates <- vector('list',length(WinterFiles))

for(i in 1:length(WinterFiles)){
lightImage(WinterData[[i]],
           offset = 19)
WinterCalibDates[[i]] <- as.POSIXct(locator(n = 2)$x, origin = "1970-01-01",tz = "GMT")
dev.off()
}

#[[Nicaragua K649]]
#[1] "2015-10-27 21:44:35 GMT" "2016-05-10 10:16:43 GMT"

#[[Guatemala Q327]]
#[1] "2015-10-14 19:27:30 GMT" "2016-05-16 19:31:43 GMT"

#[[Ecuador Q349]]
#[1] "2015-10-17 22:26:40 GMT" "2016-05-04 14:06:43 GMT"

#[[Ecuador Q353]]
#[1] "2015-10-19 16:06:40 GMT" "2016-01-20 23:00:01 GMT"

for(i in 1:length(WinterFiles)){
WinterData[[i]] <- subset(WinterData[[i]],Date > WinterCalibDates[[i]][1] &
                                          Date < WinterCalibDates[[i]][2])
}

WinterDeploy <- array(NA,c(length(WinterFiles),2))
WinterDeploy[1,] <- cbind(-86.511, 13.233) #K649 Nicaragua
WinterDeploy[2,] <- cbind(-91.508, 14.864) #Q327 Guatemala 
WinterDeploy[3,] <- cbind(-77.740, -0.683) #Q349 Ecuador
WinterDeploy[4,] <- cbind(-78.891, 0.022) #Q353 Ecuador

for(i in 1:length(WinterData)){
WinterData[[i]]$Light <- log(WinterData[[i]]$Light)
}

seed <- as.POSIXct("2016-01-01 04:00:00", format = "%Y-%m-%d %H:%M:%S",tz = "GMT")

twlWinter <- vector('list',length(WinterFiles))

for(i in 1:length(WinterFiles)){
twlWinter[[i]] <- findTwilights(tagdata = WinterData[[i]],
                           include = seed,
                           threshold = 1)
}

for(i in 1:length(WinterFiles)){

twlWinter[[i]] <- twilightEdit(twilights = twlWinter[[i]], 
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)
}
for(i in 1:4){
twlWinter[[i]] <-  subset(twlWinter[[i]],Deleted == FALSE)
}

lightImage(WinterData[[4]],offset = 16, zlim = c(0,10))
## ------------------------------------------------------------------------
# Create a vector with the dates known to be at deployment #

WinterSun<-WinterZ<-WinterZenith <-WinterTWL<-WinterDev<-WinterAlpha<-Winterfitml<-vector("list",4)

# loop through each of the two individuals #
for(i in 1:4){
  
  # Calculate solar time from calibration data 
  WinterSun[[i]]  <- solar(twlWinter[[i]][,1])
  
  # Adjust the solar zenith angle for atmospheric refraction
  WinterZ[[i]]<- refracted(zenith(sun = WinterSun[[i]],
                            lon = WinterDeploy[i,1],
                            lat = WinterDeploy[i,2]))
  
  WinterTWL[[i]]   <- twilight(tm = twlWinter[[i]][,1], 
                            lon = WinterDeploy[i,1],
                            lat = WinterDeploy[i,2],
                           rise = twlWinter[[i]][,2],
                           zenith = quantile(WinterZ[[i]],probs=0.5))
  
  # Determine the difference in minutes from when the sun rose and the geolocator said it rose 
  WinterDev[[i]] <- ifelse(twlWinter[[i]]$Rise, as.numeric(difftime(twlWinter[[i]][,1], WinterTWL[[i]], units = "mins")),
                         as.numeric(difftime(WinterTWL[[i]], twlWinter[[i]][,1], units = "mins")))
  WinterDev[[i]][WinterDev[[i]]<0]<-0.0001

  
  # Describe the distribution of the error 
  Winterfitml[[i]] <- fitdistr(WinterDev[[i]], "log-Normal")
  # save the Twilight model parameters
  WinterAlpha[[i]] <- c(Winterfitml[[i]]$estimate[1], Winterfitml[[i]]$estimate[2]) 
}

b<-unlist(WinterDev)
cols<-c("red","blue","green","yellow","orange","purple","brown","gray","black","pink")
seq <- seq(0,60, length = 100)
par(mfrow=c(1,2),mar=c(4,4,0,0))
hist(b, freq = F,
     yaxt="n",
     ylim = c(0, 0.15),
     xlim = c(0, 60),
     breaks=30,
     col="gray",
     main = "",
     xlab = "Twilight error (mins)")
axis(2,las=2)
for(i in 1:4){
lines(seq, dlnorm(seq, WinterAlpha[[i]][1], WinterAlpha[[i]][2]), col = cols[i], lwd = 3, lty = 2)
}

#Zenith angle plot
par(bty="l")
plot(median(WinterZ[[1]],na.rm=TRUE),xlim=c(1,10),ylim=c(85,100),pch=19,ylab="Zenith Angle",xlab="OSFL",col=cols[1])
segments(1,quantile(WinterZ[[1]],probs=0.025),1,quantile(WinterZ[[1]],probs=0.975),col=cols[1])
for(i in 2:4){
  par(new = TRUE)
  plot(median(WinterZ[[i]],na.rm=TRUE)~i,xlim=c(1,10),ylim=c(85,100),pch=19,yaxt="n",xaxt="n",ylab="",xlab="",col=cols[i])
  segments(i,quantile(WinterZ[[i]],probs=0.025),i,quantile(WinterZ[[i]],probs=0.975),col=cols[i])
}

WinterDeploy<-SpatialPoints(WinterDeploy)

qWZ<-lapply(WinterZ,quantile,probs = 0.5)

AEA <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

###################################################################################################
## Supplemental Plot ##
#crs(WinterDeploy)<-WGS84
#tiff(filename = "WinterZ.tiff",
#     width = 2400, height = 2400, units = "px", pointsize = 12,res = 600,
#     compression = "lzw",bg = "white")
#par(mar = c(0,0,0,0))
#plot(spTransform(WinterDeploy,CRS = AEA),xlim =c(201815.8,2718975),ylim =c(-4583537,-2534242))
#plot(spTransform(Americas,CRS = AEA),add = TRUE,col = "gray88",lwd = 0.25)
#for(i in 1:4){
#legend(spTransform(WinterDeploy,CRS = AEA)@coords[i,1],spTransform(WinterDeploy,CRS = AEA)@coords[i,2],
#       legend = round(qWZ[[i]],dig = 3),yjust = 0.5,x.intersp = 0, y.intersp = 0,cex = 0.7,
#       bty = 'n')
#}
#plot(spTransform(WinterDeploy,CRS = AEA),pch = 19, add = TRUE)

#par(fig=c(0,0.5,0,0.5), new=TRUE)
#plot(spTransform(OSFLDist,CRS = AEA),bg = "white",lwd = 0.25)
#plot(spTransform(Americas,CRS = AEA), add = TRUE,col = "gray88",lwd = 0.25)
#plot(spTransform(States,CRS = AEA), add = TRUE, lwd = 0.25)
#plot(spTransform(Canada,CRS = AEA), add = TRUE, lwd = 0.25)
#polygon(c(231815.8,231815.8,2518975,2518975),c(-4583537,-2534242,-2534242,-4583537),lwd = 2)
#box()
#dev.off()
#################################################################################################

## ------------------------------------------------------------------------
AncFiles <- list.files("Data/Anchorage_recoveries",pattern = ".lux",full.names=TRUE)
FairFiles <- list.files("Data/Fairbanks_recoveries",pattern = "driftadj.lux",full.names=TRUE)
TetlinFiles <- list.files("Data/Tetlin_recoveries",pattern = ".lux", full.names = TRUE)

# Read just the file names for Bird ID
AncNames <- list.files("Data/Anchorage_recoveries", pattern=".lux")
FairNames <- list.files("Data/Fairbanks_recoveries", pattern="driftadj.lux")
TetlinNames <- list.files("Data/Tetlin_recoveries", pattern = ".lux")

#Combine the files
OSFLFiles<-c(AncFiles,FairFiles,TetlinFiles)

# Combine the bird ID
BirdId<-c(AncNames,FairNames,TetlinNames)

# Determine the number of birds
nBirds<-length(BirdId)

## ------------------------------------------------------------------------

# Loop through all the files and read them in as LUX files #

OSFLdata <- lapply(OSFLFiles,readMTlux) 

## ------------------------------------------------------------------------
# Set the capture coordinates for each bird #
CapLocs<-array(NA,c(nBirds,2))

CapLocs[1,] <- c(-149.60425,61.38084)    # F202
CapLocs[2,] <- c(-149.74919,61.28071)    # F203
CapLocs[3,] <- c(-149.80309,61.29737)    # F311
CapLocs[4,] <- c(-146.832960,65.335090)  # K288 These are assumed
CapLocs[5,] <- c(-147.1002,64.71122)     # K648
CapLocs[6,] <- c(-146.800030,65.33139)   # K654
CapLocs[7,] <- c(-146.943180,65.34078)   # K657
CapLocs[8,] <- c(-142.30559,63.13592)    # K659
CapLocs[9,] <- c(-146.832960,65.335090)  # K660
CapLocs[10,] <- c(-147.1002,64.71122)    # K665
CapLocs[11,] <- c(-147.758300,64.936900) # K670
CapLocs[12,] <- c(-146.832960,65.335090) # S050 These are assumed
CapLocs[13,] <- c(-146.832960,65.335090) # S052 These are assumed
CapLocs[14,] <- c(-141.833333,62.666667) # Q330
CapLocs[15,] <- c(-141.833333,62.666667) # Q334

## ---- echo=FALSE,fig.cap="**Figure 2** Capture locations within Alaska"----
states<-raster::shapefile("Spatial_Layers/st99_d00.shp")
crs(states)<-WGS84
AK<-subset(states,NAME=="Alaska")

#plot(spTransform(AK,CRS = AEA),col="lightgray",border="black")
#plot(spTransform(Americas,CRS = AEA),add=TRUE,col="lightgray",border="gray")
#plot(spTransform(sp::SpatialPoints(cbind(CapLocs),proj= CRS(WGS84)),CRS =AEA),pch=19,add=TRUE)
#plot(spTransform(AK,CRS = AEA),border="black",add = TRUE)
#par(fig=c(0,0.5,0,0.5), new=TRUE)
#plot(spTransform(OSFLDist,CRS = AEA),bg = "white",lwd = 0.25)
#plot(spTransform(Americas,CRS = AEA), add = TRUE,col = "gray88",lwd = 0.25)
#plot(spTransform(States,CRS = AEA), add = TRUE, lwd = 0.25)
#plot(spTransform(Canada,CRS = AEA), add = TRUE, lwd = 0.25)
#box()
#points(spTransform(sp::SpatialPoints(CapLocs,proj = CRS(WGS84)),CRS = AEA),pch = 19, col = "red")

for(i in 1:nBirds){
OSFLdata[[i]]$Light<-log(OSFLdata[[i]]$Light)
}

twl <- vector('list',nBirds)

seed <- as.POSIXct(c(rep("2014-01-01 04:00:00",3), 
                     rep("2015-01-01 04:00:00",8),
                     rep("2016-01-01 04:00:00",4)),format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

for(i in 1:nBirds){
twl[[i]] <- findTwilights(tagdata = OSFLdata[[i]],
                        threshold = 1,
                        include = seed[i],
                        dark.min = 0) # 0 hours minimum dark period
}

twl[[4]] <- findTwilights(tagdata = OSFLdata[[4]],
                        threshold = 1,
                        include = c(seed[4],seed[nBirds]),
                        dark.min = 0) # 0 hours minimum dark period

for(i in 1:nBirds){
twl[[i]] <- preprocessLight(tagdata = OSFLdata[[i]],
                            offset = 19,
                            twilights = twl[[i]],
                            threshold = 1,
                            zlim = c(0,3))
}


### ---- Edit twl ---- ###
twlEdit <- vector('list',nBirds+1)
str(twlEdit)
for(i in c(1:3,5:nBirds)){
twlEdit[[i]] <- twilightEdit(twilights = twl[[i]], 
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)
}

twlEdit[[4]] <- twilightEdit(twilights = twl[[4]][1:594,], 
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)
twlEdit[[16]]<- twilightEdit(twilights = twl[[4]][595:nrow(twl[[4]]),], 
                    window = 4,           # two days before and two days after
                    outlier.mins = 45,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)



## ----eval = FALSE--------------------------------------------------------
 for(i in 1:nBirds){
   twlEdit[[i]]<-twilightAdjust(twilights=twlEdit[[i]], interval=300)
 }

twlEdit <- lapply(twlEdit,subset,Deleted == FALSE)

str(twlEdit[[4]])
## ----eval = FALSE--------------------------------------------------------
## for(i in 1:nBirds){
##   saveRDS(twl[[i]],paste0(BirdId[[i]],".rds")) # Saves rds file with the BirdId.rds as the name
## }

twlEdit <- readRDS("OSFL_twlEdit.rds")
## ------------------------------------------------------------------------
# Create a vector with the dates known to be at deployment #
calibration.dates <- vector('list',nBirds)

for(i in 1:3){
calibration.dates[[i]] <- c(OSFLdata[[i]][1,1],as.POSIXct("2013-07-30",tz="GMT"))
}
for(i in 4:11){
calibration.dates[[i]] <- c(OSFLdata[[i]][1,1],as.POSIXct("2014-07-10",tz="GMT"))
}
for(i in 12:nBirds){
calibration.dates[[i]] <- c(OSFLdata[[i]][1,1],as.POSIXct("2015-07-10",tz="GMT"))
}

## ------------------------------------------------------------------------
# Extract twilight data during calibration period
calibration.data<-vector('list',nBirds)
str(calibration.data)
for(i in 1:3){
  calibration.data[[i]]<-subset(twlEdit[[i]],twlEdit[[i]]$Twilight>=calibration.dates[[i]][1] & 
                                             twlEdit[[i]]$Twilight<=calibration.dates[[i]][2])
}

## ------------------------------------------------------------------------
# create empty vectors to store data #
sun<-z<-zenith0<-zenith1<-twl_t<-twl_deviation<-alpha<-fitml<-vector("list",3)

# loop through each of the two individuals #
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
  
  # Determine the difference in minutes from when the sun rose and the geolocator said it rose 
  twl_deviation[[i]] <- ifelse(calibration.data[[i]]$Rise, as.numeric(difftime(calibration.data[[i]][,1], twl_t[[i]], units = "mins")),
                         as.numeric(difftime(twl_t[[i]], calibration.data[[i]][,1], units = "mins")))
  
  twl_deviation[[i]]<-subset(twl_deviation[[i]],twl_deviation[[i]]>=0)
  
  # Describe the distribution of the error 
  fitml[[i]] <- fitdistr(twl_deviation[[i]], "log-Normal")
  # save the Twilight model parameters
  alpha[[i]] <- c(fitml[[i]]$estimate[1], fitml[[i]]$estimate[2]) 
}

## ----echo=FALSE----------------------------------------------------------
meanAlpha1<-mean(c(alpha[[1]][1],alpha[[2]][1],mean(alpha[[3]][1])))
meanAlpha2<-mean(c(alpha[[1]][2],alpha[[2]][2],mean(alpha[[3]][2])))

ALPHA<-alpha[[1]]
ALPHA[1]<-meanAlpha1
ALPHA[2]<-meanAlpha2

## ---- echo = FALSE, fig.cap = ("**Figure 4** *Left* The deviation in twilights from the true twilight at the capture site. *Right* The mean (point estimate) and 95% CI for the zenith angle at the capture location.")----

b<-unlist(twl_deviation)
cols<-c("red","blue","green","yellow","orange","purple","brown","gray","black","pink")
seq <- seq(0,60, length = 100)
par(mfrow=c(1,2),mar=c(4,4,0,0))
hist(b, freq = F,
     yaxt="n",
     ylim = c(0, 0.15),
     xlim = c(0, 60),
     breaks=15,
     col="gray",
     main = "",
     xlab = "Twilight error (mins)")
axis(2,las=2)
for(i in 1:nBirds){
lines(seq, dlnorm(seq, alpha[[i]][1], alpha[[i]][2]), col = cols[i], lwd = 3, lty = 2)
}

#Zenith angle plot
par(bty="l")
plot(median(z[[1]],na.rm=TRUE),xlim=c(1,10),ylim=c(85,100),pch=19,ylab="Zenith Angle",xlab="OSFL",col=cols[1])
segments(1,quantile(z[[1]],probs=0.025),1,quantile(z[[1]],probs=0.975),col=cols[1])
for(i in 2:nBirds){
  par(new = TRUE)
  plot(median(z[[i]],na.rm=TRUE)~i,xlim=c(1,10),ylim=c(85,100),pch=19,yaxt="n",xaxt="n",ylab="",xlab="",col=cols[i])
  segments(i,quantile(z[[i]],probs=0.025),i,quantile(z[[i]],probs=0.975),col=cols[i])
}

## ------------------------------------------------------------------------
# Create empty vectors to store objects #
d.twl<-path<-vector('list',nBirds)

zenith0<-zenith1<-rep(NA,nBirds)

# loop through the birds #
for(i in 1:nBirds){
  # Store the zenith (sun-elevation angle)
  zenith0[i] <-quantile(z[[i]],prob=0.5)
  zenith1[i]<-quantile(z[[i]],prob=0.95)
}

## ----echo = FALSE--------------------------------------------------------
zenith0

## ------------------------------------------------------------------------
# use mean Anchorage sun angle for Fairbanks birds 
AncZ<-quantile(c(z[[1]],z[[2]],z[[3]]),probs=0.95)

## ----echo=FALSE----------------------------------------------------------
AncZ
endBreed <- startBreed <- vector('list',nBirds)

for(i in 1:3){
endBreed[[i]] <- as.POSIXct("2013-11-01",tz="GMT")
startBreed[[i]] <- as.POSIXct("2014-03-01",tz="GMT")
}
for(i in 5:11){
endBreed[[i]] <- as.POSIXct("2014-11-01",tz="GMT")
startBreed[[i]] <- as.POSIXct("2015-03-01",tz="GMT")
}
for(i in 12:nBirds){
endBreed[[i]] <- as.POSIXct("2015-11-01",tz="GMT")
startBreed[[i]] <- as.POSIXct("2016-03-01",tz="GMT")
}

zenithAngles <- Dates <- vector('list',nBirds)
for(i in 1:nBirds){
Dates[[i]]<-twlEdit[[i]][,1]
zenithAngles[[i]]<-rep(AncZ,nrow(twlEdit[[i]]))
zenithAngles[[i]][which(Dates[[i]] > endBreed[[i]] & Dates[[i]] < startBreed[[i]])] <- quantile(unlist(WinterZ),probs = 0.90)
if(i == 4) {
zenithAngles[[i]][which(Dates[[i]] > as.POSIXct("2015-11-01",tz="GMT") & Dates[[i]] < as.POSIXct("2016-03-01",tz="GMT"))] <- quantile(unlist(WinterZ),probs = 0.90)
}
}

## ------------------------------------------------------------------------
path <- vector('list',nBirds)
for(i in 1:nBirds){  
   path[[i]] <- thresholdPath(twilight = twlEdit[[i]]$Twilight,
                             rise = twlEdit[[i]]$Rise,
                             zenith = AncZ,
                             tol = c(0,0.12))
}
str(path)

K228path1 <- thresholdPath(twilight = twlEdit[[4]]$Twilight[1:612],
                           rise = twlEdit[[4]]$Rise[1:612],
                           zenith = AncZ,
                           tol = c(0,0.12))

K228path2 <- thresholdPath(twilight = twlEdit[[4]]$Twilight[613:length(twlEdit[[4]]$Twilight)],
                           rise = twlEdit[[4]]$Rise[613:length(twlEdit[[4]]$Twilight)],
                           zenith = AncZ,
                           tol = c(0,0.12))

## ----echo = FALSE, fig.cap="**Figure 5** The initial annual cycle path of Olive-sided Flycatchers captured breeding in Alaska - *blue* = Fall, *green* = Spring, *red vertical lines* spring and fall equniox"----

for(i in 1:nBirds){
  print(BirdId[[i]])
  layout(matrix(c(1,3,
                  2,3), 2, 2, byrow = TRUE))
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '',xaxt="n")
  abline(h = CapLocs[i,2])
  abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 1], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '')
  abline(h = CapLocs[i,1])
  abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
  
  
  plot(Americas, col = "grey95",xlim = c(-170,-60),ylim=c(0,65))
  box()
  lines(path[[i]]$x, col = "blue")
  points(path[[i]]$x, pch = 16, cex = 0.5, col = "blue")
Sys.sleep(3)
}


  layout(matrix(c(1,2,4,5,
                  3,3,6,6), 2, 4, byrow = TRUE))
  par(mar=c(2,4,2,0))
  plot(K228path1$time, K228path1$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '',xaxt="n")
  abline(h = CapLocs[i,2])
  abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
  par(mar=c(2,4,2,0))
  plot(K228path1$time, K228path1$x[, 1], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '')
  abline(h = CapLocs[i,1])
  abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
  
  
  plot(Americas, col = "grey95",xlim = c(-170,-60),ylim=c(0,65))
  box()
  lines(K228path1$x, col = "blue")
  points(K228path1$x, pch = 16, cex = 0.5, col = "blue")
 
  plot(K228path2$time, K228path2$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '',xaxt="n")
  abline(h = CapLocs[4,2])
  abline(v = as.POSIXct("2015-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2016-03-20"),col="red",lty=2,lwd=1.5)
  par(mar=c(2,4,2,0))
  plot(K228path2$time, K228path2$x[, 1], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '')
  abline(h = CapLocs[4,1])
  abline(v = as.POSIXct("2015-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2016-03-20"),col="red",lty=2,lwd=1.5)
  
  plot(Americas, col = "grey95",xlim = c(-170,-60),ylim=c(0,65))
  box()
  lines(K228path2$x, col = "red")
  points(K228path2$x, pch = 16, cex = 0.5, col = "red")

## ------------------------------------------------------------------------
x0 <- z0 <- vector('list',nBirds)

for(i in c(1:3,5:nBirds)){
  # Take the location estimates created above
x0[[i]]<- path[[i]]$x

  # the model also needs the mid-points - generate those here
z0[[i]]<- trackMidpts(x0[[i]])
}

## ------------------------------------------------------------------------
beta <- c(0.7, 0.08)

## ------------------------------------------------------------------------
# This function sets the known location of the bird - capture location and recapture location
fixedx<-vector('list',nBirds)

for(i in 1:c(1:3,5:nBirds)){
  fixedx[[i]]<- rep(F, nrow(x0[[i]]))
  fixedx[[i]][1:5] <- T
  #fixedx[[i]][(nrow(x0[[i]])-5):nrow(x0[[i]])] <-T
  
  x0[[i]][fixedx[[i]], 1] <- CapLocs[i,1]
  x0[[i]][fixedx[[i]], 2] <- CapLocs[i,2]
  
  z0[[i]] <- trackMidpts(x0[[i]]) # update z0 positions
}

## ------------------------------------------------------------------------
# set xlim and ylim values need to span the range of your dataset
xlim <- c(-170, -60)
ylim <- c(-89, 90)

## ------------------------------------------------------------------------
## Function to construct a land/sea mask
distribution.mask <- function(xlim, ylim, n = 4, land = TRUE, shape) {
    r <- raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1], 
        xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(shape))
    r <- cover(rasterize(shape, shift = c(-360, 0), r, 1, silent = TRUE), 
        rasterize(shape, r, 1, silent = TRUE), rasterize(elide(shape, 
            shift = c(360, 0)), r, 1, silent = TRUE))
    r <- as.matrix(is.na(r))[nrow(r):1, ]
    if (land) 
        r <- !r
    xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
    ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)

    function(p) {
        r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
    }
}

## ------------------------------------------------------------------------

crs(Americas)<-WGS84
crs(OSFLDist)<-WGS84

Land <- subset(Americas, SUBREGION == 5 | SUBREGION == 21 | SUBREGION == 13)
E<-disaggregate(subset(Land, NAME == "Ecuador"))
Ecuador <- E[14,]
Land <- (subset(Land, NAME != "Ecuador"))
Land<-(gUnion(Land,Ecuador))
## ------------------------------------------------------------------------
## Define mask for Ovenbird distribution
is.dist <- distribution.mask(shape=Land,
                             xlim = xlim,
                             ylim = ylim,
                             n = 4, # 0.25 x 0.25 degree resolution
                             land = TRUE)

## ------------------------------------------------------------------------
# Define the log prior for x and z
log.prior <- function(p) {
    f <- is.dist(p)
    ifelse(f | is.na(f), 0, -10)
}

## ------------------------------------------------------------------------
# Define the threshold model - slimilar to above #
model <-  vector('list', nBirds)

for(i in 1:nBirds){
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
                            zenith = AncZ)
}

## ------------------------------------------------------------------------
for(i in c(1:3,5:nBirds)){
  model[[i]]$fixedx<-c(model[[i]]$fixedx,rep(FALSE,(dim(model[[i]]$x0)[1]-length(model[[i]]$fixedx))))
  model[[i]]$fixedx[(length(model[[i]]$fixedx)-5):length(model[[i]]$fixedx)]<-TRUE
}

## ------------------------------------------------------------------------
# This defines the error distribution around each location #
proposal.x <- proposal.z <- vector('list',nBirds)

for(i in c(1:3,5:nBirds)){
proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0[[i]]))
proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0[[i]]))
}

## ----eval = FALSE--------------------------------------------------------
 fit <- xsum <- zsum <- vector('list', nBirds)
 
 for(i in c(1:3,5:nBirds)){
 
 fit[[i]] <- estelleMetropolis(model = model[[i]],
                               proposal.x = proposal.x[[i]],
                               proposal.z = proposal.z[[i]],
                               iters = 2000, # This value sets the number of iterations to run
                               thin = 10,
                               chains = 3)

 xsum[[i]] <- locationSummary(fit[[i]]$x,collapse = TRUE)
 zsum[[i]] <- locationSummary(fit[[i]]$z,collapse = TRUE)

proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%')))
proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%')))

 fit[[i]] <- estelleMetropolis(model = model[[i]],
                              proposal.x = proposal.x[[i]],
                               proposal.z = proposal.z[[i]],
                               x0 = cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                               z0 = cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                               iters=2000, # This value sets the number of iterations to run
                               thin=10,
                               chains=3)
 
# Final Run
  xsum[[i]] <- locationSummary(fit[[i]]$x,collapse = TRUE)
  zsum[[i]] <- locationSummary(fit[[i]]$z,collapse = TRUE)

 proposal.x[[i]] <- mvnorm(chainCov(fit[[i]]$x),s=0.1)
 proposal.z[[i]] <- mvnorm(chainCov(fit[[i]]$z),s=0.1)
 
 # Note the increase in number of interations - this takes a bit longer to run
fit[[i]] <- estelleMetropolis(model = model[[i]],
                               proposal.x = proposal.x[[i]],
                               proposal.z = proposal.z[[i]],
                               x0=cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                               z0=cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                               iters=5000,  # This value sets the number of iterations to run
                               thin=10,
                               chains=3)
}

k228yr <- k228model <- vector('list',2)

k228model[[1]]<- thresholdModel(twilight = twlEdit[[4]]$Twilight[1:612],
                            rise = twlEdit[[4]]$Rise[1:612],
                            twilight.model = "ModifiedLogNormal",
                            alpha = ALPHA,
                            beta = beta,
                            # Here is where we set the constraints for land
                            logp.x = log.prior, 
                            logp.z = log.prior, 
                            x0 = K228path1$x,
                            z0 = trackMidpts(K228path1$x),
                            zenith = AncZ)


proposal.x[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(K228path1$x))
proposal.z[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(trackMidpts(K228path1$x)))

k228yr[[1]] <- estelleMetropolis(model = k228model[[1]],
                                 proposal.x = proposal.x[[4]],
                                 proposal.z = proposal.z[[4]],
                                 iters = 2000, # This value sets the number of iterations to run
                                 thin = 10,
                                 chains = 3)

 xsum[[4]] <- locationSummary(k228yr[[1]]$x,collapse = TRUE)
 zsum[[4]] <- locationSummary(k228yr[[1]]$z,collapse = TRUE)

proposal.x[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(xsum[[4]]$'Lon.50%',xsum[[4]]$'Lat.50%')))
proposal.z[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(zsum[[4]]$'Lon.50%',zsum[[4]]$'Lat.50%')))

k228yr[[1]] <- estelleMetropolis(model = k228model[[1]],
                                 proposal.x = proposal.x[[4]],
                                 proposal.z = proposal.z[[4]],
                                 x0 = cbind(xsum[[4]]$'Lon.50%',xsum[[4]]$'Lat.50%'),
                                 z0 = cbind(zsum[[4]]$'Lon.50%',zsum[[4]]$'Lat.50%'),
                                 iters=2000, # This value sets the number of iterations to run
                                 thin=10,
                                 chains=3)
 
# Final Run
  xsum[[4]] <- locationSummary(k228yr[[1]]$x,collapse = TRUE)
  zsum[[4]] <- locationSummary(k228yr[[1]]$z,collapse = TRUE)

 proposal.x[[4]] <- mvnorm(chainCov(k228yr[[1]]$x),s=0.1)
 proposal.z[[4]] <- mvnorm(chainCov(k228yr[[1]]$z),s=0.1)
 
 # Note the increase in number of interations - this takes a bit longer to run
k228yr[[1]] <- estelleMetropolis(model = k228model[[1]],
                                 proposal.x = proposal.x[[4]],
                                 proposal.z = proposal.z[[4]],
                                 x0 = cbind(xsum[[4]]$'Lon.50%',xsum[[4]]$'Lat.50%'),
                                 z0 = cbind(zsum[[4]]$'Lon.50%',zsum[[4]]$'Lat.50%'),
                                 iters=5000, # This value sets the number of iterations to run
                                 thin=10,
                                 chains=3)

fit[[4]]<-k228yr[[1]]

k228model[[2]]<- thresholdModel(twilight = twlEdit[[4]]$Twilight[613:nrow(twlEdit[[4]])],
                            rise = twlEdit[[4]]$Rise[613:nrow(twlEdit[[4]])],
                            twilight.model = "ModifiedLogNormal",
                            alpha = ALPHA,
                            beta = beta,
                            # Here is where we set the constraints for land
                            logp.x = log.prior, 
                            logp.z = log.prior, 
                            x0 = K228path2$x,
                            z0 = trackMidpts(K228path2$x),
                            zenith = AncZ)


proposal.x[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(K228path2$x))
proposal.z[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(trackMidpts(K228path2$x)))

k228yr[[2]] <- estelleMetropolis(model = k228model[[2]],
                                 proposal.x = proposal.x[[4]],
                                 proposal.z = proposal.z[[4]],
                                 iters = 2000, # This value sets the number of iterations to run
                                 thin = 10,
                                 chains = 3)

 xsum[[4]] <- locationSummary(k228yr[[2]]$x,collapse = TRUE)
 zsum[[4]] <- locationSummary(k228yr[[2]]$z,collapse = TRUE)

proposal.x[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(xsum[[4]]$'Lon.50%',xsum[[4]]$'Lat.50%')))
proposal.z[[4]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(zsum[[4]]$'Lon.50%',zsum[[4]]$'Lat.50%')))

k228yr[[2]] <- estelleMetropolis(model = k228model[[2]],
                                 proposal.x = proposal.x[[4]],
                                 proposal.z = proposal.z[[4]],
                                 x0 = cbind(xsum[[4]]$'Lon.50%',xsum[[4]]$'Lat.50%'),
                                 z0 = cbind(zsum[[4]]$'Lon.50%',zsum[[4]]$'Lat.50%'),
                                 iters=2000, # This value sets the number of iterations to run
                                 thin=10,
                                 chains=3)
 
# Final Run
  xsum[[4]] <- locationSummary(k228yr[[2]]$x,collapse = TRUE)
  zsum[[4]] <- locationSummary(k228yr[[2]]$z,collapse = TRUE)

 proposal.x[[4]] <- mvnorm(chainCov(k228yr[[2]]$x),s=0.1)
 proposal.z[[4]] <- mvnorm(chainCov(k228yr[[2]]$z),s=0.1)
 
 # Note the increase in number of interations - this takes a bit longer to run
k228yr[[2]] <- estelleMetropolis(model = k228model[[2]],
                                 proposal.x = proposal.x[[4]],
                                 proposal.z = proposal.z[[4]],
                                 x0 = cbind(xsum[[4]]$'Lon.50%',xsum[[4]]$'Lat.50%'),
                                 z0 = cbind(zsum[[4]]$'Lon.50%',zsum[[4]]$'Lat.50%'),
                                 iters=5000, # This value sets the number of iterations to run
                                 thin=10,
                                 chains=3)

#for(i in 1:nBirds){
#saveRDS(fit[[i]],paste0(BirdId[i],"_","fit.rds"))
#}
#saveRDS(k228yr[[2]],paste0(BirdId[4],"_","fit_year2_.rds"))

fit <- lapply(list.files("Data/MCMC_fit",full.names = TRUE),readRDS)
nBirds <- length(fit)


## ------------------------------------------------------------------------
# This step makes an empty raster #
r <- raster(nrows=4*diff(ylim),
            ncols=4*diff(xlim),
            xmn=xlim[1],
            xmx=xlim[2],
            ymn=ylim[1],
            ymx=ylim[2])

## ------------------------------------------------------------------------
S <- Sp <- vector('list',nBirds)

for(i in 1:nBirds){
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

Nov1<-Mar1<-rep(NA,nBirds)
for(i in 1:3){
Nov1[i]<-"2013-11-01"
Mar1[i]<-"2014-03-01"
}
for(i in c(4,6:12)){
Nov1[i]<-"2014-11-01"
Mar1[i]<-"2015-03-01"
}
for(i in c(5,13:16)){
Nov1[i]<-"2015-11-01"
Mar1[i]<-"2016-03-01"
}

source("MigSchedules.R")

schedules <- vector('list',nBirds)
for(i in 13:nBirds){
schedules[[i]]<- MigSchedule(S[[i]],prob = 0.95, known.stationary = c(Nov1[i],Mar1[i]),rm.lat.equinox = TRUE, days.omit = 2)
}

str(schedules,1)
## ------------------------------------------------------------------------
# Create a vector of the Dates in the file - here we set rise=true because the data has both rise and sunset
DATES<-vector('list',nBirds)

for(i in 1:nBirds){
DATES[[i]]<-S[[i]]$mcmc[[1]]$time[which(S[[i]]$mcmc[[1]]$rise==TRUE)]
}
DATES[[12]][1]
# Non-breeding locations # Oct 5th to April 1st
# Here we specify our dates of interest. 

Nov1<-Mar1<-July<-rep(NA,nBirds)
for(i in 1:3){
July[i]<-as.POSIXct("2013-07-30", tz = "GMT")
Nov1[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2013-11-01"))
Mar1[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2014-03-01"))
}
for(i in 4:11){
July[i]<-as.POSIXct("2014-07-30", tz = "GMT")
Nov1[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2014-11-01"))
Mar1[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2015-03-01"))
}
for(i in 12:15){
July[i]<-as.POSIXct("2015-07-30", tz = "GMT")
Nov1[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2015-11-01"))
Mar1[i]<-which(strptime(DATES[[i]],format="%Y-%m-%d")==as.POSIXct("2016-03-01"))
}

July <- as.POSIXct(July,origin = "1970-01-01",tz = "GMT")

# Breeding 
# Get time when on the breeding grounds 
NBtime<-NB<-NBprob<-breeddata<-lon<-cal<-lat<-vector('list',nBirds)

source("Lisovski_functions.R")

for(i in 1:nBirds){
NBtime[[i]]<-SGAT::sliceInterval(S[[i]],k=c(Nov1[i]:Mar1[i]))

# "Slice" the data and save all dates between Release date and July 31 2011, and June 1 2012 until capture.
NB[[i]]<-SGAT::slice(S[[i]], k=c(Nov1[i]:Mar1[i]))

breeddata[[i]]<-subset(OSFLdata[[i]],OSFLdata[[i]][,1]<July[i])

lon[[i]] <- estim.lon(datetime = breeddata[[i]][,1], light = log(breeddata[[i]][,2]))

cal[[i]] <- calibration(datetime = breeddata[[i]][,1], light = log(breeddata[[i]][,2]), crds = CapLocs[i,], max.light = 5,aggr = 5)

lat[[i]] <- estim.lat(datetime = breeddata[[i]][,1], light = log(breeddata[[i]][,2]), calib = cal[[i]], lonlim =c(-160,-130), latlim = c(50,75))
}

# Define breeding areas # 
for(i in 1:nBirds){
lat[[i]]$raster[lat[[i]]$raster<0.95]<-NA
e<-extent(x = c(as.numeric(lon[[i]]$pos[2]),as.numeric(lon[[i]]$pos[5])),
          y = c(50,70))
lat[[i]]$raster <- crop(lat[[i]]$raster,e)
lat[[i]]$raster <- mask(lat[[i]]$raster,OSFLDist)
}

library(ks)
#plot(OSFLDist)
#for(i in 1:649){
#plot(raster(kde(cbind(S[[1]]$mcmc$z[[1]][i,1,],S[[1]]$mcmc$z[[1]][i,2,]), h = 1)),add = TRUE)
#}


## ----echo=FALSE, fig.cap = "**Figure 6** Non-breeding location of each individual Olive-sided Flycather"----
for(i in 1:nBirds){
  print(BirdId[i])
  plot(Americas,col="gray74",
       ylim=c(-20,10),
       xlim=c(-115,-60))
  #plot(OSFLDist,
   #    border="gray88",
   #   col="gray88",
    #   add=TRUE)
  plot(NB[[i]],
       ylim=c(-20,10),
       xlim=c(-115,-60),
       axes=FALSE,
       add=TRUE,
       legend=FALSE,
       col=c(rep("transparent",10),rev(bpy.colors(100))),
       cex.axis=0.7)
  plot(Americas,border="black",add=TRUE) 
  }

## ----echo=FALSE, fig.cap="**Figure 7** Joint non-breeding location for all Olive-sided Flycatchers"----
for(i in 1:nBirds){
NBprob[[i]]<-NB[[i]]
NBprob[[i]][NBprob[[i]]<quantile(NBprob[[i]], probs = 0.95)]<-NA
NBprob[[i]][!is.na(NBprob[[i]])]<- 1
}

NBprob<-stack(NBprob)
plot(NBprob[[1]])
sumNBbirds<-sum(NBprob,na.rm = TRUE)
plot(sumNBbirds)
set.breaks=seq(0,15,1)
plot(Americas,col="gray74",
       ylim=c(-20,10),xlim=c(-115,-60))
  #plot(OSFLdist,border="gray88",,col="gray88",add=TRUE)
  plot(sumNBbirds,ylim=c(-20,10),xlim=c(-115,-60),
       axes=FALSE, add=TRUE,
       legend=FALSE,
       breaks=set.breaks,
        col=c("transparent",rev(bpy.colors(6))),
       cex.axis=0.7)
raster::scalebar(d=1000,lonlat=TRUE,xy=c(-98,-10),label="1000 km")
  plot(Americas,border="black",add=TRUE)  
plot(sumNBbirds,legend.only=TRUE,horizontal=TRUE,col=c("transparent",rev(bpy.colors(6))),
     smallplot=c(0.2,0.7,0.25,0.26),
     legend.width=0.25,
     legend.shrink=0.5,
     axis.args=list(at=seq(0,15,1),
                    labels=seq(0,15,1),
                    cex.axis=1,
                    mgp=c(5,0.3,0)),
     legend.args=list(text='Number of Birds', side=3, font=2, line=0.5, cex=1.25),add=TRUE)

year <- SUMMARY <- winStart <- winEnd <- endBreed <- startBreed <- vector('list',nBirds)
for(i in 1:nBirds){
year[[i]]<-slice(S[[i]], k=c(1:length(S[[i]]$weights)), chains = 3)
year[[i]][year[[i]] > 5] <- 5

plot(OSFLDist) 
plot(Americas, add = TRUE,col = "gray88")
plot(year[[i]], breaks = seq(0,5,0.01), col = c("transparent",rev(bpy.colors(250))),horizontal = TRUE,add=TRUE)
plot(Americas, add = TRUE)
SUMMARY[[i]]<-locationSummary(fit[[i]]$z)
points(cbind(yrat1$'Lon.50%',yrat1$'Lat.50%'),type = "o",pch = ".")
points(CapLocs[i,1],CapLocs[i,2],pch = 19,cex = 2)

endBreed[[i]] <- twlEdit[[i]][min(which(is.na(extract(lat[[i]]$raster,cbind(SUMMARY[[i]]$'Lon.50%',SUMMARY[[i]]$'Lat.50%'))))),1]
startBreed[[i]] <- twlEdit[[i]][max(which(is.na(extract(lat[[i]]$raster,cbind(SUMMARY[[i]]$'Lon.50%',SUMMARY[[i]]$'Lat.50%'))))),1]

winStart[[i]] <- twlEdit[[i]][min(which(!is.na(extract(NB[[i]],cbind(SUMMARY[[i]]$'Lon.50%',SUMMARY[[i]]$'Lat.50%'))))),1]
winEnd[[i]] <- twlEdit[[i]][max(which(!is.na(extract(NB[[i]],cbind(SUMMARY[[i]]$'Lon.50%',SUMMARY[[i]]$'Lat.50%'))))),1]
}

difftime(as.POSIXct(unlist(winStart),origin = "1970-01-01",tz = "GMT"), as.POSIXct(unlist(endBreed),origin = "1970-01-01",tz= "GMT"),units = "days")

source("Lisovski_functions.R")
which(OSFLdata[[4]][,1] < twlEdit[[4]][1,1])

osfl4<-subset(OSFLdata[[4]],OSFLdata[[4]][,1] > as.POSIXct("2015-05-01",tz = "GMT") & OSFLdata[[4]][,1] < as.POSIXct("2015-07-10",tz = "GMT"))

lon <- estim.lon(datetime = osfl4[,1], light = osfl4[,2])

cal <- calibration(datetime = osfl4[1:16837,1], light = osfl4[1:16837,2], crds = CapLocs[4,], max.light = 5,aggr = 5)

lat <- estim.lat(datetime = osfl4[,1], light = osfl4[,2], calib = cal, lonlim =c(-160,-130), latlim = c(50,75))

points(CapLocs[4,1],CapLocs[4,2],pch = 19, cex = 2)
abline(v = -142.1151)

lat$raster[lat$raster<0.99]<-NA
plot(lat$raster)
plot(OSFLDist,add=TRUE)
plot(States, add = TRUE)


















## ---- echo=FALSE, warning=FALSE,error=FALSE------------------------------
S<-Dates<-AnnualCycle<-AnnualCycle_stack<-AnnualCycle_proj<-SummaryPath<-PathProj<-NB_proj<-vector("list",nBirds)
stops<-rep(NA,nBirds)

xlim <- c(-170, -60)
ylim <- c(-89, 90)

r <- raster(nrows=4*diff(ylim),ncols=4*diff(xlim),xmn=xlim[1],xmx=xlim[2],ymn=ylim[1],ymx=ylim[2])

polyconic<-"+proj=poly +lat_0=0 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

Americas_proj<-spTransform(Americas,CRS=polyconic)

for(i in 1:nBirds){
crs(NB[[i]])<-WGS84
NB_proj[[i]]<-projectRaster(NB[[i]],crs=polyconic)

SummaryPath[[i]]<-SpatialLines(list(Lines(Line(cbind(zm[[i]][,"Lon.mean"],zm[[i]][,"Lat.mean"])),ID="a")))

crs(SummaryPath[[i]])<-WGS84
}
PathProj[[i]]<-spTransform(SummaryPath[[i]],CRS=polyconic)
plot(PathProj[[1]])
plot(Americas_proj,add=TRUE)
AnnualCycle[[i]]<-vector("list",stops[i])
S[[i]] <- slices(type="primary", breaks="day", mcmc=fit[[i]],grid=r)
Dates[[i]]<-S[[i]]$mcmc[[1]]$time[which(S[[i]]$mcmc[[1]]$rise==TRUE)]

for(t in 1:stops[i]){
if(t<stops[i]){
AnnualCycle[[i]][[t]]<-slice(S[[i]],k=c(which(strptime(Dates[[i]],format="%Y-%m-%d")==as.POSIXct(strptime(cL[[i]]$migTable[t,2],format="%Y-%m-%d"))):
                                         which(strptime(Dates[[i]],format="%Y-%m-%d")==as.POSIXct(strptime(cL[[i]]$migTable[t,3],format="%Y-%m-%d")))))
} else{
AnnualCycle[[i]][[t]]<-slice(S[[i]],k=c(which(strptime(Dates[[i]],format="%Y-%m-%d")==as.POSIXct(strptime(cL[[i]]$migTable[t,2],format="%Y-%m-%d"))):
                                        length(Dates[[i]])))
}
  
AnnualCycle[[i]][[t]]<-AnnualCycle[[i]][[t]]/maxValue(AnnualCycle[[i]][[t]])
AnnualCycle[[i]][[t]][AnnualCycle[[i]][[t]]<quantile(AnnualCycle[[i]][[t]],probs=0.90)]<-NA
  } # t

AnnualCycle_stack[[i]]<-stack(AnnualCycle[[i]])
AnnualCycle_proj[[i]]<-projectRaster(AnnualCycle_stack[[i]],crs=polyconic)


par(mar=c(4,0,4,0))
set.breaks<-seq(0.5,1,0.05)
plot(Americas_proj,
     col="gray75",
     border="gray88",
     main=BirdId[i],
     xlim=c(-3413915,4402162),
     ylim=c(-2462665,8710642))
plot(PathProj[[i]],add=TRUE,lty=2,col="black")
for(t in 1:stops[i]){
plot(AnnualCycle_proj[[i]][[t]],
     add=TRUE,
     breaks=set.breaks,
     col=rev(heat.colors(length(set.breaks))),
     legend=FALSE)
} # t
plot(NB_proj[[i]],
     add=TRUE,
     legend=FALSE,
     col=c(rep("transparent",10),rev(bpy.colors(100))))
for(t in 1:stops[i]){
plot(AnnualCycle[[i]][[t]],breaks=set.breaks,col=rev(heat.colors(length(set.breaks))),xlim=xlim,ylim=ylim,
     main=print(as.POSIXct(strptime(cL[[i]]$migTable[t,2],format="%Y-%m-%d"))))
plot(Americas,add=TRUE,border="gray")
}
} # i

## ----echo=FALSE, fig.cap="**Figure:** Locations where birds overlap in space during spring (left) and fall (right) migration."----
DepartBreed<-ArriveWinter<-DepartWinter<-ArriveBreed<-rep(NA,nBirds)
dpbrd<-c("2013-08-01", #F202
         "2013-08-13", #F203
         "2013-08-04", #F311
         "2014-07-24", #648
         "2014-08-08", #654
         "2014-07-31", #657
         "2014-08-05", #659
         "2014-07-31", #660
         "2014-07-15", #665
         "2014-08-15") #670

arwin<-c("2013-10-05", #202
         "2013-10-13", #203
         "2013-09-23", #311
         "2014-09-24", #648
         "2014-10-17", #654
         "2014-09-26", #657
         "2014-09-28", #659
         "2014-09-25", #660
         "2014-10-08", #665
         "2014-10-25") #670

dpwin<-c("2014-04-01", #202
         "2014-04-04", #203
         "2014-03-20", #311
         "2015-04-01", #648
         "2015-04-10", #654
         "2015-04-04", #657
         "2015-03-30", #659
         "2015-04-06", #660
         "2015-03-08", #665
         "2015-04-01") #670

arbrd<-c("2014-06-03", #202
         "2014-05-30", #203
         "2014-05-15", #311
         "2015-05-20", #648
         "2015-06-02", #654
         "2015-05-15", #657
         "2015-05-18", #659
         "2015-05-31", #660
         "2015-05-09", #665
         "2015-05-16") #670

for(i in 1:nBirds){
DepartBreed[i]<-which(strptime(Dates[[i]],format="%Y-%m-%d")==as.POSIXct(dpbrd[i]))
ArriveWinter[i]<-which(strptime(Dates[[i]],format="%Y-%m-%d")==as.POSIXct(arwin[i]))
DepartWinter[i]<-which(strptime(Dates[[i]],format="%Y-%m-%d")==as.POSIXct(dpwin[i]))
ArriveBreed[i]<-which(strptime(Dates[[i]],format="%Y-%m-%d")==as.POSIXct(arbrd[i]))
}

ylim<-c(-40,80)

##### Spring Connectivity Map #######
# This step makes an empty raster #

SpringMig<-vector("list",nBirds)
# Here for each day we will fill the raster with the locations #
# Note - you can change breaks to "month" - "week" - etc. 
for( i in 1:nBirds){
# "Slice" the data and save all dates between Release date and July 31 2011, and June 1 2012 until capture.
SpringMig[[i]]<-slice(S[[i]], k=c(DepartWinter[i]:ArriveBreed[i]))
SpringMig[[i]][is.na(SpringMig[[i]])]<-0
SpringMig[[i]]<-SpringMig[[i]]/maxValue(SpringMig[[i]])
SpringMig[[i]][SpringMig[[i]]>=quantile(SpringMig[[i]],probs=0.975)]<-1
}

SumSpring<-sum(stack(SpringMig))
SumSpring[SumSpring<1]<-NA

cols<-rev(bpy.colors(10))
par(bty="n",mar=c(3,0,0,0),mfrow=c(1,2))
par(bty="n",mar=c(3,0,0,0))
plot(SumSpring,axes=FALSE,legend=FALSE,col=cols,breaks=c(1:10),horizontal=TRUE,xlim=xlim,ylim=ylim)
plot(SumSpring,legend.only=TRUE,horizontal=TRUE,col=cols,breaks=c(1:10),
     smallplot=c(0.2,0.5,0.25,0.27),
     legend.width=0.25,
     legend.shrink=0.5,
     axis.args=list(at=c(1,5.5,9.9),
                    labels=c(1,5,10),
                    cex.axis=1,
                    mgp=c(5,0.3,0)),
     legend.args=list(text="Spring # of Birds", side=3, font=2, line=0.5, cex=1),add=TRUE)

plot(Americas,add=TRUE,border="gray")

FallMig<-vector("list",nBirds)
# Here for each day we will fill the raster with the locations #
# Note - you can change breaks to "month" - "week" - etc. 
for( i in 1:nBirds){
FallMig[[i]]<-slice(S[[i]], k=c(DepartBreed[i]:ArriveWinter[i]))
FallMig[[i]][is.na(FallMig[[i]])]<-0
FallMig[[i]]<-FallMig[[i]]/maxValue(FallMig[[i]])
FallMig[[i]][FallMig[[i]]>=quantile(FallMig[[i]],probs=0.975)]<-1
}

SumFall<-sum(stack(FallMig))
SumFall[SumFall<1]<-NA

for(i in 1:nBirds){
FallMig[[i]][FallMig[[i]]>50 | FallMig[[i]]<0.5]<-NA
}
FallMig[[10]][FallMig[[10]]>20]<-NA
par(bty="n",mar=c(3,0,0,0))
plot(SumFall,axes=FALSE,legend=FALSE,col=cols,breaks=c(1:10),horizontal=TRUE,xlim=xlim,ylim=ylim)
plot(SumFall,legend.only=TRUE,horizontal=TRUE,col=cols,breaks=c(1:10),
     smallplot=c(0.2,0.5,0.25,0.27),
     legend.width=0.25,
     legend.shrink=0.5,
     axis.args=list(at=c(1,5.5,9.9),
                    labels=c(1,5,10),
                    cex.axis=1,
                    mgp=c(5,0.3,0)),
     legend.args=list(text="Fall # of Birds", side=3, font=2, line=0.5, cex=1),add=TRUE)
plot(Americas,add=TRUE,border="gray")

## ---- echo = FALSE,warning=FALSE,error=FALSE,message=FALSE---------------
library(leaflet) 
pal <- leaflet::colorNumeric(c("#FFFF60FF","#FFDB24FF","#FFA857FF","#FF758AFF","#EF42BDFF","#9F0FF0FF","#5000FFFF","#0000FFFF","#000099FF","#000033FF"),                         values(SumSpring),
                    na.color = "transparent")

crs(SumSpring)<-WGS84

leaflet() %>% addTiles() %>%
  addRasterImage(SumSpring, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(SumSpring),
            title = "Number of Birds")

## ----echo=FALSE----------------------------------------------------------
plot(Americas,xlim=xlim,ylim=ylim,col="gray88")
plot(subset(Americas,NAME=="Guatemala" | NAME=="Ecuador" | NAME=="Nicaragua"),add=TRUE,col="red")

