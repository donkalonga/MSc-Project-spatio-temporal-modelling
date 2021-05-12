setwd("~/Don/MSc in Biostatistics/Reading Materials/THESIS/Spatio-Temporal Modelling/Books - Spatio-Temporal Data Analysis")
rm(list = ls())

library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggmap)
library(rgdal)
library(rgeos)
library(maptools)
library(dplyr)
library(plyr)
library(tidyr)
library(tmap)
library(matrixStats)
library(psych)
library(data.table)
library(dataMaid)
library(PrevMap)
library(sf)
library(leaflet)
library(maps)
library(RgoogleMaps)
library(raster)
library(rasterVis)
library(maptools)
library(sp)
library(rgdal)
library(spData)
library(spatstat)
library(lgcp)
library(rpanel)

dat<-read.csv("Merged_Final_DK_13072020.csv") 
dim(dat) # 549 cases with 15 variables
dat<-dat[!is.na(dat$sample_collected),] #removing cases without date
dim(dat) # 540 cases with 15 variables
dat$sample_collected<-as.Date(dat$sample_collected, origin="2015-03-28") # changing object class from Character to Date
#dat$sample_collected <- as.integer(dat$sample_collected) # changing object class from Date to Integer
head(dat) 
dat <- dplyr::select(dat, hh_LAT, hh_LNG, sample_collected, lineage)
#dat$sample_collected
#dat$lineage <- as.factor(dat$lineage)
#dat$lineage <- as.numeric(dat$lineage)
#dat$lineage
#dat<-dat[!is.na(dat$lineage),] #removing records without lineage data
#dim(dat)
#head(dat)

write.csv(dat, file = "ST_All-Lineages_Data.csv")

#Reading in the file for the population density of Blantyre 

BTShapeF <- st_read("MWI_BlantyreCity.shp")
extent(BTShapeF) # extent in UTM - in meters

# merging multiple raster files into one for offset use

pop_f_0 <- raster("mwi_f_0_2016.tif") # density population for female children 0 to 12 months
pop_f_0 <- as.data.frame(pop_f_0,xy=TRUE)
colnames(pop_f_0)[3] <- "age"
pop_f_0<-pop_f_0[!is.na(pop_f_0$age),]
#dim(pop_f_0)
#length(unique(pop_f_0$age))
#write.csv(dat, file = "pop_f_0.csv")
#pop_f_0
#head(pop_f_0)
#extent(pop_f_0)
#BTShapeF <- st_transform(BTShapeF, crs(pop_f_0))
#pop_f_0 <- crop(pop_f_0, BTShapeF)
#pop_f_0 <- mask(pop_f_0, BTShapeF)
#pop_f_0 <- is.owin(pop_f_0)
#class(pop_f_0)
pop_f_1 <- raster("mwi_f_1_2016.tif") # density population for female children between 1-4 years
pop_f_1 <- as.data.frame(pop_f_1,xy=TRUE)
colnames(pop_f_1)[3] <- "age"
pop_f_1<-pop_f_1[!is.na(pop_f_1$age),]
#dim(pop_f_1)
#length(pop_f_1$age)
#write.csv(dat, file = "pop_f_1.csv")
#head(pop_f_1)
#dim(pop_f_1)
#x <- rbind(pop_f_0,pop_f_1)
#x1 <- full_join(pop_f_0,pop_f_1, by = c("x" = "y"))
#head(x)
#dim(x)
#pop_f_1 <- as.owin(pop_f_1)
pop_f_2 <- raster("mwi_f_5_2016.tif") # density population for female children between 5-9 years
pop_f_2 <- as.data.frame(pop_f_2,xy=TRUE)
colnames(pop_f_2)[3] <- "age"
pop_f_2<-pop_f_2[!is.na(pop_f_2$age),]
#pop_f_2 <- as.owin(pop_f_2)

pop_m_0 <- raster("mwi_m_0_2016.tif") # density population for male children 0 to 12 months
pop_m_0 <- as.data.frame(pop_m_0,xy=TRUE)
colnames(pop_m_0)[3] <- "age"
pop_m_0<-pop_m_0[!is.na(pop_m_0$age),]
#extent(pop_f_0)
#pop_m_0 <- as.owin(pop_m_0)
pop_m_1 <- raster("mwi_m_1_2016.tif") # density population for male children between 1-4 years
pop_m_1 <- as.data.frame(pop_m_1,xy=TRUE)
colnames(pop_m_1)[3] <- "age"
pop_m_1<-pop_m_1[!is.na(pop_m_1$age),]
#pop_m_1 <- as.owin(pop_m_1)
pop_m_2 <- raster("mwi_m_5_2016.tif") # density population for male children between 5-9 years
pop_m_2 <- as.data.frame(pop_m_2,xy=TRUE)
colnames(pop_m_2)[3] <- "age"
pop_m_2<-pop_m_2[!is.na(pop_m_2$age),]
#pop_m_2 <- as.owin(pop_m_2)

#pop <- rasterFromXYZ(pop_f_0, pop_f_1, pop_f_2, pop_m_0,pop_m_1, pop_m_2)
#pop <- stack(pop_f_0, pop_f_1, pop_f_2, pop_m_0, pop_m_1, pop_m_2)
pop_f05 <- left_join(pop_f_0, pop_f_1, by = c("x" = "x"))
pop_f10 <- left_join(pop_f05, pop_f_2, by = c("x" = "x"))
pop_m05 <- left_join(pop_m_0, pop_m_1, by = c("x" = "x"))
pop_m10 <- left_join(pop_m05, pop_m_2, by = c("x" = "x"))
pop_10 <- left_join(pop_f10, pop_m10, by = c("x" = "x"))

plot(pop)


extent(pop) # extent in geographic - in degrees
BTShapeF <- st_transform(BTShapeF, crs(pop))
pop2 <- crop(pop, BTShapeF)
pop3 <- mask(pop2, BTShapeF)
extent(pop3)
plot(pop3)

# Converting dataset into a space-time planar point pattern object

Dat_ <- dplyr::select(dat, hh_LAT, hh_LNG, sample_collected, lineage)
dim(Dat_)
x <- Dat_$hh_LNG
y <- Dat_$hh_LAT
t <- Dat_$sample_collected
#t <- integerise(t)
marks <- Dat_$lineage
marks <- factor(marks)
class(marks)
marks <- as.numeric(marks)
class(marks)
marks

# converting coordinates to utm
xy <- SpatialPoints(cbind(x,y), proj4string=CRS("+proj=longlat"))
xy_utm <-spTransform(xy, CRS("+proj=utm +zone=36S"))
xy_utm

Da <- cbind(xy_utm$x,xy_utm$y,t,marks)
Da1 <- cbind(xy_utm$x,xy_utm$y,t)
summary(xy_utm)$is.projected # checking if the coordinates are projected
coordinates(xy_utm)
plot(xy_utm)
Time_lim <- c(16522,17165) # count of cases at time limit?? How about min(t) and max(t)
dim(Time_lim)
BTcity <- st_read("MWI_BlantyreCity.shp")
BTcity <- st_transform(BTcity, crs(xy_utm))
BTcity <- st_union(BTcity)
class(BTcity)
BTcity
BTcity_owin <- as.owin(BTcity)
class(BTcity_owin)
plot(BTcity_owin)
BTcity_owin
inside.owin(xy_utm$x,xy_utm$y,BTcity_owin) # testing if the points are inside the owin object
xyt <- stppp(list(data=Da, tlim = Time_lim, window = BTcity_owin))
xyt # only 239 
attr(xyt, "rejects")
class(xyt)
plot(xyt)

# computing approximate values of the process Y 
# on an appropriately transformed scale
# these approximates are used to inform size of 
# computational grid so as to capture the dependence
# properties of Y
# It is also used to inform the proposal kernel for the MCMC algorithm

minimum.contrast(xyt, model = "exponential", method = "g", intens = density(xyt), transform = log)
chooseCellwidth(xyt,cwinit = 175) # cell width is 175 metres

# Estimating spatial and temporal component

denty <- lambdaEst(xyt, axes=TRUE)
plot(denty)

atrisk  <- spatialAtRisk(denty)
atrisk

# Non-parametric estimate of mu (t)

mu_t <- muEst(xyt)
mu_t
mu_t1 <- muEst(xyt, f=2/3) # f is the lowess argument for smoothing
plot(mu_t1)
#mu_t1 <- constantInTime(xyt) # for constant time-trend
#plot(mu_t1)

# Estimating covariance parameters (sigma and phi) of the model

kin <- KinhomAverage(xyt, spatial.intensity = atrisk, temporal.intensity = mu_t1) # using the inhomogeneous K function
kin
plot(kin)

# Estimating the parameters

sigmaphi1 <- spatialparsEst(gin,sigma.range = c(0,10), phi.range = c(0,10))
sigmaphi1
sigmaphi2 <- spatialparsEst(kin,sigma.range = c(0,10), phi.range = c(0,10))
sigmaphi2

