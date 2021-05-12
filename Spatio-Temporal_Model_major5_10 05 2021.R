############
## SET-UP ##
############
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
library(tmaptools)
library(sp)
library(rgdal)
library(spData)
library(spatstat)
library(lgcp)
#library(OpenStreetMap)
library(rpanel)
library(fields)
library(miscFuncs)
library(caret)
library(stars)
library(units)
library(udunits2)
library(fasterize)

dat<-read.csv("ST_Lineage_5_Data.csv") %>% # 132 cases
  dplyr::select(hh_LAT, hh_LNG, sample_collected) %>%
  filter(!is.na(hh_LAT)) #removing records with missing GPS coordinates
#dat <- read.csv("ST_All-Linege_Data.csv")
dim(dat) # 132 cases

# Converting dataset into a space-time planar point pattern object
x <- dat$hh_LNG
y <- dat$hh_LAT
t <- dat$sample_collected
t <- as.Date(t,origin="2015-03-27")
tm <- as.integer(t - min(t))
tm <- tm+1
# Converting coordinates to utm
xy <- SpatialPoints(cbind(x,y), proj4string=CRS("+proj=longlat"))
xy_utm <-spTransform(xy, CRS("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84"))
xy_utm
plot(xy_utm)

# Defining the time object which will be used for modelling

Time_lim <- as.integer(c(min(tm),max(tm))) # first and last value for object t
#Time_lim <- as.integer(c(16522,16832))

## population density for Blantyre (2018)
pop18 <- st_read("Blantyre_City.shp") # from NSO; shape file
pop18 <- st_transform(pop18, "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84")
BT_sp <- pop18
EA_2018 <- read.csv("popEA_NSO_2018_BTCity_final.csv") # census data

pop18$EA_NUMBER <- as.integer(pop18$EA_NUMBER)
pop18 <- left_join(pop18, EA_2018, by = c("ADM_NAME" = "ADM_NAME", "EA_NUMBER" = "EA_NUMBER"))

# Removing NAs

pop18$Total_0.10[is.na(pop18$Total_0.10) & pop18$ADM_TYPE == "PA"] <- 0
#pop18$den18[is.na(pop18$den18) & pop18$ADM_TYPE == "PA"] <- 0

# Calculating population density
pop18 <- pop18 %>%
  mutate(area18=units::set_units(st_area(geometry),value="km^2")) %>%
  mutate(den18=Total_0.10/area18)

# Plotting the shapefile

plot(pop18["Total_0.10"], breaks=c(0,10,50,100,150,200,250,300,400,600,800,1000))
plot(pop18["den18"], breaks=c(0,10,50,100,500,1000,2000,30000,4000,5000,1e4,2e4))

# Owin object preparation for Blantyre City
BTShapeF <- st_read("MWI_BlantyreCity.shp")
BTShapeF <- st_transform(BTShapeF, "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84")

BTShapeF
BTShapeF <- st_union(BTShapeF)
class(BTShapeF)
BTShapeF
BTShapeF_owin <- as.owin(BTShapeF)
class(BTShapeF_owin)
plot(BTShapeF_owin)
BTShapeF_owin

# Checking if all points are inside the owin object
inside.owin(xy_utm$x,xy_utm$y,BTShapeF_owin) # 16 cases are outside the owin boundary

#Setting up a spatio-temporal object
STdat <- cbind(xy_utm$x,xy_utm$y,tm)
xyt <- stppp(list(data=STdat, tlim = Time_lim, window = BTShapeF_owin))
xyt # only 294 cases included 
xyt$t <- as.integer(xyt$t)
plot(xyt)

# computing approximate values of the process Y 
# on an appropriately transformed scale
# these approximates are used to inform size of 
# computational grid so as to capture the dependence
# properties of Y
# It is also used to inform the proposal kernel for the MCMC algorithm

# Plotting the spatial density function
denty <- lambdaEst(xyt, axes=T)
plot(denty)

# Coercing the density object to object of class spatialAtRisk for parameter estimation later

Spatrisk  <- spatialAtRisk(denty)
Spatrisk

# Plotting temporal intensity mu(t)
mu_t <- muEst(xyt, f=1/20) # f is the lowess argument for smoothing
mu_t # temporalAtRisk object
plot(mu_t)

###########################################
## SIZE OF COMP GRID & OFFSET OBJECT DEFINITION ##
###########################################

CellWidth <- 1000 # Change this value to change computational grid size
EXT <- 2 # to be used during polygon overlay
minimum.contrast(xyt, model = "exponential", method = "g", intens = density(xyt), transform = log)
chooseCellwidth(xyt,cwinit = CellWidth) # cell width is 175 metres

#pop18 <- pop18 %>% mutate(pop.count=(den18*CellWidth^2)/1e6)

# convert to SpatialPolygonsDataFrame
pop18 <- as_Spatial(pop18)
BT_sp <- as_Spatial(BT_sp)

# perform polygon overlay operations and compute computational grid

polyolay_off <- getpolyol(data = xyt, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT)
polyolay <- getpolyol(data = xyt, regionalcovariates = BT_sp, cellwidth = CellWidth, ext = EXT)

## FIRST MODEL FORMULAE

FORM <- X ~ Total_0.10 + moty
FORM_Spatial <- X ~ Total_0.10
FORM_Temporal <- t ~ moty -1

## SECOND MODEL FORMULAE

FORM_2 <- X ~ moty
FORM_Spatial_2 <- X ~ 1
FORM_Temporal_2 <- t ~ moty -1

# set the interpolation type for each variable

pop18@data <- guessinterp(pop18@data)
pop18@data <- assigninterp(df = pop18@data, vars =  c( "Total_0.10"), value = "ArealWeightedSum")
class(pop18@data$Total_0.10)

#pop3 <- setClass(pop3, slots = c(Total_0.10 = "numeric", geometry = "numeric"))
#popDen18_ <- as_Sp
## DESIGN MATRIXatial(popDen18_) # changing class type of the object to S4
#xyt <- as_Spatial(xyt)
#Zmat <- getZmat(formula = FORM, data = xyt, regionalcovariates = popDen18, cellwidth = CellWidth, ext = EXT, overl = polyolay)

Zmat <- getZmat(formula = FORM_Spatial_2, data = xyt, regionalcovariates = NULL, cellwidth = CellWidth, ext = EXT, overl = polyolay_off)
Zmat_off <- getZmat(formula = FORM_Spatial, data = xyt, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT, overl = polyolay)
#plot(Zmat)

# Specifying log of population counts to comply with Poisson model
# so that number of cases to be proportional to population at risk
# and not the exponential of population

#Zmat_off[, "Total_0.10"] <- log(Zmat_off[,"Total_0.10"])
#Zmat_off[, "Total_0.10"][is.infinite(Zmat_off[, "Total_0.10"])] <- min(Zmat_off[, "Total_0.10"][!is.infinite(Zmat_off[, "Total_0.10"])]) 
#plot(Zmat_off)

# DEFINING THE OFFSET

mm <- length(attr(Zmat_off, "mcens"))
nn <- length(attr(Zmat_off, "ncens"))

Pop.offset <- list(spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))))

# plot the spatial interpolated covariates
plot(Zmat, ask = F)
plot(Zmat_off, ask = F)

# construct dummy temporal data frame
#days <- c("Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
#tvec <- xyt$tlim[1]:xyt$tlim[2]
#da <- rep(days, length.out = length(tvec))
#tdata <- data.frame(t = tvec, dotw = da)

# construct dummy temporal data frame
#months.of.year <- lubridate::month(ymd(t), label = TRUE)
months.of.year <- c("Mar", "Apr" ,"May","Jun", "Jul", "Aug", "Sep" ,"Oct" ,"Nov", "Dec","Jan" ,"Feb")
tvec <- xyt$tlim[1]:xyt$tlim[2]
#tvec <- seq(min(ymd(t)),max(ymd(t)),by=1)
#mo <- lubridate::month(ymd(tvec), label = TRUE)
mo <- rep(months.of.year, length.out = length(tvec))
#tvec <- as.integer(tvec - min(tvec))
#tvec <- tvec + 1

tdata <- data.frame(t = tvec, moty = mo)

# choose last time point and number of proceeding time-points to include

ti <- max(tm)
ti <- as.integer(ti)
LAGLENGTH <- 25
LAGLENGTH <- as.integer(LAGLENGTH)

# bolt on the temporal covariates
ZmatList <- addTemporalCovariates(temporal.formula = FORM_Temporal_2, T = ti, laglength = LAGLENGTH, tdata = tdata, Zmat = Zmat)

#for(j in 1:length(ZmatList))
#{
#  ZmatList[[j]]<-ZmatList[[j]][,-which(colnames(ZmatList[[j]])=="Total_0.10")]
#}

## DEFINING PRIORS

# Defining priors
EtaP <- PriorSpec(LogGaussianPrior(mean = log(c(1,2000,1)), variance = diag(0.2,3)))
BetaP <- PriorSpec(GaussianPrior(mean = rep(0,12), variance = diag(1e+06, 12)))
priors <- lgcpPrior(etaprior = EtaP, betaprior = BetaP)

# set initial values for the algorithm
INITS <- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), betainit = NULL)

# choose the covariance function
CF <- CovFunction(exponentialCovFct)

# bolt on the temporal covariates
#ZmatList <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULA, T = T, laglength = LAGLENGTH, tdata = tdata, Zmat = Zmat)

# run the MCMC algorithm

DIRNAME <- getwd()

SpatioTemporal_major5 <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_2, xyt = xyt, T = ti, laglength = LAGLENGTH, ZmatList = ZmatList, model.priors = priors, 
                             model.inits = INITS, spatial.covmodel = CF, cellwidth = CellWidth, poisson.offset =  Pop.offset, 
                             mcmc.control = mcmcpars(mala.length = 20000, burnin = 1000, retain = 10, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                             targetacceptance = 0.574)), output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME,"ST_Model_01"), 
                             forceSave = TRUE)), ext = EXT) 

 save(list = ls(), file = file.path(DIRNAME, "ST_Models", "ST_Model_major5_output.RData")) 
 
# MODEL DIAGNOSTIC CHECKS by plotting the log target to check if the chain has converged to a posterior mode

plot(SpatioTemporal_major5)
plot(ltar(SpatioTemporal_major5), type = "s", xlab = "Iteration/900", ylab = "log target")

#plot(hvals(SpatioTemporal_Model_01)[2000:5000], type = "l", xlab = "Iteration", ylab = "h")

# compute and plot autocorrelations in the latent field
lagch <- c(1, 5, 15)
Sacf <- autocorr(SpatioTemporal_major5, lagch, inWindow = NULL)
library("fields")
for (i in 1:3) {
  image.plot(xvals(SpatioTemporal_major5), yvals(SpatioTemporal_major5), Sacf[, , i], zlim = c(-1, 1), axes = FALSE, 
             xlab = "", ylab = "", asp = 1, sub = paste("Lag:", lagch[i]), ask = FALSE)
  plot(xyt$window, add = TRUE, ask = FALSE)
  scalebar(5000, label = "5 km")
}

# produce traceplots of beta and eta
traceplots(SpatioTemporal_major5, ask = FALSE)

# produce autocorrelation plots of beta and eta
parautocorr(SpatioTemporal_major5, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(SpatioTemporal_major5)
# the above converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# a text summary of the model parameters
textsummary(SpatioTemporal_major5, digits = 4)

# a plot of the prior and posterior densities
priorpost(SpatioTemporal_major5, ask = FALSE)

# the posterior covariance function
postcov(SpatioTemporal_major5, ask = FALSE)

# exceedance and lower-tail exceedance probabilities
ep <- exceedProbs(c(1.5, 3))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_major5, ep)

plotExceed(ex[[1]], "ep", SpatioTemporal_major5, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
           sub = "", ask = FALSE)
scalebar(5000, label = "5 km")

 # conditional probabilities
#cp <- condProbs(SpatioTemporal_Model_01)

# segregation probabilities
#sr <- segProbs(SpatioTemporal_Model_01, domprob = 0.8)

# Inhomogeneous K-Function

kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
kin
plot(kin)

# Estimating the spatial correlation parameters (sigma and phi) of the model

#sigma_phi <- spatialparsEst(kin,sigma.range = c(0,12), phi.range = c(0,12), spatial.covmodel == "matern")
#sigma_phi
#sigma_phi$sigma
#sigma_phi$phi

#Estimating temporal correlation parameter theta

#theta <- thetaEst(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, sigma = 2, phi = 3) # the values for sigma and phi are not estimates...just supplied for now because sigma_phi above is giving an error

