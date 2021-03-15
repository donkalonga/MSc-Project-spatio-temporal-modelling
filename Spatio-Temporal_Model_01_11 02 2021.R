############
## SET-UP ##
############
setwd("~/Don/MSc in Biostatistics/Reading Materials/THESIS/Spatio-Temporal Modelling/Books - Spatio-Temporal Data Analysis")
rm(list = ls())

library(tidyverse)
library(lubridate)
#library(ggplot2)
library(ggmap)
library(rgdal)
library(rgeos)
library(maptools)
#library(dplyr)
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

dat<-read.csv("Model_01_Data.csv") %>% # 540 cases
  dplyr::select(hh_LAT, hh_LNG, sample_collected) %>%
  filter(!is.na(hh_LAT)) #removing records with missing GPS coordinates

dim(dat) # 310 cases

# Converting dataset into a space-time planar point pattern object
x <- dat$hh_LNG
y <- dat$hh_LAT
t <- dat$sample_collected

# Converting coordinates to utm
xy <- SpatialPoints(cbind(x,y), proj4string=CRS("+proj=longlat"))
xy_utm <-spTransform(xy, CRS("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84"))
xy_utm
plot(xy_utm)

# Defining the time object which will be used for modelling
Time_lim <- as.integer(c(16522,17165)) # first and last value for object t
#Time_lim <- as.integer(Time_lim)

## population density for Blantyre (2018)
pop18 <- st_read("Blantyre_City.shp") # from NSO; shape file
pop18<-st_transform(pop18, "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84")

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

###########################################
## SIZE OF COMP GRID & OFFSET OBJECT DEFINITION ##
###########################################

minimum.contrast(xyt, model = "exponential", method = "g", intens = density(xyt), transform = log)
chooseCellwidth(xyt,cwinit = CellWidth) # cell width is 175 metres

CellWidth <- 1000 # Change this value to change computational grid size
EXT <- 3 # to be used during polygon overlay
pop18 <- pop18 %>% mutate(pop.count=(den18*CellWidth^2)/1e6)

# convert to SpatialPolygonsDataFrame
pop18<-as_Spatial(pop18)

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
STdat <- cbind(xy_utm$x,xy_utm$y,t)
xyt <- stppp(list(data=STdat, tlim = Time_lim, window = BTShapeF_owin))
xyt # only 294 cases included 
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
xyt$t <- as.integer(xyt$t)
mu_t <- muEst(xyt, f=1/20) # f is the lowess argument for smoothing
mu_t # temporalAtRisk object
plot(mu_t)

# perform polygon overlay operations and compute computational grid

polyolay <- getpolyol(data = xyt, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT)
#polyolay <- getpolyol(data = xyt, regionalcovariates = popDen18$den18, cellwidth = CellWidth, ext = EXT)

## MODEL FORMULAE

FORM <- X ~ moty
FORM_Spatial <- X ~ pop.count -1
FORM_Temporal <- t ~ moty -1

# set the interpolation type for each variable

pop18@data <- guessinterp(pop18@data)
pop18@data <- assigninterp(df = pop18@data, vars =  c( "EA_NUMBER","den18","Total_0.10","area18", "pop.count"), value = "ArealWeightedMean")
class(pop18@data$pop.count)

## DESIGN MATRIX
#pop3 <- setClass(pop3, slots = c(Total_0.10 = "numeric", geometry = "numeric"))
#popDen18_ <- as_Spatial(popDen18_) # changing class type of the object to S4
#xyt <- as_Spatial(xyt)
#Zmat <- getZmat(formula = FORM, data = xyt, regionalcovariates = popDen18, cellwidth = CellWidth, ext = EXT, overl = polyolay)
Zmat <- getZmat(formula = FORM_Spatial, data = xyt, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT, overl = polyolay)
plot(Zmat)
#Zmat <- getZmat(formula = FORM.spatial, data = xyt, cellwidth = NULL, regionalcovariates = pop3, ext = NULL, overl = NULL)

# Specifying log of population counts to comply with Poisson model
# so that number of cases to be proportional to population at risk
# and not the exponential of population

#Zmat[, "den18"] <- log(Zmat[,"den18"])
#Zmat[, "den18"][is.infinite(Zmat[, "den18"])] <- min(Zmat[, "den18"][!is.infinite(Zmat[, "den18"])]) 
#plot(Zmat)

# DEFINING THE OFFSET

mm <- length(attr(Zmat, "mcens"))
nn <- length(attr(Zmat, "ncens"))
Pop.offset <- spatialAtRisk(list(X = attr(Zmat, "mcens"), Y = attr(Zmat, "ncens"), Zm = matrix(Zmat, mm, nn)))

# plot the spatial interpolated covariates

plot(Zmat, ask = F) 

# construct dummy temporal data frame
#days <- c("Saturday", "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
#tvec <- xyt$tlim[1]:xyt$tlim[2]
#da <- rep(days, length.out = length(tvec))
#tdata <- data.frame(t = tvec, dotw = da)

months.of.year <- c("March", "April", "May", "June", "July", "August", "September", "October", "November", "December", "January", "February")
tvec <- xyt$tlim[1]:xyt$tlim[2]
ma <- rep(months.of.year, length.out = length(tvec+30))
tdata <- data.frame(t = tvec, moty = ma)

# choose last time point and number of proceeding time-points to include

tim <- 16844
tim <- as.integer(tim)
LAGLENGTH <- 8
LAGLENGTH <- as.integer(LAGLENGTH)

# bolt on the temporal covariates
ZmatList <- addTemporalCovariates(temporal.formula = FORM_Temporal, T = tim, laglength = LAGLENGTH, tdata = tdata, Zmat = Zmat)

## DEFINING PRIORS

# Defining priors
EtaP <- PriorSpec(LogGaussianPrior(mean = log(c(1,2000,1)), variance = diag(0.2,3)))
BetaP <- PriorSpec(GaussianPrior(mean = rep(0,13), variance = diag(1e+06, 13)))
priors <- lgcpPrior(etaprior = EtaP, betaprior = BetaP)

# set initial values for the algorithm
INITS <- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), betainit = NULL)

# choose the covariance function
CF <- CovFunction(exponentialCovFct)

# bolt on the temporal covariates
#ZmatList <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULA, T = T, laglength = LAGLENGTH, tdata = tdata, Zmat = Zmat)

# run the MCMC algorithm

DIRNAME <- getwd()

SpatioTemporal_Model_01 <- lgcpPredictSpatioTemporalPlusPars(formula = FORM, xyt = xyt, T = tim, laglength = LAGLENGTH, ZmatList = ZmatList, model.priors = priors, 
                             model.inits = INITS, spatial.covmodel = CF, cellwidth = CellWidth, poisson.offset =  Pop.offset, 
                             mcmc.control = mcmcpars(mala.length = 2000, burnin = 100, retain = 10, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                             targetacceptance = 0.574)), output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME,"ST_Model_01"), 
                             forceSave = TRUE)), gradtrunc = Inf, ext = EXT) 

 save(list = ls(), file = file.path(DIRNAME, "ST_Model_01", "ST_Model_01_output.RData")) 
 
# MODEL DIAGNOSTIC CHECKS by plotting the log target to check if the chain has converged to a posterior mode

plot(SpatioTemporal_Model_01)
plot(ltar(SpatioTemporal_Model_01), type = "s", xlab = "Iteration/900", ylab = "log target")

#plot(hvals(SpatioTemporal_Model_01)[2000:5000], type = "l", xlab = "Iteration", ylab = "h")

# compute and plot autocorrelations in the latent field
lagch <- c(1, 5, 15)
Sacf <- autocorr(SpatioTemporal_Model_01, lagch, inWindow = NULL)
library("fields")
for (i in 1:3) {
  image.plot(xvals(SpatioTemporal_Model_01), yvals(SpatioTemporal_Model_01), Sacf[, , i], zlim = c(-1, 1), axes = FALSE, 
             xlab = "", ylab = "", asp = 1, sub = paste("Lag:", lagch[i]), ask = FALSE)
  plot(xyt$window, add = TRUE, ask = FALSE)
  scalebar(5000, label = "5 km")
}

# produce traceplots of beta and eta
traceplots(SpatioTemporal_Model_01, ask = FALSE)

# produce autocorrelation plots of beta and eta
parautocorr(SpatioTemporal_Model_01, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(SpatioTemporal_Model_01)
# the above converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# a text summary of the model parameters
textsummary(SpatioTemporal_Model_01, digits = 4)

# a plot of the prior and posterior densities
priorpost(SpatioTemporal_Model_01, ask = FALSE)

# the posterior covariance function
postcov(SpatioTemporal_Model_01, ask = FALSE)

# exceedance and lower-tail exceedance probabilities
ep <- exceedProbs(c(1.5, 3))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_Model_01, ep)

plotExceed(ex[[1]], "ep", SpatioTemporal_Model_01, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
           sub = "", ask = FALSE)
scalebar(5000, label = "5 km")

 # conditional probabilities
cp <- condProbs(SpatioTemporal_Model_01)

# segregation probabilities
sr <- segProbs(SpatioTemporal_Model_01, domprob = 0.8)

# Inhomogeneous K-Function

kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
kin
plot(kin)

# Estimating the spatial correlation parameters (sigma and phi) of the model

#sigma_phi <- spatialparsEst(kin,sigma.range = c(0,10), phi.range = c(0,10), spatial.covmodel == "matern")
#sigma_phi
#sigma_phi$sigma
#sigma_phi$phi

#Estimating temporal correlation parameter theta

#theta <- thetaEst(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, sigma = 2, phi = 3) # the values for sigma and phi are not estimates...just supplied for now because sigma_phi above is giving an error
