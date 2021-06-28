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
source("lgcpPredictMultitypeSpatialPlusParsDonKalonga.R")
source("MALAlgcpMultitypeSpatial.PlusParsDonKalonga.R")

dat<-read.csv("Multitype_Model_Data.csv") %>% dplyr::select(hh_LAT, hh_LNG, lineage)

dat$lineage[dat$lineage=="major4"] <- "major99"
dat$lineage[dat$lineage=="major5"] <- "major99" # combining major4 and major5 into a new lineage
dim(dat) # 255 cases

# Converting dataset into a space-time planar point pattern object
x <- dat$hh_LNG
y <- dat$hh_LAT
maks <- dat$lineage
maks <- as.factor(maks)
#t <- dat$sample_collected
#t <- as.Date(t,origin="2015-03-27")
#tm <- as.integer(t - min(t))
#tm <- tm+1
# Converting coordinates to utm
xy <- SpatialPoints(cbind(x,y), proj4string=CRS("+proj=longlat"))
xy_utm <-spTransform(xy, CRS("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84"))
xy_utm
plot(xy_utm)

# Owin object preparation for Blantyre City

BTShapeF <- st_read("Blantyre City Boundary-edited.shp")
#BTShapeF <- st_bbox(710385,8241972,727770,8263471) 
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
#inside.owin(xy_utm$x,xy_utm$y,BTShapeF_owin) # 16 cases are outside the owin boundary

#Setting up a multitype spatial object
#MTSdat <- cbind(xy_utm$x,xy_utm$y,maks)
MTSdat <- cbind(xy_utm$x,xy_utm$y,maks)
xym <- ppp(xy_utm$x, xy_utm$y, marks = maks, window = BTShapeF_owin)
xym # only 239 cases included 
#xyt$t <- as.integer(xyt$t)
plot(xym)

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

# computing approximate values of the process Y 
# on an appropriately transformed scale
# these approximates are used to inform size of 
# computational grid so as to capture the dependence
# properties of Y
# It is also used to inform the proposal kernel for the MCMC algorithm

# Plotting the spatial density function
denty <- lambdaEst(xym, axes=T)
plot(denty)

# minimum contrast estimation

ttx <- xym$x[xym$marks == "major0" | xym$marks == "major1" | xym$marks == 
                   "major2" | xym$marks == "major3" | xym$marks == "major6" | xym$marks == "major99"]
tty <- xym$y[xym$marks == "major0" | xym$marks == "major1" | xym$marks == 
                   "major2" | xym$marks == "major3" | xym$marks == "major6" | xym$marks == "major99"]
ttm <- as.character(xym$marks[xym$marks == "major0" | xym$marks == "major1" | xym$marks == 
                                    "major2" | xym$marks == "major3" | xym$marks == "major6" | xym$marks == "major99"])
tempppp <- ppp(x = ttx, y = tty, window = xym$window, marks = as.factor(ttm))

denls <- lapply(as.character(levels(tempppp$marks)), function(x) {
  density.ppp(tempppp[tempppp$marks == x])
})

mn <- minimum.contrast(tempppp, model = "exponential", method = "g", intens = denls, transform = log)

###########################################
## SIZE OF COMP GRID & OFFSET OBJECT DEFINITION ##
###########################################

CellWidth <- 400 # Change this value to change computational grid size
EXT <- 2 # to be used during polygon overlay
chooseCellwidth(xym,cwinit = CellWidth) 

#pop18 <- pop18 %>% mutate(pop.count=(den18*CellWidth^2)/1e6)

# convert to SpatialPolygonsDataFrame
pop18 <- as_Spatial(pop18)
BT_sp <- as_Spatial(BT_sp)

# perform polygon overlay operations and compute computational grid

polyolay_off <- getpolyol(data = xym, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT)
polyolay <- getpolyol(data = xym, regionalcovariates = BT_sp, cellwidth = CellWidth, ext = EXT)

# the statistical model for the main effets, $\beta$
#form_list <- formulaList(list(major0 ~ Total_0.10, major1 ~ Total_0.10, major2 ~ Total_0.10, major3 ~ Total_0.10, major4 ~ Total_0.10, major5 ~ Total_0.10, major6 ~ Total_0.10))
form_list <- formulaList(list(major0 ~ 1, major1 ~ 1, major2 ~ 1, major3 ~ 1, major6 ~ 1, major99 ~ 1))

## MODEL FORMULAE FOR INTERPOLATION

FORM <- X ~ Total_0.10
FORM_1 <- X ~ 1
# set the interpolation type for each variable

pop18@data <- guessinterp(pop18@data)
pop18@data <- assigninterp(df = pop18@data, vars =  c( "Total_0.10"), value = "ArealWeightedSum")
class(pop18@data$Total_0.10)

#pop3 <- setClass(pop3, slots = c(Total_0.10 = "numeric", geometry = "numeric"))
#popDen18_ <- as_Sp
## DESIGN MATRIXatial(popDen18_) # changing class type of the object to S4
#xyt <- as_Spatial(xyt)
#Zmat <- getZmat(formula = FORM, data = xyt, regionalcovariates = popDen18, cellwidth = CellWidth, ext = EXT, overl = polyolay)

Zmat_off <- getZmat(formula = FORM, data = xym, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT, overl = polyolay)
Zmat <- getZmat(formula = FORM_1, data = xym, regionalcovariates = NULL, cellwidth = CellWidth, ext = EXT, overl = polyolay_off)
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

#Pop.offset <- spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn)))

Pop.offset <- list(spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))),
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))))

# plot the spatial interpolated covariates
plot(Zmat, ask = F)
plot(Zmat_off, ask = F)


## DEFINING PRIORS

# define the priors
pr.mn <- log(c(1, 1000))
pr.vr <- c(0.2, 0.05)
priors <- list()
for (i in 1:6) {
  priors[[i]] <- lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = pr.mn, 
                  variance = diag(pr.vr))), betaprior = PriorSpec(GaussianPrior(mean = rep(0,7), variance = diag(10^6, 7))))
}
priors[[7]] <- lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = pr.mn, variance = diag(pr.vr))), betaprior = NULL)
#priors <- lgcpPrior(etaprior = PriorSpec(LogGaussianPrior(mean = pr.mn, variance = diag(pr.vr))), betaprior = NULL)

# set initial values for the algorithm
#INITS <- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), betainit = NULL)

# choose the covariance function

cfs <- list()
for (i in 1:6) {
  cfs[[i]] <- CovFunction(exponentialCovFct)
}
cfs[[7]] <- CovFunction(RandomFieldsCovFct(model = "matern", additionalparameters = 1))

#cfs <- CovFunction(exponentialCovFct)
# run the MCMC algorithm

DIRNAME <- getwd()

load("lgcpPredictMultitypeSpatialPlusParsDonKalonga.R")
load("MALAlgcpMultitypeSpatial.PlusParsDonKalonga.R")

Multitype_Spatial_Model <- lgcpPredictMultitypeSpatialPlusPars(formulaList = form_list, sd = xym, Zmat = Zmat, model.priorsList = priors, 
                             spatial.covmodelList = cfs, cellwidth = CellWidth, poisson.offset =  Pop.offset, 
                             mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                             targetacceptance = 0.574)), output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME,"ST_Models"), 
                             forceSave = TRUE)), ext = EXT) 

 save(list = ls(), file = file.path(DIRNAME, "ST_Models", "Multitype_Spatial_Model_output.RData")) 
 
# MODEL DIAGNOSTIC CHECKS by plotting the log target to check if the chain has converged to a posterior mode

plot(Multitype_Spatial_Model)
plot(ltar(Multitype_Spatial_Model), type = "s", xlab = "Iteration/900", ylab = "log target")

#plot(hvals(SpatioTemporal_Model_01)[2000:5000], type = "l", xlab = "Iteration", ylab = "h")

# compute and plot autocorrelations in the latent field
for (i in 0:6) {
  Y_i <- autocorrMultitype(lg, c(1, 5, 15), i, inWindow = NULL)
  plot(Y_i, zlim = c(-1, 1), xlab = "", ylab = "", ask = FALSE)
}
# produce traceplots of beta and eta
traceplots(Multitype_Spatial_Model, ask = FALSE)

# produce autocorrelation plots of beta and eta
parautocorr(Multitype_Spatial_Model, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(Multitype_Spatial_Model)
# the above converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# a text summary of the model parameters
textsummary(Multitype_Spatial_Model, digits = 4)

# a plot of the prior and posterior densities
priorpost(Multitype_Spatial_Model, ask = FALSE)

# the posterior covariance function
postcov(Multitype_Spatial_Model, ask = FALSE)

# conditional probabilities
cp <- condProbs(Multitype_Spatial_Model)

# segregation probabilities
sr <- segProbs(Multitype_Spatial_Model, domprob = 0.8)

# exceedance and lower-tail exceedance probabilities
ep <- exceedProbs(c(1.5, 3))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_Model_All_Cases, ep)

plotExceed(ex[[1]], "ep", SpatioTemporal_Model_All_Cases, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
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

