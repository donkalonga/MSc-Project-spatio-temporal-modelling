############
## SET-UP ##
############
setwd("~/Spatio-Temporal Modeling")
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
#library(psych)
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
library(rpanel)
library(fields)
library(miscFuncs)
library(caret)
library(stars)
library(units)
library(udunits2)
library(fasterize)
library(RandomFields)
library(RandomFieldsUtils)
library(splines)
library(graphics)

dat<-read.csv("ST_All-Lineages_Data.csv") %>% 
  dplyr::select(hh_LAT, hh_LNG, sample_collected) %>%
  filter(!is.na(hh_LAT)) #removing records with missing GPS coordinates

# Converting dataset into a space-time planar point pattern object
x <- dat$hh_LNG
y <- dat$hh_LAT
r <- dat$sample_collected
r <- as.Date(r,origin="2015-03-27")
rr <- as.numeric(r, origin="2015-03-27")
rr <- rr - min(rr) + 1
#tt <- as.Date(t,origin="2015-03-27")
#tm <- month(t)
#dat <- dat %>% add_column(r = 'r')
#dat <- dat %>% add_column(rr = 'rr')

tm <- as.integer(factor(format(as.Date(r), "%Y-%m")))
#tm <- unique(tm)
#dat <- dat %>% add_column(tm = as.integer(factor(format(as.Date(r), "%Y-%m"))))
#dat <- dat %>% add_count(sample_collected)

# number of cases per month

#lt <- table(format(as.Date(r), "%Y-%m"))
#lt <- as.vector(lt)
#dat2 <- data.frame(tm1,lt)
#datTempSum <- dat %>% table(cut(r,"month"))
#datTempSum <- datTempSum  %>% as.vector(lt)

# defining spline for long term trend
#dat <- dat %>% add_column(tm = as.integer(factor(format(as.Date(r), "%Y-%m"))))
#cases <- as.vector(table(dat$sample_collected), length = 310)
#dat <- dat %>% group_by(sample_collected) 
#dat <- dat %>% add_column(cases = count(rr))
#dat2 <-dat %>% group_by(sample_collected) %>% dplyr::summarise(nCases=n())
#dat <- cbind(dat2$nCases)
#datTempSum <-dat %>% format(as.Date(r), "%Y-%m")
#datTempSum
#datTempSum %>% ggplot(mapping=aes(x=ymd(sample_collected),y=nCases)) + geom_point()
#datTempSum <- datTempSum %>% mutate(dateNum=as.integer(ymd(sample_collected))) %>% mutate(dateNum2=dateNum-min(dateNum)+1)
#datTempSum

#splineTerms<-bs(datTempSum$dateNum2)
LTspline<-bs(tm)

#head(splineTerms)
#datTempSum<-datTempSum %>% mutate(bsTerm1=splineTerms[,1],bsTerm2=splineTerms[,2],bsTerm3=splineTerms[,3])
#head(datTempSum)
#modTemp<-lm(nCases~bsTerm1+bsTerm2+bsTerm3,data=datTempSum)
#modTemp<-lm(nCases~bs(datTempSum$dateNum2),data=datTempSum)
#modTemp <- lm(n~bs(tm), data = dat)
#newDat<-seq(1,max(datTempSum$dateNum2),length=1000)
# pred<-predict(modTemp)
# head(pred)
#datTempSumModPlot<-data.frame(dateNum=datTempSum$dateNum2,predCases=predCases)
#datTempSumModPlot<-data.frame(dateNum=datTempSum$dateNum2,predCases=pred)
#datTempSumModPlot %>% ggplot(mapping=aes(x=dateNum,y=predCases)) + geom_line()
#datTempSumModPlot %>% ggplot(mapping=aes(x=dateNum,y=predCases)) + geom_line(col="orange",lwd=1.5) + geom_point(data=datTempSum,mapping=aes(x=dateNum2,y=nCases))
#head(datTempSum)
#datTempSumModPlot %>% ggplot(mapping=aes(x=dateNum,y=predCases)) + geom_line()
#datTempSumModPlot %>% ggplot(mapping=aes(x=dateNum,y=predCases)) + geom_line(col="orange",lwd=1.5) + geom_point(data=datTempSum,mapping=aes(x=dateNum2,y=nCases))
#modTemp2<-lm(nCases~bs(dateNum2),data=datTempSum)
#summary(modTemp)

#plot(ltt, xlab = "Months", ylab = "Frequency", type = "l")
#ltt <- bs(tm)
#ltt <- bs(tm,df=20) #df 0f 60 is a stragght line 
#plot(ltt, type = 'l')
#abline(v=c(7,12,16,20),lty=2,col="darkgreen")

#lgcpExtract
#plot(ltt2,type='l')
#lines(x,y)
#tm <- unique(tm)
#tm <- rank(tm)
#tm <- tm+1

# Converting coordinates to utm
xy <- SpatialPoints(cbind(x,y), proj4string=CRS("+proj=longlat"))
xy_utm <-spTransform(xy, CRS("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84"))
xy_utm
plot(xy_utm)

# Defining the time object which will be used for modelling
Time_lim <- range(tm) # first and last value for object t
#Time_lim <- c(1,310)
## population density for Blantyre (2018)
pop18 <- st_read("Blantyre_City.shp") # shapefile from National Statistical Office (NSO)
pop18 <- st_transform(pop18, "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84")
BT_sp <- pop18
EA_2018 <- read.csv("popEA_NSO_2018_BTCity_final.csv") # census data from NSO

pop18$EA_NUMBER <- as.integer(pop18$EA_NUMBER)
pop18 <- left_join(pop18, EA_2018, by = c("ADM_NAME" = "ADM_NAME", "EA_NUMBER" = "EA_NUMBER"))

# Removing NAs
pop18$Total_0.10[is.na(pop18$Total_0.10) & pop18$ADM_TYPE == "PA"] <- 0

# Calculating population density
pop18 <- pop18 %>%
  mutate(area18=units::set_units(st_area(geometry),value="km^2")) %>%
  mutate(den18=Total_0.10/area18)

# Plotting the population of children between 0 and 10 yrs and its density

plot(pop18["Total_0.10"], breaks=c(0,10,50,100,150,200,250,300,400,600,800,1000))
plot(pop18["den18"], breaks=c(0,10,50,100,500,1000,2000,30000,4000,5000,1e4,2e4))

# Owin object preparation for Blantyre City
BTShapeF <- st_read("Blantyre City Boundary-edited.shp")
#BTShapeF <- st_read("MWI_BlantyreCity.shp")
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
table(inside.owin(xy_utm$x,xy_utm$y,BTShapeF_owin)) # 16 cases are outside the owin boundary

#Setting up a spatio-temporal object
STdat <- data.frame(xy_utm$x,xy_utm$y,tm)
xyt <- stppp(list(data=STdat, tlim = Time_lim, window = BTShapeF_owin))
xyt # only 294 cases included 
plot(xyt)

# computing approximate values of the process Y 
# on an appropriately transformed scale
# these approximates are used to inform size of 

# Plotting the spatial density function
denty <- lambdaEst(xyt, axes=T)
plot(denty)

# Coercing the density object to object of class spatialAtRisk for parameter estimation later

Spatrisk  <- spatialAtRisk(denty)
Spatrisk

# Plotting temporal intensity mu(t)
mu_t <- muEst(xyt, f=1/20) # f is the lowess argument for smoothing
plot(mu_t)

##################################################
## SIZE OF COMP GRID & OFFSET OBJECT DEFINITION ##
##################################################

CellWidth <- 350 # Change this value to change computational grid size
EXT <- 2 # to be used during polygon overlay
minimum.contrast(xyt, model  = "exponential", method = "g", intens = density(xyt), transform = log)
chooseCellwidth(xyt,cwinit = CellWidth) # cell width is 175 metres

# convert to SpatialPolygonsDataFrame
pop18 <- as_Spatial(pop18) # with population data
BT_sp <- as_Spatial(BT_sp) # without pop data

# perform polygon overlay operations and compute computational grid
polyolay_off <- getpolyol(data = xyt, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT)
polyolay <- getpolyol(data = xyt, regionalcovariates = BT_sp, cellwidth = CellWidth, ext = EXT)

## FIRST MODEL FORMULAE
FORM <- X ~ Total_0.10 + moty + monthNum
FORM_Spatial <- X ~ Total_0.10
FORM_Temporal <- t ~ moty + monthNum -1

## SECOND MODEL FORMULAE
FORM_2 <- X ~ moty + monthNum
FORM_Spatial_2 <- X ~ 1
FORM_Temporal_2 <- t ~ moty + monthNum -1

# set the interpolation type for each variable
                  
pop18@data <- guessinterp(pop18@data)
pop18@data <- assigninterp(df = pop18@data, vars =  c( "Total_0.10"), value = "ArealWeightedSum")
class(pop18@data$Total_0.10)

# Interpolation of the pop offset onto the computational grid
Zmat <- getZmat(formula = FORM_Spatial_2, data = xyt, regionalcovariates = NULL, cellwidth = CellWidth, ext = EXT, overl = polyolay_off)
Zmat_off <- getZmat(formula = FORM_Spatial, data = xyt, regionalcovariates = pop18, cellwidth = CellWidth, ext = EXT, overl = polyolay)

# Checking which cells in the computational grid have had covariate values imputed

plot(Zmat_off, misscol = 'yellow', obswin = xyt$window)

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
                   spatialAtRisk(list(X = attr(Zmat_off, "mcens"), Y = attr(Zmat_off, "ncens"), Zm = matrix(Zmat_off, mm, nn))))

# plot the spatial interpolated covariates
plot(Zmat, ask = F)
plot(Zmat_off, ask = F)

# construct dummy temporal data frame
#months.of.year <- c("Mar", "Apr" ,"May","Jun", "Jul", "Aug", "Sep" ,"Oct" ,"Nov", "Dec","Jan" ,"Feb")
#tvec <- xyt$tlim[1]:xyt$tlim[2]
#mo <- rep(months.of.year, length.out = length(tvec))
#tdata <- data.frame(t = tvec, moty = mo)


# construct dummy temporal data frame
#months.of.year <- c("Mar", "Apr" ,"May","Jun", "Jul", "Aug", "Sep" ,"Oct" ,"Nov", "Dec","Jan" ,"Feb")
monthsOfTheYear <- c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
#mo <- rep(monthsOfTheYear, length.out = length(dat$sample_collected))
#tvec<-seq(min(ymd(t)),max(ymd(t)),by="days")
#tvec <- month(tvec)
#tvec<-seq(min(ymd(t)),max(ymd(t)),by="month")
tvecNum <- xyt$tlim[1]:xyt$tlim[2]
#tvecNum <- min(tvec):max(tvec)
#mo <- lubridate::month(t, label = T, abbr = T)
#mo <- monthsOfTheYear[week(tvec)]
#tdata <- data.frame(t=tvec, tDate=tvec, motyNum=month(tvec), moty= mo) %>% dplyr::select(t,moty)
#mo <- factor(tdata$moty,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
#mo <- as.character(tdata$moty)
tdata <- data.frame(t=tvecNum, moty= monthsOfTheYear, monthNum = tvecNum)
#tdata$moty <- as.character(tdata$moty)

# choose last time point and number of proceeding time-points to include
ti <- max(tm)
ti <- as.integer(ti)
LAGLENGTH <- 21
LAGLENGTH <- as.integer(LAGLENGTH)

# bolt on the temporal covariates
ZmatList <- addTemporalCovariates(temporal.formula = FORM_Temporal, T = ti, laglength = LAGLENGTH, tdata = tdata, Zmat = Zmat)
#image plot of the ZmatList 
par(mar=c(2.1,1.1,1.1,1.1))
par(mfrow = c(4,4))
for (i in 1:22) {image.plot(ZmatList[[i]])}

# defining eta and beta priors
EtaP <- PriorSpec(LogGaussianPrior(mean = log(c(1,1000,1)), variance = diag(0.2,3)))
BetaP <- PriorSpec(GaussianPrior(mean = rep(0,13), variance = diag(1e+06,13)))
priors <- lgcpPrior(etaprior = EtaP, betaprior = BetaP)

# set initial values for the algorithm
INITS <- lgcpInits(etainit = log(c(sqrt(2.4), 475,1)), betainit = NULL)

# choose the covariance function
CF <- CovFunction(exponentialCovFct)

################################
# RUN THE MCMC ALGORITHM
################################

DIRNAME <- getwd()

SpatioTemporal_Model_All_Cases <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_2, xyt = xyt, T = ti, laglength = LAGLENGTH, ZmatList = ZmatList, model.priors = priors, 
                             model.inits = INITS, spatial.covmodel = CF, cellwidth = CellWidth, poisson.offset =  Pop.offset, 
                             mcmc.control = mcmcpars(mala.length = 5000, burnin = 1000, retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                             targetacceptance = 0.574)), output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME,"ST_Models"), 
                             forceSave = TRUE)), ext = EXT)

 save(list = ls(), file = file.path(DIRNAME, "ST_Models", "ST_Model_All_Cases_output.RData")) 

 
 #####################
 #####################################################
 # LOADING MODEL DATA INSTEAD OF RUNNING THE MCMC AGAIN AND AGAIN
 ##########################################################################
 
 load("ST_Model_All_Cases_output.RData")
 
################################
# MODEL DIAGNOSTIC CHECKS 
################################

# plotting the model
plot(SpatioTemporal_Model_All_Cases, type = "relrisk") # plotting relative risk of the model
plot(SpatioTemporal_Model_All_Cases, col=c(0,0.2,0.4,0.6,0.8,1,1.5,2,3))
plot(SpatioTemporal_Model_All_Cases, type = "serr" ) # plotting standard error of the relative risk of the model
#plot(SpatioTemporal_Model_All_Cases, type = "intensity" ) # intensity plot of the relative risk of the model

# plotting the log target to check if the chain has converged to a posterior mode
plot(ltar(SpatioTemporal_Model_All_Cases), type = "s", xlab = "Iteration/900", ylab = "log target")

# compute and plot autocorrelations in the latent field
lagch <- c(1, 5, 15)
Sacf <- autocorr(SpatioTemporal_Model_All_Cases, lagch, inWindow = NULL)
library("fields")
for (i in 1:3) {
  image.plot(xvals(SpatioTemporal_Model_All_Cases), yvals(SpatioTemporal_Model_All_Cases), Sacf[, , i], zlim = c(-1, 1), axes = FALSE, 
             xlab = "", ylab = "", asp = 1, sub = paste("Lag:", lagch[i]), ask = FALSE)
  plot(xyt$window, add = TRUE, ask = FALSE)
  scalebar(5000, label = "5 km")
}

# traceplots of beta and eta
traceplots(SpatioTemporal_Model_All_Cases, ask = FALSE)

# autocorrelation plots of beta and eta
parautocorr(SpatioTemporal_Model_All_Cases, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(SpatioTemporal_Model_All_Cases)

# the above summary converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# a text summary of the model parameters
textsummary(SpatioTemporal_Model_All_Cases, digits = 4)

# a plot of the prior and posterior densities
priorpost(SpatioTemporal_Model_All_Cases, ask = FALSE)

# the posterior covariance function
postcov(SpatioTemporal_Model_All_Cases, ask = FALSE)

# exceedance and lower-tail exceedance probabilities
ep <- exceedProbs(c(1,1.5,3))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_Model_All_Cases, ep)

plotExceed(ex[[1]], "ep", SpatioTemporal_Model_All_Cases, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
           sub = "", ask = FALSE)
scalebar(5000, label = "5 km")

# Inhomogeneous K-Function

kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
plot(kin)
