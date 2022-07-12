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
#library(psych)
library(data.table)
library(dataMaid)
library(PrevMap)
library(sf)
library(leaflet)
library(maps)
library(broom)
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
library(knitr)
library(kableExtra)
library(gridExtra)
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
library(ggpubr)

# LOADING MERGED DATA

Merg <- read.csv("Merged_Final_DK_13072020.csv")

# Plotting Bar Graph for Sub-Lineages

#ggplot(data = Merg, aes(x=as.factor(lineage)) + geom_bar(stat = "identity", fill=lineage)
colnames(Merg)[12] <- 'Lineage'
Merg$Lineage[Merg$Lineage == "major0"] <- "Clade 0"
Merg$Lineage[Merg$Lineage == "major1"] <- "Clade 1"
Merg$Lineage[Merg$Lineage == "major2"] <- "Clade 2"
Merg$Lineage[Merg$Lineage == "major3"] <- "Clade 3"
Merg$Lineage[Merg$Lineage == "major4"] <- "Clade 4"
Merg$Lineage[Merg$Lineage == "major5"] <- "Clade 5"
Merg$Lineage[Merg$Lineage == "major6"] <- "Clade 6"

###########################################################################
# SUMMARY STATISTICS
###########################################################################

calcSumStat<-function(Merg,var,lvls=NULL,is.id.var=F,is.date.var=F){
  # dat = data frame containing the data for which summary statistics are to be calculated
  # var = character string indicating the name of the variable for which summaries are to be computed
  # lvls = only to be specificed for categorical variables; gives the levels / unique values of the variable
  # is.id.var = (logical); if T the variable is an ID variable and only the number of non missing observations and the corresponding percentage will be computed
  # is.date.var = (logical); if T the variable is a date variable and no mean, SD or median are computed, only the range (min, max)
  
  # total number of observations
  N<-nrow(Merg)
  
  # numeric / non categorical variables
  if(is.null(lvls)){
    n<-sum(!is.na(Merg[,var]))
    n_perc<-n/N
    n_miss<-N-n
    n_miss_perc<-1-n_perc
    
    # numeric variables
    if(!is.id.var & !is.date.var){
      mu<-mean(Merg[,var],na.rm=T)
      sd<-sd(Merg[,var],na.rm=T)
      medianMinMax<-quantile(Merg[,var],probs=c(0.5,0,1),na.rm=T) 
      res<-c(n,n_perc,n_miss,n_miss_perc,mu,sd,medianMinMax)
      
      # date variables
    }else if(is.date.var){
      minMax<-range(ymd(Merg[,var]),na.rm=T)
      res<-c(n,n_perc,n_miss,n_miss_perc,rep("-",3),as.character(minMax))
      
      # ID variables
    }else{
      res<-c(n,n_perc,n_miss,n_miss_perc,rep("-",5)) 
    }
    
    # categorical variables
  }else{
    n<-sum(!is.na(Merg[,var]))
    res<-data.frame(Level=c("-",lvls),
                    n=c(n,rep(NA,length=length(lvls))),
                    n_perc=c(n/N,rep(NA,length=length(lvls))),
                    n_miss=c(N-n,rep(NA,length=length(lvls))),
                    n_miss_perc=c(1-n/N,rep(NA,length=length(lvls))),
                    mu=rep("-",length=1+length(lvls)),
                    sd=rep("-",length=1+length(lvls)),
                    Median=rep("-",length=1+length(lvls)),
                    Min=rep("-",length=1+length(lvls)),
                    Max=rep("-",length=1+length(lvls)),
                    stringsAsFactors = F)
    for(lvl in lvls){
      nLvl<-sum(!is.na(Merg[,var]) & Merg[,var]==lvl)
      res[res$Level==lvl,c("n","n_perc")]<-c(nLvl,nLvl/n)
    }
  }
  
  # output the results
  return(res)
}
# Read in data frame
#Merged_Data <- read.csv("Merged_Final_DK_13072020.csv", stringsAsFactors=F)

Merg$bcdate<-mdy(gsub(Merg$bcdate,pattern=" [0-9]:[0-9][0-9]",replacement="")) # reformatting the bcdate variable into a proper date

# initiate enoty data frame
mVarSum<-data.frame(variable=c("pid_tccb", "pid_epal", "bcnumber", "bcdate", "sample_collected","hh_LAT", "hh_LNG", "Lineage",rep("",7)),
                    Level=c(rep("-",8),paste(sep="","Clade ",0:6)),
                    n=NA,
                    n_perc=NA,
                    Missing=NA,
                    miss_perc=NA,
                    Mean=NA,
                    SD=NA,
                    Median=NA,
                    Minimum=NA,
                    Maximum=NA)

# compute summary statistics for each variable
mVarSum[mVarSum$variable=="pid_tccb",-c(1:2)]<-calcSumStat(Merg,"pid_tccb",is.id.var=T) # explain the use of -c(1:2)
mVarSum[mVarSum$variable=="pid_epal",-c(1:2)]<-calcSumStat(Merg,"pid_epal",is.id.var=T)
mVarSum[mVarSum$variable=="bcnumber",-c(1:2)]<-calcSumStat(Merg,"bcnumber",is.id.var=T)
mVarSum[mVarSum$variable=="bcdate",-c(1:2)]<-calcSumStat(Merg,"bcdate",is.date.var=T)
mVarSum[mVarSum$variable=="sample_collected",-c(1:2)]<-calcSumStat(Merg,"sample_collected",is.date.var=T)
# why is it that the last argument (is.date.var=T) was not used for the two statement below
mVarSum[mVarSum$variable=="hh_LAT",-c(1:2)]<-calcSumStat(Merg,var="hh_LAT")
mVarSum[mVarSum$variable=="hh_LNG",-c(1:2)]<-calcSumStat(Merg,var="hh_LNG")

# explain the code below
mVarSum[is.element(el=mVarSum$variable,set=c("Lineage",rep("",6))),-1]<-calcSumStat(Merg,var="Lineage",lvls=sort(levels(factor(Merg$Lineage))))

# reformatting for printing the table
mVarSum$n_perc<-paste(sep="",format(nsmall=1,round(digits=1,100*as.numeric(mVarSum$n_perc))),"%")
mVarSum$miss_perc<-paste(sep="",format(nsmall=1,round(digits=1,100*as.numeric(mVarSum$miss_perc))),"%")
mVarSum$Mean<-format(nsmall=2,round(digits=2,as.numeric(mVarSum$Mean)))
mVarSum$SD<-format(nsmall=2,round(digits=2,as.numeric(mVarSum$SD)))
mVarSum$Median<-format(nsmall=2,round(digits=2,as.numeric(mVarSum$Median)))
mVarSum$Minimum<-format(nsmall=2,round(digits=2,as.numeric(mVarSum$Minimum)))
mVarSum$Maximum<-format(nsmall=2,round(digits=2,as.numeric(mVarSum$Maximum)))

for(j in 5:ncol(mVarSum)){mVarSum[grepl(mVarSum[,j],pattern="NA") | is.na(mVarSum[,j]),j]<-"-"}

cnames<-colnames(mVarSum)
cnames[is.element(el=cnames,set=c("n_perc","miss_perc"))]<-"%"

# printing the table
kable(mVarSum, "latex", booktabs = T, col.names = cnames, align = c("l","c",rep("r",ncol(mVarSum)-2)), linesep="") %>%
  kable_styling(latex_options = c("stripped","scale_down"), position = "center")

###########################################################################
# BAR PLOT FOR SUB-LINEAGES
##########################################################################

png(filename="Sub-lineage.png", width = 6, height = 3, units = "in", res = 1800)
Merg2 <- Merg[!is.na(Merg$Lineage),] #%>% na.omit(lineage)
SubLinPlot <- ggplot(Merg2, aes(x=Lineage, fill = Lineage, cex.names = 0.5)) +
  geom_bar() + xlab("H58 Sub-lineages") + ylab("Count") + theme_classic() + theme(legend.position="none")
#barplot(table(Merg$lineage), col = "red", ylab = "Frequency",cex.lab=1, cex.axis=1, cex.main=1, cex.names=1)
print(SubLinPlot)
dev.off()

###############################################################################
# Graph of Cumulative Cases Over Time Per Sub-Lineage
###############################################################################

dfTotal<-Merg %>% group_by(ymd(sample_collected)) %>% tally()
colnames(dfTotal)[1]<-"date"

g1<-ggplot(dfTotal,mapping=aes(x=ymd(date), y=cumsum(n), lty = 'Typhoid fever cases')) + theme_classic() + 
  geom_line(lwd=1.5,color="red") +
  xlab("Date") +
  ylab("Cumulative Cases")

g1

# Graph of Cumulative Cases Over Time Per Sub-Lineage

dfLineage<-list()
#Merg %>% ymd(Merg$sample_collected)
for(lvl in sort(levels(factor(Merg$Lineage)))){
  dfLineage[[lvl]]<-Merg %>%
    dplyr::mutate() %>%
    dplyr::group_by(Lineage, sample_collected) %>%
    dplyr::filter(Lineage==lvl) %>%
    #dplyr::arrange(ymd(sample_collected)) %>%
    dplyr::tally() %>%
    #dplyr::summarise(n=n(), .groups = "drop") %>%
    dplyr::mutate(cs=cumsum(n)) %>% 
    dplyr::ungroup()
  
  colnames(dfLineage[[lvl]])[2]<-"date"
}

dfLineage<-bind_rows(dfLineage)

#Merg <- Merg %>% group_by(Lineage) %>% tally() %>% dplyr::mutate(cs = cumsum(n)) %>% ungroup

g2<-ggplot(dfLineage,mapping=aes(x=ymd(date), y=cs, color=Lineage)) +
  scale_color_manual(values=c("salmon","mediumorchid","orange","greenyellow","burlywood4","steelblue4","darkolivegreen4"), name="H58 Sub-Lineages") + 
  geom_line(lwd=1) + theme_classic() +
  xlab("Date") +
  ylab("Cumulative Cases") #+
  #ggtitle("Cumulative case counts over time per sub-lineage.") +
  #theme(legend.position = "right")

g2
######################################################################################
# Cumulative Cases Over Time By Sub-Lineage
#######################################################################################

png(filename="Barplot for Sub-lineages.png", width = 6, height = 3, units = "in", res = 1800)
print(SubLinPlot)
dev.off()

#################################################################################
# PLOTTING CUMULATIVE CASES OVER TIME ADJACENT TO GENOMIC SUB-LINEAGES
#################################################################################
png(filename="Cumulative Cases Over Time (All Cases) and Cumulative Cases Over Time By Sub-Lineage.png", width = 6, height = 4, units = "in", res = 1800)
#grid.arrange(g1,g2,ncol=2)
ggarrange(g1, g2, ncol = 1, nrow = 2)
dev.off()
########################################################################################
# Spatial plots by sub-lineages

BTShapeF <- st_read("Blantyre City Boundary-edited.shp")
BTShapeF <- st_transform(BTShapeF, "+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84")
BTShapeF_fortified <- tidy(BTShapeF)

g3 <- ggplot(Merg, mapping = aes(x=hh_LAT, y=hh_LNG, color = factor(Lineage))) + 
            geom_point() + stat_ellipse() + geom_sf(data=BTShapeF$geometry)
g3


register_google(key = "AIzaSyDjuNQO7WVjfHrInkqE7ZQ06MZyoTB9bS8")
#register_google(key="AIzaSyAZvE95bIaBBAHcorKibY1P8QcsNKqahKo")
#BT <- get_map(location = "blantyre")
BT <- get_map(location = c(lon = 35.041378331851824, lat = -15.772877620991228), zoom = "auto", scale = "auto", maptype = "roadmap", source = "google", force = ifelse(source == "google", TRUE, FALSE))
Lineage <- Merged_Data$lineage
qmplot(hh_LNG, hh_LAT, data = Merged_Data, col = Lineage, na.rm = T,xlab = "Lon (deg)", ylab = "Lat (deg)", title = "Spatial Plot by Lineage")

######################################################################
# SPATIO-TEMPORAL MODEL WITH ALL CASES
######################################################################

setwd("~/Don/MSc in Biostatistics/Reading Materials/THESIS/Spatio-Temporal Modelling/Books - Spatio-Temporal Data Analysis")
rm(list = ls())

load("ST_Model_All_Cases_output.RData")

# Relative Risk

png(filename="Relative Risk Plot - SP - All Cases.png", width = 6, height = 5, units = "in", res = 1800)
par(mar=c(2.1,2.1,1.1,4.1))
par(mfrow = c(6,4))
plot(SpatioTemporal_Model_All_Cases, type = "relrisk", cex.axis = 0.5, cex.lab = 0.5, xlab = "", ylab = "") # plotting relative risk of the model for every time point
dev.off()
dev.off()

# Standard Error of the relative risk

png(filename="Standard Errors Plot - SP - All Cases.png", width = 6, height = 3, units = "in", res = 1800)
plot(SpatioTemporal_Model_All_Cases, type = "serr" ) # plotting standard error of the relative risk of the model for every time point
dev.off()

#plot(SpatioTemporal_Model_All_Cases, type = "intensity" ) # plotting the mean Poisson intensity

# plotting the log target to check if the chain has converged to a posterior mode

png(filename = "Log Target Plot - ST - All Cases.png", width = 5, height = 3, units = "in", res = 1800)
plot(ltar(SpatioTemporal_Model_All_Cases), type = "s", xlab = "Iteration/90", ylab = "log target")
dev.off()

#plot(ltar(SpatioTemporal_Model_All_Cases), type = "s", xlab = "Iteration/90", ylab = "log target")

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

png(filename = "Traceplots for Beta and Eta - All Cases.png", width = 6, height = 4, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(3,3))
traceplots(SpatioTemporal_Model_All_Cases, ask = FALSE)
dev.off()

# autocorrelation plots of beta and eta

png(filename = "Autocorrelation of Beta and Eta - All Cases.png", width = 6, height = 5, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(4,3))
parautocorr(SpatioTemporal_Model_All_Cases, ask = FALSE)
dev.off()

# a summary table of beta and eta
parsum <- parsummary(SpatioTemporal_Model_All_Cases)

# the above summary converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# plotting the long term trend of the model

tt <- seq(1,22, length=1000)
DataF <- cbind(tt, bs(tt))
#yy <- area(BTShapeF) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
yy <- sum((polyolay$gridobj$cellarea) * (polyolay$gridobj$cellInside)) * exp(sum(Pop.offset[[1]]$Zm)) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
plotDf <- data.frame(t=tt,Y=yy)
plotDf %>% ggplot(mapping=aes(x=t,y=Y)) + geom_line(col="orange",lwd=2) + theme_light() + xlab("Month of the Year (March 2015 = 1)") + ylab("number of typhoid cases per month") + geom_point(data=ltDF,mapping=aes(x=tmu,y=lt),size=2)

# a text summary of the model parameters
textsummary(SpatioTemporal_Model_All_Cases, digits = 4)

# a plot of the prior and posterior densities
priorpost(SpatioTemporal_Model_All_Cases, ask = FALSE)

# the posterior covariance function
postcov(SpatioTemporal_Model_All_Cases, ask = FALSE)

# exceedance and lower-tail exceedance probabilities

png(filename = "Exceedance Probabilities - All Cases.png", width = 5, height = 3, units = "in", res = 1800)
par(mar=c(4.1,1.1,1.1,5.1))
par(mfrow=c(1,2))
ep <- exceedProbs(c(2,4))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_Model_All_Cases, ep)
plotExceed(ex[[1]], "ep", SpatioTemporal_Model_All_Cases, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
           sub = "", ask = FALSE)
scalebar(5000, label = "5 km")
dev.off()

# Inhomogeneous K-Function

kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
plot(kin)

# Extracting mean and variance of the latent field
meanfield(SpatioTemporal_Model_All_Cases)
varfield(SpatioTemporal_Model_All_Cases)

# relative risk (mean of exp{Y}) and standard error of the relative risk
rr(SpatioTemporal_Model_All_Cases)
serr(SpatioTemporal_Model_All_Cases)

######################################################################
# SPATIO-TEMPORAL MODEL FOR MAJOR 0 SUB-LINEAGE
######################################################################

setwd("~/Don/MSc in Biostatistics/Reading Materials/THESIS/Spatio-Temporal Modelling/Books - Spatio-Temporal Data Analysis")
rm(list = ls())

load("ST_Model_major0_output.RData")

# plotting the model

par(mar=c(4.1,4.1,1.1,4.1))
par(mfrow = c(1,3))

#plot(SpatioTemporal_Model_major0) # plotting relative risk of the model for every time point
plot(SpatioTemporal_Model_major0, type = "relrisk") # plotting relative risk of the model for every time point
plot(SpatioTemporal_Model_major0, type = "serr" ) # plotting standard error of the relative risk of the model for every time point


# plotting the log target to check if the chain has converged to a posterior mode
#plot(ltar(SpatioTemporal_Model_major0), type = "s", xlab = "Iteration/90", ylab = "log target")

png(filename = "Log Target Plot - ST - Major 0.png", width = 6, height = 4, units = "in", res = 1800)
plot(ltar(SpatioTemporal_Model_major0), type = "s", xlab = "Iteration/90", ylab = "log target")
dev.off()

# compute and plot autocorrelations in the latent field
lagch <- c(1, 5, 15)
Sacf <- autocorr(SpatioTemporal_Model_major0, lagch, inWindow = NULL)
library("fields")
for (i in 1:3) {
  image.plot(xvals(SpatioTemporal_Model_major0), yvals(SpatioTemporal_Model_major0), Sacf[, , i], zlim = c(-1, 1), axes = FALSE, 
             xlab = "", ylab = "", asp = 1, sub = paste("Lag:", lagch[i]), ask = FALSE)
  plot(xyt$window, add = TRUE, ask = FALSE)
  scalebar(5000, label = "5 km")
}

# traceplots of beta and eta
traceplots(SpatioTemporal_Model_major0, ask = FALSE)

# autocorrelation plots of beta and eta
parautocorr(SpatioTemporal_Model_major0, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(SpatioTemporal_Model_major0)

# the above summary converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# plotting the long term trend of the model

tt <- seq(1,19, length=1000)
DataF <- cbind(tt, bs(tt))
#yy <- area(BTShapeF) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
yy <- sum((polyolay$gridobj$cellarea) * (polyolay$gridobj$cellInside)) * exp(sum(Pop.offset[[1]]$Zm)) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
plotDf <- data.frame(t=tt,Y=yy)
plotDf %>% ggplot(mapping=aes(x=t,y=Y)) + geom_line(col="orange",lwd=2) + theme_light() + xlab("Month of the Year (April 2015 = 1)") + ylab("number of typhoid cases per month") + geom_point(data=ltDF,mapping=aes(x=tmu,y=lt),size=2)


# a text summary of the model parameters
textsummary(SpatioTemporal_Model_major0, digits = 4)

# a plot of the prior and posterior densities
priorpost(SpatioTemporal_Model_major0, ask = FALSE)

# the posterior covariance function
postcov(SpatioTemporal_Model_major0, ask = FALSE)

# exceedance and lower-tail exceedance probabilities

png(filename = "Exceedance Probabilities - Major 0.png", width = 5, height = 3, units = "in", res = 1800)
par(mar=c(4.1,1.1,1.1,5.1))
par(mfrow=c(1,2))
ep <- exceedProbs(c(2,4))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_Model_major0, ep)
plotExceed(ex[[1]], "ep", SpatioTemporal_Model_major0, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
           sub = "", ask = FALSE)
scalebar(5000, label = "5 km")
dev.off()

# Inhomogeneous K-Function

kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
plot(kin)


######################################################################
# SPATIO-TEMPORAL MODEL FOR MAJOR 2 SUB-LINEAGE
######################################################################

setwd("~/Don/MSc in Biostatistics/Reading Materials/THESIS/Spatio-Temporal Modelling/Books - Spatio-Temporal Data Analysis")
rm(list = ls())

load("ST_Model_major2_output.RData")

# plotting the model

par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow = c(1,2))

#plot(SpatioTemporal_Model_major0) # plotting relative risk of the model for every time point
plot(SpatioTemporal_Model_major2, type = "relrisk") # plotting relative risk of the model for every time point
plot(SpatioTemporal_Model_major2, type = "serr" ) # plotting standard error of the relative risk of the model for every time point

# plotting the log target to check if the chain has converged to a posterior mode
#plot(ltar(SpatioTemporal_Model_major2), type = "s", xlab = "Iteration/90", ylab = "log target")

png(filename = "Log Target Plot - ST - Major 2.png", width = 6, height = 4, units = "in", res = 1800)
plot(ltar(SpatioTemporal_Model_major2), type = "s", xlab = "Iteration/90", ylab = "log target")
dev.off()

# compute and plot autocorrelations in the latent field
lagch <- c(1, 5, 15)
Sacf <- autocorr(SpatioTemporal_Model_major2, lagch, inWindow = NULL)
library("fields")
for (i in 1:3) {
  image.plot(xvals(SpatioTemporal_Model_major2), yvals(SpatioTemporal_Model_major2), Sacf[, , i], zlim = c(-1, 1), axes = FALSE, 
             xlab = "", ylab = "", asp = 1, sub = paste("Lag:", lagch[i]), ask = FALSE)
  plot(xyt$window, add = TRUE, ask = FALSE)
  scalebar(5000, label = "5 km")
}

# traceplots of beta and eta
traceplots(SpatioTemporal_Model_major2, ask = FALSE)

# autocorrelation plots of beta and eta
parautocorr(SpatioTemporal_Model_major2, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(SpatioTemporal_Model_major2)

# the above summary converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# plotting the long term trend of the model

tt <- seq(1,15, length=1000)
DataF <- cbind(tt, bs(tt))
#yy <- area(BTShapeF) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
yy <- sum((polyolay$gridobj$cellarea) * (polyolay$gridobj$cellInside)) * exp(sum(Pop.offset[[1]]$Zm)) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
plotDf <- data.frame(t=tt,Y=yy)
plotDf %>% ggplot(mapping=aes(x=t,y=Y)) + geom_line(col="orange",lwd=2) + theme_light() + xlab("Month of the Year (May 2015 = 1)") + ylab("number of typhoid cases per month") + geom_point(data=ltDF,mapping=aes(x=tmu,y=lt),size=2)


# a text summary of the model parameters
textsummary(SpatioTemporal_Model_major2, digits = 4)

# a plot of the prior and posterior densities
priorpost(SpatioTemporal_Model_major2, ask = FALSE)

# the posterior covariance function
postcov(SpatioTemporal_Model_major2, ask = FALSE)

# exceedance and lower-tail exceedance probabilities

png(filename = "Exceedance Probabilities - Major 2.png", width = 5, height = 3, units = "in", res = 1800)
par(mar=c(4.1,1.1,1.1,5.1))
par(mfrow=c(1,2))
ep <- exceedProbs(c(2,4))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_Model_major2, ep)
plotExceed(ex[[1]], "ep", SpatioTemporal_Model_major2, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
           sub = "", ask = FALSE)
scalebar(5000, label = "5 km")
dev.off()

# Inhomogeneous K-Function

kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
plot(kin)

######################################################################
# SPATIO-TEMPORAL MODEL FOR MAJOR 1,3,4,5&6 SUB-LINEAGE
######################################################################

setwd("~/Don/MSc in Biostatistics/Reading Materials/THESIS/Spatio-Temporal Modelling/Books - Spatio-Temporal Data Analysis")
rm(list = ls())

load("ST_Model_major13456_output.RData")

# plotting the model

png(filename = "Relative Risk Plot - SP - Major 13456.png", width = 6, height = 5, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(5,4))
plot(SpatioTemporal_Model_major13456 , type = "relrisk", xlab = "Eastings", ylab = "Northings")
dev.off()
dev.off()


png(filename = "Standard Errors Plot - SP - Major 13456.png", width = 6, height = 5, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(5,4))
plot(SpatioTemporal_Model_major13456 , type = "serr", xlab = "Eastings", ylab = "Northings")
dev.off()
dev.off()

#plot(SpatioTemporal_Model_major13456, type = "relrisk") # plotting relative risk of the model for every time point
#plot(SpatioTemporal_Model_All_Cases, col=c(0,0.2,0.4,0.6,0.8,1,1.5,2,3))
#plot(SpatioTemporal_Model_major13456, type = "serr" ) # plotting standard error of the relative risk of the model for every time point
#plot(SpatioTemporal_Model_All_Cases, type = "intensity" ) # plotting the mean Poisson intensity

# plotting the log target to check if the chain has converged to a posterior mode
#plot(ltar(SpatioTemporal_Model_major13456), type = "s", xlab = "Iteration/90", ylab = "log target")

png(filename = "Log Target Plot - ST - Major 13456.png", width = 6, height = 4, units = "in", res = 1800)
plot(ltar(SpatioTemporal_Model_major13456), type = "s", xlab = "Iteration/90", ylab = "log target")
dev.off()


# compute and plot autocorrelations in the latent field

png(filename = "Autocorrelations in the Latent Field - Major 13456.png", width = 6, height = 4, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,5.1))
par(mfrow=c(1,3))
lagch <- c(1, 5, 15)
Sacf <- autocorr(SpatioTemporal_Model_major13456, lagch, inWindow = NULL)
library("fields")
for (i in 1:3) {
  image.plot(xvals(SpatioTemporal_Model_major13456), yvals(SpatioTemporal_Model_major13456), Sacf[, , i], zlim = c(-1, 1), axes = FALSE, 
             xlab = "", ylab = "", asp = 1, sub = paste("Lag:", lagch[i]), ask = FALSE)
  plot(xyt$window, add = TRUE, ask = FALSE)
  scalebar(5000, label = "5 km")
}
dev.off()

# traceplots of beta and eta
traceplots(SpatioTemporal_Model_major13456, ask = FALSE)

# autocorrelation plots of beta and eta
parautocorr(SpatioTemporal_Model_major13456, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(SpatioTemporal_Model_major13456)

# the above summary converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# plotting the long term trend of the model

tt <- seq(1,17, length=1000)
DataF <- cbind(tt, bs(tt))
#yy <- area(BTShapeF) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
yy <- sum((polyolay$gridobj$cellarea) * (polyolay$gridobj$cellInside)) * exp(sum(Pop.offset[[1]]$Zm)) * exp(log(parsum[4,1])+log(parsum[5,1])*DataF[,2] + log(parsum[6,1])*DataF[,3] + log(parsum[7,1])*DataF[,4])
plotDf <- data.frame(t=tt,Y=yy)
plotDf %>% ggplot(mapping=aes(x=t,y=Y)) + geom_line(col="orange",lwd=2) + theme_light() + xlab("Month of the Year (March 2015 = 1)") + ylab("number of typhoid cases per month") + geom_point(data=ltDF,mapping=aes(x=tmu,y=lt),size=2)

# a text summary of the model parameters
textsummary(SpatioTemporal_Model_major13456, digits = 4)

# a plot of the prior and posterior densities
priorpost(SpatioTemporal_Model_major13456, ask = FALSE)

# the posterior covariance function
postcov(SpatioTemporal_Model_major13456, ask = FALSE)

# exceedance and lower-tail exceedance probabilities

png(filename = "Exceedance Probabilities - Major 13456.png", width = 5, height = 3, units = "in", res = 1800)
par(mar=c(4.1,1.1,1.1,5.1))
par(mfrow=c(1,2))
ep <- exceedProbs(c(2,4))
ex <- lgcp:::expectation.lgcpPredict(SpatioTemporal_Model_major13456, ep)
plotExceed(ex[[1]], "ep", SpatioTemporal_Model_major13456, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
           sub = "", ask = FALSE)
scalebar(5000, label = "5 km")
dev.off()

# Inhomogeneous K-Function

kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
plot(kin)

# Extracting mean and variance of the latent field
meanfield(SpatioTemporal_Model_major13456)
varfield(SpatioTemporal_Model_major13456)

# relative risk (mean of exp{Y}) and standard error of the relative risk
rr(SpatioTemporal_Model_major13456)
serr(SpatioTemporal_Model_major13456)
#intens(SpatioTemporal_Model_All_Cases) # estimated cell-wise mean Poisson intensity

##################################################################
# MULTI-TYPE SPATIAL MODEL
##################################################################

setwd("~/Don/MSc in Biostatistics/Reading Materials/THESIS/Spatio-Temporal Modelling/Books - Spatio-Temporal Data Analysis")
rm(list = ls())

load("Multitype_Spatial_Model_output.RData")

# MODEL DIAGNOSTIC CHECKS by plotting the log target to check if the chain has converged to a posterior mode

png(filename = "Relative Risk - Multi-type.png", width = 6, height = 1.5, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(1,3))
plot(Multitype_Spatial_Model, xlab = "Eastings", ylab = "Northings")
dev.off()
dev.off()

png(filename = "Log Target - Multi-type.png", width = 5, height = 3, units = "in", res = 1800)
plot(ltar(Multitype_Spatial_Model), type = "s", xlab = "Iteration/18000", ylab = "log target")
dev.off()

#plot(hvals(SpatioTemporal_Model_01)[2000:5000], type = "l", xlab = "Iteration", ylab = "h")

# compute and plot autocorrelations in the latent field

png(filename = "Autocorrelation in the latent field - Multi-type.png", width = 6, height = 6, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,5.1))
par(mfrow=c(3,3))
for (i in 1:3) {
  Y_i <- autocorrMultitype(Multitype_Spatial_Model, c(1, 5, 15), i, inWindow = NULL)
  plot(Y_i, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "", ask = F)
}
dev.off()

# produce traceplots of beta and eta

png(filename = "Traceplot - Multi-type.png", width = 6, height = 4, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(3,4))
traceplots(Multitype_Spatial_Model, ask = FALSE)
dev.off()

# produce autocorrelation plots of beta and eta

png(filename = "Autocorrelation of Beta and Eta - Multi-type.png", width = 6, height = 5, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(4,3))
parautocorr(Multitype_Spatial_Model, ask = FALSE)
dev.off()

# a summary table of beta and etam
parsum <- parsummary(Multitype_Spatial_Model)
parsum
# the above converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# a text summary of the model parameters
textsummary(Multitype_Spatial_Model, digits = 4)

# a plot of the prior and posterior densities
priorpost(Multitype_Spatial_Model, ask = FALSE)

# the posterior covariance function

png(filename = "Posterior covariance function - Multi-type.png", width = 6, height = 6, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,1.1))
par(mfrow=c(2,2))
postcov(Multitype_Spatial_Model, ask = FALSE)
dev.off()

# conditional probabilities

png(filename = "Conditional probability - Multi-type.png", width = 6, height = 1.5, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,5.1))
par(mfrow=c(1,3))
cp <- condProbs(Multitype_Spatial_Model)
plot(cp, xlab = "", ylab = "")
dev.off()
dev.off()

# segregation probabilities
png(filename = "Segregation probability - Multi-type.png", width = 6, height = 1.5, units = "in", res = 1800)
par(mar=c(4.1,4.1,1.1,5.1))
par(mfrow=c(1,3))
sp <- segProbs(Multitype_Spatial_Model, domprob = 0.8)
plot(sp, xlab = "", ylab = "")
dev.off()
dev.off()

# exceedance and lower-tail exceedance probabilities
#ep <- exceedProbs(c(1.5, 3))
#ex <- lgcp:::expectation.lgcpPredict(Multitype_Spatial_Model, ep)

#plotExceed(ex[[1]], "ep", SpatioTemporal_Model_All_Cases, zlim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "", 
#           sub = "", ask = FALSE)
#scalebar(5000, label = "5 km")

# conditional probabilities
#cp <- condProbs(SpatioTemporal_Model_01)

# segregation probabilities
#sr <- segProbs(SpatioTemporal_Model_01, domprob = 0.8)

# Inhomogeneous K-Function

#kin <- KinhomAverage(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, time.window = xyt$tlim, rvals = NULL, correction = "iso", suppresswarnings = F) # using the inhomogeneous K function
#kin
#plot(kin)

# Estimating the spatial correlation parameters (sigma and phi) of the model

#sigma_phi <- spatialparsEst(kin,sigma.range = c(0,12), phi.range = c(0,12), spatial.covmodel == "matern")
#sigma_phi
#sigma_phi$sigma
#sigma_phi$phi

#Estimating temporal correlation parameter theta

#theta <- thetaEst(xyt, spatial.intensity = Spatrisk, temporal.intensity = mu_t, sigma = 2, phi = 3) # the values for sigma and phi are not estimates...just supplied for now because sigma_phi above is giving an error