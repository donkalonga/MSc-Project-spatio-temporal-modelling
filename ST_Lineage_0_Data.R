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
dat<-dat[which(dat$lineage =="major0"),]
dat<-dat[!is.na(dat$sample_collected),] #removing cases without date
dim(dat) # 540 cases with 15 variables
dat$sample_collected<-as.Date(dat$sample_collected,origin="2015-03-28") # changing object class from Character to Date
#dat$sample_collected <- as.integer(dat$sample_collected) # changing object class from Date to Integer
head(dat) 
dat <- dplyr::select(dat, hh_LAT, hh_LNG, sample_collected)

write.csv(dat, file = "ST_Lineage_0_Data.csv")
