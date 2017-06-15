# required libraries
#to install:
#install.packages(c("sp","gstat","RandomFields","maptools","rgeos","colorspace",
#                  "xtable","deldir","mmand","permute","plotrix","fields","optimx"))

#to install "rgeos" : sudo apt-get install libgeos-dev
#install.packages(c("rgdal","raster","ggplots2","profvis"))

#install.packages("knitr")
#install.packages("tydir")
#install.packages("dplyr")
#to install rgdal : sudo apt-get install libgdal-dev
#install.packages("rgdal")

library(sp)
library(gstat)
library(RandomFields)
library(maptools)
library(rgeos)
library(colorspace)
library(xtable)
library(deldir)
library(mmand)
library(permute)
library(plotrix)
library(fields)
library(optimx)
require("rgdal") # libgdal-dev pour debian
library(raster)
library(ggplot2)
library(profvis)

#optional - call C function
#dyn.load("./librairieC/distance/calDistance.so")

# global parameters
source("initParam.R")

# generate map or read data
source("randKmap.R")
source("randKmapGrid.R")
source("calRMmodel.R")
source("genMap.R")
source("genMapR.R")
source("genMapRNoK.R")
source("funcRandKmap.R")
source("plotMap.R")
source("readS.R")

#compute neighborhood
source("calNei.R")
source("funcCalNei.R")
#compute criteria and cost
source("calCrit.R")
source("funcCalCrit.R")
#compute distances
source("calDistance.R")
source("funcCalDistance.R")
# compute Moran and Geary indices
source("calIndices.R")
source("funcCalIndices.R")

#
# utilities
source("util.R")
source("utilZ.R")
source("visudispZ.R")
#
# generate zoning
source("funcZoning.R")
source("initialZoning.R")
source("saveZoningFromSimu.R")
source("findN.R")
source("getNs.R")
#source("contourManu.R")
#
# correct zoning
source("correctionTree.R")
source("detZoneCE.R")
source("saveZK.R")
source("updateZK.R")
source("lastPass.R")
source("sortCrit.R")
source("optiRG.R")
source("optiGrow.R")
source("zoneQ.R")
source("modNIso.R")
source("calFrame.R")
source("funcCleaning.R")


#
# test zoning
source("TESTFUNCTIONS/testZ.R")
source("TESTFUNCTIONS/testGrow.R")
source("TESTFUNCTIONS/testFusion.R")
# optim
source("loopQ1.R")
source("loopQ2.R")
source("loopQ3.R")
source("loopQ4.R")
source("loopQ5.R")
# eval perf
source("perf.R")
source("studyCriteria.R")
source("studyCriteria2.R")
source("studyCriteria3.R")
source("studyCriteria4.R")
source("studyCriteria5.R")
source("figCritN.R")
source("selMaps.R")
source("optCost.R")
source("optCostL.R")
source("optCrit.R")
source("optRank.R")

# smoothing
source("smoothingZone.R")
source("touch.border.R")
source("zone.extended.R")
source("correctBoundaryMap.R")
source("smoothingMap.R")
source("cal.max.width.R")

