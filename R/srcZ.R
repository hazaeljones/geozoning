# required libraries
#to install:
#install.packages(c("sp","gstat","RandomFields","maptools","rgeos","colorspace","xtable","deldir","mmand","permute","plotrix","fields","optimx"))

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
#source("genMapRNoK.R")
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
# compute Moran and Geary indices
source("funcCalIndices.R")

#
# utilities
source("util.R")
source("utilZ.R")
source("visudispZ.R")
source("gridXY.R")
#
# generate zoning
source("funcZoning.R")
source("initialZoning.R")
#source("saveZoningFromSimu.R")
source("findN.R")
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
#source("perf.R")
source("TESTFUNCTIONS/studyCriteria.R")
source("TESTFUNCTIONS/studyCriteria2.R")
source("TESTFUNCTIONS/studyCriteria3.R")
source("TESTFUNCTIONS/studyCriteria4.R")
source("TESTFUNCTIONS/studyCriteria5.R")
source("figCritN.R")
source("CASESTUDYFUNCTIONS/selMaps.R")
#source("optCost.R")
#source("optCostL.R")
#source("optCrit.R")
#source("optRank.R")
## Phi's functions
#source("calCost.R")
#source("correctBoundaryMap.R")
#source("smoothingMap.R")
#source("smoothingZone.R")
#source("plotM.R")
# extra functions for vignettes
source("EXTRAS/calIndices.R")
source("EXTRAS/listSeeds.R")
source("EXTRAS/plotmdist.R")
source("EXTRAS/plotmat.R")
source("EXTRAS/getBestMLoop.R")