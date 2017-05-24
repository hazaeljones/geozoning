rm(list=ls())

source("srcZ.R") # source libraries and functions - params are in initParam.R

##########################################################################
# simulate map
map=genMap(DataObj=NULL,seed=30,disp=FALSE)
# genMap calls geneData and other functions from func-mapAlea.R
##########################################################################
# generate zoning
ZK=initialZoning(qProb=c(0.1,0.4,0.7), map,disp=FALSE)
##########################################################################
# initialZoning calls zoneGeneration, calNei, calDistance, calCrit
# zoneGeneration calls contourAuto
# contourAuto calls contourLines (native R function) and extensionLigne
# get results
##########################################################################
K=ZK$resZ
Z=K$zonePolygone
qProb=K$qProb
valRef=quantile(map$krigGrid,na.rm=TRUE,prob=qProb) #data values for quantiles
#
# add contours corresponding to a given quantile value - for instance +0.1
pr=qProb+0.1
pr=pmax(pr,0.001)
pr=pmin(pr,0.999)
val = quantile(map$krigGrid,na.rm=TRUE,prob=pr)
plotZ(Z)
##########################################################################
# addContour calls contourAuto
# contourAuto calls extensionLigne
#
##########################################################################
addContour(map,val[1],col="blue")
addContour(map,val[2],col="red")
addContour(map,val[3],col="green")
##########################################################################
#
pr=qProb-0.1
pr=pmax(pr,0.001)
pr=pmin(pr,0.999)
val = quantile(map$krigGrid,na.rm=TRUE,prob=pr)
plotZ(Z)
addContour(map,val[1],col="blue")
addContour(map,val[2],col="red")
addContour(map,val[3],col="green")

# graphics.off()
