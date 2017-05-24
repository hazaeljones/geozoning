rm(list=ls())

source("srcZ.R") # source libraries and functions - params are in initParam.R

##########################################################################
# simulate map
seed=3
map=genMap(DataObj=NULL,seed=seed,disp=FALSE)
# genMap calls geneData and other functions from func-carteAlea.R
##########################################################################
# generate zoning
# zoneGeneration calls contourAuto
# contourAuto calls contourLines (native R function) and extensionLigne
# get results
##########################################################################

qProb=c(0.25,0.75)
# GridData=FALSE: case of simulated data
# GridData=TRUE: case of real data which already constitute a regular grid
Z=zoneGeneration(map,qProb=qProb,GridData=FALSE)
##########################################################################
# polygons have no comment - add comment for holes
Z = crComment(Z)
# add ids to zones
Z = setIds(Z)
plotZ(Z)
##########################################################################
# 26 zones - seed=3 - zones 5 to 8 have 1 point
# calNei removes 17 zones with 0 or 1 data pt, if remove=TRUE (default)
# calls zoneAssign and calculListZonePoint
###########################################################################
K=calNei(Z,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol)
##########################################################################
# 9 zones
# reset ids
Z=K$zonePolygone
Z = setIds(Z)
##########################################################################
# label zones
##########################################################################
K= labZone(K,qProb,map$krigGrid)
##########################################################################
Z=K$zonePolygone
plotZ(Z)
