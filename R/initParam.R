
###################### simulation parameters  ######################
#number of data pts to generate before and after kriging
nPoints=450
nPointsK=2000

# Gaussian model for data simulation and kriging
typeMod="Gau"
#
Vpsill=5
# nugget
Vnugget=0.5
#scale
Vrange=0.2
#mean
Vmean=8
#no anisotropy 
alphavario=NULL
Vanis=c(1)
Vang=c(0)

#execution ou non de la simulation conditionnelle, plante pour la carte simulée 6mai2014
nSimuCond=50

#field boundary (default=unit square for simulation)
boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0))
manualBoundary=FALSE

####Parameters for distance and criterion calculation   ###############################
#optiCrit=4 # min(min(dij^2/dii*djj))
#optiCrit=6 # min(mean(dij^2/dii*djj))
optiCrit=2 # min(min(2*dij^2/(dii^2+djj^2)))
#optiCrit=7 # mean(mean(2*dij^2/(dii^2+djj^2)))

pErr=0.9 # tolerance for distance calculations
####Parameters for small zone correction and no grow
minSize = 0.012 # valid zone surface threshold - initial step
minSizeNG= 1e-3 #threshold for both no grow and zone grow acceptance (checkContour)
# LASTPASS - if TRUE remove zones that are still too small at the last level of the correction tree
LASTPASS=TRUE

# distance to other zone >=distIsoZ -> isolated zone
distIsoZ=0.075
LEQ = 5 # length of quantile sequence to grow isolated zone
MAXP =0.1 # quantile sequence maximum shift

#used for zone growing in zoneQ and optiAgr 
tol=0.02

# tolerance for polygon simplification 
simplitol=1e-3
#
