###########################################################################################
### This Script is about testing the method for smoothing all the map #####################
###########################################################################################



rm(list=ls())
source("srcZ.R") # source libraries and functions - params are in initParam.R


seed= 1
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
plotMap(map)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone

plotZ(Z)

# zone 6 is included in zone 2
nbZ = length(Z)
for (i in 1:nbZ){
  zone = paste
}



lab = ZK$resZ$lab
lab

nbLab = length(levels(as.factor(lab)))

# on cherche les voisins des chaques zones
diagram = list()
for (i in 1:nbZ){
  nei = c()
  for (j in 1:nbZ){
    if(j!=i){
      if(is.null(gIntersection(Z[[i]],Z[[j]])) == FALSE){
        nei = c(nei, j)
      }
    }
  }
  print(nei)
  neigh = paste("neighbours",i,sep="")
  assign(neigh,nei)
}


plot(Z[[1]])
plot(Z[[2]],add =TRUE)
plot(gIntersection(Z[[1]],Z[[2]]), col ="yellow",add=TRUE)














































