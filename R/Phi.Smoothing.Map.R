###########################################################################################
### This Script is about testing the method for smoothing all the map #####################
# based on Patrice idea :)
###########################################################################################



rm(list=ls())
source("srcZ.R") # source libraries and functions - params are in initParam.R


seed= 1
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
plotMap(map)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone

plotZ(Z)

nbZ = length(Z) # nombre de zone

lab = ZK$resZ$lab # label des zones
lab

nbLab = length(levels(as.factor(lab))) # nb label

# on cherche les voisins de chaque zone

for (i in 1:nbZ){
  nei = c()
  for (j in 1:nbZ){
    if(j!=i){
      if(gDistance(Z[[i]],Z[[j]]) < 0.001){ # si la distance entre 2 zones est très petite, elles sont voisinnes
        nei = c(nei, j)
      }
    }
  }
  neigh = paste("zone",i,sep="")
  assign(neigh,nei)
}
# arbre contient la liste des voisins pour chaques zones
arbre = list(mget(paste("zone",1:nbZ,sep="")))

# on cherche le label est le moins utilisé
label.center = which.min(as.vector(table(lab)))



# tentative du lissage de la carte

# case  1 :
if (label.center==1){
  for (i in nbLab:2){
    for (j in which(lab==i)){

    }
  }
}
# case 2 :
if (label.center==nbLab){
  for(i in 1:nbLab-1){
    for (j in which(lab==i)){
      print(arbre[[j]])
    }
  }
}

# case 3 : label le moins utilisé est au centre
width = 0.01

if (label.center!=1 & label.center!=nbLab){
  for(i in nbLab:label.center+1){
    for (j in which(lab==i)){
      if (j == nbLab){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        for(k in which(lab[arbre[[1]] [[j+1]]]==1)){
          newZ = gUnion(newZ, get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }

  for(i in 1:label.center-1){
    for (j in which(lab==i)){
      if (j == 1){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        for(k in which(lab[arbre[[1]] [[j-1]]]==1)){
          newZ = gUnion(newZ,get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }
}
border = g3=readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
Z1 = gDifference(border,gUnion(Z2,Z3,Z4,Z6,Z7,Z8,Z9,Z10,Z11,Z12,Z13))


plotZ(Z)

border = g3=readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
plot(border)
plot(Z1,add=TRUE)
plot(Z2,add=TRUE)
plot(Z3,add=TRUE)
plot(Z4,add=TRUE)
plot(Z5,add=TRUE)
plot(Z6,add=TRUE)
plot(Z7,add=TRUE)
plot(Z8,add=TRUE)
plot(Z9,add=TRUE)
plot(Z10,add=TRUE)
plot(Z11,add=TRUE)
plot(Z12,add=TRUE)
plot(Z13,add=TRUE)

# TODO :
# remplir les cas 1 et 2 (cas particulier de cas 3)
# remplir la détermination des zones par les opérations de différence: gDifference
# label central est utilisé pour 2 zones, demander à Patrice pour ce cas-là



