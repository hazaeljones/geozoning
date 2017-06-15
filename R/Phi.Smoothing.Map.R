###########################################################################################
### This Script is about testing the method for smoothing all the map #####################
# based on Patrice idea :)
###########################################################################################



rm(list=ls())
source("srcZ.R") # source libraries and functions - params are in initParam.R


seed= 2
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
plotMap(map)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone

plotZ(Z)

nbZ = length(Z) # nombre de zone

lab = ZK$resZ$lab # label des zones
lab

tab = cbind(lab,1:nbZ)  # concatenation de numero des zones et leur label

nbLab = length(levels(as.factor(lab))) # nb label
nbLab
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
arbre
# on cherche le label est le moins utilisé
label.central = which.min(as.vector(table(lab)))
label.central


# tentative du lissage de la carte
#################################################################################################################################
### FUNCTION : NUM.ZONE #########################################################################################################
#################################################################################################################################
num.zone = function(listZ,tab,i){
  # description : liste des zones de label i liées à la zone Z
  # listZ : liste des zones voisines d'une zone Z dans l'arbre
  # tab : table de 2 colonnes: 1 contenant les labels , 1 contenant les numéros des zones dans l'ordre
  # i : label des zones voisinnes qu'on veut traiter
  newtab = tab[listZ,]
  newtab = newtab[which(newtab[,1]==i),]
  return(newtab[,2])
}
# exemple : liste des zones de label 1 liées à la zone 2
num.zone(arbre[[1]][[2]],tab,1)


#################################################################################################################################
### CAS 1 #######################################################################################################################
#################################################################################################################################
if (label.central==1){
  # OPERATION LISSAGE-FUSION-LISSAGE
  for(i in nbLab:2){
    for (j in which(lab==i)){
      if (i == nbLab){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        index = num.zone(arbre[[1]] [[j]], tab, i+1)
        for(k in index){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
          # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
          newZ = gUnion(newZ, get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }
  border = readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
  # OPERATION LISSAGE DES ZONES DE label.central
  # on ne lisse pas la dernière zone de label.central
  zone.central = which(lab == label.central)
  nb.zone.central = length(zone.central)

  # cas 1.1 :# si on a plusieurs zones centrales
  if(nb.zone.central>1){
    for(j in 1:nb.zone.central-1){
      newZ = Z[[zone.central[j]]]
      index = num.zone(arbre[[1]] [[zone.central[j]]], tab, 2)
      for(k in index){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
        # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
        newZ = gUnion(newZ, get(paste("Z",k,sep="")))
      }
      assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
    }
    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale, puis on peut faire l'union sur elle-même
    # union avec les zones centrales sauf la dernière
    for(j in 1:nb.zone.central-1){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }
    # union avec les zones liées à la dernière zone centrale
    index = num.zone(arbre[[1]] [[zone.central[nb.zone.central]]], tab, 2)
    for(j in index){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }

    assign(paste("Z",zone.central[nb.zone.central],sep=""), gDifference(border, zone.Union))

    # tout d'abord, on soustrait les autres zones centrales
    for(j in 1:nb.zone.central-1){
      for (k in arbre[[1]] [[zone.central[j]]]){
        assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
      }
    }
    # puis, on soustrait le reste : si on est à l'extrémité de l'arbre , on ne le fait pas, car c'est déjà fait
    for(i in 2:nbLab-1){
      for (j in which(lab==i)){
        index = num.zone(arbre[[1]] [[j]], tab, i+1)
        for(k in index){
          assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
        }
      }
    }
  }

  # cas 1.2 : si on a une seule zone centrale
  if(nb.zone.central==1){
    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale, puis on peut faire l'union sur elle-même
    # union avec les zones liées à la zone centrale
    for(j in arbre[[1]][[zone.central[1]]] ){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }
    assign(paste("Z",zone.central[1],sep=""), gDifference(border, zone.Union))
    # puis, on soustrait le reste : si on est à l'extrémité de l'arbre , on ne le fait pas
    for(i in 2:nbLab-1){
      for (j in which(lab==i)){
        for(k in which(lab[ arbre[[1]] [[j]] ]==i+1)){
          assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
        }
      }
    }
  }

}

#################################################################################################################################
### CAS 2 #######################################################################################################################
#################################################################################################################################
if (label.central==nbLab){
  # OPERATION LISSAGE-FUSION-LISSAGE
  for(i in 1:nbLab-1){
    for (j in which(lab==i)){
      if (i == 1){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        index = num.zone(arbre[[1]] [[j]], tab, i-1)
        for(k in index){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
          # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
          newZ = gUnion(newZ, get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }
  border =readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
  # OPERATION LISSAGE DES ZONES DE label.central
  # on ne lisse pas la dernière zone de label.central
  zone.central = which(lab == label.central)
  nb.zone.central = length(zone.central)

  # cas 1.1 :# si on a plusieurs zones centrales
  if(nb.zone.central>1){
    for(j in 1:nb.zone.central-1){
      newZ = Z[[zone.central[j]]]
      index = num.zone(arbre[[1]] [[zone.central[j]]], tab, nbLab-1)
      for(k in index){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
        # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
        newZ = gUnion(newZ, get(paste("Z",k,sep="")))
      }
      assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
    }
    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale, puis on peut faire l'union sur elle-même
    # union avec les zones centrales sauf la dernière
    for(j in 1:nb.zone.central-1){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }
    # union avec les zones liées à la dernière zone centrale
    index = num.zone(arbre[[1]] [[zone.central[nb.zone.central]]], tab, nbLab-1)
    for(j in index){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }

    assign(paste("Z",zone.central[nb.zone.central],sep=""), gDifference(border, zone.Union))

    # tout d'abord, on soustrait les autres zones centrales
    for(j in 1:nb.zone.central-1){
      for (k in arbre[[1]] [[zone.central[j]]]){
        assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
      }
    }
    # puis, on soustrait le reste : si on est à l'extrémité de l'arbre , on ne le fait pas, car c'est déjà fait
    for(i in nbLab-1:2){
      for (j in which(lab==i)){
        index = num.zone(arbre[[1]] [[j]], tab, i-1)
        for(k in index){
          assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
        }
      }
    }
  }

  # cas 1.2 : si on a une seule zone centrale
  if(nb.zone.central==1){ # si on a une seule zone centrale
    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale, puis on peut faire l'union sur elle-même
    # union avec les zones liées à la zone centrale
    for(j in arbre[[1]][[zone.central[1]]] ){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }
    assign(paste("Z",zone.central[1],sep=""), gDifference(border, zone.Union))
    # puis, on soustrait le reste : si on est à l'extrémité de l'arbre , on ne le fait pas
    for(i in nbLab-1:2){
      for (j in which(lab==i)){
        for(k in which(lab[ arbre[[1]] [[j]] ]==i-1)){
          assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
        }
      }
    }
  }

}
#################################################################################################################################
### CAS 3 #######################################################################################################################
#################################################################################################################################

width = 0.02

if (label.central!=1 & label.central!=nbLab){
  # OPERATION LISSAGE-FUSION-LISSAGE : de label 1 à label.central-1
  for(i in 1:(label.central-1)){
    for (j in which(lab==i)){
      if (i == 1){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        index = num.zone(arbre[[1]] [[j]], tab, i-1)
        for(k in index){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
          # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
          newZ = gUnion(newZ, get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }
  # OPERATION LISSAGE-FUSION-LISSAGE: de label nbLab à label.central+1
  for(i in nbLab:(label.central+1)){
    for (j in which(lab==i)){
      if (i == nbLab){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        index = num.zone(arbre[[1]] [[j]], tab, i+1)
        for(k in index){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
          # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
          newZ = gUnion(newZ, get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }
  border =readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
  # OPERATION LISSAGE DES ZONES DE label.central
  # on ne lisse pas la dernière zone de label.central
  zone.central = which(lab == label.central)
  nb.zone.central = length(zone.central)

  # cas : plusieurs zones centrales
  if(nb.zone.central>1){
    for(j in 1:(nb.zone.central-1)){
      newZ = Z[[zone.central[j]]]
      # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
      # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
      for(k in arbre[[1]][[zone.central[j]]]){
        newZ = gUnion(newZ, get(paste("Z",k,sep="")))
      }
      assign(paste("Z",zone.central[j],sep = ""), smoothingZone(newZ,width = width))
    }
    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = get(paste("Z",zone.central[1],sep="")) # initiliser l'union par la première zone.centrale
    # union avec les zones centrales sauf la dernière
    for(j in 1:nb.zone.central-1){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }
    # union avec les zones liées à la dernière zone centrale
    for(j in arbre[[1]][[zone.central[nb.zone.central]]]){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }

    assign(paste("Z",zone.central[nb.zone.central],sep=""), gDifference(border, zone.Union))
    # OPERATION gDifference
    # tout d'abord, on soustrait les autres zones centrales
    for(j in 1:nb.zone.central-1){
      for (k in arbre[[1]] [[zone.central[j]]]){
        assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
      }
    }
    # puis, on soustrait le reste : si on est à l'extrémité de l'arbre , on ne le fait pas, car c'est déjà fait
    for(i in label.central:2){
      if(i!= 2){
        for (j in which(lab==i)){
          index = num.zone(arbre[[1]] [[j]], tab, i-1)
          for(k in index){
            assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
          }
        }
      }
    }
    for(i in label.central:nbLab-1){
      if(i!= nbLab-1){
        for (j in which(lab==i)){
          index = num.zone(arbre[[1]] [[j]], tab, i+1)
          for(k in index){
            assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
          }
        }
      }
    }
  }
  # cas : une seule zone centrale
  if(nb.zone.central==1){

    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale, puis on peut faire l'union sur elle-même
    assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))

    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale
    # union avec les zones liées à la dernière zone centrale
    for(j in arbre[[1]][[zone.central[1]]] ){
      zone.Union = gUnion(zone.Union, get(paste("Z",j,sep="")))
    }
    assign(paste("Z",zone.central[1],sep=""), gDifference(border, zone.Union))

    # OPERATION gDifference
    for(i in label.central:2){
      if(i!= 2){
        for (j in which(lab==i)){
          index = num.zone(arbre[[1]] [[j]], tab, i-1)
          for(k in index){
            assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
          }
        }
      }
    }
    for(i in label.central:nbLab-1){
      if(i!= nbLab-1){
        for (j in which(lab==i)){
          index = num.zone(arbre[[1]] [[j]], tab, i+1)
          for(k in index){
            assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
          }
        }
      }
    }
  }

}


##########################################################################################
### function correction border between 2 neigbours #######################################
##########################################################################################

correction.Nei = fucntion()


##########################################################################################
### fucntion generates tree ##############################################################
##########################################################################################

tree.Generate = function(tree,label){
  nbLab = length(levels(as.factor(label))) # nb label
  nbZ = length(label) # nb Zone

  len = rep(0,nbLab)
  for (i in 1:nbLab){
    len[i] = length(tree[[1]][[i]])
  }
  lenMax = max(len)


  for (i in 1:nbLab){
    n = length(which(label==i))
    for(j in which(label==i)){
      x = c((1:n)/(n+1))
      y = rep(i-0.2,n)
      r = rep(1,n)
      symbols(x,y,circles = r,xlim=c(0,1),ylim=c(-1,nbLab+1), inches = 0.2)
    }
  }
}

text(x, y, labels)
ax<-c(20,20,20,80,80,80)
by<-c(80,20,20,80,20,20)

#taille
s<-c(1,1,2,1,1,2)
#couleur
couleur<-c("black","black","grey","blue","green","grey")

#on trace les cercles sur un graphique
symbols(ax,by,circles=s,inches=1,fg=couleur,xlim=c(0,100),ylim=c(0,100),lwd=2.5)



##################################################################################################
# REFAIRE ########################################################################################
##################################################################################################


################################################################################################
## CORRECTION DES FRONTIERES ###################################################################
################################################################################################


seed= 81
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone

plotZ(Z)

nbZ = length(Z) # nombre de zone

lab = ZK$resZ$lab # label des zones
lab

tab = cbind(lab,1:nbZ)  # concatenation de numero des zones et leur label

nbLab = length(levels(as.factor(lab))) # nb label
nbLab
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
arbre
# on cherche le label est le moins utilisé
label.central = which.min(as.vector(table(lab)))
label.central


# Correction de la frontière de 2 zones voisinnes

Z7 = Z[[7]]
Z8 = Z[[8]]

Z8 = correctBoundary(Z8,Z7)
plot(gUnion(Z7,Z8))

# Correction des frontières de la carte

zN = ZK$resZ$zoneN # matrice des voisins

for (i in 1:(nbZ-1)){
  for (j in (i+1):nbZ){
    if (zN[i,j]==TRUE){
      print(paste(i,j))
      Z[[i]] = correctBoundary(Z[[i]],Z[[j]])
    }
  }
}

plotZ(Z)

plot(gUnion(Z[[1]],Z[[2]]))
plot(gUnion(Z[[2]],Z[[3]]))
plot(gUnion(Z[[2]],Z[[4]]))
plot(gUnion(Z[[2]],Z[[5]]))
plot(gUnion(Z[[2]],Z[[6]]))
plot(gUnion(Z[[2]],Z[[7]]))
plot(gUnion(Z[[2]],Z[[11]]))
plot(gUnion(Z[[2]],Z[[12]]))
plot(gUnion(Z[[7]],Z[[8]]))
plot(gUnion(Z[[7]],Z[[9]]))
plot(gUnion(Z[[9]],Z[[10]]))

Z = correctBoundaryMap(Z,zN)


# OK


################################################################################################
## CORRECTION DE LA CARTE ######################################################################
################################################################################################
seed = 14

map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
width = 0.02

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
zN = ZK$resZ$zoneN # matrice des voisins

Z = correctBoundaryMap(Z,zN)
z = Z[[1]]
for ( i in 2:nbZ){
  z = gUnion(z,Z[[i]])
}
plot(z)

nbZ = length(Z)
for (i in 1:(nbZ-1)){
  for (j in (i+1):nbZ){
    if (zN[i,j]==TRUE){
      print(c(i,j))
      print(gIntersection(Z[[i]],Z[[j]]))
    }
  }
}

plotZ(Z)

Z2 = Z # clone of Z (Z2 will contain smoothed united zones)
Z3 = Z # clone of Z (Z3 will contain unsmoothed united zones)

nbZ = length(Z)
nbZ

# tree of relation
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
arbre

tabNei = rep(0,nbZ)
# tabNei = Number of Neighbours (== 0 means zone is smoothed)
for (i in 1:nbZ){
  tabNei[i] = length(arbre[[1]][[i]])
}
tabNei


nbZ.special = 0 # nb of zones which have at least 2 unsmoothed neighbours
for (i in which(tabNei>=2)){
  count.Nei.unsmoothed = 0
  for (j in arbre[[1]][[i]]){
    if (tabNei[j] > 0){
      count.Nei.unsmoothed = count.Nei.unsmoothed + 1
    }
  }
  if (count.Nei.unsmoothed >= 2){
    nbZ.special = nbZ.special + 1
  }
}
nbZ.special

iter = 0

# C'EST PARTI #################################################################################################################################"
border = readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
if(nbZ.special<=2){

  if (nbZ.special == 1){
    for (i in which(tabNei==1)){
      Z[[i]] = smoothingZone(Z[[i]],width= width)
      tabNei[i] = 0
    }
    # for the last zone : one need to remove, from the map, the neighbour zones which have been smoothed
    index = which(tabNei > 1)
    Z[[index]] = border
    for (j in arbre[[1]][[index]]){
      Z[[index]] = gDifference(Z[[index]],Z[[j]])
    }
  }

  if (nbZ.special == 2){
    for (i in which(tabNei==1)){
      Z[[i]] = smoothingZone(Z[[i]],width= width)
      Z2[[i]] = Z[[i]]
      tabNei[i] = 0
    }
    for (i in which(tabNei>1)){
      for (j in arbre[[1]][[i]]){
        if (tabNei[j]==0){
          Z3[[i]] = gUnion(Z3[[i]],Z3[[j]])
          tabNei[i] = tabNei[i]-1
        }
      }
    }
    # for the 2 last zones : one need to smooth 1 zone, then remove from the map all neighbours of the other zone
    i = which(tabNei == 1)[1]
    Z[[i]] = smoothingZone(Z3[[i]],width= width)
    Z2[[i]] = Z[[i]]
    # one need to remove the neighbour zones which have been smoothed
    neighbourSmoothed = c()
    for(k in arbre[[1]][[i]]){
      if (tabNei[k] == 0) {
        neighbourSmoothed = c(neighbourSmoothed, k)
      }
    }
    for (j in neighbourSmoothed){
      Z[[i]] = gDifference(Z[[i]],Z2[[j]])
      tabNei[i] = 0
    }
    # last zone
    index = which(tabNei == 1)
    Z[[index]] = border
    for (j in arbre[[1]][[index]]){
      Z[[index]] = gDifference(Z[[index]],Z2[[j]])
    }
  }

}else{

  while(nbZ.special > 2){

    # SMOOTHING - DIFFERENCE
    if (iter == 0){
      for (i in which(tabNei==1)){ #  smooth zones which have 1 neighbour
        Z[[i]] = smoothingZone(Z[[i]],width = width)
        Z2[[i]] = Z[[i]]
        tabNei[i] = 0
      }
    }else{
      for (i in which(tabNei==1)){ #  smooth zones which have 1 neighbour
        Z[[i]] = smoothingZone(Z3[[i]],width = width)
        Z2[[i]] = Z[[i]]
        # one need to remove the neighbour zones which have been smoothed
        neighbourSmoothed = c()
        for(k in arbre[[1]][[i]]){
          if (tabNei[k] <= 0) {
            neighbourSmoothed = c(neighbourSmoothed, k)
          }
        }
        for (j in neighbourSmoothed){
          Z[[i]] = gDifference(Z[[i]],Z2[[j]])
        }
        tabNei[i] = 0
      }
    }

    # UNION
    for (i in which(tabNei>1)){
      for (j in arbre[[1]][[i]]){
        if (tabNei[j]==0){
          Z3[[i]] = gUnion(Z3[[i]],Z3[[j]])
          tabNei[i] = tabNei[i]-1
          tabNei[j] = -1
        }
      }
    }

    # UPDATE PARAM OF LOOP
    nbZ.special = 0 # nb of zones which have at least 2 neighbours unsmoothed
    for (i in which(tabNei>=2)){
      count.Nei.unsmoothed = 0
      for (j in arbre[[1]][[i]]){
        if (tabNei[j] > 0){
          count.Nei.unsmoothed = count.Nei.unsmoothed + 1
        }
      }
      if (count.Nei.unsmoothed >= 2){
        nbZ.special = nbZ.special + 1
      }
    }
    nbZ.special
    iter = iter + 1

  }

  if (nbZ.special == 1){
    for (i in which(tabNei==1)){
      Z[[i]] = smoothingZone(Z3[[i]],width= width)
      Z2[[i]] = Z[[i]]
      # one need to remove the neighbour zones which have been smoothed
      neighbourSmoothed = c()
      for(k in arbre[[1]][[i]]){
        if (tabNei[k] == -1) {
          neighbourSmoothed = c(neighbourSmoothed, k)
        }
      }
      for (j in neighbourSmoothed){
        Z[[i]] = gDifference(Z[[i]],Z2[[j]])
      }
      tabNei[i] = 0
    }
    # for the last zone : one need to remove, from the map, the neighbour zones which have been smoothed
    index = which(tabNei > 1)
    Z[[index]] = border
    for (j in arbre[[1]][[index]]){
      Z[[index]] = gDifference(Z[[index]],Z2[[j]])
    }
  }

  if (nbZ.special == 2){
    for (i in which(tabNei==1)){
      Z[[i]] = smoothingZone(Z3[[i]],width= width)
      Z2[[i]] = Z[[i]]
      # one need to remove the neighbour zones which have been smoothed
      neighbourSmoothed = c()
      for(k in arbre[[1]][[i]]){
        if (tabNei[k] == -1) {
          neighbourSmoothed = c(neighbourSmoothed, k)
        }
      }
      for (j in neighbourSmoothed){
        Z[[i]] = gDifference(Z[[i]],Z2[[j]])
        tabNei[i] = 0
      }
    }
    # union of zones smoothed with the neighbours that aren't smoothed
    for (i in which(tabNei>1)){
      for (j in arbre[[1]][[i]]){
        if (tabNei[j]==0){
          Z3[[i]] = gUnion(Z3[[i]],Z3[[j]])
          tabNei[i] = tabNei[i]-1
          tabNei[j] = -1
        }
      }
    }

    # for the 2 last zones : one need to smooth 1 zone, then remove from the map all neighbours of the other zone
    i = which(tabNei == 1)[1]
    Z[[i]] = smoothingZone(Z3[[i]],width= width)
    Z2[[i]] = Z[[i]]
    # one need to remove the neighbour zones which have been smoothed
    neighbourSmoothed = c()
    for(k in arbre[[1]][[i]]){
      if (tabNei[k] == -1) {
        neighbourSmoothed = c(neighbourSmoothed, k)
      }
    }
    for (j in neighbourSmoothed){
      Z[[i]] = gDifference(Z[[i]],Z2[[j]])
      tabNei[i] = 0
    }
    # last zone
    index = which(tabNei == 1)
    Z[[index]] = border
    for (j in arbre[[1]][[index]]){
      Z[[index]] = gDifference(Z[[index]],Z2[[j]])
    }
  }
}

plotZ(Z)


###### TEST FUCNTION smoothingMap.R ########################################################################################################
seed = 5

map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
zN = ZK$resZ$zoneN # matrice des voisins

Z = correctBoundaryMap(Z,zN)

plotZ(Z)

width = 0.019

Z = smoothingMap(Z,width = width)

plotZ(Z)

# OKKKKKKKKKKKKKK


############################################################################################################################################
################# LISSAGE SUR UNE ZONE CORRIGÉE #############################################################################
############################################################################################################################################



qProb=c(0.1,0.2)
criti=correctionTree(qProb,map,pErr=0.9,optiCrit=2,minSize=0.012,minSizeNG=1e-3,distIsoZ=0.075,simplitol=1e-3,
                     LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=0,SAVE=TRUE,ONE=FALSE)

ZK=criti$zk
Z=ZK[[1]][[1]]$zonePolygone # zone

zN = matrix(rep(0,length(Z)^2),ncol = length(Z))
for (i in 1:length(Z)){
  for (j in 1:length(Z)){
    if (gDistance(Z[[i]],Z[[j]])<10^-3){
      zN[i,j] = TRUE
    }else{
      zN[i,j]=FALSE
    }
  }
}

Z = correctBoundaryMap(Z,zN)

plotZ(Z)

width = 0.04

Z = smoothingMap(Z,width = width)

plotZ(Z)




############################################################################################################################################
################# Test TRY CATCH ###########################################################################################################
############################################################################################################################################

readUrl <- function(x) {
  out <- tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      message("This is the 'try' part")

      readLines(con=url, warn=FALSE)
      # The return value of `readLines()` is the actual value
      # that will be returned in case there is no condition
      # (e.g. warning or error).
      # You don't need to state the return value via `return()` as code
      # in the "try" part is not wrapped insided a function (unlike that
      # for the condition handlers for warnings and error below)
    },
    error=function(cond) {
      message(paste("URL does not seem to exist:", url))
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("URL caused a warning:", url))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      return(NULL)
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>'
      message(paste("Processed URL:", url))
      message("Some other message at the end")
    }
  )
  return(out)
}


### Données réelles ###############################################################################################
donnees = load("yieldMap+Z.Rdata")
plotMap(map = map)

plotZ(Z1)

plotZ(Z4)

Z = smoothingMap(Z5,width = 0.01)

cal.max.width.Map(Z5,errMax = 0.001)




### Corriger la fonction zone extended and correctBoundary ########################################################
border2 =readWKT("POLYGON((-1 -1,2 -1, 2 2, -1 2,-1 -1))")
border =readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
p = readWKT("POLYGON((0.5 0,1 0,1 1,0.5 1,0.5 0))")
p1 = gBuffer(p,width = 0.1)
plot(border2)
plot(border,add=TRUE)
plot(p,add=TRUE)
plot(p1,add=TRUE)

x = gIntersection(p1,border)
x = gDifference(x,p)
p = gDifference(p1,x)
plot(p,add=TRUE,col = "yellow")


p = gBuffer(p,width = 0.05)
p = gBuffer(p,width = -0.05)
p = gBuffer(p,width = -0.05)
p = gBuffer(p,width = 0.05)


Z = gUnion(gUnion(Z1[[1]],Z1[[2]]),Z1[[3]])


# Corriger correctBondaryMap

nbZ = length(Z1)

# matrice de zones voisines
zN = matrix(rep(0,n*n),ncol=n)

for(i in 1:nbZ){
  for(j in 1:nbZ){
    if (gDistance(Z1[[i]],Z1[[j]])<10^-3){
      zN[i,j] = TRUE
    }else{
      zN[i,j] = FALSE
    }
  }
}






################## Adapter correctBoundary() aux cartes réelles ########################################################
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)

p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

plot(gDifference(boundary,Z2[[2]]),col="yellow")


plotZ(Z2)

z = Z2[[2]]
z.df = geom(z)
boundaryLine = gBoundary(boundary)


for (i in 1:nrow(z.df)){
  pz = readWKT(paste("POINT(",z.df[i,5],z.df[i,6],")"))
  dMin = gDistance(pz,boundaryLine)
  pointProjection = gNearestPoints(pz,boundaryLine)
  if (dMin<10^-3){
    z.df[indMin,5] = pointProjection@coords[2,1]
    z.df[indMin,6] = pointProjection@coords[2,2]
  }
}

newz = geomToPoly(z.df)

newZ.df = geom(gDifference(boundary,z))








################## Changer Boundary de la carte ########################################################################

# méthode utilisé les points proches de la bordure pour définir un nouvelle frontière##################################
load("yieldMap+Z.Rdata")

x = round(map$boundary$x, digit = 6)
y = round(map$boundary$y,digit = 6)
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

Z = Z5



result = correctBoundaryMap(Z,boundary)


Z = result[[1]]
boundary = result[[2]]

plotZ(Z)
plot(gDifference(boundary,Z[[2]]))

carte = Z[[1]]
for (i in 2:length(Z)){
  carte = gUnion(carte,Z[[i]])
}
plot(carte)
carte = gBuffer(carte,width = 0.001)
carte = gBuffer(carte,width = -0.002)
plot(carte)


newBoundary = carte

# méthode utilisé la projection ###############################"
for (i in nrow(z)){
  point = readWKT(paste("POINT(",z[i,1],z[i,2],")"))
  pointProjection = gNearestPoints(point, carte)
  z[i,1] = pointProjection@coords[2,1]
  z[i,2] = pointProjection@coords[2,2]

}
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))










# Méthode utilisé les buffers extérieurs #####################################################################

load("yieldMap+Z.Rdata")
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))



Zc = Z5
Z = Zc

# tree of relation
nbZ = length(Z)
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
# arbre contain lists of neighbours for each zone
arbre = list(mget(paste("zone",1:nbZ,sep="")))

# matrix of neighbourhood
nbZ = length(Z)
zN = matrix(rep(0,nbZ*nbZ),ncol=nbZ)
for(i in 1:nbZ){
  for(j in 1:nbZ){
    if (gDistance(Z[[i]],Z[[j]])<10^-3){
      zN[i,j] = TRUE
      if(i<j)
        print(c(i,j))
    }else{
      zN[i,j] = FALSE
    }
  }
}

## c'est parti correction boundary ##############################
width = 0.0001
res = FALSE

while (res == FALSE){
  Z = Zc
  for (i in 1:length(Z)){
    Z[[i]] = gBuffer(Z[[i]],width = width)
    Z[[i]] = gIntersection(Z[[i]],boundary)
    for (j in arbre[[1]][[i]]){
      Z[[i]] = tryCatch(
        gDifference(Z[[i]],Z[[j]]),
        error = function(e){
          return(Z[[i]])
        }
      )
    }
  }

  if(length(Z)!= nbZ){
    res = FALSE
  }
  if(length(Z)== nbZ){
    isValid  = rep(0,nbZ)
    for (i in 1:nbZ){
      isValid[i] = gIsValid(Z[[i]])
    }
    if(length(which(isValid==FALSE))>0){
      res = FALSE
    }else{
      ClassIntersection = c()
      for(i in 1:(nbZ-1)){
        for(j in (i+1):nbZ){
          if(zN[i,j]==TRUE){
            c = class(gIntersection(Z[[i]],Z[[j]]))[[1]]
            ClassIntersection = c(ClassIntersection, c)
          }
        }
      }
      if(length(which(ClassIntersection !="SpatialLines"))==0 ){
        res = TRUE
      }
    }
    #print(res)
    print(width)
  }
  width = width+0.0001
}

# test validité 1
carte = Z[[1]]
for (i in 2:length(Z)){
  carte = gUnion(carte,Z[[i]])
}
plot(carte)

# test validité 2
for(i in 1:(nbZ-1)){
  for(j in (i+1):nbZ){
    if(zN[i,j]==TRUE){
      print(paste(i,j))
      print(class(gIntersection(Z[[i]],Z[[j]]))[[1]])
    }
  }
}

# test validité 3
for (i in 1:length(Z)){
  print(paste("zone",i,":", gIsValid(Z[[i]])))
}

plotZ(Z)
plotZ(Zc)

############### TEST correctBoundaryMap.R ####################################################################

# données simulées ###########################################################################################
seed = 14

map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

X = correctBoundaryMap(Z,boundary)
Z = X[[1]]
plotZ(Z)

boundary = Z[[1]]
for (i in 2:length(Z)){
  boundary = gUnion(boundary,Z[[i]])
}

Z = smoothingMap(Z,width = 0.05, boundary = boundary)
plotZ(Z)


# test validité 1
carte = Z[[1]]
for (i in 2:length(Z)){
  carte = gUnion(carte,Z[[i]])
}
plot(carte)

# test validité 2
nbZ = length(Z)
zN = matrix(rep(0,nbZ*nbZ),ncol=nbZ)
for(i in 1:nbZ){
  for(j in 1:nbZ){
    if (gDistance(Z[[i]],Z[[j]])<10^-3){
      zN[i,j] = TRUE
      if(i<j)
        print(c(i,j))
    }else{
      zN[i,j] = FALSE
    }
  }
}
for(i in 1:(nbZ-1)){
  for(j in (i+1):nbZ){
    if(zN[i,j]==TRUE){
      print(paste(i,j))
      print(class(gIntersection(Z[[i]],Z[[j]]))[[1]])
    }
  }
}


# test validité 3
for (i in 1:length(Z)){
  print(paste("zone",i,":", gIsValid(Z[[i]])))
}


# données réelles ##################################################################################""

load("yieldMap+Z.Rdata")
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

X = correctBoundaryMap(Z1,boundary)
Z = X[[1]]
#plotZ(Z)

boundary = Z[[1]]
for (i in 2:length(Z)){
  boundary = gUnion(boundary,Z[[i]])
}

# test touch.border.R  et  zone.extended.R  
touch.border(Z[[2]],boundary)
z2.extended = zone.extended(Z[[2]],boundary)


# test smoothingZone.R
z1 = Z[[1]]
z1.smoothed = smoothingZone(z1,width = 0.01,boundary = boundary)

z5 = Z[[5]]
z5.smoothed = smoothingZone(z5,width = 0.02,boundary = boundary)
plot(z5.smoothed)
plot(zone.extended(z5,boundary),add=TRUE)


# test smoothingMap.R

width = 0.05

load("yieldMap+Z.Rdata")
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

X = correctBoundaryMap(Z4,boundary)

Z = X[[1]]
plotZ(Z)
boundary = Z[[1]]
for (i in 2:length(Z)){
  boundary = gUnion(boundary,Z[[i]])
}

Z = smoothingMap(Z,width = width,boundary = boundary)
plotZ(Z) 


boundary = Z[[1]]
for (i in 2:length(Z)){
  boundary = gUnion(boundary,Z[[i]])
}
plot(boundary)
X = correctBoundaryMap(Z,boundary)



# test validité 1
carte = Z[[1]]
for (i in 2:length(Z)){
  carte = gUnion(carte,Z[[i]])
}
plot(carte)

# test validité 2
nbZ = length(Z)
zN = matrix(rep(0,nbZ*nbZ),ncol=nbZ)
for(i in 1:nbZ){
  for(j in 1:nbZ){
    if (gDistance(Z[[i]],Z[[j]])<10^-3){
      zN[i,j] = TRUE
      if(i<j)
        print(c(i,j))
    }else{
      zN[i,j] = FALSE
    }
  }
}
for(i in 1:(nbZ-1)){
  for(j in (i+1):nbZ){
    if(zN[i,j]==TRUE){
      print(paste(i,j))
      print(class(gIntersection(Z[[i]],Z[[j]]))[[1]])
    }
  }
}


# test validité 3
for (i in 1:length(Z)){
  print(paste("zone",i,":", gIsValid(Z[[i]])))
}









