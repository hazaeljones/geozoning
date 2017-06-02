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
# seed = 14 intéressant

seed = 15
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
width = 0.01

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
zN = ZK$resZ$zoneN # matrice des voisins

Z = correctBoundaryMap(Z,zN)

nbZ = length(Z)
for (i in 1:(nbZ-1)){
  for (j in (i+1):nbZ){
    if (zN[i,j]==TRUE){
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
seed = 15
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
zN = ZK$resZ$zoneN # matrice des voisins

Z = correctBoundaryMap(Z,zN)

plotZ(Z)

width = 0.01

Z = smoothingMap(Z,width = width)

plotZ(Z)

