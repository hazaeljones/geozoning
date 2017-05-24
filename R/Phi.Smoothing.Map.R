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
label.central = which.min(as.vector(table(lab)))



# tentative du lissage de la carte
num.zone = function(listZ,tab,i){
  # listZ : liste des zones voisines d'une zone dans l'arbre
  # tab : table de 2 colonnes: 1 contenant les labels , 1 contenant les numéros des zones dans l'ordre
  # i : label des zones voisinnes qu'on veut traiter
  newtab = tab[listZ,]
  newtab = newtab[which(newtab[,1]==i),]
  return(newtab[,2])
}
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
  border = g3=readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
  # OPERATION LISSAGE DES ZONES DE label.central
  # on ne lisse pas la dernière zone de label.central
  zone.central = which(lab == label.central)
  nb.zone.central = length(zone.central)

  # cas 1.1
  if(nb.zone.central>1){ # si on a plusieurs zones centrales
    for(j in 1:nb.zone.central-1){
      newZ = Z[[zone.central[j]]]
      for(k in num.zone(arbre[[1]] [[zone.central[j]]], tab, 2) ){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
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
    for(j in num.zone(arbre[[1]] [[zone.central[nb.zone.central]]], tab, 2) ){
      zone.Union = gUnion(zone.Union, get(paste("Z",k,sep="")))
    }

    assign(paste("Z",zone.central[nb.zone.central],sep=""), gDifference(border, zone.Union))

    # tout d'abord, on soutrait les zones centrales
    for(j in 1:nb.zone.central-1){
      for (k in arbre[[1]] [[zone.central[j]]]){
        assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
      }
    }
    # puis, on soustrait le reste : si on est à l'extrémité de l'arbre , on ne le fait pas
    for(i in 2:nbLab-1){
      for (j in which(lab==i)){
        index = num.zone(arbre[[1]] [[j]], tab, i+1)
        for(k in index){
          assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
        }
      }
    }
  }

  # cas 1.2
  if(nb.zone.central==1){ # si on a une seule zone centrale
    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale, puis on peut faire l'union sur elle-même
    # union avec les zones liées à la zone centrale
    for(j in which(lab[ arbre[[1]] [[zone.central[nb.zone.central]]] ]==2) ){
      zone.Union = gUnion(zone.Union, get(paste("Z",k,sep="")))
    }
    assign(paste("Z",zone.central[nb.zone.central],sep=""), gDifference(border, zone.Union))
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
        for(k in which(lab[ arbre[[1]] [[j]] ]==i-1)){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
          # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
          newZ = gUnion(newZ, get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }
  border = g3=readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
  # OPERATION LISSAGE DES ZONES DE label.central
  # on ne lisse pas la dernière zone de label.central
  zone.central = which(lab == label.central)
  nb.zone.central = length(zone.central)

  # cas 2.1
  if(nb.zone.central>1){ # si on a plusieurs zones centrales
    for(j in 1:nb.zone.central-1){
      newZ = Z[[zone.central[j]]]
      for(k in which(lab[ arbre[[1]] [[zone.central[j]]] ]==nbLab-1)){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
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
    for(j in which(lab[ arbre[[1]] [[zone.central[nb.zone.central]]] ]==nbLab-1) ){
      zone.Union = gUnion(zone.Union, get(paste("Z",k,sep="")))
    }

    assign(paste("Z",zone.central[nb.zone.central],sep=""), gDifference(border, zone.Union))

    # tout d'abord, on soutrait les zones centrales
    for(j in 1:nb.zone.central-1){
      for (k in arbre[[1]] [[zone.central[j]]]){
        assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
      }
    }
    # puis, on soustrait le reste : si on est à l'extrémité de l'arbre , on ne le fait pas
    for(i in 2:nbLab-1){
      for (j in which(lab==i)){
        for(k in which(lab[ arbre[[1]] [[j]] ]==i-1)){
          assign(paste("Z",j,sep=""), gDifference(get(paste("Z",j,sep="")), get(paste("Z",k,sep=""))) )
        }
      }
    }
  }

  # cas 2.2
  if(nb.zone.central==1){ # si on a une seule zone centrale
    # chercher la dernière zone centrale parmi les zones centrales
    zone.Union = Z[[zone.central[1]]] # initiliser l'union par la première zone.centrale, puis on peut faire l'union sur elle-même
    # union avec les zones liées à la zone centrale
    for(j in which(lab[ arbre[[1]] [[zone.central[nb.zone.central]]] ]==2) ){
      zone.Union = gUnion(zone.Union, get(paste("Z",k,sep="")))
    }
    assign(paste("Z",zone.central[nb.zone.central],sep=""), gDifference(border, zone.Union))
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
### CAS 3 #######################################################################################################################
#################################################################################################################################

width = 0.02

if (label.central!=1 & label.central!=nbLab){
  for(i in nbLab:label.central+1){
    for (j in which(lab==i)){
      if (i == nbLab){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        for(k in which(lab[ arbre[[1]] [[i+1]] ]==i+1)){ # on cherche les zones qui sont déjà lissées et sont liées directement à newZ
          # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
          newZ = gUnion(newZ, get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }

  for(i in 1:label.central-1){
    for (j in which(lab==i)){
      if (i == 1){
        assign(paste("Z",j,sep = ""), smoothingZone(Z[[j]],width = width))
      }
      else{
        newZ = Z[[j]]
        for(k in which(lab[arbre[[1]] [[i-1]]]==i-1)){# on cherche les zones qui sont déjà lissées et sont liées directement à newZ
          # puis on fusionne ces zones avec newZ pour ensuite effectue un lissage sur l'ensemble
          newZ = gUnion(newZ,get(paste("Z",k,sep="")))
        }
        assign(paste("Z",j,sep = ""), smoothingZone(newZ,width = width))
      }
    }
  }
}

border =readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
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




## #################################################################


# TODO :
# remplir la détermination des zones par les opérations de différence: gDifference
# label central est utilisé pour 2 zones, demander à Patrice pour ce cas-là



