smoothingMap = function(Z, width = 0.01)
# function that smooth a map

# INPUT :
    # Z : zone corrected by function correctBoundaryMap
    # width : smoothing parameter

# OUTPUT :
    # Z : zone modified
{
  Z2 = Z # clone of Z (Z2 will contain smoothed united zones)
  Z3 = Z # clone of Z (Z3 will contain unsmoothed united zones)
  nbZ = length(Z)

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
  # arbre contain lists of neighbours for each zone
  arbre = list(mget(paste("zone",1:nbZ,sep="")))

  # tabNei = Number of Neighbours (== 0 means zone is smoothed, ==-1 means zone is smoothed and united with a neighbour of next level)
  tabNei = rep(0,nbZ)
  for (i in 1:nbZ){
    tabNei[i] = length(arbre[[1]][[i]])
  }

  # nb of zones which have at least 2 unsmoothed neighbours
  nbZ.special = 0
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

  iter = 0
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

  return(Z)

}
