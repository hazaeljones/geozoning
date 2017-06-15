correctBoundaryMap = function (Zi, boundary)

  # description : function that fixes the problem linked to the border between two neighbour zones and between zones and the map boundary
  
  # input :
  # Zi : list of initiales zones
  # boundary : boundary of the map
  
  
  # output :
  # Z : list of new zone
  # width : gBuffer's parameter chosen to correct the map

  
{
  Z = Zi
  
  ### INITIALIZE PARAMETERS ####################################################################################
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
  # arbre contains lists of neighbours for each zone
  arbre = list(mget(paste("zone",1:nbZ,sep="")))
  
  # matrix of neighbourhood
  nbZ = length(Z)
  zN = matrix(rep(0,nbZ*nbZ),ncol=nbZ)
  for(i in 1:nbZ){
    for(j in 1:nbZ){
      if (gDistance(Z[[i]],Z[[j]])<10^-3){
        zN[i,j] = TRUE
      }else{
        zN[i,j] = FALSE
      }
    }
  }
  
  ### CORRECTION ##############################################################################################
  width = 0.0001
  res = FALSE
  
  while (res == FALSE){
    Z = Zi
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
    }
    width = width+0.0001
  }
  
  return(list(Z, width))
}
























