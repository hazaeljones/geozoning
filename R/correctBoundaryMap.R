#' correctBoundaryMap
#' @details function for post treatment of zoning that fixes the problem linked to
#'  the border between two neighbour zones and between zones and the map boundary
#'
#' @param Zi list of initiales zones
#' @param map object returned by function genMap

#' @return new list of zones with correct boundary ang the parameter "width" used for correction
#' @importFrom rgeos gIsValid
#' @importFrom rgeos gIntersection
#' @export
#' @examples
#' map=geozoning::mapTest
#' criti = correctionTree(qProb = c(0.5), map = map)
#' Z = criti$zk[[1]][[1]]$zonePolygone
#' lab = criti$zk[[1]][[1]]$lab
#' plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
#' class(rgeos::gIntersection(Z[[1]],Z[[2]])) [1]
#' class(rgeos::gIntersection(Z[[1]],Z[[5]])) [1]
#' class(rgeos::gIntersection(Z[[2]],Z[[3]])) [1]
#' class(rgeos::gIntersection(Z[[2]],Z[[4]])) [1]
#' res = correctBoundaryMap(Zi = Z, map = map)
#' Z = res$Z
#' class(rgeos::gIntersection(Z[[1]],Z[[2]])) [1]
#' class(rgeos::gIntersection(Z[[1]],Z[[5]])) [1]
#' class(rgeos::gIntersection(Z[[2]],Z[[3]])) [1]
#' class(rgeos::gIntersection(Z[[2]],Z[[4]])) [1]
#' plotM(map = map, Z = Z, lab = lab, byLab = FALSE)


correctBoundaryMap = function (Zi, map)
{
  Z = Zi

  x = map$boundary$x
  y = map$boundary$y
  z = cbind(x,y)
  p = Polygon(z)
  ps = Polygons(list(p),ID = "p")
  boundary = SpatialPolygons(list(ps))

  # INITIALIZE PARAMETERS
  # tree of relation
  nbZ = length(Z)
  for (i in 1:nbZ){
    nei = c()
    for (j in 1:nbZ){
      if(j!=i){
        if(gDistance(Z[[i]],Z[[j]]) < 0.001){ # if the distance between 2 zones is very small, they are neighbour.
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

  # CORRECTION
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

  return(list(Z = Z, width = width))
}
























