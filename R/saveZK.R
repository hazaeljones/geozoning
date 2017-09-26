#############################################################################
#' saveZK function called by correctionTree
#'
#' @details Given a map object, a list of zonings, a current and a previous zoning, adds the current zoning to the list of zonings if it has at least 2 zones,after recalculating zone neighborhood and transferring zone labels.
#' @param map object returned by function genMap
#' @param K1 previous zoning
#' @param Z2 current zoning geometry (list of SpatialPolygons)
#' @param qProb probability vector used to generate quantile values
#' @param listOfZ list of zoning objects
#' @param indCur index of new list element
#' @param crit list of criteria
#' @param cost list of costs
#' @param costL list of per label costs 
#' @param nz list of number of zones 
#' @param mdist list of distance matrices
#' @param pErr equality tolerance for distance calculations
#' @param optiCrit criterion choice
#' @param simplitol tolerance for spatial polygons geometry simplification
#'
#' @return a  list with components
#'\describe{
#' \item{listOfZ}{updated list of zoning objects, first element corresponds to initial zoning, each other element is a list with each (last if ALL=FALSE) level zoning objects}
#' \item{mdist}{list of initial distance matrix and all (last if ALL=FALSE) level distance matrices}
#' \item{crit}{list of initial criterion and all (last if ALL=FALSE) level criteria }
#' \item{cost}{list of initial cost and all (last if ALL=FALSE) level costs }
#' \item{costL}{list of initial cost per label and all (last if ALL=FALSE) level costs per label}
#' \item{nz}{list of initial number of zones and all (last if ALL=FALSE) level number of zones}
#' }
#'
#' @keywords internal
#'
#' @examples
#' data(mapTest)
# run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7
#' criti=correctionTree(c(0.4,0.7),mapTest,LASTPASS=FALSE,SAVE=TRUE) 
#' K1=criti$zk[[1]][[1]]#initial zoning
#' Z1=K1$zonePolygone 
#' printZsurf(Z1) # 8 zones with 2 small zones (7 and 8)
#' Z2 = remove1FromZ(Z1,7,K1$zoneN)
#' printZsurf(Z2) #7 zones
#' indCur=2
#' newRes=saveZK(mapTest,K1,Z2,c(0.4,0.7),criti$zk,indCur,
#'         criti$criterion,criti$cost,criti$costL,criti$nz,criti$mdist)
#' newZ=newRes$listOfZ[[2]][[1]]$zonePolygone
#' printZsurf(newZ) #6 zones
#' # not run
saveZK=function(map,K1,Z2,qProb,listOfZ, indCur,crit,cost,costL,nz,mdist,pErr=0.9,optiCrit=2,simplitol=1e-3)
######################################################
{
  # previous zoning Z1
  	Z1=K1$zonePolygone
  # clean current zoning Z2
  # yields Z2 - 0 pt or 1 pt zones may disappear
    	K2=calNei(Z2,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol)
  	Z2=K2$zonePolygone
  # transfer zone labels in K2$lab from K1$lab
  # use zone id correspondence betwwen Z1 and Z2
  # do not reassign labels

	K2=trLabZone(K1,K2,map,qProb,disp=0)
	K2$qProb=qProb
	# per label computations - compute overall cost per label
	cL=costLab(K2,map)

  if (length(Z2)>1) # at least 2 zones in Z2
    {
      resDist1=calDistance(typedist=1,map$krigData,K2$listZonePoint,K2$zoneN,map$krigSurfVoronoi,K2$meanZone,pErr)
      resCrit1 = calCrit(resDist1$matDistanceCorr,K2$zoneNModif,optiCrit)
      if (resCrit1 == Inf) resCrit1 = 0
   
      listOfZ[[indCur]]=append(listOfZ[[indCur]],list(K2))
      crit[[indCur]]=append(crit[[indCur]],resCrit1)
      normMat=normDistMat(resDist1$matDistanceCorr,optiCrit)
      mdist[[indCur]]=append(mdist[[indCur]],list(normMat))
      cost[[indCur]]=append(cost[[indCur]],resDist1$cost)
      costL[[indCur]]=append(costL[[indCur]],cL)
      nz[[indCur]]=append(nz[[indCur]],length(Z2))
      }
	# otherwise do not record K2 in list
	return(list(listOfZ=listOfZ,crit=crit,mdist=mdist,cost=cost,costL=costL,nz=nz))
	
	}

