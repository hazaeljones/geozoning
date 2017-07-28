####################################################################################################
#' updateZK called by lastPass
#'
#' @details Given a map object, a list of zonings, a current and a previous zoning, replaces a zoning in the list of zonings
#' @param map object returned by function genMap or genMapR
#' @param qProb probability vector used to generate quantile values
#' @param le index of current level in list
#' @param kk index of current zoning in level list
#' @param listOfZ list of zoning objects
#' @param crit list of criteria
#' @param cost list of costs
#' @param costL list of per label costs
#' @param nz list of number of zones 
#' @param mdist list of distance matrices
#' @param K1 zoning to be replaced
#' @param Z2 xxxx
#' @param pErr equality tolerance for distance calculations
#' @param optiCrit xxxx
#' @param simplitol xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
updateZK=function(map,qProb,le,kk,listOfZ, crit,cost,costL,nz,mdist,K1,Z2,pErr,optiCrit,simplitol)
####################################################################################################
{
# replace zoning Z1,K1 with Z2
Z1=K1$zonePolygone
# clean current zoning Z2
# yields Z2 - 0 pt or 1 pt zones may disappear
        K2=calNei(Z2,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol)
  	Z2=K2$zonePolygone
# transfer zone labels in K2$lab from K1$lab
# use zone id correspondence between Z1 and Z2
# do not reassign labels
 	K2=trLabZone(K1,K2,Z1,Z2,map,qProb,disp=0)
	K2$qProb=qProb
# compute overall cost per label
	cL=costLab(K2,map)
#
# replace zoning
        listOfZ[[le]][[kk]] = K2
	crit[[le]][[kk]] = 0
	mdist[[le]][[kk]] = 0
	cost[[le]][[kk]] = 0
	costL[[le]][[kk]] = 0
	nz[[le]][[kk]] = 0
#
        if (length(Z2)>1) # at least 2 zones in Z2
        {
            resDist1=calDistance(typedist=1,map$krigData,K2$listZonePoint,K2$zoneN,map$krigSurfVoronoi,K2$meanZone,pErr)
            resCrit1 = calCrit(resDist1$matDistanceCorr,K2$zoneNModif,optiCrit)
            if (resCrit1 == Inf) resCrit1 = 0
            crit[[le]][[kk]] = resCrit1
	    cost[[le]][[kk]] = resDist1$cost
	    costL[[le]][[kk]] = cL
	    nz[[le]][[kk]] = length(Z2)
	    normMat=normDistMat(resDist1$matDistanceCorr,optiCrit)
	    mdist[[le]][[kk]] = normMat
       	}

	return(list(listOfZ=listOfZ,crit=crit,cost=cost,costL=costL,nz=nz,mdist=mdist,nz=nz))
}
