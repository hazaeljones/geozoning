####################################################################################################
#' updateZK
#'
#' @details description, a paragraph
#' @param map xxxx
#' @param qProb xxxx
#' @param le xxxx
#' @param kk xxxx
#' @param listOfZ xxxx
#' @param crit xxxx
#' @param cost xxxx
#' @param costL xxxx
#' @param nz xxxx
#' @param mdist xxxx
#' @param K1 xxxx
#' @param Z1 xxxx
#' @param Z2 xxxx
#' @param pErr xxxx
#' @param optiCrit xxxx
#' @param simplitol xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
updateZK=function(map,qProb,le,kk,listOfZ, crit,cost,costL,nz,mdist,K1,Z1,Z2,pErr,optiCrit,simplitol)
####################################################################################################
{
# replace zoning Z1,K1 with Z2

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
