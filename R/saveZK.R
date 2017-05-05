#############################################################################
#' saveZK
#'
#' @details description, a paragraph
#' @param map xxxx
#' @param K1 xxxx
#' @param Z2 xxxx
#' @param qProb xxxx
#' @param listOfZ xxxx
#' @param counter xxxx
#' @param crit xxxx
#' @param cost xxxx
#' @param costL xxxx
#' @param nz xxxx
#' @param mdist xxxx
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
saveZK=function(map,K1,Z2,qProb,listOfZ, counter,crit,cost,costL,nz,mdist,pErr,optiCrit,simplitol)
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

	K2=trLabZone(K1,K2,Z1,Z2,map,qProb,disp=0)
	K2$qProb=qProb
	# per label computations - compute overall cost per label
	cL=costLab(K2,map)

  if (length(Z2)>1) # at least 2 zones in Z2
    {
      resDist1=calDistance(typedist=1,map$krigData,K2$listZonePoint,K2$zoneN,map$krigSurfVoronoi,K2$meanZone,pErr)
      resCrit1 = calCrit(resDist1$matDistanceCorr,K2$zoneNModif,optiCrit)
      if (resCrit1 == Inf) resCrit1 = 0
      listOfZ[[counter+1]]=append(listOfZ[[counter+1]],list(K2))
      crit[[counter+1]]=append(crit[[counter+1]],resCrit1)
	    normMat=normDistMat(resDist1$matDistanceCorr,optiCrit)
	    mdist[[counter+1]]=append(mdist[[counter+1]],list(normMat))
	    cost[[counter+1]]=append(cost[[counter+1]],resDist1$cost)
	    costL[[counter+1]]=append(costL[[counter+1]],cL)
	    nz[[counter+1]]=append(nz[[counter+1]],length(Z2))
  }
	# else do not record K2 in listOfZ
	return(list(listOfZ=listOfZ,crit=crit,mdist=mdist,cost=cost,nz=nz,costL=costL))
}
