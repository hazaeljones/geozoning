##################################################################
#' initialZoning
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param map xxxx
#' @param pErr xxxx
#' @param simplitol xxxx
#' @param optiCrit xxxx
#' @param disp xxxx
#' @param GridData xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
initialZoning=function(qProb, map, pErr=0.9,simplitol=1e-3,optiCrit=2,disp=0,GridData=F)
##################################################################
  {
    #arguments
    #qProb=vecteur de quantiles
    #map=donnees brutes et krigees
    #choix=critere dans initParam
    #calcule zonage et criteres pour 1 vecteur de quantiles
    # attention zonage pas forcement admissible
    #chercher les valeurs des quantiles
    #simple generation de zones correspondant aux contours des quantiles
    # genere zonage a partir des donnees de map$krigGrid
    #
  Z=zoneGeneration(map,qProb,GridData) #in funcZoning
  # create comments (for holes)
  Z = crComment(Z)
  # add IDs to identify zones in zoning
  Z = setIds(Z)
  # compute neighborhood
  K=calNei(Z,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol)
  # assign zone labels
  K= labZone(K,qProb,map$krigGrid)
 # per label computations - compute overall cost per label
  cL=costLab(K,map)
#
  Z=K$zonePolygone
  # reset ids
  Z = setIds(Z)
  K$zonePolygone=Z
#
#
    tabVal=map$krigData
    listZonePoint=K$listZonePoint
    zoneN=K$zoneN
    surfVoronoi=map$krigSurfVoronoi
    meanTot=K$meanTot
    meanZone=K$meanZone

    # compute distances and criteria
  resDist=calDistance(typedist=1,tabVal,listZonePoint,zoneN,surfVoronoi,meanZone,pErr)
    # compute only on neighbor zones
  resCrit=calCrit(resDist$matDistanceCorr,K$zoneNModif,optiCrit)

   if(disp==2) dispZ(map$step,map$krigGrid,zonePolygone=Z,K=K,boundary=map$boundary,nbLvl=0,mu=2)

  return(list(resCrit=resCrit,resDist=resDist,resZ=K,cL=cL,qProb=qProb))
}
