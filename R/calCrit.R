######################################################
#' calCrit
#'
#' @details wrapper function that redirects to the proper criterion calculation function
#' according to optiCrit arg value
#' @param matDistanceCorr corrected distance matrix between zones, result of call to \code{\link{calDistance}}
#' @param zoneNModif modified zone neighborhood matrix (FALSE values on diagonal), result of call to \code{\link{calNei}}
#' @param optiCrit criterion to be optimized. Possible values are
#' * 1 for min(mean(dij^2/(dii^2+dij^2)))
#' * 2 for min(2*min(dij/(dii+djj)))
#' * 3 for min(2*min(dij/(dii+djj)))
#' * 4 for min(min(dij^2/sqrt(dii^2*djj^2)))
#' * 5 for min(median(dij^2/sqrt(dii^2*djj^2)))
#' * 7 for mean(2*mean(dij/(dii+djj)))
#'
#' @return the criterion value as a real positive number indicating the zoning quality.
#'
#' @export
#'
#' @examples
#' # compute criterion on test zoning included in package
#' # load test map with simulated data
#' data(mapTest)
#' # load zoning results from test file
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,
#'        K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' crit = calCrit(resD$matDistanceCorr,K$zoneNModif,2)
#' print(crit)
calCrit=function(matDistanceCorr,zoneNModif,optiCrit=2)
######################################################
{
#global param optiCrit
  if (optiCrit== 1)         criti = calCrit1 ( matDistanceCorr , zoneNModif )
  else if (optiCrit == 2)   criti = calCrit2 ( matDistanceCorr , zoneNModif )
  else if (optiCrit == 3)   criti = calCrit3 ( matDistanceCorr , zoneNModif )
  else if (optiCrit == 4)   criti = calCrit4 ( matDistanceCorr , zoneNModif )
  else if (optiCrit == 5)   criti = calCrit5 ( matDistanceCorr , zoneNModif )
  else if (optiCrit == 7)   criti = calCrit7 ( matDistanceCorr , zoneNModif )
  else criti=0

  return (criti)
}
