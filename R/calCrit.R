######################################################
#' calCrit
#'
#' @details wrapper function that redirects to the proper criterion calculation function
#' according to optiCrit arg value
#' @param matDistanceCorr xxxx
#' @param zoneNModif xxxx
#' @param optiCrit xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
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
