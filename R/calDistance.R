#' calDistance
#'
#' @details description, a paragraph
#' @param typedist xxxx
#' @param tabVal xxxx
#' @param listZonePoint xxxx
#' @param zoneN xxxx
#' @param surfVoronoi xxxx
#' @param meanZone xxxx
#' @param pErr xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calDistance=function(typedist,tabVal,listZonePoint,zoneN,surfVoronoi,meanZone,pErr)
{
    nbPoly=length(listZonePoint)
    resmatDistance=calMatDistance(typedist,zoneN,listZonePoint,tabVal,surfVoronoi,meanZone,pErr)

    return(resmatDistance)
}



