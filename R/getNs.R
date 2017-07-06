##################################################################
#' getNs
#'
#' @details description, a paragraph
#' @param zoneNModif zone neighborhood Logical matrix  
#' @param iZ index of current zone in zoning
#'
#' @return a Logical vector of current zone neighbors
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' K=resZTest
#' Ns=getNs(K$zoneNModif,5) # find neighbors of zone 5
##################################################################
getNs=function(zoneNModif,iZ)
{
  Ns=zoneNModif[iZ,]
  return(Ns)
}
