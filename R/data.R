#' a map
#'
#' A map object for zoning, result from the genMap function
#' @format a list of SpatialPolygons
"mapTest"

#' 
#'
#' A list of initialZoning results
#' @format a list with all results from initialZoning function: criterion value, list with components matDistance, matDistanceCorr and cost, such as returned by a call to calDistance, a list with components zoneN, zoneNModif, listZonePoint, meanTot, meanZone,listSurf, critSurf, zonePolygone, such as the object returned by calNei.
"resZTest"

#' 
#'
#' A data frame with real data used for zoning
#' @format a data frame with 6415 rows and 3 variables:
#' \describe{
#'   \item{x}{x coordinate}
#'   \item{y}{y coordinate}
#'   \item{Yield}{numeric variable - phenotype}
#' }
"yield"

#' 
#'
#' A data frame with simulated data on a regular grid
#' @format a data frame containing a regular grid with 1936 rows and 3 variables
#' \describe{
#'   \item{x}{x coordinate}
#'   \item{y}{y coordinate}
#'   \item{z}{numeric variable - simulated}
#' }
"dataReg"

#' 
#'
#' An external zoning read from a shape file
#' @format a SpatialPolygons object:
"shape1"
