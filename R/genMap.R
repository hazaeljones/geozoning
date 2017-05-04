#######################################################################################################
#' wrapper for randKmap, generate 2D map
#'
#' @details wrapper for randKmap, generate 2D map
#' @param DataObj =NULL: simulated data with seed or = a data frame with real data
#' @param seed numeric,
#' @param krig numeric, 1: kriging with vgm model, 2: inverse distance kriging
#' @param Vpsill numeric,
#' @param Vrange numeric,
#' @param Vnugget numeric,
#' @param Vmean numeric,
#' @param disp numeric,
#' @param FULL logical, if TRUE the returned list is complete
#'
#' @return a map in a list
#' \describe{
#' \item{tabAlea}{raw data, SpatialPointsDataFrame}
#' \item{surfaceVoronoi}{Voronoi polygon surfaces}
#' \item{krigTabAlea}{kriged data, SpatialPointsDataFrame}
#' \item{fitVarioAlea}{variogram}
#' \item{DataObj}{DataObj}
#' }
#'
#' @export
#'
#' @examples
#' # not run
genMap=function(DataObj=NULL,seed=80,krig=1,Vpsill=5,Vrange=0.2,Vnugget=0,Vmean=8,disp=0,FULL=FALSE)
#######################################################################################################
{
  map=randKmap(DataObj,seed=seed,Vpsill=Vpsill,Vrange=Vrange,Vnugget=Vnugget,Vmean=Vmean,krig=krig,disp=disp,FULL=FULL)
  # returns NULL if boundary pb
  #view raw data
  if(disp==2) plotMap(map)

  return(map)
}
