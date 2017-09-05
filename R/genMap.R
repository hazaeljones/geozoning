#######################################################################################################
#' wrapper for randKmap, generate 2D map
#'
#' @details wrapper for randKmap, generate 2D map with 1000 kriged data points, Gaussian field
#' @param DataObj =NULL: simulated data with seed or = a data frame with real data
#' @param seed numeric,
#' @param krig numeric, 1: kriging with vgm model, 2: inverse distance kriging
#' @param Vpsill numeric parameter of the variogram model,
#' @param Vrange numeric parameter of the variogram model,
#' @param Vnugget numeric parameter of the variogram model,
#' @param Vmean numeric parameter of the variogram model,
#' @param nPointsK number of generated points after kriging
#' @param boundary list, contains x and y coordinates of map boundaries
#' @param disp numeric, 
#' @param FULL logical, if TRUE the returned list is complete
#'
#' @return a map object as a list with components
#' \describe{
#' \item{tabAlea}{raw data, SpatialPointsDataFrame}
#' \item{surfaceVoronoi}{Voronoi polygon surfaces}
#' \item{krigTabAlea}{kriged data, SpatialPointsDataFrame}
#' \item{fitVarioAlea}{variogram}
#' \item{DataObj}{DataObj}
#' }
#'
#' @export
#' @examples
#' m=genMap(seed=1,krig=2,disp=1) #generates a map and plots data
#' mean(m$krigGrid) # mean of generated kriged data
#' # not run
genMap=function(DataObj=NULL,seed=80,krig=2,Vpsill=5,Vrange=0.2,Vnugget=0.2,Vmean=8,
       nPointsK=1000,boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),disp=0,FULL=FALSE)
#######################################################################################################
{
  map=randKmap(DataObj,seed=seed,Vpsill=Vpsill,Vrange=Vrange,Vmean=Vmean,Vnugget=Vnugget,
               boundary= boundary,krig=krig,nPoints=450,nPointsK=nPointsK,disp=disp,FULL=FULL)
  # returns NULL if boundary pb
  #view raw data
  if(disp==1) plotMap(map)

  return(map)
}
