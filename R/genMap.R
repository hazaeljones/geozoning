#######################################################################################################
#' Wrapper for randKmap, generate 2D map
#'
#' @details Wrapper for randKmap, generates a 2D map with a Gaussian field, either by simulating data or by reading data in a data frame. Kriged data are normalized so that x-coordinates are between 0 and 1. y-coordinates are normalized with the same ratio used for x-coordinates. Kriging is either done with inverse distance interpolation, or with a variogram model. It creates a structure that contains the data and parameters necessary to perform a zoning. The structure is identical wether the data are simulated or not.
#' @param DataObj =NULL: simulated data with seed or a data frame with real data
#' @param seed numeric, seed used to randomly generate data points
#' @param krig numeric, 1: kriging with vgm model, 2: inverse distance kriging
#' @param Vpsill numeric parameter of the variogram model,
#' @param Vrange numeric parameter of the variogram model,
#' @param Vnugget numeric parameter of the variogram model,
#' @param Vmean numeric parameter of the variogram model,
#' @param typeMod type of variogram model (see vgm) "Gau", "Sph", "Exp"
#' @param nPointsK number of generated points after kriging
#' @param boundary list, contains x and y coordinates of map boundaries
#' @param disp numeric,
#' @param FULL logical, if TRUE the returned list is complete
#'
#' @return a map object as a list with components
#' \describe{
#' \item{rawData}{simulated or real raw data within the boundary}
#' \item{step}{grid step}
#' \item{krigData}{kriged data as a SpatialPointsDataFrame}
#' \item{krigGrid}{kriged data in form of a grid-useful for image plots.}
#' \item{krigN}{list of neighbours of each kriged data point}
#' \item{krigSurfVoronoi}{list of areas of Voronoi polygons in the tesselation of kriged data}
#' \item{modelGen}{random fields model}
#' \item{VGMmodel}{vgm model}
#' \item{boundary}{(x,y) list of boundary points}
#' \item{ratio}{ratio used to normalize x data}
#' }
#'
#' @export
#' @examples
#' m=genMap(seed=1,krig=2,disp=1) #generates a map and plots data
#' mean(m$krigGrid) # mean of generated kriged data
#' # not run
genMap=function(DataObj=NULL,seed=80,krig=2,Vpsill=5,Vrange=0.2,Vnugget=0.2,Vmean=8,typeMod="Exp",
       nPointsK=1000,boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),disp=0,FULL=FALSE)
#######################################################################################################
{
  map=randKmap(DataObj,seed=seed,Vpsill=Vpsill,Vrange=Vrange,Vmean=Vmean,Vnugget=Vnugget,typeMod=typeMod,
               boundary=boundary,krig=krig,nPoints=450,nPointsK=nPointsK,disp=disp,FULL=FULL)
  # returns NULL if boundary pb
  #view raw data
  if(disp==1) plotMap(map)

  return(map)
}
