###################################################################
#' genMapRNoK
#'
#' @details description, a paragraph
#' @param DataObj xxxx
#' @param seed xxxx
#' @param boundary xxxx
#' @param disp xxxx
#' @param krig xxxx
#' @param FULL xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
genMapRNoK=function(DataObj,seed,boundary,disp=0,krig=1,FULL=FALSE)
####################################################################
{
  #wrapper for map generation with no kriging (data are already on a grid)
  map=randKmapGrid(DataObj,seed=seed,krig=krig,disp=disp,FULL=FULL,boundary=boundary)
  #view raw data
  if(disp==2) plotMap(map)
  return(map)
}
