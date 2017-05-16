##########################################################################
#' genMapR
#'
#' @details description, a paragraph
#' @param DataObj xxxx
#' @param seed xxxx
#' @param disp xxxx
#' @param krig xxxx
#' @param boundary xxxx
#' @param nPointsK xxxx
#' @param FULL xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
genMapR=function(DataObj,seed=NULL,disp=0,krig=1,boundary=boundary,nPointsK=7000,FULL=FALSE)
##########################################################################
{
  # wrapper for real data map generation
  map=randKmap(DataObj,seed=seed,boundary=boundary,nPointsK=nPointsK,krig=krig,disp=disp,FULL=FULL)
  #view raw data
  if(disp==2) plotMap(map)
  return(map)
}
