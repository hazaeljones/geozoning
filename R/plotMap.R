##########################
#' plot a map
#'
#' @details plot 3 different graphics of a map object
#' @param map a map object, such as returned by genMap
#'
#' @return a plot
#'
#' @export
#' @importFrom grDevices heat.colors
#' @importFrom graphics contour image locator par plot points
#'
#' @examples
#' m=genMap(seed=1,krig=2,disp=0)
#' plotMap(m)
#' # not run
#'
plotMap=function(map)
##########################
{
  #view generated data
  #persp(map$krigGrid)
  oldpar=par(mfrow=c(2,2))
  image(map$krigGrid,col=rev(heat.colors(20)))
  #
  contour(map$krigGrid)
  #
  contour(map$krigGrid)
  a=data.frame(map$rawData)
  a[,]=a[order(a[,"z"]),]
  zfac=factor(a[,"z"])
  coul=heat.colors(150+length(levels(zfac)))[length(levels(zfac)):1]#pour eviter le blanc
  points(a, col = coul)
  par(oldpar)
}
