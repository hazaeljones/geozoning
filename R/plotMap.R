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
#' @importFrom graphics contour image locator par plot points title
#'
#' @examples
#' data(mapTest)
#' plotMap(mapTest)
plotMap=function(map)
##########################
{
  #view generated data
  oldpar=par(mfrow=c(2,2))
  matVal=map$krigGrid
  cn=as.numeric(colnames(matVal))
  rn=as.numeric(rownames(matVal))
  image.plot(rn,cn,matVal,col=rev(heat.colors(20)),xlab="",ylab="")
  title("Kriged data")
  #
  contour(rn,cn,matVal)
  title("Contour lines on kriged data")
  #
  contour(rn,cn,matVal)
  a=data.frame(map$rawData)
  a[,]=a[order(a[,"z"]),]
  zfac=factor(a[,"z"])
  coul=heat.colors(150+length(levels(zfac)))[length(levels(zfac)):1]#pour eviter le blanc
  points(a, col = coul)
  title("Contour lines on kriged data \nplus raw data points")
  par(oldpar)
}
