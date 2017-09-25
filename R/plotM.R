#' plotM
#'
#' @details plot the map in color with zones and details.
#'
#' @param map object returned by function genMap or genMapR.
#' @param Z list of zones, each zone is a SpatialPolygons.
#' @param lab label of each zones.
#' @param byLab boolean, if TRUE display the label of each zone, else display the zone number.
#' @param quantile probability vector used to generate "Z". This will be displayed in the title of the plot.
#' @param crit criterion value corresponding to "Z. This will be displayed in the title of the plot.
#' @param cost cost value corresponding to "Z". This will be displayed in the title of the plot.
#' @param bestCrit best criterion value. This will be displayed in the title of the plot.
#' @param bestCost best cost value. This will be displayed in the title of the plot.
#' @param newCost new cost value. This will be displayed in the title of the plot.
#' @param line position of the title. if negative, the title goes down, otherwise, goes up.
#' @param cex text size
#' @return an empty value
#'
#' @importFrom rgeos readWKT
#' @importFrom sp plot
#' @export
#' @examples
#' seed=2
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.55,0.85),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#'
plotM = function(map, Z=NULL, lab=NULL, byLab = TRUE, 
                 quantile = NULL, crit = NULL, cost = NULL, bestCrit = NULL, bestCost = NULL, newCost=NULL, line = 0, cex = 2)
{
  step = map$step
  xsize = map$xsize
  ysize = map$ysize
  x = seq(step,xsize-step,step)
  y = seq(step,ysize-step,step)
  # plot map with color
  image.plot(x = x, y = y, z = map$krigGrid, topo.colors(20), xlim = c(0,xsize), ylim = c(0, ysize), legend.shrink = 0.6, legend.width = 0.8)

  # plot zones
  if(! is.null(Z)){
    for (i in 1:length(Z)){
      sp::plot(Z[[i]],add= TRUE, lwd = 1.5)
    }
  }

  if(! is.null(lab)){
    # plot zones label or zones numbers
    listZpt = zoneAssign(tab = map$krigData, Z = Z)
    for (i in 1:length(Z)){
      boundaryZone = gBoundary(Z[[i]])
      dist = c()
      for (p in listZpt[[i]]){
        point = readWKT(paste("POINT(",map$krigData@coords[p,1],map$krigData@coords[p,2],")"))
        d = gDistance(point, boundaryZone)
        dist = c(dist,d)
      }
      index = which.max(dist)
      indexP = listZpt[[i]][index]
      if(byLab == FALSE){
        text(map$krigData@coords[indexP,1],map$krigData@coords[indexP,2], labels = i,cex = cex)
      }else{
        text(map$krigData@coords[indexP,1],map$krigData@coords[indexP,2], labels = lab[i],cex = cex)
      }
    }
  }

  Title = NULL
  if(!is.null(quantile)){
    Title = "quantile = ["
    for(i in 1:length(quantile)){
      Title = paste(Title, quantile[i])
    }
    Title = paste(Title, "] \n")
  }
  if(! is.null(crit)){
    Title = paste(Title, "crit =", crit, "  ")
  }
  if(! is.null(cost)){
    Title = paste(Title, "cost =", cost, "  ")
  }
  if(! is.null(bestCrit)){
    Title = paste(Title, "best crit =", bestCrit, "  ")
  }
  if(! is.null(bestCost)){
    Title = paste(Title, "best cost =", bestCost, "  ")
  }
  if(! is.null(newCost)){
    Title = paste(Title, "new cost =", newCost)
  }

  title(main = Title, line = 0.2)

}
