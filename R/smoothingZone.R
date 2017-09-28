###################################################################
#' smoothingZone
#' @details function that returns a new smoothed zones. Attention: this function is just a tool for a better visualisation of the map, if it
#' doesn't work properly, please choose another value of the width parameter.
#' @param z zone to be modified (SpatialPolygon)
#' @param width smoothing parameter in gBuffer if dilatation is followed by erosion
#' @param boundary union of all zones of the corrected map (result of correctBoundaryMap())
#' @param disp logical, if TRUE, display the value of "widthExt" in case of dilatation->erosion, otherwise display "widthInt" in case of erosion->dilatation
#' @importFrom rgeos gIsValid
#' @return a zone (SpatialPolygon)
#' @importFrom sp plot
#' @export
#' @examples
#' seed=1
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' criti = correctionTree(qProb = c(0.5), map = map)
#' Z = criti$zk[[1]][[1]]$zonePolygone
#' lab = criti$zk[[1]][[1]]$lab
#' # zones' correction
#' res = correctBoundaryMap(Zi = Z, map = map)
#' Z = res$Z
#' # map boundary after correction
#' boundary = Z[[1]]
#' for(i in 2:length(Z)){
#'   boundary = rgeos::gUnion(boundary, Z[[i]])
#' }
#' # plot map
#' plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
#' # smoothing
#' zone = Z[[2]]
#' newZone = smoothingZone(z = zone, width = 0.05, boundary = boundary)
#' sp::plot(zone)
#' sp::plot(newZone)


smoothingZone = function (z, width, boundary, disp = TRUE)
{
  # this function is in the script smoothingZone.R

  zone = zone.extended(z = z, boundary = boundary)
  widthExt = cal.max.width.Zone(z = zone, step = 0.001, widthMax = width, boundary = boundary, erosion = FALSE)
  if(disp){print(paste("widthExt", widthExt))}
  dilatation1 = gBuffer(zone,width = widthExt,joinStyle="ROUND",capStyle = "ROUND")
  erosion1 = gBuffer(dilatation1,width = -widthExt,joinStyle="ROUND",capStyle = "ROUND")

  widthInt = cal.max.width.Zone(z = erosion1, step = 0.001, widthMax = width, boundary = boundary, erosion = TRUE)
  if(disp){print(paste("widthInt", widthInt))}
  erosion2 = gBuffer(erosion1,width = -widthInt,joinStyle="ROUND", capStyle = "ROUND")
  newZ = gBuffer(erosion2,width = widthInt, joinStyle="ROUND", capStyle = "ROUND")



  # search the intersection between the new smoothed zone and the map
  # 2 ways to do: (we have to check if the geometry is valid)
  # 1st way : we search directly the intersection
  # 2nd way : diff = boundary - newZ , then , newZ = boundary-diff

  # 1st way
  newZ1 = gIntersection(newZ,boundary)
  newZ1 = buffToValid(newZ1)

  # 2nd way
  diff = gDifference(boundary, newZ)
  diff = buffToValid(diff)
  newZ2 = gDifference(boundary,diff)

  if(is.null(newZ1) == FALSE){
    newZ = newZ1
  }else{
    newZ = newZ2
  }
  #print(paste("smoothing valid is", gIsValid(newZ)))
  return(newZ)
}



############################################## zone.extended ################################################################################

#' zone.extended
#' @details for a zone that has commun border with the map, it will be extended at the side of commun border. We search the commun border which is
#' a spatiaLines. This spatialLines is composed of several Lines containing only 2 points. For each Lines, we project the 2 points to the
#' convexHull of the "relaxation" of the map's boundary. We have then 4 points (2 come from a Line, 2 come from the projection). with 4 points,
#' we will have a SpatialPolygone which is the extension part of the Line.

#' @param z a zone of the map
#' @param boundary union of all zones of the corrected map (result of correctBoundaryMap())
#' @importFrom rgeos gOverlaps readWKT gNearestPoints
#' @importFrom raster geom
#' @importFrom sp plot
#' @export
#' @examples
#' seed=1
#' map = genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' criti = correctionTree(qProb = c(0.5), map = map)
#' Z = criti$zk[[1]][[1]]$zonePolygone
#' lab = criti$zk[[1]][[1]]$lab
#' # zones' correction
#' res = correctBoundaryMap(Zi = Z, map = map)
#' Z = res$Z
#' # map boundary after correction
#' boundary = Z[[1]]
#' for(i in 2:length(Z)){
#'   boundary = rgeos::gUnion(boundary, Z[[i]])
#' }
#' # plot map
#' plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
#' # extend zone
#' z = geozoning:::zone.extended(z = Z[[1]], boundary = boundary)
#' sp::plot(z)
#' sp::plot(Z[[1]],add=TRUE)


zone.extended = function (z, boundary)
{
  # this function is in the script smoothingZone.R

  boundaryLineExtend = gBoundary(gConvexHull(gBuffer(boundary,width = 0.2)))

  if(touch.border(z, boundary)){
    lineInter = gIntersection(gBoundary(boundary), z) # intersection of zone and the boundary of the map
    li = geom(lineInter) # transform lineInter into a dataframe
    level = length(levels(as.factor(li[,2]))) # compute the number of pieces of lines in the spatial line
    for(i in 1:level){
      numPoint = which(li[,2]==i) # index of the point that belongs to the piece of the line i
      text = "POLYGON(("
      for (j in numPoint){
        text = paste(text,li[j,4],li[j,5],",")
        point = readWKT(paste("POINT(",li[j,4],li[j,5],")"))
        ptProject = gNearestPoints(point,boundaryLineExtend)@coords[2,]
        text = paste(text,ptProject[1],ptProject[2],",")
      }
      text = paste(text,li[numPoint[1],4],li[numPoint[1],5],"))")
      smallZ.extend = gConvexHull(readWKT(text))
      z = gUnion(z,smallZ.extend)
      z = gBuffer(z,width = 0)
    }
  }
  return(z)
}





################################################## cal.max.width.Zone ############################################################
#' cal.max.width.Zone
#' @details function that return the maximal value of the parameter "width" in function gBuffer in order not to make zone disappear
#' or not to split a zone into 2 differents zones
#' @param z spatial polygon
#' @param step the difference between 2 values of parameter width in the function gBuffer
#' @param widthMax the maximum value of the parameter width in gBuffer
#' @param boundary union of all zones of the corrected map (result of correctBoundaryMap())
#' @param erosion logical, if TRUE, compute the maximum value of width in case erosion->dilatation, otherwise in case dilatation->erosion
#' @return maximum value of parameter width in the function smoothingZone
#' @importFrom rgeos gUnion
#' @importFrom rgeos gBuffer
#' @importFrom rgeos plot
#' @export
#' @examples
#' seed=1
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' criti = correctionTree(qProb = c(0.4,0.6), map = map)
#' Z = criti$zk[[2]][[1]]$zonePolygone
#' lab = criti$zk[[2]][[1]]$lab
#' # zones' correction
#' res = correctBoundaryMap(Zi = Z, map = map)
#' Z = res$Z
#' # map boundary after correction
#' boundary = Z[[1]]
#' for(i in 2:length(Z)){
#'   boundary = rgeos::gUnion(boundary, Z[[i]])
#' }
#' # plot map
#' plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
#' widthMax = cal.max.width.Zone(z = Z[[3]], step = 0.001,
#'            widthMax = 0.05, boundary = boundary, erosion = TRUE)
#' zone = zone.extended(z = Z[[3]], boundary = boundary)
#' erosion1 = rgeos::gBuffer(zone ,width = - (widthMax + 0.002) ,joinStyle="ROUND",capStyle = "ROUND")
#' erosion2 = rgeos::gBuffer(zone ,width = - (widthMax - 0.002) ,joinStyle="ROUND",capStyle = "ROUND")
#' rgeos::plot(erosion1)
#' rgeos::plot(erosion2)

cal.max.width.Zone = function(z, step = 0.001, widthMax = 0.05, boundary, erosion = TRUE)
{
  # this function is in the script smoothingZone.R
  Width = 0.001
  Stop = FALSE

  if(erosion == TRUE){
    while(Stop == FALSE & Width <= widthMax/2){
      buff = gBuffer(z,width = -Width)
      buff = gBuffer(buff,width = Width)
      if(is.null(buff)){
        Width = Width - step
        Stop = TRUE
      }else{
        counterPoly = length(buff@polygons[[1]]@Polygons)
        if (counterPoly>1){
          nbZvalid = 0
          for(i in 1:length(buff@polygons[[1]]@Polygons)){
            if (buff@polygons[[1]]@Polygons[[i]]@area>10^-4){
              nbZvalid = nbZvalid + 1
            }
          }
          if(nbZvalid>1){
            Width = Width - step
            Stop = TRUE
          }else{
            Width = Width + step
          }
        }else{
          Width = Width +step
        }
      }
    }
  }else{
    while(Stop == FALSE & Width <= widthMax){
      buff = gBuffer(z,width = Width)
      buff = gBuffer(buff,width = -Width)
      counterPoly = length(buff@polygons[[1]]@Polygons)
      if (counterPoly>1){
        nbZvalid = 0
        for(i in 1:length(buff@polygons[[1]]@Polygons)){
          if (buff@polygons[[1]]@Polygons[[i]]@area>10^-6){
            nbZvalid = nbZvalid + 1
          }
        }
        if(nbZvalid>1){
          Width = Width - step
          Stop = TRUE
        }else{
          Width = Width + step
        }
      }else{
        Width = Width +step
      }
    }
  }
  return(Width)
}




################################################### touch.border ###################################################
#' touch.border
#' @details verify if a zone has a commun boundary with the map
#' @param z a zone (SpatialPolygon)
#' @param boundary union of all zones of the corrected map (result of correctBoundaryMap())
#' @return logical, TRUE if zone has a commun boundary with the map, FALSE otherwise
#' @importFrom rgeos gUnion
#' @importFrom rgeos gBoundary
#' @export
#' @examples
#' seed=1
#' map = genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' criti = correctionTree(qProb = c(0.5), map = map)
#' Z = criti$zk[[1]][[1]]$zonePolygone
#' lab = criti$zk[[1]][[1]]$lab
#' # zone correction
#' res = correctBoundaryMap(Zi = Z, map = map)
#' Z = res$Z
#' # map boundary after correction
#' boundary = Z[[1]]
#' for(i in 2:length(Z)){
#'   boundary = rgeos::gUnion(boundary, Z[[i]])
#' }
#' # plot map
#' plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
#' # verification
#' for(i in 1:length(Z)){
#'   print(touch.border(z = Z[[i]], boundary = boundary))
#' }
touch.border = function (z, boundary)
{
  # this function is in the script smoothingZone.R
  lineBoundary = gBoundary(boundary) # transform the polygon into a line

  if (gDistance(z, lineBoundary) > 10^-3) {
    res = FALSE
  }
  else{
    res = TRUE
  }
  return(res)
}





############################################# buffToValid #########################################################
#' buffToValid
#' @details function that check if a zone has a valid geometry, if not , makes zone valid by using gBuffer(width = 0,...)
#' @param zone a SpatialPolygon
#' @return a new valid zone
#' @importFrom rgeos gIsValid
#' @importFrom rgeos gBuffer
#' @keywords internal


buffToValid = function(zone){
  # this function is in the script smoothingZone.R
  while(gIsValid(zone) == FALSE){
    zone = gBuffer(zone,width = 0)
  }
  return(zone = zone)
}












