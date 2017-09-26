########################################## Cost_By_Laplace ####################################################################
#' Cost_By_Laplace

#' @details : function that returns the criterion COST by approximating the valeur in a point of the grid by
#'            the linear interpolation (approximate solution of Laplace's equation. For more details see help of function
#'            Transition_Zone_Near_Boundary, Transition_Zone_Far_Boundary or Extreme_Zone)

#' @param map  object returned by function genMap 
#' @param Z : an example of zoning (a list of zones)
#' @param numZ : number of the zone in which the cost will be computed
#' @param Estimation : value of linear interpolation by solving Laplace's equation

#' @return cost computed by replacing values in zone by linear interpolation
#' @importFrom rgeos readWKT
#' @export

#' @examples
#' seed=2
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.55,0.85),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 6 is a transition zone that has commun boundary with the map
#' numZ = 6
#' Estimation = Transition_Zone_Near_Boundary(map = map, Z = Z, numZ = numZ)
#' # compute the cost
#' cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
#' cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
#' print(cL$cost_Laplace)
#' print(cM$cost_Mean)
#' # zone 6 is a zone with gradient


Cost_By_Laplace = function (map, Z, numZ, Estimation)

{
  # this function is in the script "calCost.R"
  # compute the contours (spatialLines) and their corresponding quantile value
  # krigData
  tab = map$krigData
  # reassign points to zones in Z
  listZpt = zoneAssign(tab, Z)

  Value = tab@data$var1.pred[listZpt[[numZ]]]
  Surface = map$krigSurfVoronoi[listZpt[[numZ]]]

  cost_Laplace= sum((Value - Estimation) ^ 2 * Surface) / sum(Surface)

  return(list(cost_Laplace = cost_Laplace, Surf = sum(Surface)))
}



########################################## Cost_By_Mean ###########################################################################

#' Cost_By_Mean

#' @details : function that returns the criterion COST by approximating the value at a point of the grid by
#'            the mean value of the zone.

#' @param map  object returned by function genMap
#' @param Z : an example of zoning (a list of zones)
#' @param numZ : number of the zone in which the cost will be computed

#' @return the cost value as described
#' @export
#' @examples
#' seed=2
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.55,0.85),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 6 is a transition zone that has commun boundary with the map
#' numZ = 6
#' Estimation = Transition_Zone_Near_Boundary(map = map, Z = Z, numZ = numZ)
#' # compute the cost
#' cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
#' cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
#' print(cL$cost_Laplace)
#' print(cM$cost_Mean)
#' # zone 6 is a zone with gradient

Cost_By_Mean = function(map, Z, numZ)
{
  # this function is in the script "calCost.R"
  # krigData
  tab = map$krigData
  # reassign points to zones in Z
  listZpt = zoneAssign(tab, Z)

  Value = tab@data$var1.pred[listZpt[[numZ]]]
  Surface = map$krigSurfVoronoi[listZpt[[numZ]]]

  Mean = sum((Value * Surface)) / sum(Surface)

  cost_Mean = sum((Value - Mean) ^ 2 * Surface) / sum(Surface)

  return(list(cost_Mean=cost_Mean, Surf = sum(Surface)))
}


########################################## Points_Near_Boundary ###############################################################

#' Points_Near_Boundary

#' @details function that returns a list of points in a zone that are near boundary of the map

#' @param map object returned by function genMap 
#' @keywords internal
#' @examples
#' map = mapTest
#' geozoning:::Points_Near_Boundary(map = map)

Points_Near_Boundary = function(map){
  # this function is in the script "calCost.R"
  X = map$boundary$x
  Y = map$boundary$y
  l = cbind(X,Y)
  L= Line(l) # create object of class Line
  Ls = Lines(list(L), ID = "1") # create object of class Lines
  boundaryLine = SpatialLines(list(Ls)) # create object of class SpatialLines
  result = c()
  for(i in 1:nrow(map$krigData@coords)){
    x = map$krigData@coords[i,1]
    y = map$krigData@coords[i,2]
    point = readWKT(paste("POINT(",x,y,")"))
    d = gDistance(point,boundaryLine)
    if(d <= map$step+10^-6){ # we can set the threshold equal to map$step
      result = c(result, i)
    }
  }
  return(result)
}


########################################## new_krigGrid_for_visualisation ####################################################
#' new_krigGrid_for_visualisation
#' @details Elementary function that create a new krigGrid by replacing the real values by the approximation of the function
#' "Transition_Zone_Near_Boundary","Transition_Zone_Far_Boundary" or "Extreme_Zone" in order to have a look at the new iso contour

#' @param map object returned by function genMap 
#' @param Z list of zones.
#' @param numZ number of the zone whose values will be approximated.
#' @param solution the result of function "Transition_Zone_Near_Boundary" or "Transition_Zone_Far_Boundary" or "Extreme_Zone"

#' @return new krigGrid and new data
#' importFrom sp plot
#' @export
#' @examples
#' seed=2
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.55,0.85),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 6 is a transition zone that has commun boundary with the map
#' numZ = 6
#' Estimation = Transition_Zone_Near_Boundary(map = map, Z = Z, numZ = numZ)
#' result = new_krigGrid_for_visualisation(map = map, Z = Z, numZ = numZ, solution = Estimation)
#' new_krigGrid = result$new_krigGrid
#' new_data = result$new_data
#' quant1 = quantile(map$krigData@data$var1.pred,probs = 0.55)
#' quant2 = quantile(map$krigData@data$var1.pred,probs = 0.85)
#' # plot initial isocontours
#' plotM(map = map,Z = Z,lab = lab, byLab = TRUE)
#' listContours = contourBetween(map = map, krigGrid = map$krigGrid, q1 = quant1, q2 = quant2)
#' for (i in 1:length(listContours)){
#'   sp::plot(listContours[[i]]$contour,add=TRUE,col = "red")
#' }
#' # plot modified isocontours
#' plotM(map = map,Z = Z,lab = lab, byLab = TRUE)
#' listContours = contourBetween(map = map, krigGrid = new_krigGrid, q1 = quant1, q2 = quant2)
#' for (i in 1:length(listContours)){
#'   sp::plot(listContours[[i]]$contour,add=TRUE,col = "red")
#' }


new_krigGrid_for_visualisation = function(map, Z, numZ, solution){

  # this function is in the script "calCost.R"
  tab = map$krigData
  listZpt = zoneAssign(tab = tab, Z = Z)

  new_krigGrid = map$krigGrid
  new_data = map$krigData@data$var1.pred

  for( i in 1:length(listZpt[[numZ]]) ){
    numPoint = listZpt[[numZ]] [i]
    new_data[numPoint] = solution[i]

    # find location of point on the grid

    count = 1
    X = nrow(new_krigGrid)
    Y = ncol(new_krigGrid)
    x = 1
    y = 1
    while(count<numPoint){
      if(x<X){
        if(! is.na(new_krigGrid[x,y])){
          count = count+1
        }
        x = x+1
      }else{
        if(! is.na(new_krigGrid[x,y])){
          count = count+1
        }
        x = 1
        y = y+1
      }
    }
    new_krigGrid[x,y] = solution[i]
  }
  return(list(new_krigGrid=new_krigGrid, new_data=new_data))
}



############################################# Transition_Zone_Near_Boundary ##################################################
#' Transition_Zone_Near_Boundary
#' @details funtion that approximates the value in a transition zone (which has commun boundary with the map)
#' by the solution of the Laplace's equation. The numerical resolution of the Laplace's equation will be based on the discretisation
#' of the data on the grid (map$krigGrid). The domaine of study is a transition zone which have a commun border with the map.

#' @param map object returned by function genMap
#' @param Z list of zones.
#' @param numZ number of the zone whose values will be approximated.

#' @return approximated values of the zone (numZ) given as parameter.
#' @export
#' @examples
#' seed=2
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.55,0.85),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 6 is a transition zone that has commun boundary with the map
#' numZ = 6
#' Estimation = Transition_Zone_Near_Boundary(map = map, Z = Z, numZ = numZ)
#' # compute the cost
#' cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
#' cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
#' print(cL$cost_Laplace)
#' print(cM$cost_Mean)
#' # zone 6 is a zone with gradient

Transition_Zone_Near_Boundary = function(map, Z, numZ){

  # this function is in the script "calCost.R"
  # for each zone, identify its neighbours.
  listNei = list()
  for (i in 1:length(Z)){
    nei = c()
    for (j in 1:length(Z)){
      if(j!=i){
        if(gDistance(Z[[i]],Z[[j]]) < 0.001){ # if the distance between 2 zones is small, they are neighbours
          nei = c(nei, j)
        }
      }
    }
    listNei[[i]] = nei
  }

  # reassign points to zones
  tab = map$krigData
  listZpt = zoneAssign(tab = tab, Z = Z)

  # find all points in zone "numZ" that are near boundary of the map
  pointNearBoundary = Points_Near_Boundary(map = map)
  pointNearBoundary_inZone = intersect(pointNearBoundary, listZpt[[numZ]])

  # for each point in pointNearBoundary_inZone, find the closest points that is in neighbour zones and is in the list pointNearBoundary
  # in this case, there are 2 neighbour zones, so each point in pointNearBoundary_inZone will has 2 points described as above
  listPointNei = list() # each element of the list is a couple de points, each point is in a neighbour zone
  k = 1 # counter
  for(i in pointNearBoundary_inZone){
    pointNei = c()
    # coordinates of point i
    xi = map$krigData@coords[i,1]
    yi = map$krigData@coords[i,2]

    for (zoneNei in listNei[[numZ]]){
      dist = c()
      pointInZoneNei = c()
      for(j in pointNearBoundary){
        if(j %in% listZpt[[zoneNei]]){ # check if point j is in the neighbour zone
          pointInZoneNei = c(pointInZoneNei, j)
          # coordinates of point j
          xj = map$krigData@coords[j,1]
          yj = map$krigData@coords[j,2]

          d = sqrt((xi-xj)^2+(yi-yj)^2)  # compute the distance between 2 points
          dist = c(dist, d)
        }
      }
      index = which.min(dist) # find index of the smallest distance
      pointNei = c(pointNei, pointInZoneNei[index])
    }
    listPointNei[[k]] = pointNei
    k = k+1
  }

  # now, for each point in "pointNearBoundary_inZone", interpolate linearly the value of the point by the 2 points in "listPointNei"
  # this approximation will be use for resolution of Laplace's equation
  Values_Interpolated = c()
  for(i in 1:length(listPointNei)){
    # point and coordinates
    p = pointNearBoundary_inZone[i]
    x = map$krigData@coords[p,1]
    y = map$krigData@coords[p,2]
    p1 = listPointNei[[i]][1]
    x1 = map$krigData@coords[p1,1]
    y1 = map$krigData@coords[p1,2]
    p2 = listPointNei[[i]][2]
    x2 = map$krigData@coords[p2,1]
    y2 = map$krigData@coords[p2,2]

    # compute distance between p and p1, p and p2
    d1 = sqrt((x-x1)^2+(y-y1)^2)
    d2 = sqrt((x-x2)^2+(y-y2)^2)

    # values of p1 and p2
    v1 = map$krigData@data$var1.pred[p1]
    v2 = map$krigData@data$var1.pred[p2]

    # interpolate "linearly"
    d = d1+d2
    val_interpolated = v1*d2/d + v2*d1/d
    Values_Interpolated = c(Values_Interpolated, val_interpolated)
  }

  # resolution of Laplace's equation : modeling the equation by using the grid of map, we obtain a linear form AX = B

  N = length(listZpt[[numZ]])
  A = matrix(rep(0,N*N), nrow = N)
  B = rep(0,N)

  # the diagonal of A is equal to -4
  for (i in 1:N){
    A[i,i] = -4
  }
  # fill the matrix A and the vector B
  coord = coordinates(tab)
  X = coord[,1]
  Y = coord[,2]

  #step of discretisation
  step = map$step

  for(i in 1:length(listZpt[[numZ]])){
    # coordinates of point i in zone "numZ"
    x = X[listZpt[[numZ]][i]]
    y = Y[listZpt[[numZ]][i]]
    # coordinates of the neighbour point in the West
    xw = x - step
    yw = y
    # coordinates of the neighbour point in the Nord
    xn = x
    yn = y + step
    # coordinates of the neighbour point in the Est
    xe = x + step
    ye = y
    # coordinates of the neighbour point in the South
    xs = x
    ys = y - step

    # FIND INDEX OF POINT IN THE WEST
    jw = intersect(which(abs(Y-yw)<10^-6),which(abs(X-xw)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    if(length(jw)!=0){ # if neighbour point exists
      j = which(listZpt[[numZ]]==jw)
      if(length(j)!=0){ # if neighbour point is in zone "numZ"
        A[i,j] = 1
      }else{ # if neighbour point is not in zone "numZ"
        B[i] = B[i] + tab@data[jw,1]
      }
    }else{ # if neighbour point doesn't exist (neighbour point in boundary)
      B[i] = B[i] + Values_Interpolated [which(pointNearBoundary_inZone == listZpt[[numZ]][i])]
    }
    # FIND INDEX OF POINT IN THE NORD
    jn = intersect(which(abs(Y-yn)<10^-6),which(abs(X-xn)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    if(length(jn)!=0){
      j = which(listZpt[[numZ]]==jn)
      if(length(j)!=0){
        A[i,j] = 1
      }else{
        B[i] = B[i] + tab@data[jn,1]
      }
    }else{
      B[i] = B[i] + Values_Interpolated [which(pointNearBoundary_inZone == listZpt[[numZ]][i])]
    }
    # FIND INDEX OF POINT IN THE EST
    je = intersect(which(abs(Y-ye)<10^-6),which(abs(X-xe)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    if(length(je)!=0){
      j = which(listZpt[[numZ]]==je)
      if(length(j)!=0){
        A[i,j] = 1
      }else{
        B[i] = B[i] + tab@data[je,1]
      }
    }else{
      B[i] = B[i] + Values_Interpolated [which(pointNearBoundary_inZone == listZpt[[numZ]][i])]
    }
    # FIND INDEX OF POINT IN THE SOUTH
    js = intersect(which(abs(Y-ys)<10^-6),which(abs(X-xs)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    if(length(js)!=0){
      j = which(listZpt[[numZ]]==js)
      if(length(j)!=0){
        A[i,j] = 1
      }else{
        B[i] = B[i] + tab@data[js, 1]
      }
    }else{
      B[i] = B[i] + Values_Interpolated [which(pointNearBoundary_inZone == listZpt[[numZ]][i])]
    }
  }

  # resolution
  solution = -B%*% (solve(A))
  return(solution)

}



############################################# Transition_Zone_Far_Boundary #######################################################

#' Transition_Zone_Far_Boundary
#' @details funtion that approximates the value in a transition zone (which doesn't have commun boundary with the map)
#' by the solution of the Laplace's equation. The numerical resolution of the Laplace's equation will be based on the discretisation
#' of the data on the grid (map$krigGrid).
#' @usage Transition_Zone_Far_Boundary(map, Z, numZ)
#' @param map object returned by function genMap 
#' @param Z list of zones.
#' @param numZ number of the zone whose values will be approximated.

#' @return approximated values of the zone (numZ) given as parameter.
#' @export
#' @examples
#' seed=9
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.65,0.8),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 7 is a transition zone that is far from map boundary
#' numZ = 7
#' Estimation = Transition_Zone_Far_Boundary(map = map, Z = Z, numZ = numZ)
#' # compute the cost
#' cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
#' cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
#' print(cL$cost_Laplace)
#' print(cM$cost_Mean)
#' # zone 7 is a zone with gradient

Transition_Zone_Far_Boundary = function(map, Z, numZ){

  # this function is in the script "calCost.R"

  # reassign points to zones
  tab = map$krigData
  listZpt = zoneAssign(tab = tab, Z = Z)



  # resolution of Laplace's equation : modeling the equation by using the grid of map, we obtain a linear form AX = B

  N = length(listZpt[[numZ]])
  A = matrix(rep(0,N*N), nrow = N)
  B = rep(0,N)

  # the diagonal of A is equal to -4
  for (i in 1:N){
    A[i,i] = -4
  }
  # fill the matrix A and the vector B
  coord = coordinates(tab)
  X = coord[,1]
  Y = coord[,2]

  #step of discretisation
  step = map$step

  for(i in 1:length(listZpt[[numZ]])){
    # coordinates of point i in zone "numZ"
    x = X[listZpt[[numZ]][i]]
    y = Y[listZpt[[numZ]][i]]
    # coordinates of the neighbour point in the West
    xw = x - step
    yw = y
    # coordinates of the neighbour point in the Nord
    xn = x
    yn = y + step
    # coordinates of the neighbour point in the Est
    xe = x + step
    ye = y
    # coordinates of the neighbour point in the South
    xs = x
    ys = y - step

    # FIND INDEX OF POINT IN THE WEST
    jw = intersect(which(abs(Y-yw)<10^-6),which(abs(X-xw)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    j = which(listZpt[[numZ]]==jw)
    if(length(j)!=0){ # if neighbour point is in zone "numZ"
      A[i,j] = 1
    }else{ # if neighbour point is not in zone "numZ"
      B[i] = B[i] + tab@data[jw,1]
    }

    # FIND INDEX OF POINT IN THE NORD
    jn = intersect(which(abs(Y-yn)<10^-6),which(abs(X-xn)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    j = which(listZpt[[numZ]]==jn)
    if(length(j)!=0){
      A[i,j] = 1
    }else{
      B[i] = B[i] + tab@data[jn,1]
    }

    # FIND INDEX OF POINT IN THE EST
    je = intersect(which(abs(Y-ye)<10^-6),which(abs(X-xe)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    j = which(listZpt[[numZ]]==je)
    if(length(j)!=0){
      A[i,j] = 1
    }else{
      B[i] = B[i] + tab@data[je,1]
    }

    # FIND INDEX OF POINT IN THE SOUTH
    js = intersect(which(abs(Y-ys)<10^-6),which(abs(X-xs)<10^-6)) # in stead of finding equality, we use a threshold to detect the point
    j = which(listZpt[[numZ]]==js)
    if(length(j)!=0){
      A[i,j] = 1
    }else{
      B[i] = B[i] + tab@data[js,1]
    }

  }

  # resolution
  solution = -B%*% (solve(A))
  return(solution)

}



############################################# Extreme_Zone ####################################################################
#' Extreme_Zone
#' @details funtion that approximates the value in a extreme zone (zone with label maximum or minimum, zones which have only one neighbour)
#' by the solution of the Laplace's equation. The iso contours plotted on the approximate data will take the form of concentric circles as we
#' supposed the extreme value of the zone is at the zone center (furthest point from the zone boundary.)

#' @param map object returned by function genMap
#' @param Z list of zones.
#' @param numZ number of the zone whose values will be approximated.
#' @param label.is.min boolean value that is TRUE if the label of the zone is minimum and FALSE if the label is maximum

#' @return approximated values of the values in zone (numZ).
#' @importFrom rgeos readWKT
#' @export
#' @examples
#' seed=6
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.8),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 2 is a zone with maximum label
#' numZ = 2
#' Estimation = Extreme_Zone(map = map, Z = Z, numZ = numZ, label.is.min = FALSE)
#' # compute the cost
#' cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
#' cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
#' print(cL$cost_Laplace)
#' print(cM$cost_Mean)
#' # zone 2 is homogeneous

Extreme_Zone = function(map, Z, numZ, label.is.min = TRUE){

  # this function is in the script "calCost.R"

  # reassign points to zones
  tab = map$krigData
  listZpt = zoneAssign(tab = tab, Z = Z)

  # find the farthest point from the boundary of the zone "numZ", then we will interpolate its value by the extreme value of the zone
  boundaryZone = gBoundary(Z[[numZ]])
  dist = c()
  #plot(boundaryZone)
  for(i in listZpt[[numZ]]){
    # récupérer les coords du point
    x = map$krigData@coords[i,1]
    y = map$krigData@coords[i,2]
    p = readWKT(paste("POINT(",x,y,")"))
    dist = c(dist, gDistance(p, boundaryZone))
    #text(x,y,labels = i)
  }
  index = which.max(dist)
  indexP = listZpt[[numZ]] [[index]] # index of the point

  listPt = listZpt[[numZ]] [ listZpt[[numZ]] != listZpt[[numZ]] [index] ] # list contains all points of zone except the furthest from zone's boundary

  # find the extreme value of the zone
  Val = 0
  if(label.is.min == TRUE){
    Val = min(map$krigData@data$var1.pred[listZpt[[numZ]]])
  }else{
    Val = max(map$krigData@data$var1.pred[listZpt[[numZ]]])
  }

  # find the other extreme value of the zone
  # in case that zone has commun border with the boundary of the map, we will interolate the points in boundary of the map that are close
  # to the zone by this value
  Val2 = 0
  if(label.is.min == TRUE){
    Val2 = max(map$krigData@data$var1.pred[listZpt[[numZ]]])
  }else{
    Val2 = min(map$krigData@data$var1.pred[listZpt[[numZ]]])
  }

  # resolution of Laplace's equation : modeling the equation by using the grid of map, we obtain a linear form AX = B

  N = length(listPt)
  A = matrix(rep(0,N*N), nrow = N)
  B = rep(0,N)

  # the diagonal of A is equal to -4
  for (i in 1:N){
    A[i,i] = -4
  }
  # fill the matrix A and the vector B
  coord = coordinates(tab)
  X = coord[,1]
  Y = coord[,2]

  #step of discretisation
  step = map$step


  for(i in 1:length(listPt)){
    # coordinates of point i
    x = X[listPt[i]]
    y = Y[listPt[i]]
    # coordinates of the neighbour point in the West
    xw = x - step
    yw = y
    # coordinates of the neighbour point in the North
    xn = x
    yn = y + step
    # coordinates of the neighbour point in the Est
    xe = x + step
    ye = y
    # coordinates of the neighbour point in the South
    xs = x
    ys = y - step

    # FIND INDEX OF POINT IN THE WEST
    jw = intersect(which(abs(Y-yw)<10^-6),which(abs(X-xw)<10^-6))
    if(length(jw)!=0){
      if(jw == indexP){
        B[i] = B[i] + Val
      }else{
        j = which(listPt==jw)
        if(length(j)!=0){
          A[i,j] = 1
        }else{
          B[i] = B[i] + tab@data[jw,1]
        }
      }
    }else{
      B[i] = B[i] + Val2
    }
    # FIND INDEX OF POINT IN THE NORTH
    jn = intersect(which(abs(Y-yn)<10^-6),which(abs(X-xn)<10^-6))
    if(length(jn)!=0){
      if(jn == indexP){
        B[i] = B[i] + Val
      }else{
        j = which(listPt==jn)
        if(length(j)!=0){
          A[i,j] = 1
        }else{
          B[i] = B[i] + tab@data[jn,1]
        }
      }
    }else{
      B[i] = B[i] + Val2
    }
    # FIND INDEX OF POINT IN THE EST
    je = intersect(which(abs(Y-ye)<10^-6),which(abs(X-xe)<10^-6))
    if(length(je)!=0){
      if(je == indexP){
        B[i] = B[i] + Val
      }else{
        j = which(listPt==je)
        if(length(j)!=0){
          A[i,j] = 1
        }else{
          B[i] = B[i] + tab@data[je,1]
        }
      }
    }else{
      B[i] = B[i] + Val2
    }
    # FIND INDEX OF POINT IN THE SOUTH
    js = intersect(which(abs(Y-ys)<10^-6),which(abs(X-xs)<10^-6))
    if(length(js)!=0){
      if(js == indexP){
        B[i] = B[i] + Val
      }else{
        j = which(listPt==js)
        if(length(j)!=0){
          A[i,j] = 1
        }else{
          B[i] = B[i] + tab@data[js,1]
        }
      }
    }else{
      B[i] = B[i] + Val2
    }
  }

  # resolution
  solution = -B%*% (solve(A))
  solution = append(solution, Val, after = index-1)

  return(solution)

}



############################################# list_Zone_2_Neighbours #################################################################

#' list_Zone_2_Neighbours
#' @details Returns the numbers of zones that have exactly 2 neighbours with different labels. These zone are susceptible to be
#' transitions zones
#' @param Z list of Zones
#' @param lab vector labels of zones
#' @return a vector containing zone numbers
#' @keywords internal
#' @examples
#' seed=6
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.67,0.8),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 4 and 6 are transition zones and have exactly 2 neighbours with different labels.
#' list_Zone_2_Neighbours(Z = Z, lab = lab)

list_Zone_2_Neighbours = function(Z, lab){

  # this function is in the script calCost.R()
  # find the neighbours of each zone
  listNei = list()
  for (i in 1:length(Z)){
    nei = c()
    for (j in 1:length(Z)){
      if(j!=i){
        if(gDistance(Z[[i]],Z[[j]]) < 0.001){ # if the distance between 2 zones is small, they are neighbours
          nei = c(nei, j)
        }
      }
    }
    listNei[[i]] = nei
  }
  # find zones which have exactly 2 neighbours with differents labels
  zone_transition = c()
  for (i in 1:length(Z)){
    if(length(listNei[[i]])==2){
      if(lab[listNei[[i]][1] ] != lab[ listNei[[i]][2] ] ){
        zone_transition = c(zone_transition, i)
      }
    }
  }
  return(zone_transition)
}











############################################# contourBetween #############################################
#' contourBetween
#' @details  : For the given krigGrid, this funtion returns
#'                the contourLines of the map following the 2 quantiles that defined at the beginning.
#'
#' @param map : object map defined in package geozoning
#' @param krigGrid : object that can
#' @param q1,q2 : 2 quantiles that defined zone
#' @param nbContourBetween : the number of discretisation between q1 and q2

#' @return listContours : List of Spatial Lines and the value of quantile that represent the contours generated
#' @importFrom sp plot
#' @export
#' @examples
#' seed=2
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
#' ZK=initialZoning(qProb=c(0.55,0.85),map)
#' Z=ZK$resZ$zonePolygone # list of zones
#' lab = ZK$resZ$lab # label of zones
#' plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
#' # zone 6 is a transition zone that has commun boundary with the map
#' numZ = 6
#' Estimation = Transition_Zone_Near_Boundary(map = map, Z = Z, numZ = numZ)
#' result = new_krigGrid_for_visualisation(map = map, Z = Z, numZ = numZ, solution = Estimation)
#' new_krigGrid = result$new_krigGrid
#' new_data = result$new_data
#' quant1 = quantile(map$krigData@data$var1.pred,probs = 0.55)
#' quant2 = quantile(map$krigData@data$var1.pred,probs = 0.85)
#' # plot modified isocontours
#' plotM(map = map,Z = Z,lab = lab, byLab = TRUE)
#' listContours = contourBetween(map = map, krigGrid = new_krigGrid, q1 = quant1, q2 = quant2)
#' for (i in 1:length(listContours)){
#'   sp::plot(listContours[[i]]$contour,add=TRUE,col = "red")
#' }


contourBetween = function(map,krigGrid, q1, q2, nbContourBetween = 5)
{
  # this function is the script "calCost.R"
  step = map$step
  xsize = map$xsize
  ysize = map$ysize
  ctL = contourLines(seq(step, xsize-step, by=step), seq(step, ysize-step, by=step), krigGrid, levels = seq(q1,q2,(q2-q1)/(nbContourBetween+1)))
  listContours = list()
  for (i in 1:length(ctL)){
    quantile = ctL[[i]]$level
    x = ctL[[i]]$x
    y = ctL[[i]]$y
    l =cbind(x,y) # coordinates of a line
    L= Line(l) # create object of class Line
    Ls = Lines(list(L), ID = "1") # create object of class Lines
    SLs = SpatialLines(list(Ls)) # create object of class SpatialLines
    listContours[[i]] = list(quantile = quantile, contour = SLs)
  }
  return(listContours = listContours)
}
