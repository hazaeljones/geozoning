####################################################################
#' detection of narrow zones (ratio area/perimeter^2)
#'
#' @details computes for each zone of a zoning the ratio area/squared perimeter
#' @param zonePolygone zoning
#'
#' @return a numerical value
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' calcCritNarrow(resZTest$zonePolygone)
#' # not run
calcCritNarrow=function(zonePolygone)
####################################################################
{
  #surfaces of polygons
  listSurface=as.numeric(lapply(zonePolygone,gArea))

  #criterion: (area/perim^2)
  listPerim=as.numeric(lapply(zonePolygone,gLength))
  critNarrow=listSurface/(listPerim^2)

  return(critNarrow)
}


####################################################################
#' wMean
#'
#' @details computes weighted mean or squared mean of zone data
#' @param type 1-squared mean, 2-mean
#' @param listZonePoint list of data points belonging to zone
#' @param surfVoronoi areas of Voronoi polygon corresponding to data points
#' @param data SpatialPointsDataFrame
#'
#' @return a vector of mean zone values
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' wMean(1,K$listZonePoint,mapTest$krigSurfVoronoi,mapTest$krigData)
#' # not run
wMean=function(type,listZonePoint,surfVoronoi,data)
####################################################################
{
  vectMean=numeric()
  nbPoly=length(listZonePoint)
  if(type==1)
  {
    for (i in 1:nbPoly)
    {

      # squared weighted mean
      vectMean[i]=sum((data[[1]][listZonePoint[[i]]])^2 *surfVoronoi[listZonePoint[[i]]])/sum((surfVoronoi[listZonePoint[[i]]]))

    }
  }
  else if(type==2)
  {
    for(i in 1:nbPoly)
    {
      # ordinary weighted mean
      vectMean[i]=sum(data[[1]][listZonePoint[[i]]]*surfVoronoi[listZonePoint[[i]]]) /sum(surfVoronoi[listZonePoint[[i]]])
    }

  }

  return(vectMean)
}


##########################################################
#' voronoiPolygons
#'
#' @details determines the Voronoi neighborhood of data points
#' @param spdata SpatialPointsDataFrame
#' @param gridLim list of boundary coordinates 
#' @param neighBool empty point neighborhood Logical matrix 
#' @param PTJUNCTION logical value, if FALSE (default): pts are not neighbors if their Voronoi polygons only have a vertex in common
#' @param FULL logical value, if FALSE (default): do not return Vornoi polygons
#'
#' @return a list with components
#' \describe{
#' \item{surfVoronoi}{Voronoi polygons areas}
#' \item{neighBool}{Voronoi point neighborhood Logical matrix}
#' if FULL=TRUE (warning: uses a lot of memory space), also:
#' \item{voronoi}{Voronoi polygons}
#' }
#' @importFrom deldir deldir tile.list
#'
#' @seealso http://www.carsonfarmer.com/2009/09/voronoi-polygons-with-r/
#' @export
#'
#' @examples
#' data(mapTest)
#' rx=range(mapTest$krigData$x)
#' ry=range(mapTest$krigData$y)
#' nx=nrow(mapTest$krigGrid)
#' ny=ncol(mapTest$krigGrid)
#' nB=matrix(logical((nx*ny)^2),nx*ny,nx*ny) # big matrix
#' vP=voronoiPolygons(mapTest$krigData,c(rx,ry),nB)
#' length(vP$surfVoronoi) #as many as kriged data points
#' # not run
voronoiPolygons = function(spdata,gridLim=c(0,1,0,1),neighBool,PTJUNCTION=FALSE,FULL=FALSE)
##########################################################
{
#source: http://www.carsonfarmer.com/2009/09/voronoi-polygons-with-r/
#input: spatial pts and neighborhood boolean matrix with all elements = FALSE
#output: updated neighborhood boolean matrix
# based on Delaunay tesselation
# PTJUNCTION=FALSE (default): pts are not neighbors if their Voronoi polygons only have a vertex in common

  #get coordinates
  coord = spdata@coords
  #triangulation de delaunay,dans le cadre défini par rw=c(xmin xmax ymin ymax)
  z = deldir(coord[,1], coord[,2],rw=gridLim)#
  #plot(z)

  # Voronoi polygons
  w = tile.list(z)
  polysp = vector(mode='list', length=length(w))
  #for each polygon
  for (i in seq(along=polysp)) {
    #store Voronoi polygons ( 1 per point)
    polycoord = cbind(w[[i]]$x, w[[i]]$y)
    polycoord = rbind(polycoord, polycoord[1,]) #close Voronoi polygon
    polys = Polygons(list(Polygon(polycoord)), ID=as.character(i))
    polysp[[i]] = SpatialPolygons(list(polys))
  }

  #surfVoronoi=sapply(polysp,function(x) slot(x, 'area'))
  surfVoronoi=sapply(polysp,gArea)

  for (k in (1:nrow(z$delsgs)))
  {
# update pt neighborhood matrix
# except if neighbors intersect only by one single point
# and PTJUNCTION argument is FALSE
  	  numpt1=z$delsgs$ind1[k]
 	  numpt2=z$delsgs$ind2[k]
	  if(!PTJUNCTION)
	  {
	  inter=gIntersection(polysp[[numpt1]],polysp[[numpt2]])
	  if(!is.null(inter))
		{
		if(class(inter)=="SpatialLines")
			{
			neighBool[numpt1,numpt2]=TRUE
    	  		neighBool[numpt2,numpt1]=TRUE
	  		}
	  	}
	    }
	    else
	    {
	    neighBool[numpt1,numpt2]=TRUE
    	    neighBool[numpt2,numpt1]=TRUE
	    }
}#end for
  if(FULL)
	return(list(surfVoronoi=surfVoronoi,voronoi=polysp,neighBool=neighBool))
  else
	return(list(surfVoronoi=surfVoronoi,neighBool=neighBool))
}




################################################################################
#' calZoneN
#'
#' @details calculate zone neighborhood
#' @param ptN pt neighborhood Logical matrix
#' @param zoneN empty zone neighborhood Logical matrix
#' @param listZonePoint list of indices of data points within zones
#'
#' @return a list with component zoneN holding filled zone neighborhood Logical matrix
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' ptN=mapTest$krigN
#' nZ=length(K$zonePolygone)
#' zoneN=matrix(logical(nZ*nZ),nZ,nZ)
#' listZonePoint=K$listZonePoint
#' calZoneN(ptN,zoneN,listZonePoint)
#' # not run
calZoneN=function(ptN,zoneN,listZonePoint)
################################################################################
{
#initially all elements of zoneN matrix=FALSE
#ptN=list- element k = neigbours of pt k
#value: boolean matrix of zone neighborhood
nbPoly=length(listZonePoint)
if (nbPoly<=1) return(list(zoneN=zoneN))
# at least 2 zones
	for(i in (1:(nbPoly-1)))
	{
		li=listZonePoint[[i]] # pts in zone i
		u=unique(unlist(ptN[li])) # their neigbours
		for(j in ((i+1):nbPoly))
		      {
		      # zone i has at least 1 pt that is a Voronoi neighbor of pt in zone j
		      lj=listZonePoint[[j]]
		      m=any(match(u,lj,nomatch=0))
       		      zoneN[i,j]=zoneN[j,i] = m

			}
	}
  diag(zoneN)=TRUE

  return(list(zoneN=zoneN))
}
