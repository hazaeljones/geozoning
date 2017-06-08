####################################################################
#' detection of narrow zones (ratio area/perimeter^2)
#'
#' @details description, a paragraph
#' @param zonePolygone xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param type xxxx
#' @param listZonePoint xxxx
#' @param surfVoronoi xxxx
#' @param data xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param spdata xxxx
#' @param neighBool xxxx
#' @param PTJUNCTION xxxx
#' @param FULL xxxx
#'
#' @return a ?
#' @importFrom deldir deldir tile.list
#'
#' @seealso http://www.carsonfarmer.com/2009/09/voronoi-polygons-with-r/
#' @export
#'
#' @examples
#' # not run
voronoiPolygons = function(spdata,gridLim=c(0,1,0,1),neighBool,PTJUNCTION=FALSE,FULL=FALSE)
##########################################################
{
#source: http://www.carsonfarmer.com/2009/09/voronoi-polygons-with-r/
#input: spatial pts and neighborhood boolean matrix with all elements = FALSE
#output: updated neighborhood boolean matrix
# based on Delaunay tesselation
# PTJUNCTION=FALSE (default): pts are not neighbors if their Voronoi polygons only have a vertex in common

  #récupération des coordonnées
  coord = spdata@coords
  #triangulation de delaunay,dans le cadre défini par rw=c(xmin xmax ymin ymax)
  z = deldir(coord[,1], coord[,2],rw=gridLim)#
  #plot(z)

  #polygones de voronoi
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
#' @details description, a paragraph
#' @param ptN xxxx
#' @param zoneN xxxx
#' @param listZonePoint xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
