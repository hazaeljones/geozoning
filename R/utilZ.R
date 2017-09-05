##################################################################
#' getZoneId
#'
#' @details get the zone unique identifier
#' @param zone SpatialPolygons
#'
#' @return the zone identifier ( a character vector of length 1)
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
#' Z=criti$zk[[2]][[1]]$zonePolygon
#' getZoneId(Z[[4]])
#' # not run
getZoneId=function(zone)
##################################################################
{
	id=zone@polygons[[1]]@ID
	return(id)
}

###########################################################################
#' Identify
#'
#' @details get the number of a zone with a given identifier in a zoning
#' this is necessary because correction procedures may remove zones from initial#' zoning. Therefore zone numbers change, but identifiers are conserved.
#' @param id zone identifier (character vector)
#' @param Z zoning geometry (list of SpatialPolygons) 
#'
#' @return the zone number
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
#' Z=criti$zk[[2]][[1]]$zonePolygon
#' Identify(6,Z)
#' # not run
Identify=function(id,Z)
###########################################################################
{
#returns the zone number corresponding to the id given in argument
#returns 0 if no id
  for (i in (1:length(Z)))
  {
    if (Z[[i]]@polygons[[1]]@ID == id) return(i)
  }
  return(0)
}

##################################################################
#' getId
#'
#' @details get zone identifier in a zoning
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iZ zone number
#'
#' @return a character vector giving the zone identifier
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
#' Z=criti$zk[[2]][[1]]$zonePolygon
#' getId(Z,6)
#' # not run
getId=function(Z,iZ)
##################################################################
{
	id = Z[[iZ]]@polygons[[1]]@ID
	return(id)
}

##################################################################
#' getIds
#'
#' @details get zone identifiers in a zoning
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param nums zone numbers
#'
#' @return a character vector giving the zone identifiers
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
#' Z=criti$zk[[2]][[1]]$zonePolygon
#' getIds(Z)
#' # not run
getIds=function(Z,nums=NULL)
##################################################################
{
	ids=c()
	if (is.null(nums))
	   iZZ = 1:length(Z)
	else
	   iZZ = nums
	for (iZ in iZZ)
	{
		ids = c(ids,Z[[iZ]]@polygons[[1]]@ID)
	}
	return(ids)
}


##################################################################
#' setId
#'
#' @details assign zone identifier in a zoning
#' @param Z zoning geometry (list of SpatialPolygons) 
#' @param iZ zone number
#' @param id zone identifier to assign
#'
#' @return a zoning geometry 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
#' Z=criti$zk[[2]][[1]]$zonePolygon
#' Z1=setId(Z,4,"4")
#' # not run
setId=function(Z,iZ,id)
##################################################################
{
	Z[[iZ]]@polygons[[1]]@ID = as.character(id)
	return(Z)
}

##################################################################
#' setIds
#'
#' @details set all zone identifiers in a zoning by assigning zone number to each identifier.
#' @param Z zoning geometry (list of SpatialPolygons) 
#'
#' @return a zoning geometry 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
#' Z=criti$zk[[2]][[1]]$zonePolygon
#' Z1=setIds(Z)
#' # not run
setIds=function(Z)
##################################################################
{
	for (iZ in 1:length(Z))
	{
		Z[[iZ]]@polygons[[1]]@ID = as.character(iZ)
	}
	return(Z)
}

##################################################################
#' maxId
#'
#' @details get highest number corresponding to a zone identifier in a zoning
#' @param Z zoning geometry (list of SpatialPolygons) 
#'
#' @return a number
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
#' Z=criti$zk[[2]][[1]]$zonePolygon
#' maxId(Z)
#' # not run
maxId=function(Z)
##################################################################
{
	m=as.numeric(getId(Z,1))
	if (length(Z)<2) return(m)
	for (iZ in 2:length(Z))
	{
		id=as.numeric(getId(Z,iZ))
		if (id >m) m=id
	}
	return(as.numeric(m))
}

##################################################################
#' printZid
#'
#' @details print zone identifiers in a zoning
#' @param Z zoning geometry (list of SpatialPolygons)
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' printZid(Z)
#' # not run
printZid = function(Z)
##################################################################
{
	for (ii in 1:length(Z))
	{
		print(paste("ii=",ii," ID=", Z[[ii]]@polygons[[1]]@ID))
	}
}

##################################################################
#' getSurf
#'
#' @details description, a paragraph
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iZ zone number
#'
#' @return zone area
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' getSurf(Z,1)
#' # not run
getSurf=function(Z,iZ)
##################################################################
{
	return(Z[[iZ]]@polygons[[1]]@Polygons[[1]]@area)
}

##################################################################
#' printZsurf
#'
#' @details print zone surfaces
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param minSize minimum size threshold
#'
#' @return a vector of small zone indices
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' printZsurf(Z,0.03)
#' # not run
printZsurf=function(Z,minSize=0.012)
##################################################################
{
	smallZ=c()
	gas=c()
	for ( iZ in 1:length(Z))
	{
	ga=gArea(Z[[iZ]])
	if(ga>=minSize)
		print(paste("iZ=",iZ,"area=",round(gArea(Z[[iZ]]),5)))
	else
		{
		print(paste("iZ=",iZ," area=",round(gArea(Z[[iZ]]),5)," < minSize(",minSize,")",sep=""))
		smallZ=c(smallZ,iZ)
		gas=c(gas,ga)
		}
	}
if (!is.null(gas))
   return(smallZ[order(gas)])
   else
   return(NULL)
}

##################################################################
#' getZonePts
#'
#' @details get all data points within a zone
#' @param ptsp SpatialPointsDataFrame with data values
#' @param zone SpatialPolygons defining a zone
#'
#' @return a list with components
#' \describe{
#' \item{pts}{SpatialPointsDataFrame with the data points within the zone}
#' \item{mask}{Logical vector of the within zone data points indices}
#' }
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' data(mapTest)
#' ptsp=mapTest$krigData
#' res=getZonePts(ptsp,Z[[5]])
#' plotZ(Z)
#' points(res$pts,col="blue",pch=20)
#' # not run
getZonePts=function(ptsp,zone)
##################################################################
{
# ptsp=map$krigData
#
	# pts in zone
	IN=numeric()
	poly=getPolySp(zone,1)
	IN=point.in.polygon(ptsp$x,ptsp$y,poly@coords[,1],poly@coords[,2])
	# remove pts in holes
	if(nPolySp(zone)>1)
		{
      		for(k in 2:nPolySp(zone))
      		      {
		      poly=getPolySp(zone,k)
        	      IN=IN-point.in.polygon(ptsp$x,ptsp$y,poly@coords[,1],poly@coords[,2])
     		      }
  		 }
	IN=as.logical(IN)
	ptsx=ptsp$x[IN]
	ptsy=ptsp$y[IN]
	ptsd=ptsp[[1]][IN]
	ptsub=SpatialPointsDataFrame(coords=cbind(ptsx,ptsy),data=data.frame(z=ptsd))
	return(list(pts=ptsub,mask=IN))
}


############################
#' getNumZone
#'
#' @details get zone numbers to which each point in a SpatialPointsDataFrame belongs
#' @param ptsp SpatialPointsDataFrame
#' @param Z zoning geometry (list of SpatialPolygons)
#'
#' @return the zone number
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' getNumZone(mapTest$krigData,Z)
#' # not run
getNumZone=function(ptsp,Z)
############################
{
# ptsp SpatialPointsDataframe
# Z zoning (list of SpatialPolygons)
# ptsp=map$krigData
#
numzone=rep(0,length(ptsp$x))
for (k in 1:length(Z))
    {
	res=getZonePts(ptsp,Z[[k]])
	mask=res$mask
	numzone[mask]=k
    }
return(numzone)
}

##################################################################
#' getNs
#'
#' @details get zone numbers of neighbors of a given zone
#' @param zoneN zone neighborhood Logical matrix  
#' @param iZ index of current zone in zoning
#'
#' @return a Logical vector of current zone neighbors
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' K=resZTest
#' Ns=getNs(K$zoneNModif,5) # find neighbors of zone 5
##################################################################
getNs=function(zoneN,iZ)
{
  Ns=zoneN[iZ,]
  return(Ns)
}


##################################################################
#' findNptInZone
#'
#' @details find, in a given zone, neighbor points of points belonging to another zone
#' @param K zoning object, as returned by the calNei function
#' @param i1 first zone
#' @param i2 second zone, where to search for neighbors of points in first zone
#' @param map object returned by function genMap
#'
#' @return a two-column matrix, the first column contains indices of pts in first zone which have at least one neighbor in second zone, the second column contains the neighbor indices.
#
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' # not run
findNptInZone=function(K,i1,i2,map)
##################################################################
{
ptZ1=K$listZonePoint[[i1]]
ptZ2=K$listZonePoint[[i2]]
ptN=map$krigN
#
# for each pt in Z1
mat=c()
for (ptnum in ptZ1)
{
# look for neighbor pt in Z2
nei=map$krigN[[ptnum]]
k=match(nei,ptZ2)
n2=ptZ2[k[!is.na(k)]]
if(length(n2)>0)
	{
	for (n2k in n2)
	{
	mat=rbind(mat,c(ptnum,n2k))
	}
	}
}
colnames(mat)=c("ptZ1","ptNZ2")
return(mat)
}

##################################################################
#' crComment
#'
#' @details create comment corresponding to holes in a zoning
#' @param Z zoning geometry (list of SpatialPolygons)
#'
#' @return a zoning
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' Z1=crComment(Z)
#' # not run
crComment = function(Z)
##################################################################
{
# create comments for holes
      if (!is.null(Z) & length(Z) >0)
      {
	for (iZ in 1:length(Z))
	    {
		comment(Z[[iZ]]@polygons[[1]])=createPolygonsComment(Z[[iZ]]@polygons[[1]])
	     }

      }
      return(Z)
}

##################################################################
#' testInterSpe
#'
#' @details checks if 2 zones in a zoning share somme common part (using gOverlaps)
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param i1 first zone number
#' @param i2 second zone number
#'
#' @return a Logical value, TRUE if there is an intersection, FALSE if not.
#' @importFrom rgeos gOverlaps

#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.2,0.5)
#' ZK = initialZoning(qProb, mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' plotZ(Z)
#' Z58=rgeos::gConvexHull(rgeos::gUnion(Z[[8]],Z[[5]]))
#' Z[[length(Z)+1]]=Z58 # add new zone to zoning
#' plotZ(Z)
#' testInterSpe(Z,6,length(Z))
#' # not run
testInterSpe=function(Z,i1,i2)
##################################################################
# returns TRUE if zones Z[[i1]] and Z[[i2]] share some common part
# but are not within each other
{
#if(gCrosses(Z[[i1]],Z[[i2]])) #problematic sometimes
if(gOverlaps(Z[[i1]],Z[[i2]]))
	return(TRUE)
else
	return(FALSE)
}

##################################################################
#' testInterSpeZ1
#'
#' @details checks, within a zoning, if a given zone intersects with any other zone  not within it
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iZ zone number
#'
#' @return a Logical value, TRUE if there is an intersection, FALSE if not.
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.2,0.5)
#' ZK = initialZoning(qProb, mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' plotZ(Z)
#' Z58=rgeos::gConvexHull(rgeos::gUnion(Z[[8]],Z[[5]]))
#' Z[[length(Z)+1]]=Z58 # add new zone to zoning
#' plotZ(Z)
#' testInterSpe(Z,6,length(Z))
#' # not run
testInterSpeZ1=function(Z,iZ)
##################################################################
#returns TRUE if zone Z[[iZ]] intersects with any other zone  not within it
{
	inter=FALSE
  	jj=1
  	while(!inter && (jj<=length(Z)))
  	{
		inter=testInterSpe(Z,iZ,jj)
		jj=jj+1
  		}

	return(inter)
}

##################################################################
#' testInterSpeZ
#'
#' @details checks, within a zoning, if any zone intersects with any other zone  not within it and not englobing it
#' @param Z zoning geometry (list of SpatialPolygons)
#'
#' @return a Logical value, TRUE if there is any intersection, FALSE if not
#'
#' @export
#'
#' @examples
#' qProb=c(0.2,0.5)
#' ZK = initialZoning(qProb, mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' plotZ(Z)
#' Z58=rgeos::gConvexHull(rgeos::gUnion(Z[[8]],Z[[5]]))
#' Z[[length(Z)+1]]=Z58 # add new zone to zoning
#' plotZ(Z)
#' testInterSpeZ(Z)
#' # not run
testInterSpeZ=function(Z)
##################################################################
#returns TRUE if any zone intersects with any other zone  not within it
{
	inter=FALSE
	i1=1
	while(!inter && (i1<=length(Z)))
	{
		inter=testInterSpeZ1(Z,i1)
		i1=i1+1
	}
	return(inter)
}

##################################################################
#' printInterZ
#'
#' @details checks intersection of sp and each element of Z
#' @param Z list of zones, each zone is a SpatialPolygons
#' @param sp SpatialPolygons object
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z,mapTest)
#' printInterZ(Z,Z[[1]])
#' # not run
printInterZ = function(Z,sp)
##################################################################
{
	le=length(Z)
	if (le == 0) return()
	for (i in (1:le))
	{
		print(paste("i=",i,"intersec=",gIntersects(Z[[i]],sp)))
	}
	
}

##################################################################
#' ptInZone
#'
#' @details description, a paragraph
#' @param zone SpatialPolygons
#' @param pts data points
#' @param numpt data point indices
#'
#' @return 1 if point is within zone, 0 if not
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' ptInZone(Z[[1]],mapTest$krigData,c(5,500))
#' # not run
ptInZone=function(zone,pts,numpt)
##################################################################
{
	spp=SpatialPoints(pts)
	logicalPoly=point.in.polygon(spp$x,spp$y,zone@polygons[[1]]@Polygons[[1]]@coords[,1],zone@polygons[[1]]@Polygons[[1]]@coords[,2])
    if(length(zone@polygons[[1]]@Polygons)>1)
    {
      #consider pts may be inside holes
      for(k in 2:length(zone@polygons[[1]]@Polygons))
      {
        logicalPoly=logicalPoly-point.in.polygon(spp$x,spp$y,zone@polygons[[1]]@Polygons[[k]]@coords[,1],zone@polygons[[1]]@Polygons[[k]]@coords[,2])
      }
    }
    return(logicalPoly[numpt])
}

##################################################################
#' printLabZ
#'
#' @details print zoning labels for a list of zoning objects
#' @param Klist list of zoning objects, typically result of a call to correctionTree
#'
#' @return a list of zoning objects
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.1,0.2,0.4);criti=correctionTree(qProb,mapTest) # 2 zonings at last level
#' printLabZ(criti$zk[[2]])
#' # not run
printLabZ=function(Klist)
##################################################################
{
#Klist is a list of several Ks
	lk=length(Klist)
	labZ=list()
	for (k in 1:lk)
	{
		labk=unlist(Klist[[k]]$lab)
		labq=paste(length(unique(labk))-1,"q",sep="")
		labZ[[k]]=labk
		print(paste(labq,"zone labels=",labZ[k]))
	}
	return(labZ)
}

##################################################################
#' trLabZone
#'
#' @details transfer zone labels from K1 to K2
#' @param K1 zoning object (such as returned by calNei function)
#' @param K2 zoning object (such as returned by calNei function)
#' @param map object returned by genMap function
#' @param qProb probability vector used to generate quantile values
#' @param disp 0: no info, 1: detailed info
#'
#' @return a zoning object
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Ns=getNs(K$zoneNModif,5) # find neighbors of zone 5
#' zf=zoneFusion3(K,5,Ns,mapTest,disp=0) # merge zone 5 with englobing one
#' K2=calNei(zf,mapTest$krigData,mapTest$krigSurfVoronoi,mapTest$krigN)
#' K2=trLabZone(K,K2,mapTest,K$qProb)
#' # not run
trLabZone=function(K1,K2,map,qProb,disp=0)
##################################################################
# transfer zone labels from K1 to K2
# ids of zones in Z1 and Z2 used to guide transfer

{
Z1=K1$zonePolygone
Z2=K2$zonePolygone

# only K2$lab is modified
# by taking labels from K1
lab2=rep(1,length(Z2))

lab1=K1$lab
id1s=getIds(Z1)
q1 = quantile(map$krigGrid,na.rm=TRUE,prob=qProb)
rate= max(map$krigGrid)-min(map$krigGrid)
for (iZ in 1:length(Z2))
{
	id2=getId(Z2,iZ)
	numid1=which(id1s==id2)
	if (length(numid1) >1)
	{
		print(id1s)
		print(id2)
	}
	if (length(numid1)!=0)
	   lab2[iZ]=lab1[numid1]
	else #new zone id was created - find its label
	{
		for (j in 1:length(q1))
    		{
      		if (K2$meanZone[iZ]>(q1[j]+0.01*rate))
      		   {
		   lab2[iZ]=j+1
      		   }

		}
	}

}
K2$lab=lab2

if(disp)
	{
	printZid(Z1)
	print(K1$lab)
	print(K2$lab)
	}
return(K2)
}

##################################################################
#' getClosestZone
#'
#' @details get closest non neighbor zone (i.e. excluding neighbor zones and englobing zone)
#' @param iZ current zone number
#' @param Z current zone
#' @param zoneN zone neighborhood Logical matrix
#'
#' @return the closest zone index
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' getClosestZone(4,Z,resZTest$zoneNModif)
#' # not run
getClosestZone=function(iZ,Z,zoneN)
##################################################################
{
	imin=0
	ni=1:length(Z)
	# exclude current zone
	ni=ni[ni != iZ]

	# exclude neighbor zones
	Ns=getNs(zoneN,iZ)
  	listeV=grep(TRUE, Ns)
	for (i in listeV) ni=ni[ni != i]

	# exclude englobing zone
	iE = detZoneEng(iZ,Z,zoneN)
	#
	if(iE != 0) ni = ni[ni !=iE]

	# exclude included zones
	ir=NULL
	for (i in ni)
          {
	  gb = gBuffer(gConvexHull(Z[[iZ]]),byid=TRUE,width=1e-3)
	  if(gContains(gb,Z[[i]])) ir=c(ir,i)
	  }
	for (i in ir) ni =ni[ni != i]

	# compute distances to remaining zones
	d0 = 1
	for (i in ni)
          {
		d=gDistance(Z[[iZ]],Z[[i]])
		if (d<=d0)
		   {
		   imin=i
		   d0=d
		   }
	}

return(imin)

}

#####################################################################
#' MeanVarWPts
#'
#' @details computes (weighted) mean and variance of zone data
#' @param map object returned by function genMap
#' @param zone SpatialPolygons defining a zone
#' @param w weighting vector (default NULL)
#'
#' @return a list with components
#' \describe{
#' \item{mean}{(weighted) mean of the within zone attribute value}
#' \item{var}{(weighted) variance of the within zone attribute value}
#' }
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' data(mapTest)
#' MeanVarWPts(mapTest,Z[[1]])
#' # Weights are areas of the Voronoi polygons corresponding to data points
#' MeanVarWPts(mapTest,Z[[1]],mapTest$krigSurfVoronoi) #slightly different result
#' # not run
MeanVarWPts=function(map,zone,w=NULL)
#####################################################################
{
  ptsp=map$krigData

  m=numeric()
  v=numeric()
  res=getZonePts(ptsp,zone)
  mask=res$mask

  if (is.null(w))
     w=map$krigSurfVoronoi
  else
     w=rep(1,length(mask))

  m=sum(ptsp[[1]][mask]*w[mask]) /sum(w[mask])
  d=ptsp[[1]][mask]-m
  v=sum(d*d*w[mask]) /sum(w[mask])
  #plot(zone)
  #points(ptsp$x[mask],ptsp$y[mask])

  return(list(mean=m,var=v))
}

##################################################################
#' nPolyZone
#'
#' @details number of polygons in current zone
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iC current zone number within Z
#'
#' @return the number of polygons within the current zone
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' iC=1;print(paste(nPolyZone(Z,iC),"polygons in zone",iC))
#' # not run
nPolyZone=function(Z,iC)
##################################################################
{
	return(length(Z[[iC]]@polygons[[1]]@Polygons))
}

##################################################################
#' maxDistZone
#'
#' @details maximum distance within kth polygon of current zone 
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iZ current zone index
#' @param k polygon number within current zone
#'
#' @return the maximum distance within kth polygon of the current zone
#'
#' @export
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' maxDistZone(Z,5,1)
#' # not run
maxDistZone=function(Z,iZ,k)
##################################################################
{
	return(max(dist(Z[[iZ]]@polygons[[1]]@Polygons[[k]]@coords)))
}

##################################################################
#' getPoly
#'
#' @details get the kth polygon of the current zone in zoning Z
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iZ current zone index
#' @param k polygon number within current zone
#'
#' @return a polygon (object of class Polygon)
#'
#' @export
#'
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' P1=getPoly(Z,5,1)
#' P2=getPoly(Z,5,2) # second polygon is a hole
#' plot(P1@coords,type="l")
#' lines(P2@coords,type="l",col="blue")
#' # not run
getPoly = function(Z,iZ,k)
##################################################################
{
	if (k == 0) return(Z[[iZ]])
	# k can be <0 or >0
	p=Z[[iZ]]@polygons[[1]]@Polygons[[k]]
	return(p)
}

##################################################################
#' polyToSp
#'
#' @details transforms kth polygon of zone into SpatialPolygons
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iZ zone number
#' @param k polygon number
#'
#' @return a SpatialPolygons
#'
#' @export
#'
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' sph=polyToSp(Z,5,2)
#' plotZ(Z)
#' sp::plot(sph,type="l",col="blue",add=TRUE)
#' # not run
polyToSp=function(Z,iZ,k)
##################################################################
{
	if (k == 0) return(Z[[iZ]])
	# k can be <0 or >0
	p=Z[[iZ]]@polygons[[1]]@Polygons[[k]]
	sp=SpatialPolygons(list(Polygons(list(p),1:1)))

	return(sp)
}

##################################################################
#' calcDCrit
#'
#' @details computes distances and criterion value for zoning Z
#' @param Z zoning geometry (list pf SpatialPolygons)
#' @param map object returned by function genMap
#' @param optiCrit criterion choice
#' @param pErr equality tolerance for distance calculations, default 0.9
#' @param simplitol tolerance for spatial polygons geometry simplification, default 0.001
#'
#' @return a list with components
#'\describe{
#' \item{resD}{list with uncorrected and corrected distance matrix}
#' \item{resCrit}{list with criterion and cost values}
#' }
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' Z1=zoneFusion4(Z,6,2)
#' calcDCrit(Z1,mapTest)
#' # not run
calcDCrit=function(Z,map,optiCrit=2,pErr=0.9,simplitol=1e-3)
##################################################################
{
	resZ=calNei(Z,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol=simplitol)
  le = length(resZ$zonePolygone)
	if (le <2) return(list(resD=0,resCrit=0))

  resDist=calDistance(typedist=1,map$krigData,resZ$listZonePoint,resZ$zoneN,map$krigSurfVoronoi,resZ$meanZone,pErr=pErr)
  resCrit=calCrit(resDist$matDistanceCorr,resZ$zoneNModif,optiCrit)

	return(list(resD=resDist,resCrit=resCrit))
}

##################################################################
#' normZcoords
#'
#' @details description, a paragraph
#' @param Z list of SpatialPolygons
#' @param boundary list with components x and y, used to normalize polygons in zoning
#'
#' @return a list with components
#' \describe{
#' \item{Zn}{list of normalized SpatialPolygons}
#' \item{boundaryn}{normalized boundary}
#' }
#' @export
#'
#' @examples
#  Import shape1 object (was read from a shapefile)
#' shape1 = geozoning::shape1
#' p = shape1@polygons
#' P=sp::SpatialPolygons(p) #SpatialPolygons
#' Z1=list()
#' for (kk in 1:length(P)){Z1[[kk]]=P[kk]} # transform into list of SpatialPolygons
#' bd=list(x=c(7723131,7723132,7723294,7723295,7723131),y=c(3576432,3576814,3576809,3576436,3576432))
#  Z2=normZcoords(Z1,bd)
#' # not run
normZcoords=function(Z,boundary)
##################################################################
{
	NZ=length(Z)
	Z1=list()
	for (iZ in 1:NZ)
	{
	np=nPolyZone(Z,iZ)
	pnl=list()
	for (k in 1:np)
	{
		pk = getPoly(Z,iZ,k)
		resn = spnorm(pk,boundary)
		if (is.null(resn)) return(NULL)
		pnl[[k]] = resn$pn
		boundaryn = resn$boundaryn
	}

		Z1[[iZ]] = SpatialPolygons(list(Polygons(pnl,1:1)))

	}

	return(list(Zn=Z1,boundaryn=boundaryn))
}

##################################################################
#' createHoles
#'
#' @details description, a paragraph
#' @param Z list of zones, each zone is a SpatialPolygons
#'
#' @return a list of zones where holes are distinct SpatialPolygons
#'
#' @export
#'
#' @examples
#' # not run
createHoles = function(Z)
##################################################################
{
	NZ = length(Z)
	hole = matrix(NA,nrow=NZ,ncol=NZ)
	Z1 = Z
	for (iZ in 1:NZ)
	{
		for (kZ in 1:NZ)
		{
			if (iZ == kZ) next
			h=gContains(Z[[iZ]],Z[[kZ]])
			hole[iZ,kZ]=h
			if(h)
			{
				#kz is within iz, create corresp. hole in iz
				Z1[[iZ]]=gDifference(Z[[iZ]],Z[[kZ]])
			}
		} #end for kZ
	}
	#end for iZ
	print(hole)
	return(Z1)
}

##################################################################
#' moveHoles
#'
#' @details creates SpatialPolygons excluding holes
#' @param zoneMain SpatialPolygons
#' @param zoneSuppr SpatialPolygons  inside main 
#'
#' @return a new SpatialPolygons object
#'
#' @export
#'
#' @examples
#' # not run
moveHoles = function(zoneMain,zoneSuppr)
##################################################################
{
	Zone1=zoneMain

	le=nPolySp(zoneSuppr)
	for ( i in 1:le)
	{
		poly=zoneSuppr@polygons[[1]]@Polygons[[i]]
		polys=Polygons(list(poly),"id")
		polySp = SpatialPolygons(list(polys))
		if(poly@hole) Zone1=gDifference(Zone1,polySp)
	}

	return(Zone1)
}

##################################################################
#' find contour for a given quantile value, within an envelope and englobing current zone 
#'
#' @details withing a zoning, find contour for a given vRef quantile value, contour contains current zone and is included in envel (spatial Polygon)
#' @param iC zone number
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param K zoning object (such as returned by calNei function)
#' @param map object returned by genMap function
#' @param vRef quantile value
#' @param envel SpatialPolygons within which the contour must be contained
#'
#' @return a list with components
#'\describe{
#' \item{area}{area of SpatialPolygons corresponding to contour}
#' \item{contourSp}{SpatialPolygons corresponding to contour}
#' }
#' @importFrom sp plot
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.3,0.5)
#' criti = correctionTree(qProb,mapTest)
#' best = criti$zk[[2]][[8]]
#' Z=best$zonePolygone
#' plotZ(Z)
#' iC=4
#' envel=calFrame(iC,Z,best$zoneNModif)
#' sp::plot(envel,col="blue",add=TRUE)
#' vRef=quantile(mapTest$krigGrid,0.6)
#' resp=findCinZ(iC,Z,best,mapTest,vRef,envel)
#' sp::plot(resp$contourSp,col="red",add=TRUE)
#' # not run
findCinZ = function(iC,Z,K,map,vRef,envel)
##################################################################
{
  # find contour for a given vRef quantile value
  # contour contains Z[[ic]]
  # contour included in envel (spatial Polygon)
  # init
  area=0
  #
  listContour=list()
  listContour = contourAuto(listContour,map$step,map$xsize,map$ysize,map$krigGrid,vRef,map$boundary)

  # intersect contour with boundary
  # and transform contours into sps
  sp=list()
  k=0
  for (cont in listContour)
  {
     k=k+1
     ps = interCB(cont,map$step,envel=envel)
     # returns NULL is contour is degenerate (single point)
     if(!is.null(ps)) sp[[k]]=ps else k=k-1
  }
  if (length(sp)==0) return(NULL)
  
  # check which one contains current zone
  spb = sapply(sp,gBuffer,width=1e-3)
  gc = sapply(spb,gContains,Z[[iC]])
  # test if contour is included in envel
  gci = sapply(spb, gWithin, envel)

  numc=1:length(sp) # number of non degenerated contours
  ga = sapply(sp,gArea)
  ind = numc[gc & gci]
  # ind may be empty - otherwise take the biggest contour area
  if(length(ind)==0) return(NULL)

  imax = which(ga[ind]==max(ga[ind]))
  imax = imax[1]
  im = numc[ind[imax]]
  area = max(ga[ind])
  contourSp = sp[[im]]

  return(list(area=area,contourSp=contourSp))
}

######################
#' modlm
#'
#' @details description, a paragraph
#' @param ptsp SpatialPointsDataFrame with data values
#' @param Z zoning (list of SpatialPolygons)
#'
#' @return the result of a call to lm (anova model with zone number as factor)
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' data(mapTest)
#' ptsp=mapTest$krigData
#' modlm(ptsp,Z)
#' # not run
modlm=function(ptsp,Z)
######################
{
  # ptssp SpatialPointsDataframe
  # Z zoning (list of SpatialPolygons)
  # ptsp=map$krigData

  numz=getNumZone(ptsp,Z)
  # add numz to SpatialPointsDataframe
  ptsp[[2]]=numz
  colnames(ptsp@data)=c("z","numz")
  reslm=lm(z~numz,data=ptsp@data)
  return(reslm)
}

##################################################################
#' getClosePt
#'
#' @details description, a paragraph
#' @param Z zoning (list of SpatialPolygons)
#' @param iC current zone indes
#' @param iZC close zone index
#' @param disp information level (FALSE-no info)
#'
#' @return a SpatialPoints of length 1
#'
#' @export
#' @importFrom sp SpatialPoints SpatialPointsDataFrame
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' getClosePt(Z,1,3)
#' plotZ(Z)
#' points( getClosePt(Z,1,3),col="blue",pch=20)
#' # not run
getClosePt=function(Z,iC,iZC,disp=FALSE)
##################################################################
{
#pt coordinates

  a=SpatialPoints(Z[[iC]]@polygons[[1]]@Polygons[[1]]@coords)
  b=SpatialPoints(Z[[iZC]]@polygons[[1]]@Polygons[[1]]@coords)

  #all distances between pairs of points - lengthy
  Fdist= list()
  for (i in 1:length(a))
  {
    Fdist[[i]]<-gDistance(a[i,],b,byid=TRUE)
  }

  #get the point in zone iZC which is the closest to zone iC
  min.dist <- unlist(lapply(Fdist, FUN=function(x) which(x == min(x))[1]))
  PolyDist <- unlist(lapply(Fdist, FUN=function(x) min(x)[1]))

  pZC = min.dist[which.min(PolyDist)]
  if (disp)
  {
  print(b[pZC])
  }
  return (b[pZC])
}
########################
#' findZCenter
#'
#' @details find point within zone for pretty labelling
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param num zone number
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a SpatialPoints
#'
#' @export
#'
#' @examples
#' data(mapTest)
# run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7,
# saving initial zoning and last level zonings
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE) 
#' Z=criti$zk[[2]][[1]]$zonePolygone
#' findZCenter(Z)
#' # not run
findZCenter=function(Z,num=NULL)
################################
{
if(is.null(num)) num=1:length(Z)
ptz=NULL
for (jj in num)
{
	erosion = Z[[jj]]
    	width = 0
    	while(gArea(erosion)>10^-3){
      	width = width + 0.001
      	erosion = gBuffer(Z[[jj]],width =-width)
	}
ptz=rbind(ptz,coordinates(erosion))
}
rownames(ptz)=NULL
colnames(ptz)=c("x","y")
return(ptz)
}
########################
#' findZCenterpt
#'
#' @details find point within zone for pretty labelling
#' @param data SpatialPointsDataFrame
#' @param K zoning object, as returned by the calNei function
#' @param num zone number or NULL for all zones
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a matrix of x and y coordinates for chosen points with as many rows as zones
#'
#' @export
#'
#' @examples
#' data(mapTest)
# run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7,
# saving initial zoning and last level zonings
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE) 
#' K=criti$zk[[2]][[1]]
#' data=mapTest$krigData
#' findZCenterpt(data,K)
#' # not run
findZCenterpt=function(data,K,num=NULL)
#######################################
{
Z=K$zonePolygone
if(is.null(num)) num=1:length(Z)
ptz=NULL
for (jj in num)
{
ipt=K$listZonePoint[[jj]] # indices of within zone pts
ptc=data[ipt,]
gd=rep(NA,length(ptc))
for (k in 1:length(ptc)) {gd[k]=gDistance(ptc[k,],gBoundary(Z[[jj]]))} #smallest distance from each point to zone boundary
ind=which(gd==max(gd)) # maximum distance to zone boundary
ind=ind[1] # if ties
ptz=rbind(ptz,coordinates(ptc[ind,]))
}
rownames(ptz)=NULL
colnames(ptz)=c("x","y")
return(ptz)
}
########################
getZsize=function(Z)
#' getZsize
#'
#' @details compute maximum x and y values of zoning Z
#' @param Z zoning geometry (list of SpatialPolygons)
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a vector with x and y maximum values
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' getZsize(Z)
#' # not run
########################
{
ms=sapply(Z,function(x)return(x@bbox[,"max"]))
vs=apply(ms,1,max)
return(vs)
}
##################################################################
#' interZoneC
#'
#' @details finds two intersection points of a circle with a zone. The circle radius is chosen so that it will intersect both zones given as arguments.
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iC zone number
#' @param iZC other zone number
#' @param closePt SpatialPoints object in other zone used as circle center
#'
#' @return a list with components
#'\describe{
#' \item{spi}{Two SpatialPoints to be used for the junction of the two zones}
#' \item{ord}{Order in which to use the points}
#' }
#' @importFrom rgeos gBoundary gCrosses gEnvelope gIntersection gLength
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.2,0.5)
#' ZK = initialZoning(qProb, mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' plotZ(K$zonePolygone) # zoning
#' closePt = getClosePt(Z,6,8)
#' points(closePt,col="red")
#' res  = interZoneC(Z,6,8,closePt)
#' points(res$spi,col="red")
#' # not run
interZoneC=function(Z,iC, iZC, closePt)
##################################################################
{
# returns spatialpoints = 2 points for junction of zones iC and iZC
#  iC (current zone) and iZC (close zone)
if (is.null(closePt)) return(NULL) #protection

# coordinates of close zone
b2=gBoundary(Z[[iZC]])
co1 = Z[[iC]]@polygons[[1]]@Polygons[[1]]@coords
co2sp = SpatialPoints(Z[[iZC]]@polygons[[1]]@Polygons[[1]]@coords)
end = TRUE
# draw circle with closePt as center and radius chosen so that it will intersect

  x = closePt$x
  y = closePt$y
  n = 100

  r= 0.75*max (dist(co1))
  pts = seq(0, 2 * pi, length.out = n)
  xy = cbind(x + r * sin(pts), y + r * cos(pts))
  circle= SpatialLines(list(Lines(list(Line(xy)), "line")))
  #lines(circle)

# intersect circle with close zone
  pInter = gIntersection(circle,b2)

  spi = NULL
  ord=NULL
  if (!is.null(pInter) & length(pInter)>=2)
  {

	# find  intersection pts
	i1=which.min(gDistance(pInter[1],co2sp,byid=TRUE))
  	i2=which.min(gDistance(pInter[2],co2sp,byid=TRUE))

  	#transform into Spatial Points
  	spi= co2sp[c(i1,i2),]
	# order of points
	if (i1 <i2){
    	ord = (i1+1):(i2-1)
    	end = FALSE
	}
 	 else
  	 {
    	 ord = (i2+1):(i1-1)
  	 }
	 #
  nbcoord = length(co2sp)
  if (length(ord) >nbcoord/2)
  {
    if (!end)
    {
      ord = c((i2+1):(nbcoord),1:(i1-1))
    }
    else
    {
      ord = c((i1+1):(nbcoord),1:(i2-1))
    }
  }
  }
return(list(spi=spi,ord=ord))
}
