
##################################################################
#' interZoneC
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param iC xxxx
#' @param iZC xxxx
#' @param closePt xxxx
#'
#' @return a ?
#' @importFrom rgeos gBoundary gCrosses gEnvelope gIntersection gLength
#'
#' @export
#'
#' @examples
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
##################################################################
#' getZoneId
#'
#' @details description, a paragraph
#' @param zone xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param id xxxx
#' @param Z xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param iZ xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
getId=function(Z,iZ)
##################################################################
{
	id = Z[[iZ]]@polygons[[1]]@ID
	return(id)
}
##################################################################
#' findNumZ
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param id xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
findNumZ=function(Z,id)
##################################################################
{
	ids=getIds(Z)
	num=which(ids==id)
	return(num[1])
}

##################################################################
#' getIds
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param nums xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param iZ xxxx
#' @param id xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param Z xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param Z xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' getSurf
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param iZ xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @param Z list of SpatialPloygons
#' @param minSize alarm size threshold
#'
#' @return a vector of small zone indices
#'
#' @export
#'
#' @examples

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

############################
#' getNumZone
#'
#' @details description, a paragraph
#' @param ptsp xxxx
#' @param Z xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
getNumZone=function(ptsp,Z)
############################
{
# ptssp SpatialPointsDataframe
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
#' printZid
#'
#' @details description, a paragraph
#' @param Z xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' crComment
#'
#' @details description, a paragraph
#' @param Z xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param i1 xxxx
#' @param i2 xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
testInterSpe=function(Z,i1,i2)
##################################################################
# returns TRUE if zones Z[[i1]] and Z[[i2]] share some common part
# but are not within each other
{
#if(gOverlaps(Z[[i1]],Z[[i2]]))
if(gCrosses(Z[[i1]],Z[[i2]])) #problematic sometimes
	return(TRUE)
else
	return(FALSE)
}

##################################################################
#' testInterSpeZ1
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param i1 xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
testInterSpeZ1=function(Z,i1)
##################################################################
#returns TRUE if zone Z[[i1]] intersects with any other zone  not within it
{
	inter=FALSE
  	jj=1
  	while(!inter && (jj<=length(Z)))
  	{
		inter=testInterSpe(Z,i1,jj)
		jj=jj+1
  		}

	return(inter)
}

##################################################################
#' testInterSpeZ
#'
#' @details description, a paragraph
#' @param Z xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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
#' findNptInZone
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param K xxxx
#' @param i1 xxxx
#' @param i2 xxxx
#' @param map xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
findNptInZone=function(Z,K,i1,i2,map)
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
if(length(n2)>0) mat=rbind(mat,c(ptnum,n2))
}
colnames(mat)=c("ptZ1","ptNZ2")
return(mat)
}
##################################################################
#' getNs
#'
#' @details get zone numbers of neighbors of a given zone
#' @param zoneNModif zone neighborhood Logical matrix  
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
getNs=function(zoneNModif,iZ)
{
  Ns=zoneNModif[iZ,]
  return(Ns)
}
