##################################################################
#' getZonePts
#'
#' @details 
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


#####################################################################
#' MeanVarWPts
#'
#' @details description, a paragraph
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
#' Weights are areas of the Voronoi polygons corresponding to data points
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
####################
#' r2
#'
#' @details adjusted R2
#' @param reslm result of a call to lm 
#'
#' @return the adjusted r-square of the lm model
#'
#' @export
#' @importFrom stats anova dist lm quantile sd
#'
#' @examples
#' # not run
r2=function(reslm)
####################
{
  s2T <- sum(anova(reslm)[[2]]) / sum(anova(reslm)[[1]])
  MSE <- anova(reslm)[[3]][2]
  adj.R2 <- (s2T - MSE) / s2T
  return(adj.R2)
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

###################################################################
#' getCoords
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons
#' @param k polygon number
#'
#' @return some coordinates
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' getCoords(Z[[1]])
#' # not run
getCoords=function(sp,k=1)
##################################################################
{
	coords=sp@polygons[[1]]@Polygons[[k]]@coords
	return(coords)
}

##################################################################
#' spToSL
#'
#' @details tranform SpatialPolygons into SpatialLines
#' @param sp SpatialPolygons
#'
#' @return a SpatialLines
#'
#' @export
#' @importFrom sp SpatialLines Lines
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' spToSL(Z[[5]])
#' # not run
spToSL=function(sp)
##################################################################
{
	co = getCoords(sp)
	li = Lines(list(Line(co)),ID="1")
	lis = SpatialLines(list(li))
	return(lis)
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

##################################################################
#' contourArea
#'
#' @details area corresponding to closed contour line
#' @param co contour line
#'
#' @return the area within the contour line
#'
#' @export
#' @importFrom sp SpatialPolygons SpatialPointsDataFrame Polygons Polygon
#' @importFrom maptools ContourLines2SLDF
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' contourArea(cL[[8]])
#' # not run
contourArea=function(co)
##################################################################
{
	contour = ContourLines2SLDF(list(co))
      	coords = contour@lines[[1]]@Lines[[1]]@coords
      	poly=Polygon(coords)
	polys=Polygons(list(poly),"id")
	contourSp = SpatialPolygons(list(polys))
	surface = gArea(contourSp)

	return(surface)
}

##################################################################
#' listContourArea
#'
#' @details area of all contour lines in list
#' @param cL list of contour lines
#'
#' @return a list of areas
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' listContourArea(cL)
#' # not run
listContourArea=function(cL)
##################################################################
{
	return(sapply(cL,contourArea))
}

##################################################################
#' contourToSpp
#'
#' @details transform contour line into SpatialPolygons
#' @param co contour line (list with contour level and x,y coordinates
#' @param step grid resolution
#'
#' @return a list with components
#'\describe{
#'\item{sp}{SpatialPolygons corresponding to contour line}
#'\item{contour}{SpatialLines corresponding to contour line}
#'\item{polyBuff}{SpatialPolygons corresponding to buffer around contour line}
#'\item{surface}{SpatialPolygons area}
#'}
#'
#' @export
#' @importFrom rgeos gBuffer gArea gCentroid gContains gConvexHull gDifference gDistance gIntersects gWithin
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' contourToSpp(cL[[8]],0.1)
#' # not run
contourToSpp=function(co,step)
##################################################################
{
	contour = ContourLines2SLDF(list(co))
      	coords = contour@lines[[1]]@Lines[[1]]@coords
	# attention gBuffer renvoie 2 polygones
      	polyBuff = gBuffer(contour,width=0.0001*(1/step))
	poly=Polygon(coords)
	polys=Polygons(list(poly),"id")
	contourSp = SpatialPolygons(list(polys))
	surface = gArea(contourSp)
	return(list(sp=contourSp,contour=contour,polyBuff=polyBuff,surface=surface))
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
#' nPolySp
#' @details number of polygons in SpatialPolygons
#' @param sp SpatialPolygons
#' @return the number of polygons within the current zone
#'
#' @export
#'
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' print(paste(nPolySp(Z[[2]]),"polygons"))
#' @details 
#' # not run
nPolySp =function(sp)
##################################################################
{
	return(length(sp@polygons[[1]]@Polygons))
}

##################################################################
#' holeSp
#'
#' @details number of holes in SpatialPolygons
#' @param sp SpatialPolygons
#'
#' @return the number of holes within sp
#'
#' @export
#'
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' holeSp(Z[[5]]) #zone 5 has 1 hole
#' # not run
holeSp = function(sp)
##################################################################
{
	le=nPolySp(sp)
	nh=0
	if (le == 0) return(0)
	for (i in 1:le)
	{
		poly=sp@polygons[[1]]@Polygons[[i]]
		if(poly@hole) nh=nh+1

	}
	return (nh)
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
#' maxDistSP
#'
#' @details maximum distance within kth polygon of current zone 
#' @param sp SpatialPolygons
#' @return the maximum distance within sp
#'
#' @export
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' maxDistSp(Z[[5]])
#' # not run
maxDistSP=function(sp)
##################################################################
{
	return(max(dist(sp@polygons[[1]]@Polygons[[1]]@coords)))
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
#' getPolySp
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
#' sp=Z[[5]]
#' P1=getPolySp(sp,1)
#' P2=getPolySp(sp,2) # second polygon is a hole
#' plot(P1@coords,type="l")
#' lines(P2@coords,type="l",col="blue")
#' # not run
getPolySp = function(sp,k=1)
##################################################################
{
	p=sp@polygons[[1]]@Polygons[[k]]
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
#' lines(sph,type="l",col="blue")
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
#' polyToSp2
#'
#' @details transforms polygon into SpatialPolygons
#' @param p polygon 
#'
#' @return a SpatialPolygons
#'
#' @export
#'
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' sp=Z[[5]]
#' P1=getPolySp(sp,1)
#' sph=polyToSp2(P1)
#' plotZ(Z)
#' lines(sph,col="blue",lwd=2)
#' # not run
polyToSp2=function(p)
##################################################################
{
	sp=SpatialPolygons(list(Polygons(list(p),1:1)))
	return(sp)
}

##################################################################
#' lineToSp
#'
#' @details description, a paragraph
#' @param lin xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#  lin=data.frame(x=cL[[8]]$x,y=cL[[8]]$y)
#' sp=lineToSp(lin)
#' # not run
lineToSp=function(lin)
##################################################################
{
	sp=SpatialPolygons(list(Polygons(list(Polygon(lin,hole = FALSE)), "1")))
	return(sp)
}

##################################################################
#' normZcoords
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param boundary xxxx
#'
#' @return a list
#'
#' @export
#'
#' @examples
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
#' spnorm
#'
#' @details description, a paragraph
#' @param p xxxx
#' @param boundary xxxx
#'
#' @return a list
#'
#' @export
#'
#' @examples
#' # not run
spnorm = function(p, boundary)
##################################################################
{
	xmin=min(boundary$x)
	xmax=max(boundary$x)
	ymin=min(boundary$y)
	ymax=max(boundary$y)
	if (abs(xmax-xmin) < 1e-6) return(NULL)
	if (abs(ymax-ymin) < 1e-6) return(NULL)

	x= p@coords[,1]
	y = p@coords[,2]
	x = (x-xmin)/(xmax-xmin)
	y = (y-ymin)/(ymax-ymin)
	pn = Polygon(cbind(x,y))

	bx=boundary$x
	by=boundary$y
	bx = (bx-xmin)/(xmax-xmin)
	by = (by-ymin)/(ymax-ymin)
	bn = list(x=bx,y=by)
	return(list(pn=pn, boundaryn=bn))
}

##################################################################
#' normalize data coords and border
#'
#' @details normalize boundary between 0 and 1 and data coordinates accordingly
#' @param data data frame with x and y components
#' @param bd boundary (list with x and y components)
#'
#' @return a list with components
#'\describe{
#' \item{dataN}{normalized data}
#' \item{boundaryN}{normalized boundary}
#' \item{xmin}{minimum vaue of x within boundary}
#' \item{xmax}{maximum vaue of x within boundary}
#' \item{ymin}{minimum vaue of y within boundary}
#' \item{ymax}{maximum vaue of y within boundary}
#' }
#'
#'
#' @export
#'
#' @examples
#' # not run
datanorm = function(data, bd)
##################################################################
{
  xmin=min(bd$x)
  xmax=max(bd$x)
  ymin=min(bd$y)
  ymax=max(bd$y)
  if (abs(xmax-xmin) < 1e-6) return(NULL)
  if (abs(ymax-ymin) < 1e-6) return(NULL)

  x = data$x
  y = data$y
  x = (x-xmin)/(xmax-xmin)
  y = (y-ymin)/(ymax-ymin)
  data$x = x
  data$y = y

  #normalize bder
  bx=bd$x
  by=bd$y
  bx = (bx-xmin)/(xmax-xmin)
  by = (by-ymin)/(ymax-ymin)
  bdn = list(x=bx,y=by)
  return(list(dataN=data, boundaryN=bdn,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
}

##################################################################
#' normalize data coords
#'
#' @details normalize data coordinates between 0 and 1
#' @param data frame with x and y components
#'
#' @return a normalized data frame
#'
#' @export
#'
#' @examples
#' # not run
datanormXY = function(data)
##################################################################
{
  xmin=min(data$x)
  xmax=max(data$x)
  ymin=min(data$y)
  ymax=max(data$y)
  if (abs(xmax-xmin) < 1e-6) return(NULL)
  if (abs(ymax-ymin) < 1e-6) return(NULL)

  x = data$x
  y = data$y
  x = (x-xmin)/(xmax-xmin)
  y = (y-ymin)/(ymax-ymin)
  data$x = x
  data$y = y


  return(data)
}

##################################################################
#' normalize data coords with same ratio (for non square field)
#'
#' @details normalize x between 0 and 1, y and boundary with same ratio
#' @param data frame with x and y components
#' @param bd list with x and y components
#'
#' @return a list with components
#'\describe{
#' \item{dataN}{normalized data}
#' \item{boundaryN}{normalized boundary}
#' \item{ratio}{normalizing ratio}
#' \item{xmin}{minimum value of x within boundary}
#' \item{xmax}{maximum value of x within boundary}
#' \item{ymin}{minimum value of y within boundary}
#' \item{ymax}{maximum value of y within boundary}
#' }
#'
#' @export
#'
#' @examples
#' # not run
datanormX = function(data,bd)
##################################################################
{
  xmin=min(data$x)
  xmax=max(data$x)
  ymin=min(data$y)
  ymax=max(data$y)
  
  if (abs(xmax-xmin) < 1e-6) return(NULL)
  
  x = data$x
  y = data$y
  ratio=xmax-xmin
  x = (x-xmin)/ratio
  y = (y-ymin)/ratio
  data$x = x
  data$y = y
  #normalize boarder
  bx=bd$x
  by=bd$y
  bx = (bx-xmin)/ratio
  by = (by-ymin)/ratio
  bdn = list(x=bx,y=by)

  return(list(dataN=data, boundaryN=bdn,ratio=ratio,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
}

##################################################################
#' normSize
#'
#' @details description, a paragraph
#' @param boundaryN xxxx
#' @param minSize xxxx
#' @param minSizeNG xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
normSize=function(boundaryN,minSize,minSizeNG)
##################################################################
{
# normalize threshold for small zone detection
#
coords=boundaryN
poly=Polygon(coords)
polys=Polygons(list(poly),"id")
boundarySp = SpatialPolygons(list(polys))
boundaryArea=gArea(boundarySp)
minSize=minSize/boundaryArea
minSizeNG=minSizeNG/boundaryArea
print(paste("after standardization minSize=",minSize,"minSizeNG=",minSizeNG))

return(list(minSize=minSize,minSizeNG=minSizeNG))
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
#' plotListC
#'
#' @details add contour lines to a plot
#' @param cL list of contour lines
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' plot(mapTest$boundary,type="l")
#' plotListC(cL)
#' # not run
plotListC = function(cL,col="red")
##################################################################
{
	le=length(cL)
	if (le >0)
	   for (i in 1:le)
	   {
		ci=cL[[i]]
		lines(ci,col=col)
	   }

}



##################################################################
#' genQseq
#'
#' @details description, a paragraph
#' @param qProb probability vector used to generate quantile values
#' @param K zoning object, as returned by the calNei function
#' @param map object returned by function genMap or genMapR
#' @param i1 current zone index
#' @param i2 englobing zone index
#' @param LEQ length of quantile sequence 
#' @param MAXP maximum shift from center for quantile sequence 
#' @param disp 0: no info, 1: some info
#'
#' @return a plot
#'
#' @export
#' data(mapTest)
#' @examples
#' qProb=c(0.4,0.7)
#' ZK=initialZoning(qProb,mapTest)
#' K=ZK$resZ
#' print(K$lab)
#' genQseq(qProb,K,mapTest,1,2) # from label 3 to label 2
#' # not run
genQseq = function(qProb,K,map,i1,i2,LEQ=5,MAXP=0.1,disp=0)
##################################################################
{
# i1 = current zone index, i2= englobing zone index
#
  lab1 = K$lab[i1]
  lab2 = K$lab[i2]
  nq = length(qProb)

  inc = TRUE # increase quantile value
  if (lab1 >= lab2) inc = FALSE # decrease quantile value
  valRef= quantile(map$krigGrid,na.rm=TRUE,prob=qProb)
  # find quantile value corresponding to zone i1
  q1min = max(1,lab1-1)
  q1max = min(lab1,nq)

  if(inc)
	q1=q1max
  else
	q1=q1min

  prob1 = qProb[q1]
  q1= quantile(map$krigGrid,na.rm=TRUE,prob=prob1)
  # increase or decrease quantile value
  if (inc)
  {
  q2 = quantile(map$krigGrid,na.rm=TRUE,prob= min(prob1 + MAXP, 0.99))
  Qseq =seq(q1,q2,length.out=LEQ+1)
  }
  else
  {
  q3 = quantile(map$krigGrid,na.rm=TRUE,prob= max(prob1 - MAXP, 0.01))
  Qseq = seq(q1,q3,length.out=LEQ+1)
   }
  Qseq=Qseq[-1]
  Qseq=sort(unique(Qseq))
  if (disp>0)
     {
	print("Qseq=")
  	print(Qseq)
	}
  return(Qseq)
}

##################################################################
#' pointsSp
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons object
#' @param k polygon number
#' @param col color
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' pointsSp(Z[[1]])
#' # not run
pointsSp = function(sp,k=1,col="red")
##################################################################
{
	p = sp@polygons[[1]]@Polygons[[k]]
	points(p@coords,col=col)
}

##################################################################
#' linesSp
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons object
#' @param k polygon number
#' @param lty line type
#' @param col color
#' @param lwd line width
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' linesSp(Z[[4]])
#' # not run
linesSp = function(sp,k=1,lty=1,col="red",lwd=1)
##################################################################
{
	p = sp@polygons[[1]]@Polygons[[k]]
	lines(p@coords,lty=lty,col=col)
}

##################################################################
#' plotSp
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons object
#' @param k polygon number
#' @param xlim x range
#' @param ylim y range
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotSp(Z[[1]],xlim=c(0,1),ylim=c(0,1))
#' # not run
plotSp = function(sp,k=1,xlim,ylim)
##################################################################
{

	p = sp@polygons[[1]]@Polygons[[k]]
	cp=p@coords
	if (missing(xlim)) xlim=range(cp[,"x"])
	if (missing(ylim)) ylim=range(cp[,"y"])
	plot(cp,type="b",col="blue",xlim=xlim,ylim=ylim)
}

##################################################################
#' plotZ
#'
#' @details description, a paragraph
#' @param Z zoning geometry (list pf SpatialPolygons)
#' @param map map object returned by function genMap
#' @param id logical value, if TRUE display zone ids, if FALSE display zone numbers
#' @param noXY logical value, if TRUE do not display x and y axes
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z,mapTest)
#' # not run
plotZ = function(Z,map=NULL,id=FALSE,noXY=FALSE)
##################################################################
{

	if (!is.null(map))
   	   dispZ(map$step,matVal=map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0,noXY=noXY)
	else
          {
	  dispZ(map$step,matVal=NULL, nbLvl=0, zonePolygone=Z,id=id)
	  }
	 
}

##################################################################
#' checkContour
#'
#' @details check admissibility for contour line: surface >minSizeNG and refPoint close enough
#' @param contourSp SpatialPolygons
#'
#' @export
#'
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' sph=polyToSp(Z,5,2)
#' plotZ(Z)
#' lines(sph,type="l",col="blue")
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
#' polyToSp2
#'
#' @details transforms polygon into SpatialPolygons
#' @param p polygon 
#'
#' @return a SpatialPolygons
#'
#' @export
#'
#' @examples
#' ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' sp=Z[[5]]
#' P1=getPolySp(sp,1)
#' sph=polyToSp2(P1)
#' plotZ(Z)
#' lines(sph,col="blue",lwd=2)
#' # not run
polyToSp2=function(p)
##################################################################
{
	sp=SpatialPolygons(list(Polygons(list(p),1:1)))
	return(sp)
}

##################################################################
#' lineToSp
#'
#' @details description, a paragraph
#' @param lin xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#  lin=data.frame(x=cL[[8]]$x,y=cL[[8]]$y)
#' sp=lineToSp(lin)
#' # not run
lineToSp=function(lin)
##################################################################
{
	sp=SpatialPolygons(list(Polygons(list(Polygon(lin,hole = FALSE)), "1")))
	return(sp)
}

##################################################################
#' normZcoords
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param boundary xxxx
#'
#' @return a list
#'
#' @export
#'
#' @examples
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
#' spnorm
#'
#' @details description, a paragraph
#' @param p xxxx
#' @param boundary xxxx
#'
#' @return a list
#'
#' @export
#'
#' @examples
#' # not run
spnorm = function(p, boundary)
##################################################################
{
	xmin=min(boundary$x)
	xmax=max(boundary$x)
	ymin=min(boundary$y)
	ymax=max(boundary$y)
	if (abs(xmax-xmin) < 1e-6) return(NULL)
	if (abs(ymax-ymin) < 1e-6) return(NULL)

	x= p@coords[,1]
	y = p@coords[,2]
	x = (x-xmin)/(xmax-xmin)
	y = (y-ymin)/(ymax-ymin)
	pn = Polygon(cbind(x,y))

	bx=boundary$x
	by=boundary$y
	bx = (bx-xmin)/(xmax-xmin)
	by = (by-ymin)/(ymax-ymin)
	bn = list(x=bx,y=by)
	return(list(pn=pn, boundaryn=bn))
}

##################################################################
#' normalize data coords and border
#'
#' @details normalize boundary between 0 and 1 and data coordinates accordingly
#' @param data data frame with x and y components
#' @param bd boundary (list with x and y components)
#'
#' @return a list with components
#'\describe{
#' \item{dataN}{normalized data}
#' \item{boundaryN}{normalized boundary}
#' \item{xmin}{minimum vaue of x within boundary}
#' \item{xmax}{maximum vaue of x within boundary}
#' \item{ymin}{minimum vaue of y within boundary}
#' \item{ymax}{maximum vaue of y within boundary}
#' }
#'
#'
#' @export
#'
#' @examples
#' # not run
datanorm = function(data, bd)
##################################################################
{
  xmin=min(bd$x)
  xmax=max(bd$x)
  ymin=min(bd$y)
  ymax=max(bd$y)
  if (abs(xmax-xmin) < 1e-6) return(NULL)
  if (abs(ymax-ymin) < 1e-6) return(NULL)

  x = data$x
  y = data$y
  x = (x-xmin)/(xmax-xmin)
  y = (y-ymin)/(ymax-ymin)
  data$x = x
  data$y = y

  #normalize bder
  bx=bd$x
  by=bd$y
  bx = (bx-xmin)/(xmax-xmin)
  by = (by-ymin)/(ymax-ymin)
  bdn = list(x=bx,y=by)
  return(list(dataN=data, boundaryN=bdn,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
}

##################################################################
#' normalize data coords
#'
#' @details normalize data coordinates between 0 and 1
#' @param data frame with x and y components
#'
#' @return a normalized data frame
#'
#' @export
#'
#' @examples
#' # not run
datanormXY = function(data)
##################################################################
{
  xmin=min(data$x)
  xmax=max(data$x)
  ymin=min(data$y)
  ymax=max(data$y)
  if (abs(xmax-xmin) < 1e-6) return(NULL)
  if (abs(ymax-ymin) < 1e-6) return(NULL)

  x = data$x
  y = data$y
  x = (x-xmin)/(xmax-xmin)
  y = (y-ymin)/(ymax-ymin)
  data$x = x
  data$y = y


  return(data)
}

##################################################################
#' normalize data coords with same ratio (for non square field)
#'
#' @details normalize x between 0 and 1, y and boundary with same ratio
#' @param data frame with x and y components
#' @param bd list with x and y components
#'
#' @return a list with components
#'\describe{
#' \item{dataN}{normalized data}
#' \item{boundaryN}{normalized boundary}
#' \item{ratio}{normalizing ratio}
#' \item{xmin}{minimum value of x within boundary}
#' \item{xmax}{maximum value of x within boundary}
#' \item{ymin}{minimum value of y within boundary}
#' \item{ymax}{maximum value of y within boundary}
#' }
#'
#' @export
#'
#' @examples
#' # not run
datanormX = function(data,bd)
##################################################################
{
  xmin=min(data$x)
  xmax=max(data$x)
  ymin=min(data$y)
  ymax=max(data$y)
  
  if (abs(xmax-xmin) < 1e-6) return(NULL)
  
  x = data$x
  y = data$y
  ratio=xmax-xmin
  x = (x-xmin)/ratio
  y = (y-ymin)/ratio
  data$x = x
  data$y = y
  #normalize boarder
  bx=bd$x
  by=bd$y
  bx = (bx-xmin)/ratio
  by = (by-ymin)/ratio
  bdn = list(x=bx,y=by)

  return(list(dataN=data, boundaryN=bdn,ratio=ratio,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
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
#' plotListC
#'
#' @details add contour lines to a plot
#' @param cL list of contour lines
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' plot(mapTest$boundary,type="l")
#' plotListC(cL)
#' # not run
plotListC = function(cL,col="red")
##################################################################
{
	le=length(cL)
	if (le >0)
	   for (i in 1:le)
	   {
		ci=cL[[i]]
		lines(ci,col=col)
	   }

}



##################################################################
#' genQseq
#'
#' @details description, a paragraph
#' @param qProb probability vector used to generate quantile values
#' @param K zoning object, as returned by the calNei function
#' @param map object returned by function genMap or genMapR
#' @param i1 current zone index
#' @param i2 englobing zone index
#' @param LEQ length of quantile sequence 
#' @param MAXP maximum shift from center for quantile sequence 
#' @param disp 0: no info, 1: some info
#'
#' @return a plot
#'
#' @export
#' data(mapTest)
#' @examples
#' qProb=c(0.4,0.7)
#' ZK=initialZoning(qProb,mapTest)
#' K=ZK$resZ
#' print(K$lab)
#' genQseq(qProb,K,mapTest,1,2) # from label 3 to label 2
#' # not run
genQseq = function(qProb,K,map,i1,i2,LEQ=5,MAXP=0.1,disp=0)
##################################################################
{
# i1 = current zone index, i2= englobing zone index
#
  lab1 = K$lab[i1]
  lab2 = K$lab[i2]
  nq = length(qProb)

  inc = TRUE # increase quantile value
  if (lab1 >= lab2) inc = FALSE # decrease quantile value
  valRef= quantile(map$krigGrid,na.rm=TRUE,prob=qProb)
  # find quantile value corresponding to zone i1
  q1min = max(1,lab1-1)
  q1max = min(lab1,nq)

  if(inc)
	q1=q1max
  else
	q1=q1min

  prob1 = qProb[q1]
  q1= quantile(map$krigGrid,na.rm=TRUE,prob=prob1)
  # increase or decrease quantile value
  if (inc)
  {
  q2 = quantile(map$krigGrid,na.rm=TRUE,prob= min(prob1 + MAXP, 0.99))
  Qseq =seq(q1,q2,length.out=LEQ+1)
  }
  else
  {
  q3 = quantile(map$krigGrid,na.rm=TRUE,prob= max(prob1 - MAXP, 0.01))
  Qseq = seq(q1,q3,length.out=LEQ+1)
   }
  Qseq=Qseq[-1]
  Qseq=sort(unique(Qseq))
  if (disp>0)
     {
	print("Qseq=")
  	print(Qseq)
	}
  return(Qseq)
}

##################################################################
#' pointsSp
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons object
#' @param k polygon number
#' @param col color
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' pointsSp(Z[[1]])
#' # not run
pointsSp = function(sp,k=1,col="red")
##################################################################
{
	p = sp@polygons[[1]]@Polygons[[k]]
	points(p@coords,col=col)
}

##################################################################
#' linesSp
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons object
#' @param k polygon number
#' @param lty line type
#' @param col color
#' @param lwd line width
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' linesSp(Z[[4]])
#' # not run
linesSp = function(sp,k=1,lty=1,col="red",lwd=1)
##################################################################
{
	p = sp@polygons[[1]]@Polygons[[k]]
	lines(p@coords,lty=lty,col=col)
}

##################################################################
#' plotSp
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons object
#' @param k polygon number
#' @param xlim x range
#' @param ylim y range
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotSp(Z[[1]],xlim=c(0,1),ylim=c(0,1))
#' # not run
plotSp = function(sp,k=1,xlim,ylim)
##################################################################
{

	p = sp@polygons[[1]]@Polygons[[k]]
	cp=p@coords
	if (missing(xlim)) xlim=range(cp[,"x"])
	if (missing(ylim)) ylim=range(cp[,"y"])
	plot(cp,type="b",col="blue",xlim=xlim,ylim=ylim)
}

##################################################################
#' plotZ
#'
#' @details description, a paragraph
#' @param Z zoning geometry (list pf SpatialPolygons)
#' @param map map object returned by function genMap
#' @param id logical value, if TRUE display zone ids, if FALSE display zone numbers
#' @param noXY logical value, if TRUE do not display x and y axes
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z,mapTest)
#' # not run
plotZ = function(Z,map=NULL,id=FALSE,noXY=FALSE)
##################################################################
{

	if (!is.null(map))
   	   dispZ(map$step,matVal=map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0,noXY=noXY)
	else
          {
	  dispZ(map$step,matVal=NULL, nbLvl=0, zonePolygone=Z,id=id)
	  }
	 
}

##################################################################
#' checkContour
#'
#' @details check admissibility for contour line: surface >minSizeNG and refPoint close enough
#' @param contourSp SpatialPolygons corresponding to closed contou line
#' @param step grid resolution
#' @param refPoint referene point
#' @param minSizeNG zone area threshold under which a zone is not admissible
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
checkContour = function(contourSp,step,refPoint,minSizeNG=1e-3)
##################################################################
{
	#  contourSp  spatial object

	polyBuff = gBuffer(contourSp,width=0.0001*(1/step)) #spatialPolygon

	surface = gArea(contourSp)
 	condi = surface <=minSizeNG || gDistance(gCentroid(polyBuff),refPoint)>=0.1
	if (condi)
	   return (NULL)
	else
	   return (list(contourSp=contourSp,polyBuff=polyBuff))
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
#' createHoles
#'
#' @details description, a paragraph
#' @param Z list of zones, each zone is a SpatialPolygons
#'
#' @return a list of zones where holes are SpatialPolygons
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
#' @details creates zone excluding holes
#' @param zoneMain main zone
#' @param zoneSuppr zone inside main zone
#'
#' @return a new zone
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
#' cleanSp
#'
#' @details removes from sp polygons that are too small (artefacts of gDifference)
#' @param sp SpatialPolygons
#' @param tol mnimum area for removal
#' @return a SpatialPolygons
#'
#' @export
#'
#' @examples
#' # not run
cleanSp = function(sp,tol=1e-5)
##################################################################
{
  # number of polygons
  n=nPolySp(sp)
  polyl=list()
  k=0
  # remove polygons that are too small (<1e-6)
  # (artefacts of gDifference)
  for (i in 1:n)
  {
	  poly=sp@polygons[[1]]@Polygons[[i]]
	  area=poly@area
	  if (area >= tol)
	   {
	   k=k+1
	   polyl[[k]]=poly
	   }
  }

  if (length(polyl)==0) return(NULL)

  polys=Polygons(polyl,"id")
  spc=  SpatialPolygons(list(polys))
  return(spc)
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
#' ptInSp
#'
#' @details finds data points in sp
#' @param sp SpatialPolygons
#' @param pts data points
#' @param hole if TRUE also consider points in holes
#'
#' @return a data frame with data points within sp
#'
#' @export
#'
#' @examples
#' @examples
#' data(mapTest)
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' ptsInSp(Z[[5]],mapTest$krigData) # 5 data points within zone 5
#'
#' # not run
ptsInSp=function(sp,pts,hole=FALSE)
##################################################################
{
	PointsSpatiaux=SpatialPoints(pts)
	logicalPoly=point.in.polygon(PointsSpatiaux$x,PointsSpatiaux$y,sp@polygons[[1]]@Polygons[[1]]@coords[,1],sp@polygons[[1]]@Polygons[[1]]@coords[,2])
    if(hole && (length(sp@polygons[[1]]@Polygons)>1))
    {
      #pts in holes
      for(k in 2:length(sp@polygons[[1]]@Polygons))
      {
        logicalPoly=logicalPoly-point.in.polygon(PointsSpatiaux$x,PointsSpatiaux$y,sp@polygons[[1]]@Polygons[[k]]@coords[,1],sp@polygons[[1]]@Polygons[[k]]@coords[,2])
      }
    }

    return(pts[logicalPoly!=0,])
}

##################################################################
#' searchNODcrit
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param le xxxx
#' @param zk xxxx
#' @param critere xxxx
#' @param cost xxxx
#' @param costL xxxx
#' @param nz xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
searchNODcrit=function(qProb,le,zk,critere,cost,costL,nz)
##################################################################
{
# zk: list of zonings
# le: list index
# crit: list of criteria
   critf=unlist(critere[[le]])
   costf=unlist(cost[[le]])
   costfL=unlist(costL[[le]])
   nzf=unlist(nz[[le]])
#
 # check for degenerated zonings
  # number of zone labels < number of quantiles+1
  # if labels start at 1
    nq0=length(qProb)+1
  # compute number of labs for each solution
    nqf=c()
    ind=list()
    bestcrit=list()
    bestcost=list()
    bestcostL=list()
    bestnz=list()
    nq=list()

    lk=length(zk[[le]])
    kk=1:lk

    for (ilab in 1:lk)
    {
	u=unique(zk[[le]][[ilab]]$lab)
	lu=length(u)
	nqf=c(nqf,lu)
    }

    while (nq0>1)
    {
    maskNOD=(nqf==nq0)
    critq=critf[maskNOD]

    if(!is.na(critq) && (length(critq)>0))
	{
	ii=which(critq == max(critq)) #search from the best non degenerated ones
	critmax=critq[ii]
	# find original index in critf vector
	k0=kk[maskNOD]
	jj=k0[ii]
	labq=paste("q",nq0-1,sep="")
    	ind[[labq]]=jj
	bestcrit[[labq]]=critf[jj]
	bestcost[[labq]]=costf[jj]
	bestcostL[[labq]]=costfL[jj]
	bestnz[[labq]]=nzf[jj]
	nq[[labq]]=nq0-1

 	}
     nq0=nq0-1

    }#end while
  return(list(ind=ind,critList=bestcrit,costList=bestcost,costLList=bestcostL,nzList=bestnz,nq=nq))
}
##################################################################
searchNODcrit1=function(qProb,crit)
#' searchNODcrit1
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param crit xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
##################################################################
{
# zk: list of zonings
# le: list index

zk=crit$zk
cr=crit$criterion
cost=crit$cost
costL=crit$costL
nz=crit$nz
le=length(zk)

# crit: list of criteria
   critf=unlist(cr[[le]])
   costf=unlist(cost[[le]])
   costfL=unlist(costL[[le]])
   nzf=unlist(nz[[le]])
#
 # check for degenerated zonings
  # number of zone labels < number of quantiles+1
  # if labels start at 1
    nq0=length(qProb)+1
  # compute number of labs for each solution
    nqf=c()
    ind=list()
    bestcrit=list()
    bestcost=list()
    bestcostL=list()
    bestnz=list()
    nq=list()

    lk=length(zk[[le]])
    kk=1:lk
    
    for (ilab in 1:lk)
    {
	u=unique(zk[[le]][[ilab]]$lab)
	lu=length(u)
	nqf=c(nqf,lu)
    }
  
    while (nq0>1)
    {
    maskNOD=(nqf==nq0)
    critq=critf[maskNOD]
    
    if(!is.na(critq) && (length(critq)>0))
	{
	ii=which(critq == max(critq)) #search from the best non degenerated ones
	critmax=critq[ii]
	# find original index in critf vector
	k0=kk[maskNOD]
	jj=k0[ii]
	labq=paste("q",nq0-1,sep="")
    	ind[[labq]]=jj
	bestcrit[[labq]]=critf[jj]
	bestcost[[labq]]=costf[jj]
	bestcostL[[labq]]=costfL[jj]
	bestnz[[labq]]=nzf[jj]
	nq[[labq]]=nq0-1
    
 	}
     nq0=nq0-1
    	
    }#end while
 
  
    
    return(list(ind=ind,critList=bestcrit,costList=bestcost,costLList=bestcostL,nzList=bestnz,nq=nq))
}

##################################################################
#' printLabZ
#'
#' @details description, a paragraph
#' @param zkl xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
printLabZ=function(zkl)
##################################################################
{
#zkl is a list of several Ks
	lk=length(zkl)
	labZ=list()
	for (k in 1:lk)
	{
		labk=unlist(zkl[[k]]$lab)
		labq=paste(length(unique(labk))-1,"q",sep="")
		labZ[[k]]=labk
		print(paste(labq,"zone labels=",labZ[k]))
	}
	return(labZ)
}

##################################################################
#' normDistMat
#'
#' @details normalize distance matrix so that diagonal is equal to 1
#' @param matDistanceCorr corrected distance matrix using pErr, result of calDistance
#' @param optiCrit criterion choice
#'
#' @return a normalized distance matrix
#'
#' @export
#'
#' @examples
#' # not run
normDistMat=function(matDistanceCorr,optiCrit)
##################################################################
{
# other values of optiCrit not managed
#
	if(optiCrit==4||optiCrit==6)
		normMD=distanceNormalisationSqrt(matDistanceCorr)
	if(optiCrit==2)
		normMD=distanceNormalisationSum(matDistanceCorr)

return(normMD)

}

##################################################################
#' trLabZone
#'
#' @details description, a paragraph
#' @param K1 xxxx
#' @param K2 xxxx
#' @param Z1 xxxx
#' @param Z2 xxxx
#' @param map xxxx
#' @param qProb xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
trLabZone=function(K1,K2,Z1,Z2,map,qProb,disp=0)
##################################################################
#transfer zone labels from K1 to K2
# ids of zones in Z1 and Z2 used to guide transfer

{
# only K2$lab is modified
# by taking labels from K1
lab2=rep(1,length(Z2))

lab1=K1$lab
id1s=getIds(Z1)
q1 = quantile(map$krigGrid,na.rm=TRUE,prob=qProb)

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
      		if (K2$meanZone[iZ]>(q1[j]+0.1))
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
#' @details description, a paragraph
#' @param iZ xxxx
#' @param zoneNModif xxxx
#' @param K xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
getClosestZone=function(iZ,Z,zoneNModif)
##################################################################
{
	imin=0
	ni=1:length(Z)
	# exclude current zone
	ni=ni[ni != iZ]

	# exclude neighbor zones
	Ns=getNs(zoneNModif,iZ)
  	listeV=grep(TRUE, Ns)
	for (i in listeV) ni=ni[ni != i]

	# exclude englobing zone
	iE = detZoneEng(iZ,Z,zoneNModif)
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

##################################################################
#' find contour for a given vRef quantile value
#'
#' @details description, a paragraph
#' @param iC xxxx
#' @param Z xxxx
#' @param K xxxx
#' @param map xxxx
#' @param vRef xxxx
#' @param envel xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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

##################################################################
#' linesC
#'
#' @details add contour Lines to plot
#' @param list of contour lines
#' @param col line color
#'
#' @return an empty value
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' plot(mapTest$boundary)
#' linesC(cL,col="black")
#' # not run
linesC = function(listContour,col="blue")
##################################################################
{
for (i in 1:length(listContour))
{
	lines(listContour[[i]],col=col)
	}
return()
}

##################################################################
#' interCB
#'
#' @details description, a paragraph
#' @param co xxxx
#' @param step xxxx
#' @param bd xxxx
#' @param envel xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
interCB = function(co,step,bd=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),envel,disp=0)
##################################################################
{
  #returns spatial polygon corresponding to intersection of contour with boundary
  #
	polygoneGlobal=SpatialPolygons(list(Polygons(list(Polygon(bd, hole = FALSE)), "1")))
	contourL = ContourLines2SLDF(list(co))
	polyBuff=gBuffer(contourL,width=0.0001*step)
	polyDiff=gDifference(polygoneGlobal,polyBuff)
  recupPoly=separationPoly(polyDiff)

	ler=length(recupPoly)
	if(ler<2) return(NULL) # intersection=cadre -> degenerate contour

	sp1=recupPoly[[1]] #
	sp2=recupPoly[[2]] #
	# keep the smallest one that is within the envelope
	sp=sp2
	if (gContains(envel,sp1))
	{
		sp=sp1
		if (gContains(envel,sp2) & (gArea(sp2)<gArea(sp1))) sp=sp2
	}
	if(disp) linesSp(sp)
	return(sp)
}

##################################################################
#' getNq
#'
#' @details description, a paragraph
#' @param critList xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
getNq=function(critList)
##################################################################
{
	n=names(critList)
	for (qq in 1:length(critList))
    	{
		nq=sapply(strsplit(n[qq],"q"),function(x){return(x[2])})
		nq=as.numeric(nq)
		print(paste("criterion=",round(critList[[qq]][1],2),"nq=",nq))
	}
	return(nq)
}

##################################################################
#' addContour
#'
#' @details description, a paragraph
#' @param map xxxx
#' @param val xxxx
#' @param col xxxx
#' @param super xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
addContour=function(map,val,col="blue",super=T)
##################################################################
{
# IS 19/05/2017 change boundary by map$boundary!
for ( v in val)
  {
    listeContour=list()
    listeContour = contourAuto(listeContour,map$step,map$xsize,map$ysize,map$krigGrid,v,map$boundary)
    lc = length(listeContour)
    # intersect contour with boundary
    # and transform contours into sps
    if (lc >0)
    {
	    sp=list()
	    if (!super) plot(map$boundary,type="l") # new plot
    	for (k in 1:lc)
    	   lines(listeContour[[k]],col=col) # add contour to existing plot
  	}
  }
  return()
}

##########################################################################
#' extractionPoly
#'
#' @details description, a paragraph
#' @param polyTot xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
extractionPoly=function(polyTot)
###########################################################################
{
#fonction permettant d'extraire chaque sous-élément d'un spatial polygone dans un Spatial Polygone qui lui est propre
#(gestion identique des trous et des polygones pleins)
#entrée:SpatialPolygons
#sortie:liste de SpatialPolygons contenant les sous-éléments du précédent

  contenuPoly=polyTot@polygons[[1]]@Polygons
  listePolyExtract=list()
  for(i in (1:length(contenuPoly)))
  {
    listePolyExtract[[i]]=SpatialPolygons(list(Polygons( contenuPoly[[i]] , "1")))

  }


  return(listPolyExtract)
}

###################################
#' plotVario
#'
#' @details description, a paragraph
#' @param map xxxx
#' @param ylim xxxx
#'
#' @return a ?
#'
#' @export
#' @importFrom RandomFields RFempiricalvariogram
#'
#' @examples
#' # not run
plotVario=function(map,ylim=NULL)
###################################
{
  modelGen=map$modelGen
  data=map$rawData #raw data
  dataK=map$krigData #kriged data
  empvario=RFempiricalvariogram(data=data)
  empvarioK=RFempiricalvariogram(data=dataK)

  if(!is.null(modelGen))
	{
	  plot(empvario,model=modelGen,ylim=ylim)
	  #kriged variogram
    # IS 19/05/2017: add comment for x11
    #x11()
	  plot(empvarioK,model=modelGen,ylim=ylim,col="blue")
	}
  else	plot(empvario)
}

#########################
#' costLab
#'
#' @details description, a paragraph
#' @param K xxxx
#' @param map xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
costLab=function(K,map)
#########################
{
#number of labels
	uni=unique(K$lab)
	nL=length(uni)
	nZ=length(K$lab)
	listZonePoint=K$listZonePoint
	tabVal=map$krigData
	surfVoronoi=map$krigSurfVoronoi
#find which zones are assigned to each label
	zlab=list()
	vZ=1:nZ
	for (ilab in uni)
       	{
	mask=which(K$lab==ilab)
	zlab[[ilab]]=vZ[mask]
       	}
# compute costL (per label)
       res = SigmaL2(zlab,listZonePoint,tabVal,surfVoronoi)
       return(res$cL)
}

################################################################
#' SigmaL2
#'
#' @details description, a paragraph
#' @param zlab xxxx
#' @param listZonePoint xxxx
#' @param tabVal xxxx
#' @param surfVoronoi xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
SigmaL2=function(zlab,listZonePoint,tabVal,surfVoronoi)
################################################################
{
# zlab: list with zone numbers for each zone label
# compute overall mean and std of all zones for each label
# and overall cost per label
# remove empty labels from zlab
  zlabNULL=sapply(zlab,is.null)
  zlab=zlab[!zlabNULL]
  nL=length(zlab)
  mL=c()
  SL=c()
  vLab=list()
  voroLab=list()
  #first compute mean
  for (k in 1:nL)
      {
      zlabK=zlab[[k]]
      vLabK=c()
      voroLabK=c()
      for (j in zlabK)
      {
	vLabK=c(vLabK,tabVal@data[listZonePoint[[j]],1])	#all data values for zones with label k
      	voroLabK=c(voroLabK,surfVoronoi[listZonePoint[[j]]])  # corresp. voronoi surfaces
      }
      vLab[[k]]=vLabK
      voroLab[[k]]=voroLabK
      SL[k]= sum(voroLabK)
      mL[k]=sum(vLabK*voroLabK)/SL[k]
      }
  # then compute sd
  sigmaL2=rep(0,nL)
  for (k in 1:nL)
      {
      sigmaL2[k]=sum(((vLab[[k]]-mL[k])^2)*voroLab[[k]])/SL[k]
      }

  cL=sum(sigmaL2*SL)/sum(SL)

  return(list(cL=cL,sigmaL2=sigmaL2,SL=SL,mL=mL,voroLab=voroLab))
}

################################################################
#' meanL
#'
#' @details description, a paragraph
#' @param zlab xxxx
#' @param listZonePoint xxxx
#' @param tabVal xxxx
#' @param surfVoronoi xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
meanL=function(zlab,listZonePoint,tabVal,surfVoronoi)
################################################################
{
# zlab: list with zone numbers for each zone label

  nL=length(zlab)
  mL=c()
  SL=c()
  for (k in 1:nL)
      {
      zlabK=zlab[[k]]
      vLabK=c()
      voroLabK=c()
      for (j in zlabK)
      {
	vLabK=c(vLabK,tabVal@data[listZonePoint[[j]],1])	#all data values for zones with label k
      	voroLabK=c(voroLabK,surfVoronoi[listZonePoint[[j]]])  # corresp. voronoi surfaces
      }
      SL[k]= sum(voroLabK)
      mL[k]=sum(vLabK*voroLabK)/SL[k]
      }

      return(list(mL=mL,SL=SL))
}

######################################################################
#' listSeeds
#'
#' @details description, a paragraph
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
listSeeds=function()
#######################################################################
{
a = system("ls res-simuseed*-4q*.csv",intern=TRUE)
vseed=c()
for (f in a)
{
	f1=strsplit(f,split="res-simuseed")[[1]][2]
	f2=as.numeric(strsplit(f1,split="-")[[1]][1])
	vseed=c(vseed,f2)
}
return(vseed)
}

###########################
#' plotmdist
#'
#' @details description, a paragraph
#' @param md xxxx
#' @param pdf xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
plotmdist=function(md,pdf)
###########################
{
	system("cp MOD.matDistance1.tex mdist.tex")
	 write("\\tiny",file="mdist.tex",append=TRUE)
 	write(print(xtable(md,table.placement=NULL,latex.environment="",floating.environment="center")),file="mdist.tex",append=TRUE)
	write("\\end{document}",file="mdist.tex",append=TRUE)
	system("pdflatex mdist")
	system(paste("mv mdist.pdf",pdf))
	return()
}

###########################
#' plotmat
#'
#' @details description, a paragraph
#' @param m xxxx
#' @param pdf xxxx
#'
#' @return a ?
#'
#' @export
#' @importFrom xtable xtable
#'
#' @examples
#' # not run
plotmat=function(m,pdf)
###########################
{
	system("cp MOD.mat.tex mdist.tex")
	 write("\\tiny",file="mdist.tex",append=TRUE)
 	write(print(xtable(m,table.placement=NULL,latex.environment="",floating.environment="center")),file="mdist.tex",append=TRUE)
	write("\\end{document}",file="mdist.tex",append=TRUE)
	system("pdflatex mdist")
	system(paste("mv mdist.pdf",pdf))
	return()
}

###########################
#' plotOpt
#'
#' @details description, a paragraph
#' @param seed xxxx
#' @param pdf xxxx
#' @param k xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
plotOpt=function(seed,pdf,k=5)
###########################
{
	system("cp MOD.matDistance1.tex mdist.tex")
	 write("\\tiny",file="mdist.tex",append=TRUE)
	 file0=paste("opt-seed",seed,"crit-cost.pdf",sep="")
	 system(paste("pdfcrop",file0,"toto.pdf"))
	 system(paste("mv toto.pdf",file0))
	 file1=paste("opt-seed",seed,"-critA.pdf",sep="")
	 system(paste("pdfcrop",file1,"toto.pdf"))
	 system(paste("mv toto.pdf",file1))
	 file2=paste("opt-seed",seed,"-costA.pdf",sep="")
	 system(paste("pdfcrop",file2,"toto.pdf"))
	 system(paste("mv toto.pdf",file2))
	 file3=paste("opt-seed",seed,"-m.pdf",sep="")
	 system(paste("pdfcrop",file3,"toto.pdf"))
	 system(paste("mv toto.pdf",file3))
	 # adj plots
	 gr0=paste("\\includegraphics[width=\\linewidth]{",file0,"}\\\\",sep="")
	 write(gr0,file="mdist.tex",append=TRUE)
	 # adj tables
	gr2=paste("\\begin{tabular}{cc}\\includegraphics[width=0.5\\linewidth]{",file1,"}&\\includegraphics[width=0.5\\linewidth]{",file2,"}\\end{tabular}",sep="")
 	write(gr2,file="mdist.tex",append=TRUE)
	gr3=paste("\\vspace{-5cm}\\includegraphics[scale=0.9]{",file3,"}\\\\",sep="")
 	write(gr3,file="mdist.tex",append=TRUE)
	write("\\end{document}",file="mdist.tex",append=TRUE)
	#
	system("pdflatex mdist")
	system(paste("mv mdist.pdf",pdf))

	return()
}

#####################################
#' getBestMloop
#'
#' @details description, a paragraph
#' @param seed xxxx
#' @param thr xxxx
#'
#' @return a ?
#'
#' @export
#' @importFrom utils read.table
#'
#' @examples
#' # not run
getBestMloop=function(seed,thr=0.5)
#####################################
{

if(is.null(m1)) m1=read.table(paste("res-simuseed",seed,"-1q-pE0.3.csv",sep=""))
if(is.null(m2)) m2=read.table(paste("res-simuseed",seed,"-2q-pE0.3.csv",sep=""))
if(is.null(m3)) m3=read.table(paste("res-simuseed",seed,"-3q-pE0.3.csv",sep=""))
if(is.null(m4)) m4=read.table(paste("res-simuseed",seed,"-4q-pE0.3.csv",sep=""))
if(is.null(m5)) m5=read.table(paste("res-simuseed",seed,"-5q-pE0.3.csv",sep=""))
# remove degenerate quantiles
mask2=m2[,"nq"]==2
mb2=m2[mask2,]
mask3=m3[,"nq"]==3
mb3=m3[mask3,]
mask4=m4[,"nq"]==4
mb4=m4[mask4,]
mask5=m5[,"nq"]==5
mb5=m5[mask5,]
crit1=m1[1,"crit"]
crit2=mb2[1,"crit"]
crit3=mb3[1,"crit"]
crit4=mb4[1,"crit"]
crit5=mb5[1,"crit"]

mask1=(-m1[,"crit"]+crit1)<=thr
mask2=(-m2[,"crit"]+crit2)<=thr
mask3=(-m3[,"crit"]+crit3)<=thr
mask4=(-m4[,"crit"]+crit4)<=thr
mask5=(-m5[,"crit"]+crit5)<=thr

return(list(mb1=m1[mask1,],mb2=m2[mask2,],mb3=m3[mask3,],mb4=m4[mask4,],mb5=m5[mask5,]))
}

#####################################
#' getKMloop
#'
#' @details description, a paragraph
#' @param seed xxxx
#' @param m1 xxxx
#' @param m2 xxxx
#' @param m3 xxxx
#' @param m4 xxxx
#' @param m5 xxxx
#' @param k xxxx
#' @param pE xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
getKMloop=function(seed,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,k=5,pE=0.9)
#####################################
{

if(is.null(m1)) m1=read.table(paste("res-simuseed",seed,"-1q-pE",pE,".csv",sep=""))
if(is.null(m2)) m2=read.table(paste("res-simuseed",seed,"-2q-pE",pE,".csv",sep=""))
if(is.null(m3)) m3=read.table(paste("res-simuseed",seed,"-3q-pE",pE,".csv",sep=""))
if(is.null(m4)) m4=read.table(paste("res-simuseed",seed,"-4q-pE",pE,".csv",sep=""))
if(is.null(m5)) m5=read.table(paste("res-simuseed",seed,"-5q-pE",pE,".csv",sep=""))
# remove degenerate quantiles
mask2=m2[,"nq"]==2
mb2=m2[mask2,]
mask3=m3[,"nq"]==3
mb3=m3[mask3,]
mask4=m4[,"nq"]==4
mb4=m4[mask4,]
mask5=m5[,"nq"]==5
mb5=m5[mask5,]

mb1=m1[1:k,]
mb1=cbind(mb1[,-ncol(mb1)],matrix(NA,ncol=4,nrow=k),mb1[,ncol(mb1),])
rownames(mb1)=paste("m1-",1:k,sep="")
colnames(mb1)=c("crit","cost","costL","nz",paste("q",1:5,sep=""),"nq")
mb2=m2[1:k,]
mb2=cbind(mb2[,-ncol(mb2)],matrix(NA,ncol=3,nrow=k),mb2[,ncol(mb2),])
rownames(mb2)=paste("m2-",1:k,sep="")
colnames(mb2)=c("crit","cost","costL","nz",paste("q",1:5,sep=""),"nq")
mb3=m3[1:k,]
mb3=cbind(mb3[,-ncol(mb3)],matrix(NA,ncol=2,nrow=k),mb3[,ncol(mb3),])
rownames(mb3)=paste("m3-",1:k,sep="")
colnames(mb3)=c("crit","cost","costL","nz",paste("q",1:5,sep=""),"nq")
mb4=m4[1:k,]
mb4=cbind(mb4[,-ncol(mb4)],matrix(NA,ncol=1,nrow=k),mb4[,ncol(mb4),])
rownames(mb4)=paste("m4-",1:k,sep="")
colnames(mb4)=c("crit","cost","costL","nz",paste("q",1:5,sep=""),"nq")
mb5=m5[1:k,]
rownames(mb5)=paste("m5-",1:k,sep="")
colnames(mb5)=c("crit","cost","costL","nz",paste("q",1:5,sep=""),"nq")

return(rbind(mb1,mb2,mb3,mb4,mb5))
}

######################################
#' meansdSimu
#'
#' @details description, a paragraph
#' @param vseed xxxx
#' @param krig xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
meansdSimu=function(vseed=NULL,krig=1)
######################################
{

if (is.null(vseed))
{
	vseed=listSeeds()
}
m=c()
km=m
sdd=m
sdkd=m
for (seed in vseed)
{

    map=genMap(DataObj=NULL,seed=seed,disp=0,krig=krig)
    v=map$rawData
    v=v@data[,1]
    kv=map$krigData
    kv=kv@data[,1]
    m=c(m,mean(v))
    sdd=c(sdd,sd(v))
    km=c(km,mean(kv))
    sdkd=c(sdkd,sd(kv))
    }
mat=cbind(m,km,sdd,sdkd)
return(mat)
}

######################################
#' calSurf
#'
#' @details description, a paragraph
#' @param polyL xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calSurf=function(polyL)
######################################
{
  listeSurf=0

  #list of polys
  for (i in 1:length(polyL))
  {
    #création d'une variable de type 'Polygons' à partir des sommets définissant les polygones
    tmpPoly=Polygons(list(Polygon(data.frame(polyL[[i]]$x,polyL[[i]]$y))),paste("Pol",i,sep=""))
    #then use gArea
    listSurf[i]=gArea(SpatialPolygons(list(tmpPoly)))
  }

  return(listSurf)
}

##########################################################
#' calRapPolygone
#'
#' @details description, a paragraph
#' @param Surf xxxx
#' @param polyL xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calRapPolygone=function(Surf,polyL)
##########################################################
{
	perim=0
  #compute polygon perimeter
	for (ib in 1:(length(polyL[[1]])-1))
	{
		perim=perim+as.numeric(dist(matrix( c(polyL[ib],polyL[ib+1],polyL[ib],polyL[ib+1]),2,2)))
	}
	return(Surf/(perim^2))
}


######################################
#' transfoSpPoly
#'
#' @details description, a paragraph
#' @param polyL xxxx
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
transfoSpPoly=function(polyL)
######################################
{
	polyLSp=list(0)
  #add to polyLSp a dataframe with polygon vertices
  # make it a spatial object
	for (i in 1:length(polyL))
	{
		polyLSp[[i]]=data.frame(x=polyL[[i]]$x,y=polyL[[i]]$y,col=0)
		sp::coordinates(polyLSp[[i]])=~x+y
	}
	return(polyLSp)
}

###################
valZ=function(map,K)
#' valZ
#'
#' @details description, a paragraph
#' @param map xxxx
#' @param K xxxx
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
###################
{
tab=map$krigData
pt=K$listZonePoint
mu=K$meanZone
ord=order(mu)
val=list()
k=0
    for(ii in ord)
    {
      k=k+1
      val[[k]]=tab[[1]][pt[[ii]]]
    
  }
    
  return(list(val=val,ord=ord))
}
########################
getZsize=function(Z)
#' getZsize
#'
#' @details description, a paragraph
#' @param Z xxxx
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
########################
{
ms=sapply(Z,function(x)return(x@bbox[,"max"]))
vs=apply(ms,1,max)
return(vs)
}
########################
findZCenterpt=function(data,K,num=NULL)
#' findZCenter
#'
#' @details description, a paragraph
#' @param data xxxx
#' @param K xxxx
#' @param num xxxx
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
########################
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
return(ptz)
}
########################
findZCenter=function(Z,num=NULL)
#' findZCenter
#'
#' @details description, a paragraph
#' @param data xxxx
#' @param num xxxx
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
########################
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
return(ptz)
}
#############################
#' superLines
#'
#' @details converts boundary (list of x and y pts) into Spatial Lines
#' @param boundary list, contains x and y coordinates of map boundaries
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a SpatialLines object
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' superL=superLines(mapTest$boundary)
#' plot(superL)
#' # not run
superLines=function(boundary)
#############################
{
  #transform boundary into SpatialLines
  boundary = data.frame(boundary)
  sp::coordinates(boundary)=~x+y
  bl=Line(coordinates(boundary))
  bspl=SpatialLines(list(Lines(list(bl),'1')))
  bdLines = bspl@lines[[1]]@Lines[[1]]
  listBdLines=list(Lines(list(Line(bdLines@coords[1:2,])),'1'))
  for (i in 2:(length(bdLines@coords)/2-1))
  {
    listBdLines[[i]] = Lines(list(Line(bdLines@coords[i:(i+1),])),paste(i))
  }
  SuperLines = SpatialLines(listBdLines)
  return(SuperLines)
}