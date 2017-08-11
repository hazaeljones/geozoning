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
#' getPolySp
#'
#' @details get the kth polygon of the current SpatialPolygons
#' @param k polygon number
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
#' @details transform closed line into SpatialPolygons
#' @param lin list with x and y line coordinates
#'
#' @return a SpatialPolygons
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
#' spnorm
#'
#' @details normise Polygon according to border limits
#' @param sp object of class Polygons
#' @param boundary list with x and y components, used to normalize sp
#'
#' @return a list with components
#' \describe{
#' \item{pn}{normalized Polygon}
#' \item{boundaryn}{normalized boundary}
#' }
#'
#' @export
#'
#' @examples
#' z=readS("Field_8_zones.shp",dir="../data/")
#' bb=list(x=z$sp@bbox[1,],y=z$sp@bbox[2,])
#' P1=getPolySp(z$sp,1)
#' NP1=spnorm(P1,bb)$pn
#' Nbb=spnorm(P1,bb)$boundaryn
#' plot(NP1@coords,xlim=Nbb$x,ylim=Nbb$y)
#' # not run
spnorm = function(sp, boundary)
##################################################################
{
	xmin=min(boundary$x)
	xmax=max(boundary$x)
	ymin=min(boundary$y)
	ymax=max(boundary$y)
	if (abs(xmax-xmin) < 1e-6) return(NULL)
	if (abs(ymax-ymin) < 1e-6) return(NULL)

	x= sp@coords[,1]
	y = sp@coords[,2]
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
#' normalize data coordinates and border
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
#' x=runif(100, min=0, max=1)
#' y=runif(100, min=0.2, max=1.7)
#' range(x) # not [0,1]
#' tabData=data.frame(x=x,y=y)
#' bd=list(x=c(0,0,1,1,0), y=c(0.2,1.7,1.7,0.2,0.2))
#' res=datanorm(tabData,bd) 
#' apply(res$dataN,2,range)# 
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

  #normalize border
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
#' @details normalize data coordinates between 0 and 1 with different ratios for x and y
#' @param data frame with x and y components
#'
#' @return a normalized data frame
#'
#' @export
#'
#' @examples
#' nPoints=500
#' x=runif(nPoints, min=0, max=1)
#' y=runif(nPoints, min=0, max=1)
#' range(x) # not [0,1]
#' tabData=data.frame(x=x,y=y)
#' tabData=datanormXY(tabData) # x,y ranges are now [0,1]
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
#' x=runif(100, min=0, max=1)
#' y=runif(100, min=0.2, max=1.7)
#' range(x) # not [0,1]
#' tabData=data.frame(x=x,y=y)
#' bd=list(x=c(0,0,1,1,0), y=c(0.2,1.7,1.7,0.2,0.2))
#' res=datanormX(tabData,bd) 
#' apply(res$dataN,2,range)# x range is now [0,1], not y range
#' res$ratio # normalization ratio
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
#' @details normalize thresholds for small zone detection and no grow zone, considering mapo boundary
#' @param boundaryN normalized map boundary
#' @param minSize minimum size threshold
#' @param minSizeNG no grow size threshold 
#'
#' @return a list with components
#' \describe{
#' \item{minSize}{normalized minimum size threshold}
#' \item{minSizeNG}{normalized no grow size threshold}
#' }
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' resT=normSize(mapTest$boundary,0.012,0.001)#normalize thresholds relatively to map boundary area
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
#' genQseq
#'
#' @details description, a paragraph
#' @param qProb probability vector used to generate quantile values
#' @param K zoning object, as returned by the calNei function
#' @param map object returned by function genMap
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
#' checkContour
#'
#' @details check admissibility for contour line: surface >minSizeNG and refPoint close enough
#' @param contourSp SpatialPolygons corresponding to closed contour line
#' @param step grid resolution
#' @param refPoint referene point
#' @param minSizeNG zone area threshold under which a zone is not admissible
#'
#' @return Null if contour is not admissible or a list with components
#' \describe{
#' \item{contourSp}{SpatialPolygons corresponding to admissible contour}
#' \item{}{polyBuff}{SpatialPolygons corresponding to gBuffer around admissible contour}
#' }
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=contourAuto(list(),mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' pG=polyToSp2(Polygon(mapTest$boundary)) #SpatialPolygons corresponding to map boundary
#' plot(pG)
#' sp8 = contourToSpp(cL[[8]],0.1)$sp
#' refPoint = gCentroid(sp8)
#' resp=checkContour(sp8,mapTest$step,refPoint)
#' lines(resp$contourSp,col="red")
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
#' cleanSp
#'
#' @details removes from sp polygons that are too small (artefacts of gDifference)
#' @param sp SpatialPolygons
#' @param tol minimum area for removal
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
#' @param qProb probability vector used to generate quantile values
#' @param le level index
#' @param zk list of zonings
#' @param criterion list of criteria
#' @param cost list of costs
#' @param costL list of per label costs 
#' @param nz list of numbers of zones
#'
#' @return a list with components
#'\describe{
#' \item{ind}{index of last level zoning that has the higher criterion value}
#' \item{critList}{criterion value corresponding to best last level zoning}
#' \item{costlist}{cost value corresponding to best last level zoning}
#' \item{costLlist}{cost per label value corresponding to best last level zoning}
#' \item{nzList}{number of zones of best last level zoning}
#' \item{nq}{lenght of quantile vector}
#' }
#'
#' @export
#'
#' @examples
#'qProb=c(0.1,0.2);criti=correctionTree(qProb,map)
#'res=searchNODcrit1(qProb,criti)
#' # not run
searchNODcrit=function(qProb,le,zk,criterion,cost,costL,nz)
##################################################################
{
# zk: list of zonings
# le: list index
# crit: list of criteria
   critf=unlist(criterion[[le]])
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
#' @param qProb probability vector used to generate quantile values
#' @param crit result of call to correctionTree (with SAVE=TRUE)
#'
#' @return a list with components
#'\describe{
#' \item{ind}{index of last level zoning that has the higher criterion value}
#' \item{critList}{criterion value corresponding to best last level zoning}
#' \item{costlist}{cost value corresponding to best last level zoning}
#' \item{costLlist}{cost per label value corresponding to best last level zoning}
#' \item{nzList}{number of zones of best last level zoning}
#' }
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.1,0.2,0.4);criti=correctionTree(qProb,mapTest) # 2 zonings at last level
#' res=searchNODcrit1(qProb,criti)# best one is frist element of last level
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
#' # load test map with simulated data
#' data(mapTest)
#' # load zoning results from test file
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' normDistMat(resD$matDistanceCorr,2)
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
#' @details generates SpatialPolygons object corresponding to intersection of contour with boundary, must be within SpatialPolygons given in envel argument
#' @param co contour line
#' @param step map grid resolution
#' @param bd map boundary
#' @param envel envelope
#' @param disp info level (0-no info, 1- add lines to plot)
#'
#' @return a SpatialPolygons
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' pG=polyToSp2(Polygon(mapTest$boundary)) #SpatialPolygons corresponding to map boundary
#' cL=contourAuto(list(),mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' ps = interCB(cL[[8]],mapTest$step,mapTest$boundary,pG)#envelope is the whole map
#' plot(pG)
#' lines(ps,col="red")
#' # not run
interCB = function(co,step,bd=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),envel,disp=0)
##############################################################################
{
  #returns spatial polygon corresponding to intersection of contour with boundary
  #
	polygoneGlobal=SpatialPolygons(list(Polygons(list(Polygon(bd, hole = FALSE)), "1")))
	contourL = ContourLines2SLDF(list(co))
	polyBuff=gBuffer(contourL,width=0.0001*step)
	polyDiff=gDifference(polygoneGlobal,polyBuff)
  	recupPoly=separationPoly(polyDiff)

	ler=length(recupPoly)
	if(ler<2) return(NULL) # intersection=boundary -> degenerate contour

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
#' @details determine size of quantile  in result from correctionTree
#' @param critList component critList of result from correctionTree
#'
#' @return a vector with the size of quantile vectors for each zoning corresponding to critList
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=FALSE) # run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7, save only best last level results
#' getNq(criti$critList)
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
#' @details add contour lines to plot
#' @param map object returned by function genMap
#' @param val quantile value vector
#' @param col color parameter
#' @param super if TRUE add to existing plot lines coresponding to contour, if FALSE plot boundary and add lines
#'
#' @return void
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' addContour(mapTest,c(5,7),super=F)
#' # not run
addContour=function(map,val,col="blue",super=T)
##################################################################
{
for ( v in val)
  {
    cL=list()
    cL = contourAuto(cL,map$step,map$xsize,map$ysize,map$krigGrid,v,map$boundary)
    lc = length(cL)
    if (lc >0)
    {
	    sp=list()
	    if (!super) plot(map$boundary,type="l") # new plot
    	for (k in 1:lc)
    	   lines(cL[[k]],col=col) # add contour to existing plot
  	}
  }
  return()
}

##########################################################################
#' extractionPoly
#'
#' @details extract all elements from SpatialPolygons, holes and full polygons are handled equally 
#' @param polyTot SpatialPolygons
#'
#' @return a list of SpatialPolygons
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' ZK=initialZoning(qProb=c(0.2,0.4,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' extractionPoly(Z[[5]]) # returns 2 SpatialPolygons
#' # not run
extractionPoly=function(polyTot)
###########################################################################
{
  inPoly=polyTot@polygons[[1]]@Polygons
  listPolyExtract=list()
  for(i in (1:length(inPoly)))
  {
    listPolyExtract[[i]]=SpatialPolygons(list(Polygons(list(inPoly[[i]]), "1")))
  }


  return(listPolyExtract)
}

###################################
#' plotVario
#'
#' @details plot empirical variogram for model and data in map (raw data plus kriged data)
#' @param map object returned by function genMap
#' @param ylim range of y axis 
#'
#' @return a plot
#'
#' @export
#' @importFrom RandomFields RFempiricalvariogram
#'
#' @examples
#' data(mapTest)
#' plotVario(mapTest)
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
	  plot(empvarioK,model=modelGen,ylim=ylim,col="blue")
	}
  else	plot(empvario)
}

#########################
#' costLab
#'
#' @details description, a paragraph
#' @param K zoning object, as returned by the calNei function
#' @param map object returned by genMap function
#'
#' @return the sum of per label costs
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE) # run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7, save initial zoning and last level zonings
#' K=criti$zk[[1]][[1]] # initial zoning
#' costLab(K,mapTest) #identical to criti$costL[[1]][[1]]
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
#' @details compute overall mean and variance of all zones for each label plus sum of them for all labels
#' @param zlab  list with zone numbers for each zone label
#' @param listZonePoint list of indices of data points within zones, result of call to \code{\link{calNei}}
#' @param tabVal SpatialPointsDataFrame containing data values
#' @param surfVoronoi Surfaces of the Voronoi polygons corresponding to data pts
#'
#' @return a list with components
#' \describe{
#' \item{cL}{weighted (with Voronoi surfaces) average of per label variances}
#' \item{SigmaL2}vector of per label variances}
#' \item{SL}{vector of per label Voronoi surfaces}
#' \item{mL}{vector of weighted (with Voronoi surfaces) per label average values}
#' \item{voroLab}{vector of per label data}
#' }
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE) # run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7, save initial zoning and last level zonings
#' K=criti$zk[[2]][[1]]
#' uni=unique(K$lab)
#' zlab=sapply(uni,function(x){(1:length(K$lab))[K$lab==x]})
#' sig=SigmaL2(zlab,K$listZonePoint,mapTest$krigData,mapTest$krigSurfVoronoi)
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
#' @details compute overall mean of all zones for each label 
#' @param zlab  list with zone numbers for each zone label
#' @param listZonePoint list of indices of data points within zones, result of call to \code{\link{calNei}}
#' @param tabVal SpatialPointsDataFrame containing data values
#' @param surfVoronoi Surfaces of the Voronoi polygons corresponding to data pts
#'
#' @return a list with components
#' \describe{
#' \item{mL}{vector of weighted (with Voronoi surfaces) per label average values}
#' \item{SL}{vector of per label Voronoi surfaces}
#' }
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE) # run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7, save initial zoning and last level zonings
#' K=criti$zk[[2]][[1]]
#' uni=unique(K$lab)
#' zlab=sapply(uni,function(x){(1:length(K$lab))[K$lab==x]})
#' resL=meanL(zlab,K$listZonePoint,mapTest$krigData,mapTest$krigSurfVoronoi)
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

######################################
#' meansdSimu
#'
#' @details computes mean and standard deviation of a set of map simulations
#' @param vseed list of simulation seeds
#' @param krig type of kriging (1-variogram model-based, 2-inverse distance-based)
#'
#' @return a matrix with as many rows as simulations, and 4 columns, the first two columns give mean and standard deviation of generated raw data, the last two columns give mean and standard deviation of kriged data
#'
#' @export
#'
#' @examples
#' meansdSimu(c(1,2))
#' # not run
meansdSimu=function(vseed=NULL,krig=2)
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
    map=genMap(DataObj=NULL,seed=seed,krig=krig,disp=0)
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
rownames(mat)=paste(vseed)
return(mat)
}

######################################
#' meanVarsimu
#'
#' @details computes mean and standard deviation of a set of map simulations
#' @param map object generated by genMap
#'
#' @return a vector with  4 elements, the first two  give mean and standard deviation of generated raw data, the last two  give mean and standard deviation of kriged data
#'
#' @export
#'
#' @examples
#' meansdSimu(c(1,2))
#' # not run
meanvarSimu=function(map)
######################################
{
    v=map$rawData
    v=v@data[,1]
    kv=map$krigData
    kv=kv@data[,1]
    m=mean(v)
    km=mean(kv)
    sdd=sd(v)
    sdkd=sd(kv)
    vec=c(m,km,sdd,sdkd)
    names(vec)=c("raw mean","kriged mean","raw sd","kriged sd")
    return(vec)
}

###################
valZ=function(map,K)
#' valZ
#'
#' @details description, a paragraph
#' @param map map object returned by genMap function
#' @param K zoning object (such as returned by calNei function)
#'
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#'
#' @return a list with components
#' \describe{
#' \item{val}{list with vector of data values for each zone, zones are sorted by increasing mean values}
#' \item{ord}{order of zones sorted by increasing mean values}
#' }
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE) # run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7, save initial zoning and last level zonings
#' K=criti$zk[[2]][[1]]
#' valZ(mapTest,K)
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
