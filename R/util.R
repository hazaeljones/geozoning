###################################################################
#' getPtsZone
#'
#' @details description, a paragraph
#' @param ptsp xxxx
#' @param zone xxxx,
#'
#' @return a map in a list
#' \describe{
#' \item{pts}{ptsub}
#' \item{mask}{mask}
#' }
#'
#' @export
#'
#' @examples
#' # not run
getPtsZone=function(ptsp,zone)
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
#' @param map xxxx
#' @param zone xxxx,
#' @param w default NULL
#'
#' @return a in a list
#' \describe{
#' \item{mean}{mean}
#' \item{var}{var}
#' }
#'
#' @export
#'
#' @examples
#' # not run
MeanVarWPts=function(map,zone,w=NULL)
#####################################################################
{
  ptsp=map$krigData

  m=numeric()
  v=numeric()
  res=getPtsZone(ptsp,zone)
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
#' @param reslm an lm object
#'
#' @return a numeric
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
#' @param ptsp xxxx
#' @param Z a zone?
#'
#' @return a lm object
#'
#' @export
#'
#' @examples
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
#' @param sp xxxx
#' @param k xxxxxx
#'
#' @return some coordinates
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param sp xxxx
#'
#' @return a list
#'
#' @export
#' @importFrom sp SpatialLines Lines
#'
#' @examples
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
#' @param Z xxxx
#' @param indiceCourant xxxx
#' @param indiceZP xxxx
#' @param disp xxxx
#'
#' @return a numeric?
#'
#' @export
#' @importFrom sp SpatialPoints SpatialPointsDataFrame
#'
#' @examples
#' # not run
getClosePt=function(Z,indiceCourant,indiceZP,disp=FALSE)
##################################################################
{
#coordonnées des points de la petite zone et la  zoneProche

  a=SpatialPoints(Z[[indiceCourant]]@polygons[[1]]@Polygons[[1]]@coords)
  b=SpatialPoints(Z[[indiceZP]]@polygons[[1]]@Polygons[[1]]@coords)

  #on recupere les distances entre chaque paire de points
  #on pourra peut etre utiliser un gSimplify si cela prend trop de temps
  Fdist= list()
  for (i in 1:length(a))
  {
    Fdist[[i]]<-gDistance(a[i,],b,byid=TRUE)
  }

  #on recupere le point de la zoneProche le plus proche de la petite zone
  min.dist <- unlist(lapply(Fdist, FUN=function(x) which(x == min(x))[1]))
  PolyDist <- unlist(lapply(Fdist, FUN=function(x) min(x)[1]))

  pProche = min.dist[which.min(PolyDist)]
  if (disp)
  {
  print(b[pProche])
  }
  return (b[pProche])
}

##################################################################
#' cadArea
#'
#' @details description, a paragraph
#' @param cad xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' # not run
cadArea = function(cad)
##################################################################
  {
	rn=as.numeric(rownames(cad))
	cn=as.numeric(colnames(cad))
	xrn=range(rn)
	yrn=range(cn)
	cadSurf=(max(xrn)-min(xrn))*(max(yrn)-min(yrn))
	return(cadSurf)
  }

##################################################################
#' contourArea
#'
#' @details description, a paragraph
#' @param contour1 xxxx
#'
#' @return a numeric?
#'
#' @export
#' @importFrom sp SpatialPolygons SpatialPointsDataFrame Polygons Polygon
#' @importFrom maptools ContourLines2SLDF
#'
#' @examples
#' # not run
contourArea=function(contour1)
##################################################################
{
	contour = ContourLines2SLDF(list(contour1))
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
#' @details description, a paragraph
#' @param listContour xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' # not run
listContourArea=function(listContour)
##################################################################
{
	le=length(listeContour)
	if (le ==0) return(NA)
	surface=list()
	for (k in 1:length(listeContour))
	{
		surface[[k]] = contourArea(listeContour[[k]])
	}
	return(surface)
}

##################################################################
#' contourToSpp
#'
#' @details description, a paragraph
#' @param contour1 xxxx
#' @param step xxxx
#'
#' @return a numeric?
#'
#' @export
#' @importFrom rgeos gBuffer gArea gCentroid gContains gConvexHull gDifference gDistance gIntersects gWithin
#'
#' @examples
#' # not run
contourToSpp=function(contour1,step)
##################################################################
{
	contour = ContourLines2SLDF(list(contour1))
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param ind xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' # not run
nPolyZone=function(Z,ind)
##################################################################
{
	return(length(Z[[ind]]@polygons[[1]]@Polygons))
}

##################################################################
#' nPolySp
#'
#' @details description, a paragraph
#' @param sp xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' # not run
nPolySp =function(sp)
##################################################################
{
	return(length(sp@polygons[[1]]@Polygons))
}

##################################################################
#' holeSp
#'
#' @details description, a paragraph
#' @param sp xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param iZ xxxx
#' @param k xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' # not run
maxDistZone=function(Z,iZ,k)
##################################################################
{
	return(max(dist(Z[[iZ]]@polygons[[1]]@Polygons[[k]]@coords)))
	#return(gArea(Z[[iZ]]))
}

##################################################################
#' maxDistSP
#'
#' @details description, a paragraph
#' @param sp xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
#' # not run
maxDistSP=function(sp)
##################################################################
{
	return(max(dist(sp@polygons[[1]]@Polygons[[1]]@coords)))
}

##################################################################
#' getPoly
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param iZ xxxx
#' @param k xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param sp xxxx
#' @param k xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param iZ xxxx
#' @param k xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
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
#' polyToSp
#'
#' @details description, a paragraph
#' @param p xxxx
#'
#' @return a numeric?
#'
#' @export
#'
#' @examples
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
#' @details description, a paragraph
#' @param data xxxx
#' @param bd xxxx
#'
#' @return a list
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
#' @details description, a paragraph
#' @param data xxxx
#'
#' @return a list
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
#' normalize polygons in zoning
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param xmin xxxx
#' @param xmax xxxx
#' @param ymin xxxx
#' @param ymax xxxx
#'
#' @return a dataframe?
#'
#' @export
#'
#' @examples
#' # not run
ZnormXY = function(Z,xmin,xmax,ymin,ymax)
##################################################################
{
# normalize polygons in zoning
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
#' calcDCrit
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param map xxxx
#' @param optiCrit xxxx
#' @param pErr xxxx, default 0.9
#' @param simplitol xxxx, default 0.001
#'
#' @return a list
#'
#' @export
#'
#' @examples
#' # not run
calcDCrit=function(Z,map,optiCrit,pErr=0.9,simplitol=1e-3)
##################################################################
{
	resZ=calNei(Z,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol=simplitol)
  le = length(resZ$zonePolygone)
	if (le <2) return(list(resD=0,resCrit=0))

  resDist=calDistance(typedist=1,map$krigData,resZ$listZonePoint,resZ$zoneN,map$krigSurfVoronoi,resZ$meanZone,pErr=pErr)
  resCrit=calCrit(resDist$matDistanceCorr,resZ$zoneNModif,optiCrit)

	return(list(resD=resDist,resCrit=resCrit))
}

# ##################################################################
# ATTENTION, CETTE FONCTION EST EN DOUBLE!!! cf lignes 981
# JE METS CELLE-CI EN COMMENTAIRES VU QUE DANS UN SOURCE C'4'EST LA 2e
# FONCTION QUI EST STOCKEE...
# #' linesSp
# #'
# #' @details description, a paragraph
# #' @param contourSp xxxx
# #' @param k xxxx
# #'
# #' @return a list
# #'
# #' @export
# #'
# #' @examples
# #' # not run
# linesSp=function(contourSp,k)
# ##################################################################
# {
# 	co=contourSp[[k]]
# 	lines(co@lines[[1]]@Lines[[1]]@coords,col="blue")
# }

##################################################################
#' plotCad
#'
#' @details description, a paragraph
#' @param cad xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
plotCad=function(cad)
##################################################################
{
	rn=as.numeric(rownames(cad))
	cn=as.numeric(colnames(cad))
	xmin=min(rn)
	xmax=max(rn)

	ymin=min(cn)
	ymax=max(cn)

	x=c(xmin,xmin,xmax,xmax,xmin)
	y=c(ymin,ymax,ymax,ymin,ymin)

	lines(x,y,col="blue")
}

##################################################################
#' calcCad
#'
#' @details description, a paragraph
#' @param cad xxxx
#'
#' @return a SpatialPolygons object
#'
#' @export
#'
#' @examples
#' # not run
calcCad=function(cad)
##################################################################
{
	rn=as.numeric(rownames(cad))
	cn=as.numeric(colnames(cad))
	xmin=min(rn)
	xmax=max(rn)

	ymin=min(cn)
	ymax=max(cn)

	x=c(xmin,xmin,xmax,xmax,xmin)
	y=c(ymin,ymax,ymax,ymin,ymin)
	polys=Polygon(cbind(x,y))
	lpolys=Polygons(list(polys),"id")
	return(SpatialPolygons(list(lpolys)))
}

##################################################################
#' plotListeC
#'
#' @details description, a paragraph
#' @param listeC xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
plotListeC = function(listeC)
##################################################################
{
	le=length(listeC)
	if (le >0)
	   for (i in 1:le)
	   {
		ci=listeC[[i]]
		lines(ci,col="red")
	   }

}


##################################################################
#' genQseq0
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param K xxxx
#' @param map xxxx
#' @param i1 xxxx
#' @param i2 xxxx
#' @param leq xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
genQseq0 = function(qProb,K,map,i1,i2,leq,disp=0)
##################################################################
{
  labs = c(K$lab[i1],K$lab[i2])
  #labs = unlist(labs)

  ind = unique(labs)
  le = length(qProb)
  if (le==1) ind=1
  Qseq=c()

  for (ii in 1:length(ind))
  {
    	#
	k = min(length(qProb),ind[ii])
	prob1 = qProb[k]
  	q1= quantile(map$krigGrid,na.rm=TRUE,prob=prob1)
    	q2 = quantile(map$krigGrid,na.rm=TRUE,prob= min(prob1 + 0.1, 0.99))
  	Qseq = c(Qseq,(seq(q1,q2,length.out=leq)))
  	q3 = quantile(map$krigGrid,na.rm=TRUE,prob= max(prob1 - 0.1, 0.01))
	Qseq = c(Qseq,(seq(q1,q3,length.out=leq)))

    }
  Qseq=sort(unique(Qseq))
  if (disp>0)
     {
	print("Qseq=")
  	print(Qseq)
	}
  return(Qseq)
}

##################################################################
#' genQseq
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param K xxxx
#' @param map xxxx
#' @param i1 current zone index
#' @param i2 englobing zone index
#' @param LEQ xxxx
#' @param MAXP xxxx
#' @param disp xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
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
#' @param sp xxxx
#' @param k xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
pointsSp = function(sp,k=1)
##################################################################
{
	p = sp@polygons[[1]]@Polygons[[k]]
	points(p@coords,col="red")
}

##################################################################
#' linesSp
#'
#' @details description, a paragraph
#' @param sp xxxx
#' @param k xxxx
#' @param lty xxxx
#' @param col xxxx
#' @param lwd xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
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
#' @param sp xxxx
#' @param k xxxx
#' @param xlim xxxx
#' @param ylim xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
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
#' plotSp
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param map xxxx
#' @param id xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
plotZ = function(Z,map=NULL,id=FALSE)
##################################################################
{
	# IS 16/05/2017: comments for x11 device
  #x11()
	if (!is.null(map))
   	   dispZ(map$step,matVal=map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
	else
  {
	  dispZ(map$step,matVal=NULL, nbLvl=0, zonePolygone=Z,id=id)
	}
}

##################################################################
#' plotzf
#'
#' @details description, a paragraph Attention, voir la sortie de cette fonction en package!!!
#' @param qProb xxxx
#' @param map xxxx
#' @param pdf1 xxxx
#' @param FULL xxxx
#' @param i xxxx
#' @param vecj xxxx
#' @param data xxxx
#'
#' @return a plot
#'
#' @export
#' @importFrom grDevices dev.off pdf x11
#' @importFrom graphics lines title
#'
#' @examples
#' # not run
plotzf = function(qProb,map,pdf1=NULL,FULL=F,i,vecj,data=NULL)
##################################################################
{
  # following funcCritereCNO
  valRef=quantile(map$krigGrid,na.rm=TRUE,prob=qProb)
  if (missing(i) & missing(vecj))
  {
    best=searchNODcrit(qProb,length(zk),zk,critere,cost,costL,nz) #global variables
    i=length(zf)
    lq=length(qProb)
    lev=paste("q",lq,sep="")
    vecj=best[["ind"]][[lev]]
    vecj=vecj[1]
  }

  if(!is.null(pdf1)) pdf("fig.pdf")
  # IS 19/05/2017: add comment for x11
  #else x11()

  for (j in vecj)
  {
    if(length(zf[[i]][[j]])>9) FULL=TRUE
    dispZ(map$step,map$krigGrid,zonePolygone=zf[[i]][[j]],K=zk[[i]][[j]],boundary=map$boundary,nbLvl=0,id=FALSE)
    title(paste(data, " valRef=[",toString(round(valRef,2)),"]   critere=",round(critere[[i]][[j]],2),sep=""))
  }

  if(!is.null(pdf1))
  {
	  dev.off() #close fig.pdf file
	  system("cp MOD.matDistance1.tex mdist.tex")
	  if(!FULL) write("\\begin{minipage}[l]{0.5\\linewidth}",file="mdist.tex",append=TRUE)
    #insertion du graphique à partir d'un pdf,fin de la minipage
  	write(paste("\\includegraphics[width=\\linewidth]{","fig.pdf}",sep=""),file="mdist.tex",append=TRUE)
	  if(!FULL)
		{
		  write("\\end{minipage}\\hfill",file="mdist.tex",append=TRUE)
      #nouvelle minipage avec la matrice des distances 1 à coté du graphique
  	  write("\\begin{minipage}[r]{0.5\\linewidth}",file="mdist.tex",append=TRUE)
	    write("\\tiny",file="mdist.tex",append=TRUE)
		}
	  write(print(xtable(mdist[[i]][[j]]),table.placement=NULL,latex.environment="",floating.environment="center"),
	        file="mdist.tex",append=TRUE)

  	if(!FULL) write("\\end{minipage}",file="mdist.tex",append=TRUE)
  	write("\\end{document}",file="mdist.tex",append=TRUE)

	  system("pdflatex mdist")
	  system(paste("mv mdist.pdf",pdf1))
  }
  return(c(i,j))
}

##################################################################
#' calczf
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param map xxxx
#' @param optiCrit xxxx
#' @param minSize xxxx
#' @param minSizeNG xxxx
#' @param disp xxxx
#' @param pdf1 xxxx
#' @param FULL xxxx
#' @param data xxxx
#'
#' @return a list
#'
#' @export
#'
#' @examples
#' # not run
calczf = function(qProb,map,optiCrit,minSize,minSizeNG,disp=0,pdf1=NULL,FULL=FALSE,data=NULL)
##################################################################
{
  crit2=correctionTree(qProb,map,pErr=0.9,optiCrit=2,minSize=minSize,minSizeNG=minSizeNg,distIsoZ=distIsoZ,
                       simplitol=simplitol,LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=0,SAVE=TRUE,ONE=FALSE)
  ij=plotzf(qProb,map,pdf1=pdf1,FULL=FULL,data=data)
  return(list(crit=crit2,ij=ij))
}

##################################################################
#' checkContour
#'
#' @details description, a paragraph
#' @param contourSp xxxx
#' @param step xxxx
#' @param refPoint xxxx
#' @param minSizeNG xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
checkContour = function(contourSp,step,refPoint,minSizeNG)
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param sp xxxx
#'
#' @return a message?
#'
#' @export
#'
#' @examples
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
#' normD
#'
#' @details description, a paragraph
#' @param DataObj xxxx
#' @param boundary xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
normD = function(DataObj,boundary)
##################################################################
{

	xmin=min(boundary$x)
	xmax=max(boundary$x)
	ymin=min(boundary$y)
	ymax=max(boundary$y)
	if (abs(xmax-xmin) < 1e-6) return(NULL)
	if (abs(ymax-ymin) < 1e-6) return(NULL)

	x=DataObj[,1]
	y=DataObj[,2]
	x = (x-xmin)/(xmax-xmin)
	y = (y-ymin)/(ymax-ymin)
	DataObj[,1]=x
	DataObj[,2]=y

	return(DataObj)

}

##################################################################
#' denormD
#'
#' @details description, a paragraph
#' @param DataObj xxxx
#' @param boundary xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
denormD = function(DataObj,boundary)
##################################################################
{
	xmin=min(boundary$x)
	xmax=max(boundary$x)
	ymin=min(boundary$y)
	ymax=max(boundary$y)
	x=DataObj[,1]
	y=DataObj[,2]
	x=(x*(xmax-xmin))+xmin
	y=(y*(ymax-ymin))+ymin
	DataObj[,1]=x
	DataObj[,2]=y

	return(DataObj)
}

##################################################################
#' createHoles
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
#' @details description, a paragraph
#' @param zonePrinc xxxx
#' @param zoneSuppr xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
moveHoles = function(zonePrinc,zoneSuppr)
##################################################################
{
	Zone1=zonePrinc

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
#' @details description, a paragraph
#' @param sp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
cleanSp = function(sp)
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
	  if (area >= 1e-5)
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
#' interContourZ
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param contour xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
interContourZ=function(Z,contour)
##################################################################
{
# contour vient de contourAuto (fonction contourLines)
# le transformer en spatialLines
        contourSp=ContourLines2SLDF(list(contour))
	listePolygoneTemp=list()
	polyBuff=gBuffer(contourSp,width=1e-3)
	listId=list()
	for (j in (1:length(Z)))
      {

          #On tente d'intersecter la ligne de niveau avec ce polygone
          # Si l'intersection existe bien, on renvoie les deux polygones qui en résultent
        if(gIntersects(Z[[j]],contourSp))
          {
            polyDiff=gDifference(Z[[j]],polyBuff)

            # - Fonction qui récupère dans une structure le premier polygone non trou,avec tous ses voisins,puis le deuxième,puis le troisième,etc.
            recupPoly=separationPoly(polyDiff)

            listePolygoneTemp[[2*(j)-1]]=recupPoly[[1]]
            #continuer execution mais générer un warning si jamais il y a plus de 3 polygones non trou

            if (length(recupPoly) > 1){
              listePolygoneTemp[[2*(j)]]=recupPoly[[2]]
            }

            #if(length(recupPoly)>2)
            #{
             # warning("un isocontour définit 3 zones ou plus: incohérent (dans foncMapAlea-zoneGeneration)")
            #}
          }
          else #Sinon on renvoie original et un NULL
          {
            listePolygoneTemp[[2*(j)-1]]=Z[[j]]
            listePolygoneTemp[[2*(j)]]=NULL
          }
      } # fin boucle sur polygones

      # On ne garde que les éléments non nuls
      iNonNuls=grep(FALSE,lapply(listePolygoneTemp,is.null))
      if (length(iNonNuls)==0) return(NULL)
      Z=listePolygoneTemp[iNonNuls]
}

##################################################################
#' ptInZone
#'
#' @details description, a paragraph
#' @param zone xxxx
#' @param pts xxxx
#' @param numpt xxxx
#'
#' @return a ?
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' # not run
ptInZone=function(zone,pts,numpt)
##################################################################
{
	PointsSpatiaux=SpatialPoints(pts)
	logicalPoly=point.in.polygon(PointsSpatiaux$x,PointsSpatiaux$y,zone@polygons[[1]]@Polygons[[1]]@coords[,1],zone@polygons[[1]]@Polygons[[1]]@coords[,2])
    if(length(zone@polygons[[1]]@Polygons)>1)
    {
      #on prend en compte les points qui sont dans les trous
      for(k in 2:length(zone@polygons[[1]]@Polygons))
      {
        logicalPoly=logicalPoly-point.in.polygon(PointsSpatiaux$x,PointsSpatiaux$y,zone@polygons[[1]]@Polygons[[k]]@coords[,1],zone@polygons[[1]]@Polygons[[k]]@coords[,2])
      }
    }
    return(logicalPoly[numpt])
}

##################################################################
#' ptInZone
#'
#' @details description, a paragraph
#' @param sp xxxx
#' @param pts xxxx
#' @param hole xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
ptsInSp=function(sp,pts,hole=FALSE)
##################################################################
{
	PointsSpatiaux=SpatialPoints(pts)
	logicalPoly=point.in.polygon(PointsSpatiaux$x,PointsSpatiaux$y,sp@polygons[[1]]@Polygons[[1]]@coords[,1],sp@polygons[[1]]@Polygons[[1]]@coords[,2])
    if(hole && (length(sp@polygons[[1]]@Polygons)>1))
    {
      #on prend en compte les points qui sont dans les trous
      for(k in 2:length(sp@polygons[[1]]@Polygons))
      {
        logicalPoly=logicalPoly-point.in.polygon(PointsSpatiaux$x,PointsSpatiaux$y,sp@polygons[[1]]@Polygons[[k]]@coords[,1],sp@polygons[[1]]@Polygons[[k]]@coords[,2])
      }
    }

    return(pts[logicalPoly!=0,])
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
#' @details description, a paragraph
#' @param matDistanceCorr xxxx
#' @param optiCrit xxxx
#'
#' @return a ?
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
#' @param iC xxxx
#' @param Z xxxx
#' @param K xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
getClosestZone=function(iC,Z,K)
##################################################################
{
	imin=0
	ni=1:length(Z)
	# exclude current zone
	ni=ni[ni != iC]

	# exclude neighbor zones
	Ns=getNs(K,iC)
  	listeV=grep(TRUE, Ns)
	for (i in listeV) ni=ni[ni != i]

	# exclude englobing zone
	iE = detZoneEng(iC,Z,K)
	#
	if(iE != 0) ni = ni[ni !=iE]

	# exclude included zones
	ir=NULL
	for (i in ni)
          {
	  gb = gBuffer(gConvexHull(Z[[iC]]),byid=TRUE,width=1e-3)
	  if(gContains(gb,Z[[i]])) ir=c(ir,i)
	  }
	for (i in ir) ni =ni[ni != i]

	# compute distances to remaining zones
	d0 = 1
	for (i in ni)
          {
		d=gDistance(Z[[iC]],Z[[i]])
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
  listeContour=list()
  listeContour = contourAuto(listeContour,map$step,map$krigGrid,vRef,map$boundary)

  # intersect contour with boundary
  # and transform contours into sps
  sp=list()
  k=0
  for (cont in listeContour)
  {
     k=k+1
     ps = interCB(cont,map$step,envel=envel)
     # returns NULL is contour is degenerate (single point)
     if(!is.null(ps)) sp[[k]]=ps else k=k-1
  }

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
#' @details description, a paragraph
#' @param listeContour xxxx
#' @param col xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
linesC = function(listeContour,col="blue")
##################################################################
{
for (i in 1:length(listeContour))
{
	lines(listeContour[[i]],col=col)
	}
return()
}

##################################################################
#' interCB
#'
#' @details description, a paragraph
#' @param contour1 xxxx
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
interCB = function(contour1,step,bd=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),envel,disp=0)
##################################################################
{
  #returns spatial polygon corresponding to intersection of contour with boundary
  #
	polygoneGlobal=SpatialPolygons(list(Polygons(list(Polygon(bd, hole = FALSE)), "1")))
	contourL = ContourLines2SLDF(list(contour1))
	polyBuff=gBuffer(contourL,width=0.0001*(1/step))
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
    listeContour = contourAuto(listeContour,map$step,map$krigGrid,v,map$boundary)
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
