###############################################################################
#' dispZ
#'
#' @details plots a color image representation of values and zones
#' @param step grid resolution
#' @param matVal data frame of values
#' @param nbLvl number of contour lines to generate
#' @param zonePolygone zoning geometry (list of SpatialPolygons)
#' @param K zoning object, as returned by the calNei function
#' @param colBreaks if vector of length 1 number of color breaks, or else color breaks themselves
#' @param texMain main title
#' @param boundary map boundary (list with x and y values)
#' @param id logical value (if TRUE display zone ids on plot)
#' @param valQ quantile values to use for contour lines
#' @param palCol color palette
#' @param noXY if TRUE do not draw axes
#' @param iZ index of zone to outline in red
#' @param mu mu=1-only display zone number or id, mu=2-also display mean zone value
#' @param cex text size
#' @param ptz zone id location, if NULL automatically find the best locations
#'
#' @return an empty value
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics text
#' @importFrom graphics plot
#' @importFrom fields image.plot
#' @importFrom rgeos plot
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' dispZ(mapTest$step,mapTest$krigGrid)
dispZ=function(step,matVal,nbLvl=0,zonePolygone=NULL,K=NULL,colBreaks=0,texMain="",boundary=NULL,id=FALSE,valQ=NULL,
               palCol=colorRampPalette(c("brown","yellow")),noXY=FALSE,iZ=0,mu=1,cex=1,ptz=NULL)
###############################################################################
{
  # for PA figures palCol=colorRampPalette(c("brown","yellow"))
  # phi palCol=topo.colors(20)
  # mixed palCol=colorRampPalette(topo.colors(20))
  Z=zonePolygone
  nbPoly=length(Z) #0 if Z is NULL
  levelsList=NULL

  if(!is.null(valQ)){
    #generate isocontours with given valQ
    levelsList=valQ
  } else if(nbLvl!=0) {
    #else generate  nbLvls regularly spaced isocontours (if valQ not given)
	  levelsList=seq(min(matVal,na.rm=TRUE),max(matVal,na.rm=TRUE),length=nbLvl)
  }

  # if data, plot them as image
  if (!is.null(matVal)){
	  cn=as.numeric(colnames(matVal))
	  rn=as.numeric(rownames(matVal))
    #if color breaks are given in colBreaks vector
 	  if(length(colBreaks)!=1){
 	    image.plot(main=texMain,rn, cn,matVal,breaks=colBreaks,
                 col=palCol(length(colBreaks)-1),legend.shrink=0.4,legend.width=0.6)
  	}	else {
       image.plot(rn,cn,matVal,col=palCol(20),legend.shrink=0.4,legend.width=0.6,xlab="X",ylab="Y",cex.axis=1.5,cex.lab=1.5)
    }
  } #end case where data are available
  else # if no data just plot axes
  {
	  vs=getZsize(Z)
  	xsize=vs["x"]
  	ysize=vs["y"]
	  plot(0,0,xlim=c(0,xsize),ylim=c(0,ysize),type="n") #no data, only set axes
	}

  #add contour lines if there are data and if levels were given
  if(!is.null(levelsList))  {
    if (!is.null(matVal)) {
		  contour(rn,cn,matVal,levels = levelsList,add=TRUE)
    }
  }

  if(!is.null(boundary)) # add boundary if given
  {
    lines(boundary)
  }

  if(nbPoly!=0) # if there are zones
  {
    for (j in (1:nbPoly)) {
      #display zones
      lines(Z[[j]]@polygons[[1]]@Polygons[[1]]@coords,type='l')
	    #pointLabel=Z[[j]]@polygons[[1]]@Polygons[[1]]@labpt
	    # find good location to display (zone centroid may not be right)
      if(is.null(ptz))	pointLabel=findZCenter(Z,j)
	    else pointLabel=ptz[j,]
      # use zone id instead of zone number
      if(id) jid = getZoneId(Z[[j]])
      else jid=j
	    if(mu==2 & !is.null(K)){
	      vecmu=K$meanZone
     	  text(pointLabel[1],pointLabel[2],paste(jid,"(",round(vecmu[j],1),")",sep=""),cex=cex)
	    } else if(mu==1) text(pointLabel[1],pointLabel[2],jid,cex=cex)
    }

    if (iZ !=0){
      lines(Z[[iZ]]@polygons[[1]]@Polygons[[1]]@coords,type='l',col='red')
      lines(boundary)
    }
  }
  return()
}

#########################################################
#'
#' dispZmap
#' @details plots a color representation of values and zones
#' @param map map object returned by genMap function
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param qProb quantile associated probability vector
#' @param valbp values used for boxplots
#' @param scale field scale
#' @param lev number of color levels
#' @param palCol color palette
#' @param legend.width relative width of legend
#' @param parG graphics parameters (result of call to par)
#' @param ptz zone id location, if NULL automatically find the best locations
#'
#' @return an empty value
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics text boxplot title polygon
#' @importFrom sp plot
#' @importFrom fields image.plot
#'
#' @export
#'
#' @examples
#' seed=80
#' data(mapTest)
#' ZK=initialZoning(c(0.5,0.7),mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' #order zone ids by attribute mean value
#' ord=order(K$meanZone)
#'   Z=orderZ(Z,ord)
#'   plotZ(Z,id=TRUE)
dispZmap=function(map,Z=NULL,qProb=NULL,valbp=NULL,scale=NULL,lev=20,palCol=colorRampPalette(c("brown","yellow")),legend.width=1,parG=NULL,ptz=NULL)
#########################################################
{
  step=map$step
  xsize=map$xsize
  ysize=map$ysize
  matVal=map$krigGrid
  bd=map$boundary
  xlim=c(0,xsize);ylim=c(0,ysize)
  cn=as.numeric(colnames(matVal))
  rn=as.numeric(rownames(matVal))

  if(is.null(parG))
  {
    par(fig=c(0,1,0.3,1))#c(x1,x2,y1,y2) bottom=0,0 - top=1,1
    par(mar=c(0,0,5,0))#c(bottom, left, top, right))
  }

  image.plot(rn, cn,matVal,col=palCol(lev),asp=1,xlim=xlim,ylim=ylim,xlab="",ylab="",xaxt="n",yaxt="n")
  image.plot( zlim=range(matVal,na.rm=TRUE), nlevel=20,legend.only=TRUE, vertical=FALSE,col=palCol(20))
  if (!is.null(qProb)) title(paste("Best zoning for nL=",length(qProb)+1,"\nm=[",paste(qProb,collapse=","),"]",sep=""),cex.main=1.5)

  dec=0.03
  decy=0.05
  xn=-0.2
  yn=1.8
  p1=polygon(list(x=c(xn,xn+dec,xn+dec,xn),y=c(yn,yn+decy,yn+2*decy,yn)))
  p2=polygon(list(x=c(xn+dec,xn+dec,xn+2*dec,xn+dec),y=c(yn+decy,yn+2*decy,yn,yn+decy)),density=100,angle=45,col="black")
  text(xn+dec,yn-decy,"N",font=2,cex=1)

  # add scale
  if(!is.null(scale))
  {
    met=100/scale
    ys=0.78
    xs=0.8
    lines(c(xs,xs+met),c(ys,ys))
    tau=0.02
    lines(c(xs,xs),c(ys-tau,ys+tau))
    lines(c(xs+met,xs+met),c(ys-tau,ys+tau))
    text(xs,ys-3*tau,"0")
    text(xs+met,ys-3*tau,"100m")
  }
  # add zones
  par(fig=c(0,1,0.3,1))#c(x1,x2,y1,y2) bottom=0,0 - top=1,1
  par(mar=c(0,0,5,0))#c(bottom, left, top, right))

  if(!is.null(Z))
  {
    for (k in 1:length(Z))
    {
      lines(getCoords(Z[[k]]))
      if(is.null(ptz))
      #pointLabel=Z[[k]]@polygons[[1]]@Polygons[[1]]@labpt
      pointLabel=findZCenter(Z,k)
      else
      pointLabel=ptz[k,]

      jid=getZoneId(Z[[k]])
      text(pointLabel[1],pointLabel[2],jid,cex=1.5)
    }
  }

  if(!is.null(valbp))
  {
    if(is.null(parG))
    {
      par(fig=c(0.2,0.9,0,0.3),new=TRUE)
      par(mar=c(2,0,3,0))#c(bottom, left, top, right)
    }
    boxplot(valbp,cex.axis=1,col=palCol(length(valbp)))
    title("Distribution of values within zones",cex.main=1)
  }

  return()
}

##################################################################
#' pointsSp
#'
#' @details description, a paragraph
#' @param sp SpatialPolygons object
#' @param k polygon number
#' @param col color
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' pointsSp(Z[[1]])
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
#' @return NULL
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' linesSp(Z[[4]])
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
#' @return NULL
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotSp(Z[[1]],xlim=c(0,1),ylim=c(0,1))
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
#' @details wrapper function for dispZ
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param map map object returned by genMap function
#' @param id logical value, if TRUE display zone ids, if FALSE display zone numbers
#' @param noXY logical value, if TRUE do not display x and y axes
#' @param palCol argument of colorRampPalette
#'
#' @return an empty value
#' @importFrom grDevices topo.colors
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z,mapTest)
plotZ = function(Z,map=NULL,id=FALSE,noXY=FALSE,palCol=colorRampPalette(topo.colors(20)))
##################################################################
{

	if (!is.null(map))
   	   dispZ(map$step,matVal=map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0,noXY=noXY,palCol=palCol)
	else
          {
	  dispZ(map$step,matVal=NULL, nbLvl=0, zonePolygone=Z,id=id,palCol=palCol)
	  }

}

##################################################################
#' plotListC
#'
#' @details add contour lines to a plot
#' @param cL list of contour lines
#' @param col color to use
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




