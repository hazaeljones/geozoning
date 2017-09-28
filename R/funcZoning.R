###############################################################################
#' zoneGeneration
#'
#' @details Generates zones from map data using quantile values associated to given probabilities
#' @param map object returned by function genMap
#' @param qProb probability vector used to generate quantile values
#' @param GridData logical value indicating if data are already on a regular grid (no kriging in that case)
#'
#' @return a list of zones, each zone is a SpatialPolygons
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' Z=zoneGeneration(mapTest)
#' # not run
zoneGeneration=function(map,qProb=c(0.25,0.75),GridData=FALSE)
################################################################################
# data in map argument

{
  matVal=map$krigGrid
  boundary=map$boundary
  step=map$step
   xsize=map$xsize
  ysize=map$ysize
  
  # Define field boundaries
  if(!is.null(boundary))
  {
    frame=boundary
  }else {frame=structure(list(x = c(0,0,1,1,0), y = c(0,1,1,0,0)))}

  # Global Polygon
  pG=polyToSp2(Polygon(frame))

  # Init lists
  Z=list(pG)
  cL=list()
  cLSp=list()
  # quantile values for data
  valQuant=quantile(matVal,na.rm=TRUE,prob=qProb)

  for(i in (1:length(valQuant)))
  {
     cL=contourAuto(cL,step,xsize,ysize,matVal,vRef=valQuant[i],frame,GridData)
     if(is.null(cL)) return(NULL)
  }

  # For each isocontour
  for (jContour in (1:length(cL)))
    {
      # convert to spatial lines
      cLSp[[jContour]] <- ContourLines2SLDF(list(cL[[jContour]]))
      #  level attribute conserved
      # starting from frame, intersect with each isocontour
      # and iterate

      listPolyTmp=list()

      # define buffer: contour and another line around contour
      polyBuff=gBuffer(cLSp[[jContour]],width=0.0001*step)
      # Loop on already defined zone
      for (j in (1:length(Z)))
      {
        #intersect current zone with buffer -> creates several polygons
	#because polygons may have holes inside
        if(gIntersects(Z[[j]],cLSp[[jContour]]))
          {
            polyDiff=gDifference(Z[[j]],polyBuff)
            # separate Polygons into holes and non holes
            recupPoly=separationPoly(polyDiff)
            listPolyTmp[[2*(j)-1]]=recupPoly[[1]]
            if (length(recupPoly) > 1){
              listPolyTmp[[2*(j)]]=recupPoly[[2]]
            }
          }
          else #no intersection with current zone
          {
            listPolyTmp[[2*(j)-1]]=Z[[j]]
            listPolyTmp[[2*(j)]]=NULL
          }
      } # end loop on current zones

      # keep only non null elements
      iNonNulls=grep(FALSE,lapply(listPolyTmp,is.null))
      Z=listPolyTmp[iNonNulls]
    }
  # end loop on contours
  return(Z)
}

##############################################################################
#' contourAuto
#'
#' @details builds contout Lines qith the quantile vector given in argument and closes them with the map border
#' @param cL empty or existing list of contour lines
#' @param step grid step as returned by calStep
#' @param xsize size of map along x-axis
#' @param ysize size of map along y-axis
#' @param matVal dataframe with data values organized into a grid
#' @param vRef quantile vector 
#' @param boundary list, contains x and y dy on a regular grid
#' @param GridData logical value indicating if data are already on a regular grid 
#'
#' @return a list of contour lines
#' @importFrom grDevices contourLines
#' @importFrom sp Line coordinates
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' plot(mapTest$boundary,type="l",col="red")
#' geozoning:::linesC(cL)
#' # not run
contourAuto=function(cL,step,xsize,ysize,matVal,vRef,boundary,GridData=FALSE)
##############################################################################
{
  # find isocontours (level attribute)
  if(!GridData)
	cLplus=contourLines(seq(step, xsize-step, by=step),seq(step, ysize-step, by=step),matVal, levels = vRef)
 else

	cLplus=contourLines(z=matVal, levels = vRef)
  # merge with previous isocontours (other level)
  cL=c(cL,cLplus)
  if(length(cL)==0) return(NULL) # no contour
  
  # transform boundary into spatial object
  bdSP =SpatialPoints(boundary)
  superL=superLines(boundary)
  #for each isocontour, extend it to frame
  for(jContour in (1:length(cL)))
  	{
  	#test if contour is closed or not - if not, close it
	iso=cL[[jContour]]
    	np=length(iso$x) #number of pts in contour line
     	#if(iso$x[1]!= iso$x[np] || iso$y[1] != iso$y[np])#not exact sometimes
	# correct small discrepancies
	dx=abs(iso$x[1]-iso$x[np])
	if(dx<=1e-3)
		{
		iso$x[np]=iso$x[1]
		dx=0.0
		}
	dy=abs(iso$y[1]-iso$y[np])
	if(dy<=1e-3)
		{
		iso$y[np]=iso$y[1]
		dy=0.0
		}
	if (dx >1e-3 || dy >1e-3)
    	  {
	  # if contour is not closed, add projection on frame
     	  ######################################################
	  cL[[jContour]]=extensionLine(iso,step,bdSP,superL)
	  ######################################################
    	  }
 	 }

  return(cL)
}

#########################################################################
#' zoneAssign
#'
#' @details assigns points to zones
#' @param tab data frame with data values
#' @param Z zoning object
#'
#' @importFrom sp coordinates
#'
#' @return a list of data points within each zone
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' ZK=initialZoning(qProb=c(0.4,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' listZpts=zoneAssign(mapTest$krigData,Z)
#' #identical to ZK$resZ$listZonePoint
#' listZptsRaw=zoneAssign(mapTest$rawData,Z)
#' plotZ(Z)
#' points(mapTest$rawData[listZptsRaw[[1]],],col="blue") # add raw data for zone 1
#' # not run
zoneAssign=function(tab,Z)
#########################################################################
{
  #nb zones
  nbZ=length(Z)

  #n data points
  n=nrow(tab)
  zone=rep(0,n)

  #nZone[i]=1 or more (nb of zones to which pt i belongs)
  nZone=zone
  inZone=list()
  pts=as.data.frame(coordinates(tab))

  for(i in 1:nbZ)
  {
    inZone[[i]]=0
    co=getCoords(Z[[i]])
    #point.in.polygon returns 0 (exterior),1 (interior) or 2 (vertex)
    pti=point.in.polygon(pts$x,pts$y,co[,1],co[,2])
    pti[pti>1]=1
    inZone[[i]]=pti

    np=nPolySp(Z[[i]])

    if(np>1)
	{
      	# do not consider pts within holes
      	for(k in 2:np)
      	      {
		coh=getCoords(Z[[i]],k)
		ptk=point.in.polygon(pts$x,pts$y,coh[,1],coh[,2])
		ptk[ptk>1]=1
        	inZone[[i]]=inZone[[i]]-ptk
      		}
    	}

    nZone=nZone + inZone[[i]]
    mask=as.numeric(nZone<2)   #pt is in 1 zone or more
    #in that case do not reassign pt
    zone=zone + i*inZone[[i]]*mask

   }
  ind=1:nrow(tab)

  listZpt= vector("list", nbZ)
  for (k in 1:nbZ)
      {
	listZpt[[k]]=ind[zone==k]
      }

#check that all pts belong to a zone
v=1:nrow(pts)
us=unlist(listZpt)
noZ=v[!(v %in% us)]
# correct by using distances to zones
for (ind in noZ)
    {
    kk=1
    gdref=Inf
    iZ=0
    for(ii in 1:nbZ)
    	  {
	  gd=gDistance(tab[ind,],Z[[ii]])
	  if (gd<gdref)
	     {
	     gdref=gd
	     iZ=ii
	     }
	  }
     if(iZ!=0) listZpt[[iZ]]=append(ind,listZpt[[iZ]])
    }
  return(listZpt)
}



################################################################################
#' separationPoly
#'
#' @details separates holes and non holes
#' @param polyTot a SpatialPolygons object
#'
#' @return a SpatialPolygons with holes in separate polygons
#'
#' @keywords internal
#'
#' @examples
#' data(mapTest)
#' cL=list()
#' cL=contourAuto(cL,mapTest$step,mapTest$xsize,mapTest$ysize,mapTest$krigGrid,c(5,7),mapTest$boundary)
#' plot(mapTest$boundary,type="l",col="red")
#' graphics::lines(cL[[8]])
#' pG=polyToSp2(sp::Polygon(mapTest$boundary)) # transform boundary into SpatialPolygons objects
#' cLSp=maptools::ContourLines2SLDF(list(cL[[8]])) # transform contour line into SpatialLines objects
#' polyBuff=rgeos::gBuffer(cLSp,width=0.00001) # extend geometry
#' polyDiff=rgeos::gDifference(pG,polyBuff)
#' recupPoly=geozoning:::separationPoly(polyDiff)
#' Z1=list(recupPoly[[1]],recupPoly[[2]])
#' plotZ(Z1)
#' # not run
separationPoly=function(polyTot)
################################################################################

{
  # indFull indicates if Polygon has no hole (TRUE) or has holes (FALSE)
  indFull=numeric()
  listePoly=list()

  #reach Polygons structure
  polyList=polyTot@polygons[[1]]@Polygons

  # create comments for holes
  # ex : 0 1 0 1 3
  # ->indicates that polygones 1 and 3 have no holes, 2 and 4 are holes belognigng to 1 and 5 is a hole belonging to 3
  # 0 = no hole in polygon
  polyCom=createPolygonsComment(polyTot@polygons[[1]])
  listeCom<-as.numeric(unlist(strsplit(polyCom, " ")))
  indFull=grep(TRUE,listeCom==0)

  # For each polygon which has no hole, associate holes
  for(k in (1:length(indFull)))
  {
    i=indFull[k]
    indHole=grep(TRUE,listeCom==i)
    newP=c(polyList[i],polyList[indHole])
    listePoly[[k]]=SpatialPolygons(list(Polygons( newP , "1")))
  }

  return(listePoly)
}

################################################################################
#' extensionLine
#'
#' @details closes contour lines by extending them to their interesection with the map border
#' @param contourL contour line
#' @param step grid step as returned by calStep
#' @param bdSP list, contains x and y coordinates of map boundaries
#' @param superLines object returned by superLines(bdSP)
#'
#' @return a list
#' @importFrom sp SpatialPoints
#' @keywords internal
#'
#' @examples
#' data(mapTest)
#' step=mapTest$step
#' xsize=mapTest$xsize
#' ysize=mapTest$ysize
#' cL=contourLines(seq(step, xsize-step, by=step),seq(step, ysize-step, by=step),
#'                mapTest$krigGrid, levels = c(5,7))
#' plot(mapTest$boundary,type="l",col="red")
#' lines(cL[[1]])#contour line is not closed
#' lines(geozoning:::extensionLine(cL[[1]],step,sp::SpatialPoints(mapTest$boundary),
#'      geozoning:::superLines(mapTest$boundary)),col="red") #contour line is closed
#' # not run
extensionLine=function(contourL=NULL,step=NULL,bdSP,superLines)
################################################################################

{
  #gets end pts of contourL line
  contourL2=  cbind(contourL$x,contourL$y)
  colnames(contourL2)=c("x","y")
  contourL2 = SpatialPoints(contourL2)

  # extend contour line
  p3 = contourL2[1]
  p4=contourL2[length(contourL2)]
  gDist1= gDistance(p3,superLines,byid=TRUE)
  gDist2= gDistance(p4,superLines, byid=TRUE)
  indMin1 = which.min(gDist1)
  indMin2 = which.min(gDist2)

  right1 = bdSP[indMin1:(indMin1+1)]

  p1 = right1[1,]
  p2 = right1[2,]
  u = ((p3$x - p1$x)*(p2$x - p1$x)+ (p3$y - p1$y)*(p2$y - p1$y)) / ((p2$x-p1$x)^2+(p2$y-p1$y)^2)
  x3 = p1$x + u*(p2$x-p1$x)
  y3 = p1$y + u*(p2$y-p1$y)

  right2 = bdSP[indMin2:(indMin2+1)]
  p1 = right2[1,]
  p2 = right2[2,]
  u = ((p4$x - p1$x)*(p2$x - p1$x)+ (p4$y - p1$y)*(p2$y - p1$y)) / ((p2$x-p1$x)^2+(p2$y-p1$y)^2)
  x4 = p1$x + u*(p2$x-p1$x)
  y4 = p1$y + u*(p2$y-p1$y)

  contourL$x = c(x3,contourL$x,x4)
  contourL$y = c(y3,contourL$y,y4)
  #necessary to avoid duplicate rownames
  names(contourL$x)=NULL
  names(contourL$y)=NULL

  return(contourL)
}

