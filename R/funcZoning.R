###############################################################################
#' zoneGeneration
#'
#' @details description, a paragraph
#' @param map xxxx
#' @param qProb xxxx
#' @param GridData xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
zoneGeneration=function(map,qProb=c(0.25,0.75),GridData=FALSE)
#################################################################################renvoie une liste de polygones,détection grace aux isocontours et quantiles

#entrée:1/step = intervalle de repartition des points(numeric),
#       taille en x et y de la map(numeric),vecteur des valeurs des points(numeric),
#       dataframe des coord.et valeurs des points initiaux(data.frame),matrice des valeurs des pts (krigés par exemple)(matrix),
#       boundary de la parcelle, quantiles
#sortie:liste des coordonnées des polygones(list(list(numeric)))
# data in matVal argument

{
  matVal=map$krigGrid
  boundary=map$boundary
  step=map$step
  # On définit le contour de la parcelle
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
     cL=contourAuto(cL,step,matVal,vRef=valQuant[i],frame,GridData)
  }

  # For each isocontour
  for (jContour in (1:length(cL)))
    {
      # On le passe en structure spatiale (2eme structure de données, Objet Spatial)
      cLSp[[jContour]] <- ContourLines2SLDF(list(cL[[jContour]]))
      # attribut level conserve
      # tableau temporaire qui remplacera le tableau original
      # une fois toutes les operations effectuees pour cet isocontour = intersection des polygones existants avec isocontour

      listPolyTmp=list()

      # define buffer: contour and another line around contour
      polyBuff=gBuffer(cLSp[[jContour]],width=0.0001*(1/step))
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
#' @details description, a paragraph
#' @param cL xxxx
#' @param step xxxx
#' @param matVal xxxx
#' @param vRef xxxx
#' @param boundary xxxx
#' @param GridData xxxx
#'
#' @return a ?
#' @importFrom grDevices contourLines
#' @importFrom sp Line coordinates
#'
#' @export
#'
#' @examples
#' # not run
contourAuto=function(cL,step,matVal,vRef,boundary,GridData=FALSE)
##############################################################################
{
 #---------------------------------------------------------------------------------------------------------------------------------#
#fonction qui construit un contour de zones à partir de la liste des isocontours
# donnees dans matVal
#entrée:numéro du contour pour ce polygone(numeric),dataframe des positions des points et de leurs valeurs(spDataFrame)
#   1/step=ecart entre deux points sur la grille (numeric),taille du frame en x et y(numeric),matrice des valeurs des points krigés(matrix)
#   liste des contours deja definis pour ce polygone(list(list(numeric))),liste des polygones deja définis(list(list(numeric)))
#sortie:contour (list(numeric))/à vérifier

  # find isocontours (level attribute)
  if(!GridData)
	cLplus=contourLines(seq(1/step, 1 - 1/step, by=1/step),seq(1/step, 1 - 1/step, by=1/step),matVal, levels = vRef)
 else

	cLplus=contourLines(z=matVal, levels = vRef)
  # merge with previous isocontours (other level)
  cL=c(cL,cLplus)

  boundary = data.frame(boundary)
  sp::coordinates(boundary)=~x+y
  bl=Line(coordinates(boundary))
  bSPL1=SpatialLines(list(Lines(list(bl),'1')))

  #for each isocontour, extend it to frame
  for(jContour in (1:length(cL)))
  	{
  	#test if contour is closed or not - if not, close it
	iso=cL[[jContour]]
    	np=length(iso$x) #number of pts in contour line
     	if(iso$x[1]!= iso$x[np] || iso$y[1] != iso$y[np])
    	  {
	  # if contour is not closed, add projection on frame
     	  ######################################################
	  cL[[jContour]]=extensionLine(iso,step,boundary,bSPL1)
	  ######################################################
    	  }
 	 }

  return(cL)
}

#########################################################################
#' zoneAssign
#'
#' @details description, a paragraph
#' @param tab xxxx
#' @param Z xxxx
#'
#' @importFrom sp coordinates
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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


  return(listZpt)
}



################################################################################
#' separationPoly
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
separationPoly=function(polyTot)
################################################################################
#fonction qui retourne les différents sous-éléments d'un Spatial Polygons( qui  ne sont pas des trous) normalement 2 au plus,
#en associant à chacun les bon trous
#entrée:spatialPolygon
#sortie: liste de spatialPolygons

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

  # For each polygon wich has no hole, associate holes
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
#' @details description, a paragraph
#' @param contour xxxx
#' @param step xxxx
#' @param boundary xxxx
#' @param bspl xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
extensionLine=function(contour=NULL,step=NULL,boundary,bspl)
################################################################################
################################################################################
#fonction qui complète les lignes définies sur un frame ou x et y appartiennent à l'intervalle [1/step, 1-1/step]:rajoute un point correspondant à la projection
#de leur extrémités sur le frame (0,0) (0,1) (1,1) (0,1)

#entrée:liste contenant x et y (list(numeric)),step=ecart entre les points sur la  grille(numeric),taille du frame en x et y(numeric)
#sortie:liste de coordonnées de points représentant un contour (list(numeric))

{
  #gets end pts of contour line

  boundary =SpatialPoints(boundary)
  contour2=  cbind(contour$x,contour$y)
  colnames(contour2)=c("x","y")
  contour2 = SpatialPoints(contour2)

  lignes = bspl@lines[[1]]@Lines[[1]]

  listLignes=list(Lines(list(Line(lignes@coords[1:2,])),'1'))

  for (i in 2:(length(lignes@coords)/2-1))
  {
    listLignes[[i]] = Lines(list(Line(lignes@coords[i:(i+1),])),paste(i))
  }

  p3 = contour2[1]
 # p4 = tail(contour2,1)#correction bch septembre 2015
  p4=contour2[length(contour2)]
  SuperLines = SpatialLines(listLignes)
  gDist1= gDistance(p3,SuperLines,byid=TRUE)
  gDist2= gDistance(p4,SuperLines, byid=TRUE)
  indMin1 = which.min(gDist1)
  indMin2 = which.min(gDist2)

  right1 = boundary[indMin1:(indMin1+1)]

  p1 = right1[1,]
  p2 = right1[2,]
  u = ((p3$x - p1$x)*(p2$x - p1$x)+ (p3$y - p1$y)*(p2$y - p1$y)) / ((p2$x-p1$x)^2+(p2$y-p1$y)^2)
  x3 = p1$x + u*(p2$x-p1$x)
  y3 = p1$y + u*(p2$y-p1$y)

  right2 = boundary[indMin2:(indMin2+1)]
  p1 = right2[1,]
  p2 = right2[2,]
  u = ((p4$x - p1$x)*(p2$x - p1$x)+ (p4$y - p1$y)*(p2$y - p1$y)) / ((p2$x-p1$x)^2+(p2$y-p1$y)^2)
  x4 = p1$x + u*(p2$x-p1$x)
  y4 = p1$y + u*(p2$y-p1$y)

  contour$x = c(x3,contour$x,x4)
  contour$y = c(y3,contour$y,y4)
  #necessary to avoid duplicate rownames
  names(contour$x)=NULL
  names(contour$y)=NULL

  return(contour)
}
