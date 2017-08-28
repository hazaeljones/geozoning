##############################################
#' detectSmallZones
#'
#' @details description, a paragraph
#' @param zonePolygone list of zones, each zone is a SpatialPolygons
#' @param minSize zone area threshold under which a zone is too small to be manageable
#'
#' @return a vector pf small zones indices
#' @importFrom rgeos gArea
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' ZK=initialZoning(qProb=c(0.4,0.7),mapTest)
#' Z=ZK$resZ$zonePolygone
#' iSmall=detectSmallZones(Z,minSize) # 2 small zones 
#' # not run
detectSmallZones=function(zonePolygone,minSize)
##############################################
{
  # On détecte la taille des zones (ici la taille est la plus grande distance entre deux points du polygone)
  #  et leur largeur(on tente une érosion,si elle échoue la zone etait trop étroite)
  #  et on renvoie leurs inds en vue d'une suppression
  # bch septembre 2015
  # gArea(zonePolygone[[i]] renvoie la surface
  vectSize=numeric()
  vectIndex=numeric()

  # Pour chaque zone
  for(i in (1:length(zonePolygone)))
  {
     vectSize[i] = gArea(zonePolygone[[i]])
    if (vectSize[i]<minSize)
    {
      vectIndex=append(vectIndex,i)
    }
  }

  # sort by increasing size

 surf = vectSize[vectIndex]
 mask=order(surf)

 return(list(vectIndex=vectIndex[mask]))
}


##################################################################
#' zoneFusion2 basic function for merging 2 zones
#'
#' @details merge 2 zones, called by zoneFusion3 and zoneFusion4
#' @param zoneMain zone to merge into
#' @param zoneSuppr zone to remove by merging it into main zone
#' @param simplitol tolerance for spatial polygons geometry simplification
#'
#' @return a zone
#' @importFrom rgeos createPolygonsComment gSimplify gUnion
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' plotZ(Z)
#' plot(zoneFusion2(Z[[6]],Z[[2]]),add=T,col="blue")
#' # not run
zoneFusion2 = function(zoneMain,zoneSuppr,simplitol=1e-3)
##################################################################
{
  comment(zoneMain@polygons[[1]])=createPolygonsComment(zoneMain@polygons[[1]])
  comment(zoneSuppr@polygons[[1]])=createPolygonsComment(zoneSuppr@polygons[[1]])
  zoneSupprBuff=gBuffer(zoneSuppr,width=simplitol)

  # on tente de regrouper les deux zones concernées
  zoneTot=gUnion(zoneMain,zoneSupprBuff)
  zoneTot=cleanSp(zoneTot) # remove artefacts
  comment(zoneTot@polygons[[1]])=createPolygonsComment(zoneTot@polygons[[1]])

  return(zoneTot)
}

######################################################################
#' zoneFusion3
#'
#' @details merge current zone #iC with neighbor zone in zoning. If there are several neighbor zones, the selected one is the zone whose area is greater than the admissible size threshold that has the closest average value to the current one. 
#' @param K zoning object, as returned by the calNei function
#' @param iC index of current zone in zoning
#' @param Ns zone neighborhood Boolean matrix  
#' @param map object returned by function genMap or genMapR
#' @param minSize  minimum admissible zone size
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param disp information level (0-no info, 1-print info, 2-plot)
#'
#' @return a zone obtained by merging current zone with neighbor zone
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Ns=getNs(K$zoneNModif,5) # find neighbors of zone 5
#' zoneFusion3(K,5,Ns,mapTest,disp=2) # merge and plot result of merging
zoneFusion3=function(K,iC,Ns,map,minSize=1e-2,simplitol=1e-3,disp=0)
######################################################################
{
  #########################################
  #merge zone iC with neighbor in Ns
  #########################################
	#
	Z=K$zonePolygone
  	listN =  grep( TRUE , Ns)
	indZV = findN(K,listN,iC,minSize) # find the neighbor zone with which to merge the current zone

	if (indZV == 0) return(NULL) # no neighbour found for merging

	  # could also be written as
    # slot(slot(Z[[indZV]],"polygons")[[1]],"ID")
    #########################################
	  if(disp>0) print(paste("merging zone",iC," with main zone",indZV))
	  # case when merging  ZA with ZB and ZB is within ZA
	  # ids must be swapped 
	  # keep outer polygon ID
	  IN=gContains(gEnvelope(Z[[iC]]),Z[[indZV]])
	  if (IN)
	     newId= Z[[iC]]@polygons[[1]]@ID
	  else
	     newId= Z[[indZV]]@polygons[[1]]@ID #
	  #
	  Z[[indZV]] = zoneFusion2(  Z[[indZV]], Z[[ iC ]],simplitol)
 	  # reset polygon ID
    	  Z[[indZV]]@polygons[[1]]@ID = newId
    	  #remove zone that was merged
    	  Z[[iC]]=NULL
	if (disp==2)
	{
	   # plot (bch)
	   # IS 19/05/2017: add comment for x11
	   x11()
	   #IS 19/05/2017: modify this call...
     dispZ(map$step,map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
	}
  return(Z)
}

######################################################################
#' zoneFusion4
#'
#' @details merge 2 zones from given zoning
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iSmall index of zone to remove by merging it into other zone
#' @param iBig index of zone to merge into
#' @param map map object returned by function genMap
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param disp 0: no info, 1: some info
#'
#' @return a new zoning geometry
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' zoneFusion4(Z,5,4,mapTest,disp=2)
#' # not run
zoneFusion4=function(Z,iSmall,iBig,map,simplitol=1e-3,disp=0)
######################################################################
{
#########################################
#merge zone iSmall with zone iBig
#########################################

	 if(disp>0) print(paste("merging zone",iSmall," into main zone",Z[[iBig]]@polygons[[1]]@ID))
	 newId= Z[[iBig]]@polygons[[1]]@ID # keep polygon ID
	 Z[[iBig]] = zoneFusion2(  Z[[iBig]], Z[[ iSmall ]],simplitol)
   # reset polygon ID
   Z[[iBig]]@polygons[[1]]@ID = newId
	 # hack to avoid self intersection pbs
   Z[[iBig]] = gSimplify(Z[[iBig]],tol=simplitol)

   #remove zone that was merged
   Z[[iSmall]]=NULL
	 if (disp==1)
	 {
	    # plot resulting zoning
	    # IS 19/05/2017: add comment for x11
      #x11()
      dispZ(map$step,map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
	 }

  return(Z)
}

############################################################################
#' zoneGrow
#'
#' @details either grow isolated zone or group 2 zones together
#' if isolated zone, run optimization procedure to find the new quantile
#' if zone very small (area < minSizeNG) do not grow it
#' @param K zoning object, such as returned by the calNei function
#' @param map object returned by function genMap
#' @param iC index of current zone
#' @param optiCrit criterion choice
#' @param minSizeNG zone area threshold under which a zone will be removed
#' @param distIsoZ threshold distance to next zone, above which a zone is considered to be isolated
#' @param LEQ length of quantile sequence used to grow isolated zone
#' @param MAXP quantile sequence maximum shift from center
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param disp information level (0-no info, 1-print info)
#'
#' @return a zone obtained by growing current zone
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.2,0.5)
#' ZK = initialZoning(qProb, mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' plotZ(K$zonePolygone) # plot zoning
#' kmi=zoneGrow(K,mapTest,6) # grow zone 6 by grouping it with its closest neighbor with same label
#' linesSp(kmi[[7]])
#' qProb=c(0.3,0.5)
#' criti = correctionTree(qProb,mapTest)
#' best = criti$zk[[2]][[8]]
#' Z=best$zonePolygone
#' plotZ(Z)
#' refPoint = gCentroid(Z[[4]])
#' plot(refPoint,add=T,col="blue",pch=21)
#' zg=zoneGrow(best,mapTest,4) #grow isolated zone 4 by searching for other quantile
#' plotZ(zg)
#' # not run
zoneGrow=function(K,map,iC,optiCrit=2,minSizeNG=1e-3,distIsoZ=0.075,LEQ=5,MAXP=0.1,simplitol=1e-3,disp=0)
############################################################################
{
	# either grow isolated zone or group 2 zones together
	# if isolated zone, optim procedure to find the new quantile
	# if zone very small, skip (useless) growing step
#
  Z=K$zonePolygone
  if(disp>0) print(paste("trying to grow zone id",getZoneId(Z[[iC]]), "- new number", iC))
  Ns = getNs(K$zoneNModif,iC)
  qProb=K$qProb
  if(is.null(qProb)) return(NULL)
  refSurf = gArea(Z[[iC]])
  if (refSurf < minSizeNG) return(NULL)

  resC = detZoneClose(iC,Z,K$zoneNModif,distIsoZ) # renvoie FALSE si zone trop proche dune autre, TRUE sinon
  ##############################################################
  InterZoneSpace = resC$InterZoneSpace
  zoneClose = resC$zoneClose
  step=map$step

  # keep centroid of small zone
  # used to check that it is the same zone that has grown
  # contours can be anywhere on the plot
  ##########################################################
  refPoint = gCentroid(Z[[iC]])
  ##########################################################
	Zopti=NULL
  if (InterZoneSpace)#isolated zone
  {
	  if (disp>0) print(paste("growing isolated zone: ",getZoneId(Z[[iC]])))
    ##############################################################
	  #searches for the best quantile to grow zone
    ##############################################################
    resZ = optiGrow(K,iC,qProb,refPoint, map,optiCrit,minSize,minSizeNG,distIsoZ,LEQ,MAXP,simplitol,disp)
		# returns NULL if dead end
		if (!is.null(resZ))
		{
			if(disp) print("growing successful")
			Zopti = resZ$Zopti
			# create comments for holes
			Zopti = crComment(Zopti)
		}
	 ##############################################################
  }
  else #non isolated zone
  {
	   if (length(zoneClose)==0) return(NULL) # no close zone
	   if (disp>0) print(paste("growing non isolated zone ",getZoneId(Z[[iC]]), "(close to zone",getZoneId(Z[[zoneClose[[1]]]]),")"))
     # reuse zoneClose
     Kopti = zoneModifnonIso(K,qProb,map,zoneClose,iC,simplitol,disp)
     # create comments for holes
     if (!is.null(Kopti$zonePolygone)) Zopti = crComment(Kopti$zonePolygone)
     Kopti$zonePolygone=Zopti
	}
  ##########################################################

  if (disp==2 && !is.null(Zopti))
	{
	   dispZ(map$step,map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
	 }
	return(Zopti)
}


###################################################
#' remove1FromZ
#'
#' @details description, a paragraph
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param iC current zone index
#' @param zoneN zone neighborhood Logical matrix
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param disp 0: no info, 1: some info
#'
#' @return a new zoning where current zone has been removed
#'
#' @export
#' @examples
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' plotZ(remove1FromZ(Z,2,K$zoneN))
#' # not run
remove1FromZ = function(Z,iC,zoneN,simplitol=1e-3,disp=0)
########################################################
{
  # remove zone iC from zoning
  # first find neighbor zone for merging zone iC
  # then delete zone iC
  diag(zoneN)=FALSE #exclude zone is its own neighbor

  iNP=grep(TRUE,zoneN[iC,]) # may be empty if no pt
	if (length(iNP)==0)
	{
	  indN = (1:length(Z))[-iC] #exclude current zone
	  dN=rep(1,length(Z))
	  for (k in indN)
    {
      iN = k
	    gd=gDifference(Z[[iN]],Z[[iC]])
	    dN[k] = gDistance(gd,Z[[iC]])
	  }
	  iN=which(dN==min(dN))
	  iN=iN[1]
	} # end no pt
	else
	{
    # if more than 1 neighbor, take the closest zone
    d=sapply(Z,gDistance,Z[[iC]])
    kk = which(d[iNP]==min(d[iNP]))
	  iN = iNP[kk[1]]
	}
  newId= Z[[iN]]@polygons[[1]]@ID
  Z=zoneFusion4(Z,iC,iN,map,simplitol,disp)

	return(Z)
}


##########################################################################
#' removeFromZ
#'
#' @details description, a paragraph
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param zoneN zone neighborhood Logical matrix
#' @param ptN indices of data pts neighbours
#' @param listZonePoint list of indices of data points within zones, result of call to \code{\link{calNei}}
#' @param data SpatialPointsDataFrame with data values
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param n minimal number of points below which a zone is removed from zoning
#'
#' @return a list with components
#'\describe{
#' \item{Z}{new zoning geometry (list of SpatialPolygons)} where zones with less than n points were removed
#' \item{zoneN}{new zone neighborhood Logical matrix}
#' \item{listZonePoint}{new list of indices of data points within zones}
#' }
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Z=K$zonePolygone
#' plotZ(Z)
#' remove from Z all zones with less than 10 data points
#' Z2=removeFromZ(Z,K$zoneN,K$krigN,K$listZonePoint,mapTest$krigData,n=10)
#' printZid(Z2$Z)
#' # not run
removeFromZ = function(Z,zoneN,ptN,listZonePoint,data,simplitol=1e-3,n=1)
##########################################################################
{
# remove from Z all zones with npts<=n or area<minSizeNG

  mask1 = sapply(listZonePoint,function(x){return(length(x)<=n)})
  mask2 = sapply(Z,function(x){return(gArea(x)<minSizeNG)})
  nbZ=length(Z)
  ind = 1:nbZ
  ind = ind[mask1 | mask2]
  ids=c()
  if (!is.null(ind))
     ids = getIds(Z,ind)
  for (zid in ids)
      {
      iC = Identify(zid,Z)
      Z = remove1FromZ(Z,iC,zoneN,simplitol)
      nbZ=length(Z)
      zoneN=matrix(logical(nbZ^2),nbZ,nbZ)
      #update zone assignment
      listZonePoint=zoneAssign(data,Z)
      #recreate zone neighbors
      vZ=calZoneN(ptN,zoneN,listZonePoint)
      zoneN = vZ$zoneN
      }

  return(list(Z=Z,zoneN=zoneN,listZonePoint=listZonePoint))
}
