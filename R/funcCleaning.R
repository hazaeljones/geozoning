##############################################
#' detectSmallZones
#'
#' @details description, a paragraph
#' @param zonePolygone xxxx
#' @param minSize xxxx
#'
#' @return a ?
#' @importFrom rgeos gArea
#'
#' @export
#'
#' @examples
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
#' zoneFusion2
#'
#' @details description, a paragraph
#' @param zoneMain xxxx
#' @param zoneSuppr xxxx
#' @param simplitol xxxx
#'
#' @return a ?
#' @importFrom rgeos createPolygonsComment gSimplify gUnion
#'
#' @export
#'
#' @examples
#' # not run
zoneFusion2 = function( zoneMain,zoneSuppr,simplitol)
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param K xxxx
#' @param iC xxxx
#' @param Ns xxxx
#' @param map xxxx
#' @param minSize xxxx
#' @param simplitol xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
zoneFusion3=function(Z,K,iC,Ns,map,minSize,simplitol,disp=0)
######################################################################
{
  #########################################
  #merge zone iC with neighbor in Ns
  #########################################
	#
  listN =  grep( TRUE , Ns)
	indZV = findN(Z,K,listN,iC,minSize) # chercher le voisin (parmi tous) avec lequel fusionner

	if (indZV == 0) return(NULL) # step de voisin avec lequel fusionner

	  # on pourrait ecrire aussi :
    # slot(slot(Z[[indZV]],"polygons")[[1]],"ID")
    #########################################
	  if(disp>0) print(paste("merging zone",iC," with main zone",indZV))
	  # attention, modif bch octobre 2016
	  # cas ou on fusionne ZA avec ZB et ZB est inclus dans ZA
	  # il faut alors echanger les ids, pour garder le label de ZB
	  # (contour le plus exterieur)
	  # keep outer polygon ID
	  IN=gContains(gEnvelope(Z[[iC]]),Z[[indZV]])
	  if (IN)
	     newId= Z[[iC]]@polygons[[1]]@ID
	  else
	     newId= Z[[indZV]]@polygons[[1]]@ID #
	  #
	  Z[[indZV]] = zoneFusion2(  Z[[indZV]], Z[[ iC ]],simplitol)
    #Z[[indZV]]=gUnion(  Z[[indZV]], Z[[ iC ]] )
    # reset polygon ID
    Z[[indZV]]@polygons[[1]]@ID = newId
    #remove zone that was merged
    Z[[iC]]=NULL
	if (disp==2)
	{
	   # plot (bch)
	   # IS 19/05/2017: add comment for x11
	   #x11()
	   #IS 19/05/2017: modify this call...
     dispZ(map$step,map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
	}
  return(Z)
}

######################################################################
#' zoneFusion4
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param iSmall xxxx
#' @param iBig xxxx
#' @param map xxxx
#' @param simplitol xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
zoneFusion4=function(Z,iSmall,iBig,map,simplitol,disp=0)
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
#' @details description, a paragraph
#' @param Z xxxx
#' @param K xxxx
#' @param iC xxxx
#' @param Ns xxxx
#' @param map xxxx
#' @param optiCrit xxxx
#' @param valRef xxxx
#' @param qProb xxxx
#' @param minSizeNG xxxx
#' @param distIsoZ xxxx
#' @param LEQ xxxx
#' @param MAXP xxxx
#' @param simplitol xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
zoneGrow=function(Z,K,iC,Ns,map,optiCrit,valRef,qProb,minSizeNG,distIsoZ,LEQ,MAXP,simplitol,disp=0)
############################################################################
{
	## On va lancer une procedure optimisation pour determiner la taille de lagrandissement
  ## Si espace est suffisamment grand, alors on lagrandit, en choisisst la meilleure taille possible
  ## parmi une courte liste (souci de vitesse execution).
	# modif bch juin 2016
	# if zone very small, skip (useless) growing step
	# param minSizeNG in initParam.R
	if(disp>0) print(paste("trying to grow zone",getZoneId(Z[[iC]])))

  refSurf = gArea(Z[[iC]])
	if (refSurf < minSizeNG) return(NULL)

  resC = detZoneClose(iC,Z,K) # renvoie FALSE si zone trop proche dune autre, TRUE sinon
  ##############################################################
  InterZoneSpace = resC$InterZoneSpace
  zoneClose = resC$zoneClose
  step=map$step

  # on conserve le centroide de la petite zone
  # utilise apres agrandissement pour verifier que c'est la meme zone
	##########################################################
  refPoint = gCentroid(Z[[iC]])
  ##########################################################
	Zopti=NULL
  if (InterZoneSpace)#isolated zone
  {
	  if (disp>0) print(paste("growing isolated zone: ",getZoneId(Z[[iC]])))
    ##############################################################
	  #fonction qui trouve le meilleur quantile pour agrandir la zone
    ##############################################################
    resZ = optiGrow(Z,K,iC,qProb,refPoint, map,optiCrit,minSize,minSizeNG,distIsoZ,LEQ,MAXP,simplitol,disp)
		# renvoie NULL si voie sans issue
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
     Zopti = zoneModifnonIso(Z,K,qProb,map,zoneClose,iC,simplitol,disp)
	 	 # create comments for holes
		 Zopti = crComment(Zopti)
	}
  ##########################################################

  if (disp==2 && !is.null(Zopti))
	{
	   #x11()
	   dispZ(map$step,map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
	 }
	return(Zopti)
}


###################################################
#' remove1FromZ
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param iC xxxx
#' @param zoneN xxxx
#' @param simplitol xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
remove1FromZ = function(Z,iC,zoneN,simplitol,disp=0)
###################################################
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
#' @param Z xxxx
#' @param zoneN xxxx
#' @param ptN xxxx
#' @param listZonePoint xxxx
#' @param data xxxx
#' @param simplitol xxxx
#' @param n xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
removeFromZ = function(Z,zoneN,ptN,listZonePoint,data,simplitol,n=1)
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
