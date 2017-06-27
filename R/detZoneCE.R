#################################################################
#' detZoneClose
#'
#' @details description, a paragraph
#' @param iZ xxxx
#' @param Z xxxx
#' @param K xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
detZoneClose=function(iZ,Z,zoneNModif)
##################################################################
  {
    # iZ=current zone
    # Z=list of zones
    # neighbor zone = shares some points
    # close zone = not neighbor and not too far
    # distIsoZ = maximum distance for isolated zone (has no close zones)
    # december 2016

    Ns=getNs(zoneNModif,iZ)
    notN = grep( FALSE, Ns) # list of non neighbor zones
    notN = notN[notN != iZ] # remove current zone from that list

    InterZoneSpace=TRUE # no close zone (except neighbor zones)

    D=numeric(0)
    zoneClose = list()

          for (zC in notN)
          {
	    d=gDistance(Z[[iZ]],Z[[zC]])
            if (d<distIsoZ)
            {
              InterZoneSpace=FALSE
              zoneClose =append(zoneClose,zC)
	      D=append(D,d)
            }
          }
# sort by increasing distance
  if (! InterZoneSpace) zoneClose=zoneClose[order(D)]
# check  englobing zone is the same (no crossing)
  le = length(zoneClose)
  indE=detZoneEng(iZ,Z,zoneNModif)

 	if (le>0)
	{
		for (ip in le:1)
 		{
		zip=zoneClose[[ip]]
		indEP=detZoneEng(zip,Z,zoneNModif)
		mask = (indE == 0 ) | (indEP != indE) | (indEP == 0)
		if ( mask) zoneClose[[ip]]=NULL

 		}
	}


 return(list(InterZoneSpace=InterZoneSpace,zoneClose=zoneClose))
  }

##################################################################
#' detZoneP0
#'
#' @details description, a paragraph
#' @param iZ xxxx
#' @param Z xxxx
#' @param K xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
detZoneP0=function(iZ,Z,zoneNModif)
##################################################################
  {
    # iZ=indice dans la liste de la zone à agrandir
    # Z=liste of zones
    # distIsoZ = dist max for isolated zone
    # modif bch :

    Ns=getNs(zoneNModif,iZ)

    InterZoneSpace=TRUE

     listeV=grep(TRUE, Ns)
   # zone englobing current zone
   indEG=detZoneEng(iZ,Z,zoneNModif)
   zoneClose = list()
   if (indEG ==0)
      {
	InterZoneSpace = FALSE
	return(list(InterZoneSpace=InterZoneSpace,zoneClose=zoneClose))
	}

    #zone list except englobing zone
    notN = grep( FALSE, Ns)
    notN = notN[notN!=iZ]
      #

	  D=numeric(0)

          for (zC in notN)
          {
	    d=gDistance(Z[[iZ]],Z[[zC]])
            if (d<distIsoZ)
            {
              InterZoneSpace=FALSE
              zoneClose =append(zoneClose,zC)
	      D=append(D,d)
            }
          }
# trier par distance croissante
if (! InterZoneSpace)
 {
	zoneClose=zoneClose[order(D)]

 }


 if (indEG == 0)
    zoneClose=list()
 else
 {
	#
 	# check  englobing zone is the same (no crossing)
 	if (length(zoneClose)>0)
	{
		for (ip in length(zoneClose):1)
 		{
		zip=zoneClose[[ip]]
		indEGP=detZoneEng(zip,Z,K$zoneNModif)
		if (indEGP != indEG) zoneClose[[ip]]=NULL

 		}
	}
 #eliminate close zones included in current zone
 le=length(zoneClose)
 if (le > 0)
    {
    for (ip in le:1)
     	{
	zip=zoneClose[[ip]]
	gc=gContains(gBuffer(gConvexHull(Z[[iZ]]),byid=TRUE,width=1e-3),Z[[zip]])
	if (gc) zoneClose[[zip]]=NULL
 	}
     }
}
 return(list(InterZoneSpace=InterZoneSpace,zoneClose=zoneClose))
}

##################################################################
#' detZoneEng0
#'
#' @details description, a paragraph
#' @param iZ xxxx
#' @param Z xxxx
#' @param K xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
detZoneEng0=function(iZ,Z,K)
##################################################################
{
    # iZ=indice dans la liste de la zone à agrandir
    # Z=liste des zones
    #On recupere la zone la plus proche qui englobe notre petite zone
    #

  listeVois = grep( TRUE ,K$zoneNModif[iZ,] )
  indiceZe=0
  for (iv in listeVois)
  {
	# enveloppe convexe pose pb parfois
   	gc=gContains(gBuffer(gConvexHull(Z[[iv]]),byid=TRUE,width=1e-3),Z[[iZ]])
	if (gc) indiceZe=iv
  }

  return(indiceZe)
}

##################################################################
#' detZoneEng
#'
#' @details description, a paragraph
#' @param iZ index of zone for which englobing zone is searched
#' @param Z zoning
#' @param zoneNModif modified zone neighborhood matrix (FALSE values on diagonal)
#'
#' @return an integer value (0 if no englobing zone was found, englobing zone index otherwise)
#'
#' @export
#'
#' @examples
#' # load zoning results from test file
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' zoneNModif=resZTest$zoneNModif
#' detZoneEng(3,Z,K$zoneNModif) # zone 2 englobes zone 3
#' detZoneEng(2,Z,K$zoneNModif) # no englobing zone for zone 2 
detZoneEng=function(iZ,Z,zoneNModif)
##################################################################
{
    # iZ=index of zone for which englobing zone is searched
    # Z=zoning 
    # result = closest englobing zone

  listN = grep( TRUE ,zoneNModif[iZ,] )
  indZe=NULL
  ar=NULL
  for (iN in listN)
  {
	# convex envelope (problematic sometimes)
	ge=gEnvelope(Z[[iN]])
	arv=gArea(ge)
   	gc=gContains(ge,Z[[iZ]])
	if (gc)
	   {
		ar=c(ar,arv)
		indZe=c(indZe,iN)
		}

  }

  if(is.null(indZe)) return(0)
  mask=(ar==min(ar)) # if several englobing zones
  iZe=indZe[mask]
  iZe=iZe[1]
  return(iZe)
}

##################################################################
#' detZoneEng2
#'
#' @details description, a paragraph
#' @param iZ xxxx
#' @param Z xxxx
#' @param K xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
detZoneEng2=function(iZ,Z,zoneNModif)
##################################################################
{
    # iZ=indice dans la liste de la zone à agrandir
    # Z=liste des zones
    #On recupere la zone la plus petite qui englobe notre petite zone
    # version sans enveloppe convexe
    #ne marche pas a cause des pbs de bordure...
  listeN = grep( TRUE ,zoneNModif[iZ,] ) # only consider neighbour zones
  indiceZe=NULL

  for (iv in listeN)
  {
	ge=polyToSp2(getPolySp(Z[[iv]],1)) #external polygon


   	gc=gContains(ge,Z[[iZ]])
	if (gc)
	   {

		indiceZe = c(indiceZe,iv)
		}

  }

  if(is.null(indiceZe)) return(0)

  iZe=indiceZe[1]
  return(iZe)
}

##################################################################
#' detZonePIso
#'
#' @details description, a paragraph
#' @param iZ xxxx
#' @param Z xxxx
#' @param K xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
detZonePIso = function(iZ,Z,zoneNModif)
##################################################################
{
# closest zone (excluding neighboring zones) to a given zone
# used to generate quantile sequence to grow isolated zone
  Ns=getNs(zoneNModif,iZ)
  notN = grep( FALSE, Ns)
  notN = notN[notN!=iZ]
  zP=c()
  D=numeric(0)

          for (zC in notN)
          {
	    d=gDistance(Z[[iZ]],Z[[zC]])

            zP=append(zP,zC)
	    D=append(D,d)

          }
# sort
	zP=zP[order(D)]
	return(zP[1])
 }
