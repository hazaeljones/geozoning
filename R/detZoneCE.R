#################################################################
#' detZoneClose
#'
#' @details determines zones that are close to current zone, but not neighbors (common border). Therefore embedded or englobing zones are excluded. 
#' @param iZ zone number
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param zoneN modified zone neighborhood Logical matrix (FALSE values on diagonal)
#'
#' @return a list with components
#' \describe{
#' \item{InterZoneSpace}{TRUE if zone is isolated, FALSE otherwise}
#' \item{zoneClose}indices of zones close to zone iZ, empty if zone is isolated}
#' }
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' zoneN=resZTest$zoneNModif
#' plotZ(Z)
#' detZoneClose(4,Z,zoneN) # zone 4 is close to zone 3
#' detZoneClose(6,Z,zoneN) # zone 6 is isolated (no zone at a distance smaller than 0.075).
#' # not run
detZoneClose=function(iZ,Z,zoneN,distIsoZ=0.075)
##################################################################
  {
    # iZ=current zone
    # Z=list of zones
    # neighbor zone = shares some points
    # close zone = not neighbor and not too far
    # distIsoZ = maximum distance for isolated zone (has no close zones)
    # december 2016

    Ns=getNs(zoneN,iZ)
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
  indE=detZoneEng(iZ,Z,zoneN)

 	if (le>0)
	{
		for (ip in le:1)
 		{
		zip=zoneClose[[ip]]
		indEP=detZoneEng(zip,Z,zoneN)
		mask = (indE == 0 ) | (indEP != indE) | (indEP == 0)
		if ( mask) zoneClose[[ip]]=NULL

 		}
	}


 return(list(InterZoneSpace=InterZoneSpace,zoneClose=zoneClose))
  }


##################################################################
#' detZoneEng
#'
#' @details description, a paragraph
#' @param iZ index of zone for which englobing zone is searched
#' @param Z zoning
#' @param zoneN modified zone neighborhood matrix (FALSE values on diagonal)
#'
#' @return an integer value (0 if no englobing zone was found, englobing zone index otherwise)
#'
#' @export
#'
#' @examples
#' # load zoning results from test file
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' zoneN=resZTest$zoneNModif
#' detZoneEng(3,Z,zoneN) # zone 2 englobes zone 3
#' detZoneEng(2,Z,zoneN) # no englobing zone for zone 2 
detZoneEng=function(iZ,Z,zoneN)
##################################################################
{
    # iZ=index of zone for which englobing zone is searched
    # Z=zoning 
    # result = closest englobing zone

  listN = grep( TRUE ,zoneN[iZ,] )
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
