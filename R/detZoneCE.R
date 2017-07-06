#################################################################
#' detZoneClose
#'
#' @details description, a paragraph
#' @param iZ zone number
#' @param Z zoning
#'
#' @return a list with components
#' \describe{
#' \item{InterZoneSpace}{TRUE if zone is isolated, FALSE otherwise}
#' \item{zoneClose}{numbers of zones close to zone iZ, empty if zone is isolated}
#' }
#'
#' @export
#'
#' @examples
#'
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
