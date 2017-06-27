##################################################################
#' calFrame
#'
#' @details description, a paragraph
#' @param iZ index of zone for which the envelope is searched
#' @param Z zoning
#' @param  zoneNModif modified zone neighborhood matrix (FALSE values on diagonal
#' @param distIsoZ threshold distance above which a zone is considered as isolated
#'
#' @return a apatial polygon corresponding to the frame within which grown zone must be contained 
#'
#' @export
#'
#' @examples
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' zoneNModif=resZTest$zoneNModif
#' f=calFrame(6,Z,zoneNModif)
#' plotZ(Z)
#' plot(f,add=TRUE)
calFrame = function(iZ,Z,zoneNModif,distIsoZ=0.075)
# returns spatial polygon = frame within which grown zone must be contained
##################################################################
{
	iE = detZoneEng(iZ,Z,zoneNModif)
	if (iE == 0 ) return(NULL)
	# closest zone (excluding hole)
	# eliminate from search zones included in current zone
	iP = getClosestZone(iZ,Z,zoneNModif)
	if(iP == 0) return(NULL)
	# dilate zone Z[[iZ]] so that closest zone distance = distIsoZ
	dP = gDistance(Z[[iZ]],Z[[iP]])
	pc = getPolySp(Z[[iZ]],1)
	spc = polyToSp2(pc)

	if ((dP-distIsoZ) > 1e-3)
	   Ze = gBuffer(spc,width=(dP-distIsoZ),byid=TRUE)
	else
	   return(NULL)

	# get the  polygon (in case of hole, may be several polys)
	# the first one is the non hole one
	p = getPolySp(Ze,1)
	#transform p into sp
	spN = polyToSp2(p)

	return(spN)
}



