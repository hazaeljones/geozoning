##################################################################
#' calFrame
#'
#' @details description, a paragraph
#' @param iC xxxx
#' @param Z xxxx
#' @param K xxxx
#' @param distIsoZ xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calFrame = function(iC,Z,K,distIsoZ=0.075)
# returns spatial polygon = envelope within which grown zone must be contained
##################################################################
{
	iE = detZoneEng(iC,Z,K)
	if (iE == 0 ) return(NULL)
	# closest zone (excluding hole)
	# eliminate from search zones included in current zone
	iP = getClosestZone(iC,Z,K)
	if(iP == 0) return(NULL)
	# dilate zone Z[[iC]] so that closest zone distance = distIsoZ
	dP = gDistance(Z[[iC]],Z[[iP]])
	pc = getPolySp(Z[[iC]],1)
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



