##################################################################
#' zoneQ
#'
#' @details description, a paragraph
#' @param contourSp xxxx
#' @param iC xxxx
#' @param iE xxxx
#' @param Z xxxx
#' @param K xxxx
#' @param map xxxx
#' @param simplitol xxxx
#'
#' @return a ?
#' @importFrom rgeos createSPComment
#'
#' @export
#'
#' @examples
#' # not run
zoneQ = function (contourSp,iC,iE,Z,K,map,simplitol)
##################################################################
{
	# add one contour to replace zone in existing zoning
   	if (is.null(contourSp)) return(NULL)
	le=length(Z)
	if (le == 0) return(NULL)
	# save englobing zone and current zone ids
	idE = getId(Z,iE)
	idC = getId(Z,iC) # keep current zone id from initial zoning
	# first merge current zone and englobing zone
	Zopti = zoneFusion4(Z,iC,iE,map,simplitol,disp=FALSE)

	if (is.null(Zopti)) return(NULL)
	# englobing zone number may have changed-find it using zone id
	iBig=findNumZ(Zopti,idE)
	# then intersect new contour with englobing zone
	polyDiff=gDifference(Zopti[[iBig]],contourSp)
	if(is.null(polyDiff)) return(NULL)

	recupPoly=separationPoly(polyDiff)
	# correct difference bug for small discrepancies
	# may be use cleanSp instead
	kp=1
	if (length(recupPoly) >1)
	{
	for (kk in 1:length(recupPoly))
	{
		ga=gArea(recupPoly[[kk]])
		if(ga > simplitol) kp=kk
	}
	}
	rpc1=createSPComment(recupPoly[[kp]])
	#
	Zopti[[iBig]]=gBuffer(rpc1, byid=TRUE, width=0)
	Zopti=setId(Zopti,iBig,idE)

	# find if new zone has zones within it (=holes)
	# if so remove them (create holes)
	ct=sapply(Zopti,gWithin,contourSp)
	ind=(1:(le-1))[ct]
	for (kk in ind)
	    {
	    contourSp=gDifference(contourSp,Zopti[[kk]])
	    }
	# set zones and ids
	# use gBuffer with width=0 to avoid self inter pbs
	Zopti[[le]]=cleanSp(gBuffer(contourSp, byid=TRUE, width=0))
	Zopti=setId(Zopti,le,idC)

	#print(gCentroid(Zopti[[1]])@coords[1,])

	return(Zopti)
}
