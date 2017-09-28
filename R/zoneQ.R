#############################################################################
#' zoneQ
#'
#' @details called by optiGrow,replaces the current zone by a bigger one
#' @param contourSp contour line transformed into SpatialPolygons 
#' @param iC zone to grow
#' @param iE englobing zone
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param K zoning object (such as returned by calNei function)
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @return a zoning geometry updated with the grown zone(list of SpatialPolygons)
#' @importFrom rgeos createSPComment
#' @keywords internal
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.3,0.5)
#' criti = correctionTree(qProb,mapTest)
#' K = criti$zk[[2]][[1]]
#' Z=K$zonePolygone
#' plotZ(Z)
#' iC=4
#' iE=detZoneEng(iC,Z,K$zoneNModif)
#' envel=calFrame(iC,Z,K$zoneNModif)
#' sp::plot(envel,add=TRUE,col="blue")
#' Qseq = genQseq(qProb,K,mapTest,iC,iE)
#' resi = geozoning:::findCinZ(iC,Z,K,mapTest,Qseq[5],envel)
#' Zopti=geozoning:::zoneQ(resi$contourSp,iC,iE,Z,K)
#' plotZ(Zopti)
#' # not run
zoneQ = function (contourSp,iC,iE,Z,K,simplitol=1e-3)
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
	Zopti = zoneFusion4(Z,iC,iE,simplitol,disp=FALSE)

	if (is.null(Zopti)) return(NULL)
	# englobing zone number may have changed-find it using zone id
	iBig=Identify(idE,Zopti)
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
