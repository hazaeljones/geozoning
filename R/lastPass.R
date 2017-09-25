###############################################
#' lastPass
#'
#' @details description, a paragraph
#' @param carte xxxx
#' @param qProb xxxx
#' @param listOfZ xxxx
#' @param crit xxxx
#' @param cost xxxx
#' @param costL xxxx
#' @param nz xxxx
#' @param mdist xxxx
#' @param pErr xxxx
#' @param optiCrit xxxx
#' @param simplitol xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
lastPass=function(carte,qProb,listOfZ,crit,cost,costL,nz,mdist,pErr,optiCrit,simplitol,disp=F)
###########################################################################
# simply remove zones of last level zonings that are too small and recalculate criteria
#
{
le = length(crit)
lef = length(crit[[le]])
if(lef>1)
	{
	if(disp) print(paste(lef,"zonings in last level"))
	for ( kk in 1:lef)
	    {
	    # for each zoning
	    K1 = listOfZ[[le]][[kk]]
	    Z1 = K1$zonePolygone
	    ZIF = detectSmallZones(Z1,minSizeNG) #returns small zone numbers
	    vNum = ZIF$vectIndex

	    # get small zone ids
	    vId=c()
	    for (jj in vNum) vId=c(vId,getZoneId(Z1[[jj]]))
	    #
	    for (zid in vId)
	    	{
	  	iC = Identify(zid,Z1) # get current zone number
		if(disp) print(paste("in lastPass zid=",zid,"iC=",iC))
	    	Z2 = zoneFusion3(Z1,K1,iC, getNs(K1,iC),simplitol ,carte,F)

		ZK = updateZK(carte,lvlQuant,le,kk,listOfZ, crit,cost,costL,nz,mdist,K1,Z1,Z2,pErr,optiCrit,simplitol) # updates listOfZ, crit,mdist
		#
		listOfZ = ZK$K
		mdist = ZK$mdist
		crit=ZK$crit
		cost=ZK$cost
		costL=ZK$costL
		nz=ZK$nz
		#
	    	K1 = listOfZ[[le]][[kk]]
		Z1=K1$zonePolygone
		}
	    }
	}
	return(list(listOfZ=listOfZ,crit=crit,cost=cost,costL=costL,nz=nz,mdist=mdist))
}
