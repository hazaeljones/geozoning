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
lastPass=function(map,qProb,listOfZ,crit,cost,costL,nz,mdist,pErr,optiCrit,minSize,simplitol,disp=F)
###########################################################################
# simply remove zones of last level zonings that are too small and recalculate criteria
#
{
	le = length(crit)
	lef = length(crit[[le]])
	if(lef<1) return(list(listOfZ=listOfZ,crit=crit,cost=cost,costL=costL,nz=nz,mdist=mdist)) # protection against empty result in last level

	if(disp) print(paste(lef,"zonings in last level"))
	for ( kk in 1:lef) # for each final sublevel
	    {
	    # for each zoning
	    K1 = listOfZ[[le]][[kk]]
	    Z1 = K1$zonePolygone
	    ZIF = detectSmallZones(Z1,minSize) #returns small zone numbers
	    if(disp) print(paste("lef=",lef,"kk=",kk))
	    if(disp) print(ZIF,collapse="")
	   
	    vNum = ZIF$vectIndex

	    # get small zone ids
	    vId=c()
	    for (jj in vNum) vId=c(vId,getZoneId(Z1[[jj]]))
	    #
	    K0=K1
	    Z0=Z1
	    for (zid in vId)
	    	{
	  	iC = Identify(zid,Z1) # get current zone number
		if(disp) print(paste("in lastPass zid=",zid,"iC=",iC))
	    	Z2 = zoneFusion3(K1,iC, getNs(K1$zoneNModif,iC) ,map,minSize,simplitol)
		K2=calNei(Z2,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol)
  		Z2=K2$zonePolygone
		K3=trLabZone(K1,K2,Z1,Z2,map,qProb,disp=0)
		Z3=K3$zonePolygone
		K1=K3
		Z1=Z3
		} #end loop on small zone ids
		# keep the final result
	     if(length(vNum)>0) # then update zoning kk
	        {
		ZK = updateZK(map,qProb,le,kk,listOfZ, crit,cost,costL,nz,mdist,K0,Z0,Z3,pErr,optiCrit,simplitol) # updates listOfZ, crit,mdist
		
		listOfZ= ZK$listOfZ
		mdist= ZK$mdist
		crit=ZK$crit
		cost=ZK$cost
		costL=ZK$costL
		nz=ZK$nz
		} #end update zoning kk
	    } #end loop on last level zonings

	return(list(listOfZ=listOfZ,crit=crit,cost=cost,costL=costL,nz=nz,mdist=mdist))
}
