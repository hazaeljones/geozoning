###############################################
#' lastPass
#'
#' @details description, a paragraph
#' @param map object returned by function genMap or genMapR
#' @param qProb probability vector used to generate quantile values
#' @param listOfZ list of zoning objects (such as returned by calNei function)
#' @param crit criterion value list
#' @param cost cost value list
#' @param costL cost per lable value list
#' @param nz number of zones list
#' @param mdist distance matrix list
#' @param pErr equality tolerance for distance calculations
#' @param optiCrit criterion choice
#' @param minSize zone area threshold under which a zone is too small to be manageable
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param disp 0: no info, 1: detailed info
#'
#' @return a list with components
#'\describe{
#' \item{listZ}{list of zoning objects (such as returned by calNei function)
#' \item{crit}{criterion value list}
#' \item{cost}{cost value list}
#' \item{costL}{cost per label value list}
#' \item{nz}{number of zones list}
#' \item{mdist{distance matrix list}
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' criti=correctionTree(c(0.4,0.7),mapTest,LASTPASS=FALSE) # run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7
#' Z=criti$zk[[1]][[1]]$zonePolygone #initial zoning
#' printZsurf(Z) # 8 zones with 2 small zones (7 and 8)
#' newRes=lastPass(mapTest,c(0.4,0.7),criti$zk[1],criti$criterion[1],criti$cost[1],criti$costL[1],criti$nz[1],criti$mdist[1])
#' newZ=newRes$listOfZ[[1]][[1]]$zonePolygone
#' printZsurf(newZ) # 6 zones
#' # not run
lastPass=function(map,qProb,listOfZ,crit,cost,costL,nz,mdist,pErr=0.9,optiCrit=2,minSize=0.012,simplitol=1e-3,disp=F)
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
		if(iC==0) next
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
		ZK = updateZK(map,qProb,le,kk,listOfZ, crit,cost,costL,nz,mdist,K0,Z3,pErr,optiCrit,simplitol) # updates listOfZ, crit,mdist
		
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
