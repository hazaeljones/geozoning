
#################################################################################
#' zoneModifnonIso
#'
#' @details modify non isolated zone (depends on distIsoZ parameter) so that it is joined to the closest neighbour zone with the same label.
#' @param K zoning object (such as returned by calNei function)
#' @param qProb probability vector used to generate quantile values
#' @param map object returned by function genMap or genMapR
#' @param zoneClose indices of close zones
#' @param iC current zone index
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param disp 0: no info, 1: detailed info
#'
#' @return a zoning object
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.2,0.5)
#' ZK = initialZoning(qProb, mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' plotZ(Z) 
#' resP=detZoneClose(6,Z,K$zoneNModif) # zone 6 is close to zone 5 and zone 7
#' zoneClose = resP$zoneClose
#' kmi = zoneModifnonIso(K,qProb,mapTest,zoneClose,6,disp=1)
#' plotZ(kmi$zonePolygone) # zones 6 and 7 are joined into new zone 6
#' # now it is the turn of zone 5
#' Z=kmi$zonePolygone
#' resP=detZoneClose(5,Z,kmi$zoneNModif) # zone 5 is close to zone 7 and zone 6
#' kmi2 = zoneModifnonIso(kmi,qProb,mapTest,resP$zoneClose,5,disp=1)
#' plotZ(kmi2$zonePolygone) # zones 5 and 6 are joined into new zone 5
#' # not run
zoneModifnonIso=function(K,qProb,map,zoneClose,iC,simplitol=1e-3,disp=0)
#################################################################################
#
{
	    Kopti=NULL
	    Z=K$zonePolygone
            # loop on  zones close to current one
            if (length(zoneClose) == 0)
	       {
		if(disp) print(paste("zone: ",iC,"is isolated - no possible junction"))
		return(NULL)
		}
		for (ip in 1:length(zoneClose))
		{
		 indZC = zoneClose[[ip]]
		 # swap indices if zoneClose is smaller
                 if (gArea(Z[[indZC]])<=gArea(Z[[iC]]))
                       {
		       if(disp) print("swap")
		       tmp=iC
		       iC = indZC
		       indZC = tmp
		     }
		  # check if same label
                  if (K$lab[iC] == K$lab[indZC])
                          {
                  # then group them
                            Kopti = optiRG(K,map,iC,indZC,simplitol,disp)
			    #
                      	   if (is.null(Kopti))
			       next
			    else
				{
				break
				}
                          }
                        else next
			#if different label, next iteration
	           }
         if(disp)
		{
		if (!is.null(Kopti))
		   print(paste("junction of non isolated zone: ",iC,"and zone",indZC))
		   else
		    print(paste("no junction of non isolated zone: ",iC,"- no close zone with same lab or junction intersects other zones"))
		}
         return(Kopti)
}
