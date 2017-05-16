
#################################################################################
#' zoneModifnonIso
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param K xxxx
#' @param qProb xxxx
#' @param map xxxx
#' @param zoneClose xxxx
#' @param iC xxxx
#' @param simplitol xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
zoneModifnonIso=function(Z,K,qProb,map,zoneClose,iC,simplitol,disp=0)
#################################################################################
#
{
	    Zopti=NULL
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
                            Zopti = optiRG(Z,K,map,iC,indZC,simplitol,disp)
			    #
                      	   if (is.null(Zopti))
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
		if (!is.null(Zopti))
		   print(paste("junction of non isolated zone: ",iC,"and zone",indZC))
		   else
		    print(paste("no junction of non isolated zone: ",iC,"- no close zone with same lab or junction intersects other zones"))
		}
         return(Zopti)
}
