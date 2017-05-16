#########################################
#' computeZoneInt
#'
#' @details description, a paragraph
#' @param qi xxxx
#' @param Z0 xxxx
#' @param i0 xxxx
#' @param map xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
computeZoneInt = function(qi,Z0,i0,map)
#########################################
{

	#x11()
	res=func1Zonage(qi, map,display=F)
	Z=res$zonePolygone

	NZ=length(Z)
	zn=list()
	n=vector()
	for (i1 in 1:NZ)
	{
		ok=gWithin(Z[[i1]],gBuffer(Z0[[i0]],width=0.01))
		if (ok)
		{
		zn=append(zn,Z[[i1]])
		n=c(n,i1)
		}
	}
	print(n)
	# prendre la plus grande
	area=sapply(zn,gArea)
	ind=rev(order(area))
	Z[[i0]]=zn[[ind[1]]]

	return(Z[[i0]])
}
