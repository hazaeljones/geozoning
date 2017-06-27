#########################################
# ALL IN COMMENTS 19/05/2017, fUNCTION NOT USED
# computeZoneInt
#
# @details description, a paragraph
# performs a new zoning on map with quantile probability qi
# and returns the biggest zone from this new zoning that is included in zone #i0 from initial zoning
# @param qi probability
# @param Z0 zoning
# @param i0 zone number in zoning Z0
# @param map map on which perform zoning
#
# @return a zone (SpatialPolygons object)
#
# @export
#
# @examples
# # not run
# computeZoneInt = function(qi,Z0,i0,map)
# #########################################
# {
#
# 	res=initialZoning(qi, map)
# 	Z=res$zonePolygone
#
# 	NZ=length(Z)
# 	zn=list()
# 	n=vector()
# 	for (i1 in 1:NZ)
# 	{
# 		ok=gWithin(Z[[i1]],gBuffer(Z0[[i0]],width=0.01))
# 		if (ok)
# 		{
# 		zn=append(zn,Z[[i1]])
# 		n=c(n,i1)
# 		}
# 	}
# 	print(n)
# 	# keep the biggest zone within
# 	area=sapply(zn,gArea)
# 	ind=rev(order(area))
# 	if (length(zn)>0)
#	{
#	Z[[i0]]=zn[[ind[1]]]
# 	return(Z[[i0]])
#	}
#	else return(NULL)
# }
