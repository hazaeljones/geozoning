#########################################
# ALL IN COMMENTS 19/05/2017, fUNCTION NOT USED
# computeZoneExt
#
# @details description, a paragraph
# performs a new zoning on map with quantile probability qi
# and returns the smallest zone from this new zoning that includes zone #i0 from initial zoning
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
# computeZoneExt = function(qi,Z0,i0,map)
# #########################################
# {
# 	res=initialZoning(qi, map)
# 	Z=res$zonePolygone
#
# 	NZ=length(Z)
# 	zn=list()
# 	n=vector()
# 	for (i1 in 1:NZ)
# 	{
# 		ok=gContains(gBuffer(Z[[i1]],width=0.001),Z0[[i0]])
# 		if (ok)
# 		{
# 		zn=append(zn,Z[[i1]])
# 		n=c(n,i1)
# 		}
# 	}
# 	print(n)
# 	#take the smallest one
# 	area=sapply(zn,gArea)
# 	ind=order(area)
#	if(length(zn)>0)
#	{
# 	Z[[i0]]=zn[[ind[1]]]
# 	return(Z[[i0]])
#	}
#	else return(NULL)
# }
