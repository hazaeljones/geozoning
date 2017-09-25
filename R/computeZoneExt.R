#########################################
# ALL IN COMMENTS 19/05/2017, fUNCTION NOT USED
# computeZoneExt
#
# @details description, a paragraph
# @param qi xxxx
# @param Z0 xxxx
# @param i0 xxxx
# @param map xxxx
#
# @return a plot
#
# @export
#
# @examples
# # not run
# computeZoneExt = function(qi,Z0,i0,map)
# #########################################
# {
#   # IS 19/05/2017: add comment for x11
# 	#x11()
# 	res=func1Zonage(qi, map,display=TRUE)
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
# 	Z[[i0]]=zn[[ind[1]]]
#
# 	return(Z[[i0]])
# }
