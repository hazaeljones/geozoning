###########################################################
#' optiRG join two zones close to each other
#'
#' @details Within a zoning, two zones close to each other are geometrically joined. The zoning is updated accordingly. If the zone resulting from the junction is not valid, i.e. if it crosses another zone, the function returns NULL. 
#' @param K zoning object (such as returned by calNei function)
#' @param map object returned by function genMap
#' @param iC first zone
#' @param iZC second zone
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param disp 0: no info, 1: detailed info
#'
#' @return a zoning object
#' @importFrom rgeos createSPComment
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.2,0.5)
#' ZK = initialZoning(qProb, mapTest)
#' K=ZK$resZ
#' Z=K$zonePolygone
#' plotZ(K$zonePolygone) # zoning
#' kmi=optiRG(K,mapTest,6,7,disp=1)
#' #zones 6 and 7 are joined into new zone 6
#' sp::plot(kmi$zonePolygone[[6]],col="red",add=TRUE)
#' # not run
optiRG = function(K,map,iC,iZC,simplitol=1e-3,disp=0)
###########################################################
{
#regroup (aggregate) 2 close zones iC and iZC
# swap zones -> smaller one first
#
 Z=K$zonePolygone
 qProb=K$qProb
 if (gArea(Z[[iZC]])<=gArea(Z[[iC]]))
                       {
		       tmp=iC
		       iC = iZC
		       iZC = tmp
		       }

 # zone englobing current zone
  iZE = detZoneEng(iC,Z,K$zoneNModif)
  if (iZE == 0) return(NULL)
#
# pt in close zone near current zone
  closePt = getClosePt(Z,iC,iZC)
# find 2 pts to join zones
  res  = interZoneC(Z, iC,iZC,closePt)
  if (is.null(res)) return(NULL)
  #
  spi=res$spi
  ord=res$ord # order for pts
  if(is.null(spi)) return(NULL)

  # there is an intersection
  #Union of current zone and extension
  Zopti=Z
  polyUni =gUnion(Zopti[[iC]],spi)
  #add pts of close zone that lie between intersection pts
  polyUni@pointobj@coords = rbind(polyUni@pointobj@coords,Zopti[[iZC]]@polygons[[1]]@Polygons[[1]]@coords[ord,])

  polyUni=gConvexHull(polyUni)
  #merge zone iZC with polyUni-keep only envelope
  tmpZ = gUnion(polyUni,Zopti[[iZC]])
  tmpZ = gBuffer(tmpZ,byid=TRUE,width=0)
  # tmpZ MUST NOT INTERSECT WITH OTHER ZONES EXCEPT THE 2 ZONES TO JOIN
  # PLUS THE ENGLOBING ONE
  nother=1:length(Z)
  nother=nother[-match(iC,nother)]
  nother=nother[-match(iZC,nother)]
  nother=nother[-match(iZE,nother)]
  for (kk in nother)
  {
  if(gIntersects(tmpZ,Z[[kk]])) return(NULL)
  }
  #fusion + removal + ID management
  # assign current zone id to new merged zone
  # to avoid handling new merged zone again in small zone loop
  idZC= getId(Zopti,iZC)
  ie0 = getId(Zopti,iZE)

  # replace close zone with new one including the 2 joined zones
  Zopti[[iZC]] = tmpZ
  # reassign id to merged zone
  Zopti=setId(Zopti,iZC,idZC)

  #redefine englobing zone
  Zopti[[iZE]] = gDifference(Zopti[[iZE]],Zopti[[iZC]])
  Zopti[[iZE]] = cleanSp(Zopti[[iZE]])
  Zopti[[iZE]] = gSimplify(Zopti[[iZE]],simplitol)
  # reassign id to englobing zone
  Zopti=setId(Zopti,iZE,ie0)

# manage newly created zones in englobing zone if any
# separationPoly separates polygons that contain holes and polygons that are just contours, it returns non holes polygons
  lp=separationPoly(Zopti[[iZE]])
  le=length(lp)
  for (k in 1:le)
      {
      lp[[k]]=createSPComment(lp[[k]])
      }


  if (le >2) # msg from separationPoly, save elements to test and exit
  # if more than one (non hole) polygon there is an intersection pb
  {
	if (disp >0) print("problem in optiRG - no junction")
#	Zoptipb <<-Zopti
#	Zpb <<- Z
#	indpbZE <<- iZE
#	indpbC <<- iC
#	indpbP <<- iZC
	return(NULL)
  }

  #

  if (le >1)
  {
	 Zopti[[iZE]]=lp[[1]] # new englobing zone
  	 Zopti=setId(Zopti,iZE,ie0)

	for (k in le:2)
  	{
		# add new zone
		j=length(Zopti)
		Zopti[[j+1]]=lp[[k]]
		newid=1+maxId(Zopti) #new zone
		Zopti=setId(Zopti,j+1,newid)
  		}
  }
# remove current zone NOW  !!!!
  Zopti[[iC]]=NULL

# manage intersection of new united zone and other zones
# test: if new regrouped zone intersects with others (except englobing zone, or inner zone), pb, do not keep it
# loop on all zones
# create comments for holes
  Zopti=crComment(Zopti)
  Kopti = calNei(Zopti,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol)
  Zopti = Kopti$zonePolygone
  Kopti=trLabZone(K,Kopti,map,qProb,disp=0)
  Kopti$qProb=K$qProb
# find merged zone number in new zoning
  index=Identify(idZC,Zopti)
# must not intersect with other zones except itself and included zones
# already done with tmpZ
#  inter=testInterSpeZ1(Zopti,index)
#
#  if(inter) Kopti=NULL

  return(Kopti)
}

