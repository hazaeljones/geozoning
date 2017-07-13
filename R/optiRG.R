###########################################################
#' optiRG
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param K xxxx
#' @param map xxxx
#' @param iC xxxx
#' @param iZC xxxx
#' @param simplitol xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
optiRG = function(Z,K,map,iC, iZC,simplitol,disp=0)
###########################################################
{
#regroup (aggregate) 2 close zones iC and iZC
# swap zones -> petite en premier
#
 if (gArea(Z[[iZC]])<=gArea(Z[[iC]]))#zone proche plus petite - on echange les 2 pour le regroupement
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


  #Attention il peut y avoir un bug à cause de ça. En effet lenveloppe convexe peut empieter sur une zone incluse dans la zone
  #dindice courant. Il peut en résulter un changement daffectation de point qui fait que la zone incluse na plus de points.
  #cette zone se retrouve donc eliminée par la fonction calNei au prochain appel. --> erreur
  #solution : il faut modifier foncRegroup pour faire en sorte que si on se trouve dans un tel cas de figure, on fasse quelque chose de
  #different.

  polyUni=gConvexHull(polyUni)



  #fusion + removal + ID management
  # assign current zone id to new merged zone
  # to avoid handling new merged zone again in small zone loop
  idZC= getId(Zopti,iZC)
  ie0 = getId(Zopti,iZE)

  #merge zone iZC with polyUni-keep only envelope
  Zopti[[iZC]]= gUnion(polyUni,Zopti[[iZC]])
  Zopti[[iZC]] = gBuffer(Zopti[[iZC]],byid=TRUE,width=0)

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
	if (disp >0) print("pb in optiRG")
	Zoptipb <<-Zopti
	Zpb <<- Z
	indpbZE <<- iZE
	indpbC <<- iC
	indpbP <<- iZC
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
# find merged zone number in new zoning
  index=findNumZ(Zopti,idZC)
# must not intersect with other zones except itself and included zones
  inter=testInterSpeZ1(Zopti,index)
#
  if(inter) Zopti=NULL

  return(Zopti)
}

