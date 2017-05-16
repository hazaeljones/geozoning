###########################################################################
#' calIndice
#'
#' @details description, a paragraph
#' @param tabAlea xxxx
#' @param voisin xxxx
#' @param zoneNmodif xxxx
#' @param listZonePoint xxxx
#' @param listSurf xxxx
#' @param meanZone xxxx
#' @param surfVoronoi xxxx
#' @param matDistanceMoranB xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calIndice=function(tabAlea,voisin,zoneNmodif,listZonePoint,listSurf,meanZone,surfVoronoi,matDistanceMoranB)
###########################################################################
{
  ### moran and  geary correlation indices ###

  nbPoly=length(listZonePoint)
  #indice de moran entre tous les points de toutes les zones voisines(version globale/locale)
  indiceMoranB=calMoranBTot(zoneNmodif,matDistanceMoranB,listSurf)
  vectMoranB=calMoranBLocal(zoneNmodif,matDistanceMoranB,listSurf)

  #calcul des indices moran/geany entre les zones en ne prenant que leurs moyennes
  #global/local
  vectMoranInter=calMoranLoc(zoneNmodif,meanZone,mean(meanZone),listSurf)
  vectGearyInter=calGearyLoc(zoneNmodif,meanZone,mean(meanZone),listSurf)
  indiceMoranInter=calMoranGlo(zoneNmodif,meanZone,mean(meanZone),listSurf)
  indiceGearyInter=calGearyGlo(zoneNmodif,meanZone,mean(meanZone),listSurf)

  #calcul des indices moran/geany sur tous les points à l'intérieur d'une même zone
  vectMoranIntra=list(NULL)
  vectGearyIntra=list(NULL)
  indiceMoranIntra=list(NULL)
  indiceGearyIntra=list(NULL)
  for (i in (1:nbPoly))
  {
    vectMoranIntra[[i]]=calMoranLoc(voisin[listZonePoint[[i]],listZonePoint[[i]]],data.frame(tabAlea[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
    vectGearyIntra[[i]]=calGearyLoc(voisin[listZonePoint[[i]],listZonePoint[[i]]],data.frame(tabAlea[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
    indiceMoranIntra[[i]]=calMoranGlo(voisin[listZonePoint[[i]],listZonePoint[[i]]],data.frame(tabAlea[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
    indiceGearyIntra[[i]]=calGearyGlo(voisin[listZonePoint[[i]],listZonePoint[[i]]],data.frame(tabAlea[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
  }
  ############################################################
  return(list(indiceMoranB=indiceMoranB,vectMoranB=vectMoranB,vectMoranInter=vectMoranInter,vectGearyInter=vectGearyInter,vectMoranIntra=vectMoranIntra,vectGearyIntra=vectGearyIntra,indiceMoranInter=indiceMoranInter,indiceGearyInter=indiceGearyInter,indiceMoranIntra=indiceMoranIntra,indiceGearyIntra=indiceGearyIntra))
}
