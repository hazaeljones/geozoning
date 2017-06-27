###########################################################################
#' calMoranGeary
#'
#' @details description, a paragraph
#' @param spdata xxxx
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
#' data(resZTest)
#' Z=resZTest
#' data(mapTest)
#' ptN=mapTest$krigN
#' spdata=mapTest$krigData
#' zoneNmodif=Z$zoneNmodif
#' listZonePoint=Z$listZonePoint
#' listSurf=Z$listSurf
#' meanZone=Z$meanZone
#' surfVoronoi=mapTest$surfVoronoi
#' calMoranGeary(spdata,ptN,zoneNmodif,listZonePoint,listSurf,meanZone,surfVoronoi,matDistanceMoranB)

calMoranGeary=function(spdata,ptN,zoneNmodif,listZonePoint,listSurf,meanZone,surfVoronoi,matDistanceMoranB)
###########################################################################
{
  ### moran and  geary correlation indices ###

  nbPoly=length(listZonePoint)
  #Moran index for all pts of neighboring zones (global/local version)
  indiceMoranB=calMoranBTot(zoneNmodif,matDistanceMoranB,listSurf)
  vectMoranB=calMoranBLocal(zoneNmodif,matDistanceMoranB,listSurf)

  #Moran index between zones considering aonly the average zone values
  #global/local
  vectMoranInter=calMoranLoc(zoneNmodif,meanZone,mean(meanZone),listSurf)
  vectGearyInter=calGearyLoc(zoneNmodif,meanZone,mean(meanZone),listSurf)
  indiceMoranInter=calMoranGlo(zoneNmodif,meanZone,mean(meanZone),listSurf)
  indiceGearyInter=calGearyGlo(zoneNmodif,meanZone,mean(meanZone),listSurf)

  #Moran index within zone
  vectMoranIntra=list(NULL)
  vectGearyIntra=list(NULL)
  indiceMoranIntra=list(NULL)
  indiceGearyIntra=list(NULL)
  for (i in (1:nbPoly))
  {
    vectMoranIntra[[i]]=calMoranLoc(ptN[listZonePoint[[i]],listZonePoint[[i]]],spdata.frame(spdata[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
    vectGearyIntra[[i]]=calGearyLoc(ptN[listZonePoint[[i]],listZonePoint[[i]]],spdata.frame(spdata[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
    indiceMoranIntra[[i]]=calMoranGlo(ptN[listZonePoint[[i]],listZonePoint[[i]]],spdata.frame(spdata[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
    indiceGearyIntra[[i]]=calGearyGlo(ptN[listZonePoint[[i]],listZonePoint[[i]]],spdata.frame(spdata[listZonePoint[[i]],])[,3],meanZone[i],surfVoronoi[listZonePoint[[i]]])
  }
  ############################################################
  return(list(indiceMoranB=indiceMoranB,vectMoranB=vectMoranB,vectMoranInter=vectMoranInter,vectGearyInter=vectGearyInter,vectMoranIntra=vectMoranIntra,vectGearyIntra=vectGearyIntra,indiceMoranInter=indiceMoranInter,indiceGearyInter=indiceGearyInter,indiceMoranIntra=indiceMoranIntra,indiceGearyIntra=indiceGearyIntra))
}
