## ----echo=TRUE,message=FALSE, warning=FALSE------------------------------
  library(geozoning)

## ----echo=TRUE,message=FALSE, warning=FALSE------------------------------
  #source("srcZ.R") # 
  # result of call to saveZoningFromSimu(seed=20,qProb=0.4,disp=0)
  # Z.Rdata
  # read zone polygons k zones=k polygones
  # read zoning Z, nb of zones NBZ, initial data and seed
  seed<-20
  NBZ<-6
  Vang<-c(0)
  Vanis<-c(1)

## ----echo=TRUE,message=FALSE, warning=FALSE------------------------------
  Z=list()
  Z$zonePolygone=ZPSTest
  # kriged data
  #map=genMap(DataObj,seed=seed,disp=0)
  map=genMap(DataObj=NULL,seed=seed,krig=1,Vpsill=5,Vrange=0.2,Vnugget=0,Vmean=8,disp=0,FULL=FALSE)
  #plotMap(map) 

## ----echo=TRUE,message=FALSE, warning=FALSE------------------------------
  # calNei removez zones with 0 or 1 pt
    #resZ =calNei(Z$zonePolygone,map$krigData,map$krigSurfVoronoi,map$krigN)
  #returns
  #     resZ$zoneN: matrice carree de voisinages entre zones (booleenne)
  #     resZ$zoneNModif: matrice carree de voisinages entre zones (booleenne)
  #     avec FALSE sur la diagonale (zone n'est pas sa propre voisine)
  #     resZ$listZonePoint : liste avec autant d'elements que de zones,
  #     pour chaque zone, les numeros des points kriges qu'elle contient
  #     resZ$meanTot : la moyenne de toutes les valeurs
  #      resZ$meanZone : la moyenne des valeurs de chaque zone
   #zonePolygone =resZ$zonePolygone
  # calcul indicateurs de la partition representee par Z et par les points de map$tabAlea
  # travaille sur les valeurs krigees
  # attention les surfaces de Voronoi sont celles des valeurs krigees
  # renvoie la matrice des distances dans resDistanceA$matDistanceCorr
  # type distance=1
   #resD=calDistance(typedist=1,map$krigData,resZ$listZonePoint,resZ$zoneN,map$krigSurfVoronoi,resZ$meanZone,pErr=pErr)
  # calcul du critere - choix =4
   #resCrit=calCrit(resD$matDistanceCorr,resZ$zoneNModif)
   #print(paste("Crit=",round(resCrit,3)))

## ----session,echo=FALSE,message=FALSE, warning=FALSE---------------------
  sessionInfo()

