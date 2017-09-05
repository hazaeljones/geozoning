#' calDistance
#'
#' @details calculates matrix of heterogeneities between neighbour zones.
#' max(sigmai2[i],(fxmean*pErr/100)^2) + max(sigmai2[j],(fymean*pErr/100)^2) + (fxmean-fymean)^2
#' @param typedist default value is 1, other values not implemented yet.
#' @param zoneN zone neighborhood matrix (TRUE values on diagonal), result of call to \code{\link{calNei}}
#' @param listZonePoint list of indices of data points within zones, result of call to \code{\link{calNei}}
#' @param tabVal SpatialPointsDataFrame, contains data points to be used for zoning (spatial coordinates plus attribute values)
#' result of call to \code{\link{genMap}}
#' @param surfVoronoi vector of Voronoi polygon surfaces corresponding to all data points,result of call to \code{\link{genMap}}
#' @param meanZone vector of average attribute values for all zones
#' @param pErr error percentage for correcting distances
#'
#' @return a list with components
#'\describe{
#' \item{matDistance}{matrix of real values, corresponding to heterogeneities between neighbour zones. All other values are set to 0.}
#' \item{matDistanceCorr}{corrected distance matrix using pErr}
#' \item{cost}{sum or errors obtained by replacing all data values within a zone by the zone mean value}
#'}
#'
#' @export
#'
#' @examples
#' # load test map with simulated data
#' data(mapTest)
#' # load zoning results from test file
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
calDistance=function(typedist=1,tabVal=NULL,listZonePoint=NULL,zoneN=NULL,surfVoronoi=NULL,meanZone=NULL,pErr=0.9)
{
  # bch
  # type distance=1
  # tabVal holds data
  #on initialise la matrice qui contiendra les distances
  nbPoly=length(listZonePoint)
  matDistance=matrix(0,nbPoly,nbPoly)
  matDistanceCorr=matDistance

  #keep only lower triangle (matrix is symmetric)
  zoneNInf = (lower.tri(zoneN)*zoneN)==1 #(on repasse en booléens avec ==1)
  diag(zoneNInf)=diag(zoneN)

  #on cree un tableau qui contient les indices des voisinages de zones valant "vrai",sous forme de vecteur
  iN=grep(TRUE,zoneNInf)

  #on transforme l'indice d'un voisinage (correspondant à une position dans un vecteur) en indices x et y (indice de voisinage entre deux zones dans une matrice)
  xIN=(iN-1)%%nbPoly +1
  yIN=(iN-1)%/%nbPoly +1

  if (length(iN)>=1)
  {
    if(typedist==1)
    {
       # compute intra variances sigma i ^2 and zone areas SI
       num=unique(xIN)
       sigmai2=rep(0,nbPoly)
       SI=sigmai2
       for (k in 1:nbPoly)
       {
       res = Sigmai2(k,listZonePoint,tabVal,surfVoronoi,meanZone)
       sigmai2[k]=res$sigmai2
       SI[k]=res$SI
       }
       # compute total cost (per zone)
       cost=sum(sigmai2*SI)/sum(SI)
       #

      for (i in 1:length(iN))
      {
        #on appelle la fonction de calcul de distances de zones 1 (sur valeurs krigées)
	      resDij=DIJ(xIN[i],yIN[i],sigmai2,meanZone,pErr)
        matDistance[xIN[i],yIN[i]]=matDistance[yIN[i],xIN[i]] = resDij$d
	      matDistanceCorr[xIN[i],yIN[i]]=matDistanceCorr[yIN[i],xIN[i]] = resDij$dCorr
      }
    }
     else print("errorr, typedist not implemented")

  } #end more than 1 zone
   return(list(matDistance=matDistance,matDistanceCorr=matDistanceCorr,cost=cost))
}

################################################################
#' Sigmai2
#'
#' @details computes mean (weighted) variance and Voronoi area for zone
#' @param index zone number
#' @param listZonePoint list of pts within each zone
#' @param tabVal SpatalPointsDataFrame holding data
#' @param surfaceVoronoi vector of Voronoi surfaces associated to data values
#' @param meanZone Zone mean values
#'
#' @return a list with components sigmai2 and SI
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Sigmai2(5,K$listZonePoint,mapTest$krigData,mapTest$krigSurfVoronoi,K$meanZone)
#' # not run
Sigmai2=function(index,listZonePoint,tabVal,surfaceVoronoi,meanZone)
################################################################
{
 valZoneI=tabVal@data[listZonePoint[[index]],1]	#liste des valeurs d'altitude des points de la zone
 voroZoneI=surfaceVoronoi[listZonePoint[[index]]]    #vecteur des surfaces des points de la zone indice1
  SI= sum(voroZoneI)
  fxmean= meanZone[[index]]
  sigmai2=sum(((valZoneI-fxmean)^2)*voroZoneI)/SI

  return(list(sigmai2=sigmai2,SI=SI))
}


################################################################
#' DIJ
#'
#' @details description, a paragraph
#' @param i zone index
#' @param j neighbor zone index
#' @param sigmai2 vector of zone variances
#' @param meanZone list of zone mean values 
#' @param pErr tolerance for distance correction
#'
#' @return a list with components d and dCorr
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' nz=length(K$zonePolygone)
#' si2=rep(NA,nz)
#' for (kk in 1:nz){
#' si2[kk]=Sigmai2(kk,K$listZonePoint,mapTest$krigData,
#'         mapTest$krigSurfVoronoi,K$meanZone)$sigmai2
#' }
#' d12=DIJ(1,2,si2,K$meanzone,0.9)
#' # not run
DIJ=function(i,j,sigmai2,meanZone,pErr)
{
  fxmean = meanZone[[i]]
  fymean = meanZone[[j]]
  res = sigmai2[i] + sigmai2[j] + (fxmean-fymean)^2
  resCorr = max(sigmai2[i],(fxmean*pErr/100)^2) + max(sigmai2[j],(fymean*pErr/100)^2) + (fxmean-fymean)^2

  return(list(d=res,dCorr=resCorr))
}

#################################################################################
#' distanceNormalisationSqrt
#'
#' @details normalize all MIJ terms of the distance matrix by dividing it by square root of diagonal terms MII*MJJ
#' @param matDistance distance matrix as returned by a call to calDistance
#'
#' @return a normalized distance matrix
#'
#' @export
#'
#' @examples
#' # load test map with simulated data
#' data(mapTest)
#' # load zoning results from test file
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,
#'        mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' distanceNormalisationSqrt(resD$matDistanceCorr)
#' # not run
distanceNormalisationSqrt=function(matDistance)
{
  #matrice termeNormal[i,j]=1/sqrt(dii^2*djj^2)
  nbPoly=nrow(matDistance)
  termeNormal=matrix(rep(diag(matDistance),nbPoly),nbPoly,nbPoly)
  termeNormal=1/sqrt(termeNormal*t(termeNormal))
  matDistanceNorm=matDistance*termeNormal

  return(matDistanceNorm)
}
#################################################################################
#' distanceNormalisationSum
#'
#' @details normalize all MIJ terms of the distance matrix by dividing it by sum  of squared diagonal terms sum(MII^2+MJJ^2)
#' @param matDistance distance matrix as returned by a call to calDistance
#'
#' @return a normalized distance matrix
#'
#' @export
#'
#' @examples
#' # load test map with simulated data
#' data(mapTest)
#' # load zoning results from test file
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,
#'        mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' distanceNormalisationSqrt(resD$matDistanceCorr)
#' # not run
distanceNormalisationSum=function(matDistance)
{
  nbPoly=nrow(matDistance)
  normalTerm=matrix(rep(diag(matDistance),nbPoly),nbPoly,nbPoly)
  normalTerm=1/(normalTerm+t(normalTerm))
  matDistanceNorm=2*matDistance*normalTerm

  return(matDistanceNorm)
}
