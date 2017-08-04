####################################################################
#' local Geary criteria
#'
#' @details computes local Geary indices
#' @param matN neighborhood (zone or point) matrix
#' @param vectMean vector of mean zone values
#' @param meanTot global mean
#' @param vectSurface vector of zone areas
#'
#' @return a vector of local Geary criteria
#'
#' @export
#' @examples
#' K=resZTest
#' zoneA=sapply(K$zonePolygone,gArea)
#' calGearyLoc(K$zoneNModif,K$meanZone,K$meanTot,zoneA)
#' # not run
calGearyLoc=function(matN,vectMean,meanTot,vectSurface)
{
#--------------------------------------------------------------------------------
  # number of zones or points
  nbZones=length(vectMean)
  #transform logical neighborhood into integer
  matNum=matrix(as.numeric(matN),nbZones,nbZones)

  #matrix with vectMean columns
  matMean=vectMean%*%t(rep(1,nbZones))

  # compute matWZ[i,j]=wij*(zi-zj)^2*Sj
  #                                 wij=1 si i et j sont Ns,0 sinon
  #                                 zi=meanI
  #                                 Sj=surfaceJ
  matWZ=matNum*((matMean-t(matMean))^2)*vectSurface

  #compute indiceMoranLoc[i]=nbzones*zi*SUM(wij zj)/SUM(wij)*SUM(zi-z)^2
  indiceGearyLoc=(apply(FUN=sum,matWZ,MARGIN=1))/(((vectMean-meanTot)^2)*apply(FUN=sum,MARGIN=1,matNum*t(vectSurface%*%t(1:nbZones)) ))

  return(indiceGearyLoc)
}


####################################################################
#' computes Moran criterion on whole zoning
#'
#' @details computes Moran criterion on zoning
#' @param NZone xxxx
#' @param matDistanceMoranB xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranBTot=function(NZone,matDistanceMoranB,vectSurface)
{
#################################################################################
#-------------------------------------------------------------------------------

  #nombre de zones
  n=length(diag(matDistanceMoranB))
  #nombre de coefficients non nuls sur l'extradiagonale(=somme des wij)
  m=length(grep(TRUE,NZone))
  #on applique à la matrice des distances un masque pour enlever les termes diagonaux
  matDistanceModif=matDistanceMoranB*NZone
  # numerator: si pour i et j Ns
  #                      Mij=distance(ij)*surfacei*surfacej
  #                      numMoran=somme(Mij)
  numMoran=sum(matDistanceModif * (vectSurface%*%t(vectSurface)))

  # denominator:Mi=distance(ii)*surface(i)*somme(surfaces des Ns de i)
  #                       denomMoran=somme(Mi)
  denomMoran=sum( diag(matDistanceMoranB) * vectSurface * apply(FUN=sum,MARGIN=1,NZone*t(vectSurface%*%t(1:n))) )
  iMoran=(n/m)*numMoran/denomMoran
  return(iMoran)
}

####################################################################
#' compute local Moran indices (per zone)
#'
#' @details description, a paragraph
#' @param NZone xxxx
#' @param matDistanceMoranB xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranBLocal=function(NZone,matDistanceMoranB,vectSurface)
{
#################################################################################
  nbZones=length(vectSurface)
  # numérateur: si pour i et j Ns
  #             Mij=dij*surface zone i* surface zone j,
  #             numérateur=somme(Mij) sur les lignes
  #dénominateur:si Mi=dii*somme(surfaces des Ns de i)
  #             denominateur=somme (Mi)
  vectMoran=((apply( FUN=sum,MARGIN=1 , matDistanceMoranB * NZone * t(vectSurface%*%t(1:nbZones) )))
             /(diag(matDistanceMoranB) * apply( FUN=sum,MARGIN=1 , NZone * t(vectSurface%*%t(1:nbZones)) )))
  return(vectMoran)
}

####################################################################
#' computes specific Moran criterion
#'
#' @details description, a paragraph
#' @param matNZone xxxx
#' @param vectMean xxxx
#' @param meanTot xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranGlo=function(matNZone,vectMean,meanTot,vectSurface)
{
#-------------------------------------------------------------------------------
################################################################################
#nombre de zones ou de points
  nbZones=length(vectMean)

  #transformation de matrice de booleéns en entiers
  matNum=matrix(as.numeric(matNZone),nbZones,nbZones)
  #somme de la matrice des Nages
  sommeNage=sum(matNum)

  #numérateur:si pour i et j Ns
  #             Mij=(mI-mGlob)*(mJ-mGlob)*surfaceI*surfaceJ
  #             numerateur=somme(Mij)
  #denominateur:si Mi=(mI-mGlob)^2 * surfaceI*somme(surfaces des Ns de I)
  #             denominateur=somme(Mi)

  indiceMoranGlo=(sum( ((vectMean-meanTot)%*%t(vectMean-meanTot)) * matNum * (vectSurface%*%t(vectSurface)) )
                  /sum(((vectMean-meanTot)^2)*vectSurface*apply(t(vectSurface%*%t(1:nbZones)*matNum),FUN=sum,MARGIN=1)) )
  return(indiceMoranGlo)
}

####################################################################
#' calMoranLoc
#'
#' @details description, a paragraph
#' @param matN xxxx
#' @param vectMean xxxx
#' @param meanTot xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranLoc=function(matN,vectMean,meanTot,vectSurface)
{
#--------------------------------------------------------------------------------#################################################################################
  #nombre de zones ou de points
  nbZones=length(vectMean)
  #transform logical meighborhood into integer
  matNum=matrix(as.numeric(matN),nbZones,nbZones)

  #on compute la matrice,matWZ[i,j]=wij*zj*Sj
  #                                 wij=1 si i et j sont Ns,0 sinon
  #                                 zj=mJ-mGlob
  #                                 Sj=surfaceJ
  matWZ=t(t(matNum)*((vectMean-meanTot)*vectSurface))

  #on compute le vecteur indiceMoranLoc[i]=nbzones*zi*SUM(wij zj)/SUM(wij)*SUM(zi-z)^2
  indiceMoranLoc=((vectMean-meanTot)*apply(FUN=sum,matWZ,MARGIN=1))/(((vectMean-meanTot)^2)*apply(FUN=sum,MARGIN=1,matNum*t(vectSurface%*%t(1:nbZones)) ))

  return(indiceMoranLoc)
}


####################################################################
#' calGearyGlo
#'
#' @details computes global Geary criterion
#' @param matN xxxx
#' @param vectMean xxxx
#' @param meanTot xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calGearyGlo=function(matN,vectMean,meanTot,vectSurface)
{
#--------------------------------------------------------------------------------#################################################################################
  #nombre de zones ou de points
  nbZones=length(vectMean)

  #transform logical meighborhood into integer
  matNum=matrix(as.numeric(matN),nbZones,nbZones)
  sommeNage=sum(matNum)

  #matrice ou chaque colone est vectMean
  matMean=vectMean%*%t(rep(1,nbZones))

  #numerator:
  #           Mij=(mI-mJ)^2 * surfaceI*surfaceJ
  #           numerator=sum(Mij)*(number of zones -1)
  #denominator:
  #           Mi=(mI-mGlob)*somme(surface of Ns  of I)
  #           denominatorr=sum(Mi)*(2*number of neighborhoods)

  indiceGeary=(((nbZones-1)/2*sommeNage)*sum(matNum * ((matMean-t(matMean))^2)*(vectSurface%*%t(vectSurface)) )
               /sum( ((vectMean-meanTot)^2) * apply(FUN=sum,MARGIN=1 , matNum * t(vectSurface%*%t(1:nbZones))) ))

  return(indiceGeary)
}
