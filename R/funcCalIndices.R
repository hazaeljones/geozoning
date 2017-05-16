####################################################################
#' vecteur d'indices de Geary locaux
#'
#' @details description, a paragraph
#' @param matVoisin xxxx
#' @param vectMoy xxxx
#' @param moyTot xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calGearyLoc=function(matVoisin,vectMoy,moyTot,vectSurface)
{
#---------------------------------------------------------------------------------------------------------------------------------#
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#
#entrée:matrice de voisinages(zones ou points),vecteur des moyennes par zone, moyenne totale
#sortie:vecteur d'indices de Geary locaux
#---------------------------------------------------------------------------------------------------------------------------------#
  #nombre de zones ou de points
  nbZones=length(vectMoy)
  #transformation de matrice de voisinage booleéns->entiers
  matNum=matrix(as.numeric(matVoisin),nbZones,nbZones)

  #matrice ou chaque colone est vectMoy
  matMoy=vectMoy%*%t(rep(1,nbZones))

  #on compute la matrice,matWZ[i,j]=wij*(zi-zj)^2*Sj
  #                                 wij=1 si i et j sont voisins,0 sinon
  #                                 zi=moyenneI
  #                                 Sj=surfaceJ
  matWZ=matNum*((matMoy-t(matMoy))^2)*vectSurface

  #on compute le vecteur indiceMoranLoc[i]=nbzones*zi*SUM(wij zj)/SUM(wij)*SUM(zi-z)^2
  indiceGearyLoc=(apply(FUN=sum,matWZ,MARGIN=1))/(((vectMoy-moyTot)^2)*apply(FUN=sum,MARGIN=1,matNum*t(vectSurface%*%t(1:nbZones)) ))

  return(indiceGearyLoc)
}


####################################################################
#' fonction qui compute un critère de moran sur toute la carte
#'
#' @details description, a paragraph
#' @param voisinZone xxxx
#' @param matDistanceMoranB xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranBTot=function(voisinZone,matDistanceMoranB,vectSurface)
{
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#

#fonction qui compute un critère de moran sur toute la carte
#entrée:matrice des voisinages sans la diagonale(matrix),matrice des distance entre zones(matrix) ( de Moran)
#sortie:indice de Moran sur le découpage(numeric)
#---------------------------------------------------------------------------------------------------------------------------------#
  #nombre de zones
  n=length(diag(matDistanceMoranB))
  #nombre de coefficients non nuls sur l'extradiagonale(=somme des wij)
  m=length(grep(TRUE,voisinZone))
  #on applique à la matrice des distances un masque pour enlever les termes diagonaux
  matDistanceModif=matDistanceMoranB*voisinZone
  #cal du numérateur: si pour i et j voisins
  #                      Mij=distance(ij)*surfacei*surfacej
  #                      numMoran=somme(Mij)
  numMoran=sum(matDistanceModif * (vectSurface%*%t(vectSurface)))

  #cal du denominateur:Mi=distance(ii)*surface(i)*somme(surfaces des voisins de i)
  #                       denomMoran=somme(Mi)
  denomMoran=sum( diag(matDistanceMoranB) * vectSurface * apply(FUN=sum,MARGIN=1,voisinZone*t(vectSurface%*%t(1:n))) )
  iMoran=(n/m)*numMoran/denomMoran
  return(iMoran)
}

####################################################################
#' fonction calant un vecteur d'indices de moran locaux(= par zone)
#'
#' @details description, a paragraph
#' @param voisinZone xxxx
#' @param matDistanceMoranB xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranBLocal=function(voisinZone,matDistanceMoranB,vectSurface)
{

  #---------------------------------------------------------------------------------------------------------------------------------#
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#

#fonction calant un vecteur d'indices de moran locaux(= par zone)
#entrée:matrice des voisinages,matrice des distance de Moran entre zones (matrix)
#sortie:vecteur des indices de Moran locaux sur le découpage(numeric)
#---------------------------------------------------------------------------------------------------------------------------------#
  nbZones=length(vectSurface)
  # numérateur: si pour i et j voisins
  #             Mij=dij*surface zone i* surface zone j,
  #             numérateur=somme(Mij) sur les lignes
  #dénominateur:si Mi=dii*somme(surfaces des voisins de i)
  #             denominateur=somme (Mi)
  vectMoran=((apply( FUN=sum,MARGIN=1 , matDistanceMoranB * voisinZone * t(vectSurface%*%t(1:nbZones) )))
             /(diag(matDistanceMoranB) * apply( FUN=sum,MARGIN=1 , voisinZone * t(vectSurface%*%t(1:nbZones)) )))
  return(vectMoran)
}

####################################################################
#' fonction calant un indice de Moran(adapté) entre les différentes zone,ou les différents points
#'
#' @details description, a paragraph
#' @param matVoisinZone xxxx
#' @param vectMoy xxxx
#' @param moyTot xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranGlo=function(matVoisinZone,vectMoy,moyTot,vectSurface)
{
#---------------------------------------------------------------------------------------------------------------------------------#
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#
#fonction calant un indice de Moran(adapté) entre les différentes zone,ou les différents points
#entrée:matrice des voisinages de zones OU POINTS (bien qu'il soit précisé matVoisin"Zone") (diagonale=False),
#vecteur des moyennes des zones ou des valeurs des points, moyenne totale,
#sortie:indice associé à un découpage spécifique de la carte
#---------------------------------------------------------------------------------------------------------------------------------#
#nombre de zones ou de points
  nbZones=length(vectMoy)

  #transformation de matrice de booleéns en entiers
  matNum=matrix(as.numeric(matVoisinZone),nbZones,nbZones)
  #somme de la matrice des voisinages
  sommeVoisinage=sum(matNum)

  #numérateur:si pour i et j voisins
  #             Mij=(moyenneI-moyenneTotale)*(moyenneJ-moyenneTotale)*surfaceI*surfaceJ
  #             numerateur=somme(Mij)
  #denominateur:si Mi=(moyenneI-moyenneTotale)^2 * surfaceI*somme(surfaces des voisins de I)
  #             denominateur=somme(Mi)

  indiceMoranGlo=(sum( ((vectMoy-moyTot)%*%t(vectMoy-moyTot)) * matNum * (vectSurface%*%t(vectSurface)) )
                  /sum(((vectMoy-moyTot)^2)*vectSurface*apply(t(vectSurface%*%t(1:nbZones)*matNum),FUN=sum,MARGIN=1)) )
  return(indiceMoranGlo)
}

####################################################################
#' calMoranLoc
#'
#' @details description, a paragraph
#' @param matVoisin xxxx
#' @param vectMoy xxxx
#' @param moyTot xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMoranLoc=function(matVoisin,vectMoy,moyTot,vectSurface)
{
#---------------------------------------------------------------------------------------------------------------------------------#
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#
#entrée:matrice de voisinages(zones ou points),vecteur des moyennes par zone, moyenne totale
#sortie:vecteur d'indices de Moran locaux
#---------------------------------------------------------------------------------------------------------------------------------#
  #nombre de zones ou de points
  nbZones=length(vectMoy)
  #transformation de matrice de voisinage booleéns->entiers
  matNum=matrix(as.numeric(matVoisin),nbZones,nbZones)

  #on compute la matrice,matWZ[i,j]=wij*zj*Sj
  #                                 wij=1 si i et j sont voisins,0 sinon
  #                                 zj=moyenneJ-moyenneTotale
  #                                 Sj=surfaceJ
  matWZ=t(t(matNum)*((vectMoy-moyTot)*vectSurface))

  #on compute le vecteur indiceMoranLoc[i]=nbzones*zi*SUM(wij zj)/SUM(wij)*SUM(zi-z)^2
  indiceMoranLoc=((vectMoy-moyTot)*apply(FUN=sum,matWZ,MARGIN=1))/(((vectMoy-moyTot)^2)*apply(FUN=sum,MARGIN=1,matNum*t(vectSurface%*%t(1:nbZones)) ))

  return(indiceMoranLoc)
}


####################################################################
#' calGearyGlo
#'
#' @details description, a paragraph
#' @param matVoisin xxxx
#' @param vectMoy xxxx
#' @param moyTot xxxx
#' @param vectSurface xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calGearyGlo=function(matVoisin,vectMoy,moyTot,vectSurface)
{
#---------------------------------------------------------------------------------------------------------------------------------#
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#
#entrée:matrice de vosinages de zones ou de points,vecteur des moyennes par zone(ou valeurs par points ), moyenne totale ou moyenne de la zone
#sortie:indice de geary
#---------------------------------------------------------------------------------------------------------------------------------#
  #nombre de zones ou de points
  nbZones=length(vectMoy)

  #transformation de booleéns en entiers
  matNum=matrix(as.numeric(matVoisin),nbZones,nbZones)
  sommeVoisinage=sum(matNum)

  #matrice ou chaque colone est vectMoy
  matMoy=vectMoy%*%t(rep(1,nbZones))

  #numerateur:si pour i et j voisins
  #           Mij=(moyenneI-moyenneJ)^2 * surfaceI*surfaceJ
  #           numerateur=somme(Mij)*(nombre de zones -1)
  #denominateur:si pour i et j voisins
  #           Mi=(moyenneI-moyenneTotale)*somme(surface des voisins  de I)
  #           denominateur=somme(Mi)*(2*nombre de voisinages)

  indiceGeary=(((nbZones-1)/2*sommeVoisinage)*sum(matNum * ((matMoy-t(matMoy))^2)*(vectSurface%*%t(vectSurface)) )
               /sum( ((vectMoy-moyTot)^2) * apply(FUN=sum,MARGIN=1 , matNum * t(vectSurface%*%t(1:nbZones))) ))

  return(indiceGeary)
}
