#' calMatDistance
#'
#' @details description, a paragraph
#' @param typedist xxxx
#' @param zoneN xxxx
#' @param listZonePoint xxxx
#' @param tabVal xxxx
#' @param surfVoronoi xxxx
#' @param meanZone xxxx
#' @param pErr xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calMatDistance=function(typedist=1,zoneN=NULL,listZonePoint=NULL,tabVal=NULL,surfVoronoi=NULL,meanZone=NULL,pErr)
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
#' @details description, a paragraph
#' @param index xxxx
#' @param listZonePoint xxxx
#' @param tabVal xxxx
#' @param surfaceVoronoi xxxx
#' @param meanZone2 xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
Sigmai2=function(index,listZonePoint,tabVal,surfaceVoronoi,meanZone2)
################################################################
{
 valZoneI=tabVal@data[listZonePoint[[index]],1]	#liste des valeurs d'altitude des points de la zone
 voroZoneI=surfaceVoronoi[listZonePoint[[index]]]    #vecteur des surfaces des points de la zone indice1
  SI= sum(voroZoneI)
  fxmean= meanZone2[[index]]
  sigmai2=sum(((valZoneI-fxmean)^2)*voroZoneI)/SI

  return(list(sigmai2=sigmai2,SI=SI))
}


################################################################
#' DIJ
#'
#' @details description, a paragraph
#' @param i xxxx
#' @param j xxxx
#' @param sigmai2 xxxx
#' @param meanZone2 xxxx
#' @param pErr xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
DIJ=function(i,j,sigmai2,meanZone2,pErr)
{
  fxmean = meanZone2[[i]]
  fymean = meanZone2[[j]]
  res = sigmai2[i] + sigmai2[j] + (fxmean-fymean)^2
  resCorr = max(sigmai2[i],(fxmean*pErr/100)^2) + max(sigmai2[j],(fymean*pErr/100)^2) + (fxmean-fymean)^2

  return(list(d=res,dCorr=resCorr))
}

################################################################
#---------------------------------------------------------------------------------------------------------------------------------#
#fonction qui renvoie la valeur de la "distance" intrazone ou interzone(se calcule pareil),avec des valeurs krigees
#entree:indices des deux zones a mesurer, liste des points appartenant a chaque zone,dataframe des valeurs et coordonnees de chaque point (krig?s)
#sortie:valeur entre les deux zones sélectionnées(integer)
#---------------------------------------------------------------------------------------------------------------------------------#
#' calDistZone1
#'
#' @details description, a paragraph
#' @param indice1 xxxx
#' @param indice2 xxxx
#' @param listZonePoint xxxx
#' @param tabVal xxxx
#' @param surfaceVoronoi xxxx
#' @param meanZone2 xxxx
#' @param pErr xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calDistZone1=function(indice1,indice2,listZonePoint=NULL,tabVal=NULL,surfaceVoronoi,meanZone2,pErr=0.9)
{
  # bch
  # tabVal contient les valeurs krigees ou brutes selon appel

  valZoneI=tabVal@data[listZonePoint[[indice1]],1]	#liste des valeurs d'altitude des points de la zone indice1
  valZoneJ=tabVal@data[listZonePoint[[indice2]],1]	#liste des valeurs d'altitude des points de la zone indice2

  polyZoneI=surfaceVoronoi[listZonePoint[[indice1]]]    #vecteur des surfaces des points de la zone indice1
  polyZoneJ=surfaceVoronoi[listZonePoint[[indice2]]]		#vecteur des surfaces des points de la zone indice2

  SI= sum(polyZoneI)
  SJ= sum(polyZoneJ)

  fxmean= meanZone2[[indice1]]
  fymean=meanZone2[[indice2]]

sigmai2=sum(((valZoneI-fxmean)^2)*polyZoneI)/SI
sigmaj2=sum(((valZoneJ-fymean)^2)*polyZoneJ)/SJ



  res= sigmai2 + sigmaj2 + (fxmean-fymean)^2
  resCorr= max(sigmai2,(fxmean*pErr/100)^2) + max(sigmaj2,(fymean*pErr/100)^2) + (fxmean-fymean)^2

  return(list(matDistance=res,matDistanceCorr=resCorr))
}
#####################################################################################################

#####################################################################################################
#####################################################################################################
#renvoie un vecteur contenant la valeur du variogramme max- le variogramme pour une distance donnée
#entrée:distance entre les points choisis(integer),plateau du variogramme(integer),modèle du variogramme sur la carte"variogramModel" "data.frame"
#' calculMatriceVario
#'
#' @details description, a paragraph
#' @param distij xxxx
#' @param maxVario xxxx
#' @param fitVarioAlea xxxx
#'
#' @return a ?
#' @importFrom gstat variogramLine
#'
#' @export
#'
#' @examples
#' # not run
calculMatriceVario=function(distij,maxVario,fitVarioAlea)
{
  #rajouter une matrice différences x
  #                     différences y
  return(maxVario - variogramLine(fitVarioAlea,maxdist=10,n=1,dist_vector=distij)[,2])

}
#####################################################################################################

#####################################################################################################
# IS 19/05/2017: comments on calDistZone2() calling a .C function!
#####################################################################################################
#fonction qui renvoie la valeur de la "distance" intrazone ou interzone(se calcule pareil),en multipliant
# par une fonction du variogramme
#entree:indices des deux zones a mesurer(integer), liste des points appartennant a chaque zone(list(integer)),
#       dataframe des valeurs et coordonnees de chaque point(data.frame),vecteur des surfaces de voronoi des points(numeric)
#sortie:distance entre les deux zones(integer)
#####################################################################################################
# calDistZone2
#
# details description, a paragraph
# param indice1 xxxx
# param indice2 xxxx
# param listZonePoint xxxx
# param tabVal xxxx
# param surfaceVoronoi xxxx
# param fitVarioAlea xxxx
#
# return a ?
# importFrom sp coordinates
# importMethodsFrom sp coordinates
#
# export
#
# examples
# # not run
# calDistZone2=function(indice1,indice2,listZonePoint,tabVal,surfaceVoronoi,fitVarioAlea)
# {
#
#   tabCoord= as.data.frame(coordinates(tabVal))
#   d=0
#   valeursZoneI=data.frame(tabVal[listZonePoint[[indice1]],])[,3]    #vecteur des valeurs d'altitude des points de la zone indice1
#   valeursZoneJ=data.frame(tabVal[listZonePoint[[indice2]],])[,3]		#vecteur des valeurs d'altitude des points de la zone indice2
#
#   polyZoneI=surfaceVoronoi[listZonePoint[[indice1]]]  	#vecteur des surfaces des points de la zone indice1
#   polyZoneJ=surfaceVoronoi[listZonePoint[[indice2]]]		#vecteur des surfaces des points de la zone indice2
#
#   taille1=length(valeursZoneI)
#   taille2=length(valeursZoneJ)
#
#   matValeursZoneI=matrix(rep(valeursZoneI,taille2),taille1,taille2)				#creation d'une matrice ayant sur chaque colone les valeurs des points de la zone i
#   matValeursZoneJ=t(matrix(rep(valeursZoneJ,taille1),taille2,taille1))    #creation d'une matrice ayant sur chaque ligne les valeurs des points de la zone j
#
#
#   #denom=0
#   maxVario=sum(fitVarioAlea[,2] )     #plateau du variogramme(asymptote)
#   matDist=matrix(0,taille1,taille2)
#
#   ##########################
#   # appel à une fonction C #
#   ##########################
#   vectDist=as.numeric(matDist)
#
#   resC<-.C("calDistance",as.double(tabCoord$x[listZonePoint[[indice1]]]),
#                         as.double(tabCoord$y[listZonePoint[[indice1]]]),
#                         as.double(tabCoord$x[listZonePoint[[indice2]]]),
#                         as.double(tabCoord$y[listZonePoint[[indice2]]]),
#                         as.double(vectDist),
#                         length(listZonePoint[[indice1]]),
#                         length(listZonePoint[[indice2]]))
#   vectDist<-resC[[5]]
#   matDist=matrix(vectDist,taille1,taille2)
#
#   #on génère un vecteur de longueur taille1*taille2 contenant les valeurs de variogramme pour toutes les paires
#   #de points des zones considérées,puis on transforme en matrice
#   matCoeffVario=matrix(calculMatriceVario(distij=as.vector(matDist),maxVario,fitVarioAlea),taille1,taille2)
#
#   #matrice contenant la multiplication des surfaces voronoi i par voronoi j à la position(i,j),par produit matriciel
#   matPolyZone=polyZoneI%*%t(polyZoneJ)
#
#   #==Somme sur i (Somme sur j (f(Xkrig i)-f(Xkrig j))*surface voronoi i j) / somme des multiplication des surfaces voronoi
#   d=sum(((matValeursZoneI-matValeursZoneJ)^2)*matPolyZone*matCoeffVario)/(sum(matPolyZone*matCoeffVario))
#
#   return(d)
# }
#####################################################################################################

#####################################################################################################
#####################################################################################################
#renvoie la valeur de la "distance" intrazone ou interzone(se calcule pareil),simple
#entree:indices des deux zones a mesurer(integer), liste des points appartenant a chaque zone(list(integer)),dataframe des valeurs et coordonn?es de chaque point(data.frame),vecteur des surfaces de voronoi des points(numeric)
#sortie:distance entre les deux zones(integer)
#####################################################################################################

# fonction commentee car identique à la fonction calDistZone1
# calDistZone3=function(indice1,indice2,listZonePoint,tabVal,surfaceVoronoi)
# {
#
#
#   d=0
#
#
#   valeursZoneI=data.frame(tabVal[listZonePoint[[indice1]],])[,3]		#vecteur des valeurs d'altitude des points de la zone indice1
#   valeursZoneJ=data.frame(tabVal[listZonePoint[[indice2]],])[,3]		#vecteur des valeurs d'altitude des points de la zone indice2
#
#   polyZoneI=surfaceVoronoi[listZonePoint[[indice1]]]  	#vecteur des surfaces des points de la zone indice1
#   polyZoneJ=surfaceVoronoi[listZonePoint[[indice2]]]		#vecteur des surfaces des points de la zone indice2
#
#
#
#   taille1=length(valeursZoneI)
#   taille2=length(valeursZoneJ)
#   matValeursZoneI=matrix(rep(valeursZoneI,taille2),taille1,taille2)				#creation d'un matrice ayant sur chaque colone les valeurs des points de la zone i
#   matValeursZoneJ=t(matrix(rep(valeursZoneJ,taille1),taille2,taille1))    #creation d'une matrice ayant sur chaque ligne les valeurs des points de la zone j
#
#   matPolyZone=polyZoneI%*%t(polyZoneJ)
#
#   d=sum(((matValeursZoneI-matValeursZoneJ)^2)*matPolyZone)/sum(matPolyZone)						#==Somme sur i (Somme sur j (f(Xkrig i)-f(Xkrig j))*surface voronoi i j) / somme des multiplication des surfaces voronoi
#
#   return(d)
# }
#####################################################################################################

#####################################################################################################
#####################################################################################################
#fonction calculant la distance de Moran entre deux zones
#entrée:indices des zones concernées par le calcul de distance(integer),liste des points appartennant à chaque zone(list(integer)),dataframe des points et des leurs position(data.frame),moyenne totale de lacarte(integer)
#sortie:valeur de la distance entre les deux zones concernées(integer)
#' calDistZoneMoranB
#'
#' @details description, a paragraph
#' @param indice1 xxxx
#' @param indice2 xxxx
#' @param listZonePoint xxxx
#' @param tabVal xxxx
#' @param surfaceVoronoi xxxx
#' @param moyTot xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calDistZoneMoranB=function(indice1,indice2,listZonePoint,tabVal,surfaceVoronoi,moyTot)
{
  polyZoneI=surfaceVoronoi[listZonePoint[[indice1]]]    #vecteur des surfaces des points de la zone indice1
  polyZoneJ=surfaceVoronoi[listZonePoint[[indice2]]]  	#vecteur des surfaces des points de la zone indice2

  valeursZoneI=(data.frame(tabVal[listZonePoint[[indice1]],])[,3] - moyTot)*polyZoneI    #vecteur des valeurs d'altitude des points de la zone indice1-moyenne
  valeursZoneJ=(data.frame(tabVal[listZonePoint[[indice2]],])[,3] - moyTot)*polyZoneJ  	#vecteur des valeurs d'altitude des points de la zone indice2-moyenne

  d=sum(valeursZoneI%*%t(valeursZoneJ))/(sum(polyZoneI)*sum(polyZoneJ))


  return(d)
}


#---------------------------------------------------------------------------------------------------------------------------------#
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#

#fonction normalisant chaque terme Mij de la matrice par la racine carrée de Mii*Mjj
#entr?e:matrice des distances classique,
#sortie:matrice des distance normalis?e

#---------------------------------------------------------------------------------------------------------------------------------#

#' distanceNormalisationSqrt
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
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

#' distanceNormalisationSum
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
distanceNormalisationSum=function(matDistance)
{
  #matrice termeNormal[i,j]=1/sum(dii^2+djj^2)
  nbPoly=nrow(matDistance)
  termeNormal=matrix(rep(diag(matDistance),nbPoly),nbPoly,nbPoly)
  termeNormal=1/(termeNormal+t(termeNormal))
  matDistanceNorm=2*matDistance*termeNormal

  return(matDistanceNorm)
}

#---------------------------------------------------------------------------------------------------------------------------------#
#################################################################################################################################
#---------------------------------------------------------------------------------------------------------------------------------#

#fonction qui ajoute à la matrice des distances normalisée:les termes diagonaux non normalisés et le rapport de correction effectué sur les termes diagonaux
#entr?e:matrice des distance normalisée ,classique et corrigée(matrix)
#sortie:dataframe des distances normalis?e avec une colone pour la diagonale classique et une pour la correction(data.frame)

#---------------------------------------------------------------------------------------------------------------------------------#
#' distanceExtension
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param matDistanceNorm xxxx
#' @param matDistanceCorr xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
distanceExtension=function(matDistance,matDistanceNorm,matDistanceCorr)
{
  #vecteur,boolCorrection[i] indique si la valeur du dii  à été corrigée ou non
  boolCorrection=as.numeric(diag(matDistance)!=diag(matDistanceCorr))
  rapportCorrection=boolCorrection*(diag(matDistance)/(diag(matDistanceCorr)))

  #dataframe de la matrice corrigée normalisée,avec valeurs de la diagonale corrigée et le rapport de correction
  tabExt=data.frame(matDistanceNorm,diag(matDistanceCorr),rapportCorrection)

  names(tabExt)[length(diag(matDistance))+1]="Dii"
  names(tabExt)[length(diag(matDistance))+2]="Corr"
  return( tabExt)

}


## essai en pénalisant plus les gros écarts

#' calDistZone12
#'
#' @details description, a paragraph
#' @param indice1 xxxx
#' @param indice2 xxxx
#' @param listZonePoint xxxx
#' @param tabVal xxxx
#' @param surfaceVoronoi xxxx
#' @param meanZone2 xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calDistZone12=function(indice1,indice2,listZonePoint=NULL,tabVal=NULL,surfaceVoronoi,meanZone2)
{

  d=0
  valeursZoneI=data.frame(tabVal[listZonePoint[[indice1]],])[,3]  	#liste des valeurs d'altitude des points de la zone indice1
  valeursZoneJ=data.frame(tabVal[listZonePoint[[indice2]],])[,3]		#liste des valeurs d'altitude des points de la zone indice2

  polyZoneI=surfaceVoronoi[listZonePoint[[indice1]]]    #vecteur des surfaces des points de la zone indice1
  polyZoneJ=surfaceVoronoi[listZonePoint[[indice2]]]		#vecteur des surfaces des points de la zone indice2

  matPolyZone=polyZoneI%*%t(polyZoneJ)

  taille1=length(valeursZoneI)
  taille2=length(valeursZoneJ)
  matValeursZoneI=matrix(rep(valeursZoneI,taille2),taille1,taille2)				#creation d'un matrice ayant sur chaque colone les valeurs des points de la zone i
  matValeursZoneJ=t(matrix(rep(valeursZoneJ,taille1),taille2,taille1))    #creation d'une matrice ayant sur chaque ligne les valeurs des points de la zone j


  d=sum(matPolyZone*(matValeursZoneI-matValeursZoneJ)^4)/sum(matPolyZone)					#==Somme sur i (Somme sur j (f(Xkrig i)-f(Xkrig j))?) / nbptkrig1*nbptkrig2

  return(d)
}
#---------------------------------------------------------------------------------------------------------------------------------#
