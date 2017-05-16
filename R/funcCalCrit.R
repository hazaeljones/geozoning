################################################
#' calCrit1
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calCrit1=function(matDistance,zoneNModif)
################################################
{
#
# returns min(mean(dij^2/(dii^2+dij^2)))

  nbPoly=length(diag(matDistance))
  #on fixe la valeur initiale du critere a un nombre tres grand
  val=Inf

  #zone pour laquelle la plus petite valeur du crit?re est calcul?e
  zoneVal=0

  for (i in 1:nbPoly)
  {
    nbVois=0
    temp=0
    for (j in 1:nbPoly)
    {
      #pour chaque paire de zones,si elles sont voisines
      if(zoneNModif[i,j])
      {
        #on ajoute au résultat déja calculé dij/(dii+djj)
        temp=temp+(matDistance[i,j]/(matDistance[j,j]+matDistance[i,i]))

        #print(temp)
        nbVois=nbVois+1

      }
    }
    #on divise par le nombre de voisins de la zone i	pour avoir la moyenne
    if(nbVois!=0){
      temp=temp/nbVois
    }

    #on garde en mémoire la plus petite valeur
    if(temp<val && temp!=0)
    {
      zoneVal=i
      val=temp
    }
  }
  #print("ok calcul critere 1")

  return(val)
}



##############################################
#' calCrit2
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calCrit2=function(matDistance,zoneNModif)
##############################################
{
#returns min(2*min(dij/(dii+djj)))
# with dii, djj, dij matrices of squared distances
  nbPoly=length(diag(matDistance))

  #initial criterion value
  val=Inf

  zoneVal=0

  #for each zone
  for (i in 1:nbPoly)
  {
    tmpi=Inf
    #for zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #conpute dij/(dii+djj)
        tmpj=(2*matDistance[i,j]/(matDistance[j,j]+matDistance[i,i]))

        #if current value smaller than previous one, store it
        if(tmpj<tmpi)
        {
               tmpi=tmpj
        }
      }
    }
    #
    if(tmpi<val)
    {
      #store into val
      zoneVal=i
      val=tmpi
    }
  }

  return(val)
}


#################################################
#' calCrit2bis
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calCrit2bis=function(matDistance,zoneNModif)
#################################################

{
##Juste pour voir ce que cela donne si l'on souhaite plus discriminer le manque d'homogénéité intra en divisant
##par la somme des carrés des indices intra
  #returns min(min(dij/(dii^2+dij^2)))

  nbPoly=length(diag(matDistance))

  #on fixe la valeur initiale du critere a un nombre tres grand
  val=Inf

  #zone pour laquelle la plus petite valeur du crit?re est calcul?e
  zoneVal=0

  #pour chaque zone i
  for (i in 1:nbPoly)
  {
    tmpi=Inf
    #pour chaque zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #on calcule la "distance" dij/(dii+djj)
        tmpj=(2*matDistance[i,j]/(matDistance[j,j]^2+matDistance[i,i]^2))

        #si la valeur calculée pour ij est plus petite que la plus petite valeur courante pour i
        if(tmpj<tmpi)
        {
          #on stocke la valeur calculée
          tmpi=tmpj
        }
      }
    }
    #si la plus petite valeur calculée pour ce i est plus petite que la plus petite valeur courante
    if(tmpi<val)
    {
      #on la stocke dans val
      zoneVal=i
      val=tmpi
    }
  }


  return(val)
}

#############################################
#' calCrit3
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calCrit3=function(matDistance,zoneNModif)
#############################################


{
#variant of criterion 1
# mais avec une normalisation par racine de multipl. et non une somme
#entree:la matrice des distances entre zones(matrix), ainsi que la matrice de voisinages modifi?e(une zone n'est pas sa propre voisine)(matrix)
#sortie:critère(numeric)

  #returns min(mean(dij^2/sqrt(dii^2*dij^2)))

  nbPoly=length(diag(matDistance))
  #on fixe la valeur initiale du critere a un nombre tres grand
  val=Inf

  #zone pour laquelle la plus petite valeur du critere est calculee
  zoneVal=0


  for (i in 1:nbPoly)
  {
    nbVois=0
    tmp=0
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #pour chaque paire de zones, si les zones sont voisines
        #on ajoute au resultat deja calculé dij/sqrt(dii*djj)
        if (matDistance[j,j]!='NaN' && matDistance[i,j]!= 'NaN' && matDistance[i,i]!= 'NaN'){
          tmp=tmp+(matDistance[i,j]/sqrt(matDistance[j,j]*matDistance[i,i]))
        }
        nbVois=nbVois+1
      }
    }
    #on divise par le nombre de voisins pour obtenir une moyenne
    if (nbVois !=0){
      tmp=tmp/nbVois
    }
    #on garde la plus petite valeur calculée
    if(tmp<val && tmp!=0)
    {
      zoneVal=i
      val=tmp
    }
  }

  return(val)
}

###############################################
#' calCrit4
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calCrit4=function(matDistance,zoneNModif)
###############################################
{
#renvoie critere 4 variante du critere 2 mais avec une normalisation par racine de multipl. et non une somme
  #returns min(min(dij^2/sqrt(dii^2*djj^2)))

nbPoly=length(diag(matDistance))
  #on fixe la valeur initiale du critere a un nombre tres grand
  val=Inf

  #zone pour laquelle la plus petite valeur du critere est calculee
  zoneVal=0

  #pour chaque zone i
  for (i in 1:nbPoly)
  {
    tmpi=Inf

    #pour chaque zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #si les zones sont voisines
        #on calcule dii/sqrt(dii*dij)
        tmpj=(matDistance[i,j]/sqrt(matDistance[j,j]*matDistance[i,i]))
        if(tmpj<tmpi)
        {
          #on garde la plus petite valeur calculée pour cette zone i
          tmpi=tmpj
        }
      }
    }
    if(tmpi<val)
    {
      #on garde la plus petite valeur calculée parmi toutes les zones
      zoneVal=i
      val=tmpi
    }
  }

  return(val)
}
###############################################
#' calCrit5
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#' @importFrom stats median
#'
#' @export
#'
#' @examples
#' # not run
calCrit5=function(matDistance,zoneNModif)
###############################################
{

## variante : utilisation de la mediane pour remplacer la moyenne ou le minimum. Normalisation geométrique.

  nbPoly=length(diag(matDistance))
  #on fixe la valeur initiale du critere a un nombre tres grand
  val=Inf

  #zone pour laquelle la plus petite vvaleur du crit?re est calcul?e
  zoneVal=0

  mat= as.data.frame(matrix(0,nrow=nbPoly,ncol=nbPoly))
  v=list()
  #pour chaque zone i
  for (i in 1:nbPoly)
  {
    v[[i]] = numeric()
    #pour chaque zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #si les zones sont voisines
        #on calcule dii/sqrt(dii*dij)
        v[[i]]=append(v[[i]],(matDistance[i,j]/sqrt(matDistance[j,j]*matDistance[i,i])))
      }
    }
  }
  a = numeric()
  for (i in 1:nbPoly )
  {
    a = append(a,median(v[[i]]))
  }
  #print("ok calcul critere 4")
  #on retourne: min(min(dij^2/sqrt(dii^2*djj^2)))
  return(min(a))
}

#################################################
#' calCritMinMean
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calCritMinMean=function(matDistance,zoneNModif)
#################################################
{
# variant of criterion 4 (with mean instead of min)
# min(mean(dij^2/sqrt(dii^2*djj^2)))

  nbPoly=length(diag(matDistance))
  val=Inf

  # for each zone i
  for (i in 1:nbPoly)
  {
    tmpi=0
    nb=0
    #pour chaque zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        # for neighboring zones
        # compute dii/sqrt(dii*dij)
        tmpj=(matDistance[i,j]/sqrt(matDistance[j,j]*matDistance[i,i]))
        tmpi=tmpi+tmpj #sum over j neighboring zones
	nb=nb+1
       }
    }
    if (nb>0) tmpi=tmpi/nb #mean over j neighboring zones
    if(tmpi<val)
    {
      #keep minimum value for i
      val=tmpi
    }
  }

  return(val)
}

##############################################
#' calCrit7
#'
#' @details description, a paragraph
#' @param matDistance xxxx
#' @param zoneNModif xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
calCrit7=function(matDistance,zoneNModif)
##############################################
{
#returns mean(2*mean(dij/(dii+djj)))
# with dii, djj, dij matrices of squared distances

  nbPoly=length(diag(matDistance))
  val=0

  #for each zone
  for (i in 1:nbPoly)
  {
	 tmpi=0
    	 nb=0
       for (j in 1:nbPoly) #for its neighbors
       	   {
       	   if(zoneNModif[i,j])
       	   {
           #compute dij/(dii+djj)
           tmpj=(matDistance[i,j]/(matDistance[j,j]+matDistance[i,i]))
	   tmpi=tmpi+tmpj #sum over j neighboring zones
	   nb=nb+1
    	   }
	   }
    	if (nb>0) tmpi=tmpi/nb
	val=val+tmpi
  }

 val=2*val/nbPoly
 return(val)
}
