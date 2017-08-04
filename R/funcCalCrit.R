################################################
#' calCrit1
#'
#' @details computes a quality criterion equal to min(mean(dij^2/(dii^2+dij^2)))
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCrit1(resD$matDistanceCorr,K$zoneNModif)
#' # not run
calCrit1=function(matDistance,zoneNModif)
################################################
{
#
# returns min(mean(dij^2/(dii^2+dij^2)))

  nbPoly=length(diag(matDistance))
  #set the initial value
  val=Inf

  #index of zone corresponding to minimum criterion value
  zoneVal=0

  for (i in 1:nbPoly)
  {
    nbNb=0
    tmp=0
    for (j in 1:nbPoly)
    {
      #for each pair of neighbor zones
      if(zoneNModif[i,j])
      {
        #add dij/(dii+djj)
        tmp=tmp+(matDistance[i,j]/(matDistance[j,j]+matDistance[i,i]))
        nbNb=nbNb+1

      }
    }
    #divide by number of neighbors to obtain mean
    if(nbNb!=0){
      tmp=tmp/nbNb
    }

    #keep smallest value
    if(tmp<val && tmp!=0)
    {
      zoneVal=i
      val=tmp
    }
  }

  return(val)
}



##############################################
#' calCrit2
#'
#' @details computes a quality criterion equal to min(2*min(dij/(dii+djj)))
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value equal to min(mean(dij^2/(dii^2+dij^2)))
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCrit2(resD$matDistanceCorr,K$zoneNModif)
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
        #compute dij/(dii+djj)
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
#' @details computes a quality criterion equal to  min(min(dij/(dii^2+dij^2)))
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCrit2(resD$matDistanceCorr,K$zoneNModif)
#' # not run

calCrit2bis=function(matDistance,zoneNModif)
#################################################

{
##Juste pour voir ce que cela donne si l'on souhaite plus discriminer le manque d'homogénéité intra en divisant
##par la somme des carrés des indices intra
  #returns min(min(dij/(dii^2+dij^2)))

  nbPoly=length(diag(matDistance))

  #set initial criterion value
  val=Inf

  #zone index for minimum value
  zoneVal=0

  #for each zone i
  for (i in 1:nbPoly)
  {
    tmpi=Inf
    #for each zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #compute dij/(dii+djj)
        tmpj=(2*matDistance[i,j]/(matDistance[j,j]^2+matDistance[i,i]^2))

        # keep smaller current value
        if(tmpj<tmpi)
        {
          tmpi=tmpj
        }
      }
    }
    #
    if(tmpi<val)
    {
      #store minimum in val
      zoneVal=i
      val=tmpi
    }
  }


  return(val)
}

#############################################
#' calCrit3
#' @details computes a quality criterion equal to min(mean(dij^2/sqrt(dii^2*dij^2))) 
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCrit3(resD$matDistanceCorr,K$zoneNModif)
#' # not run
calCrit3=function(matDistance,zoneNModif)
#############################################
{
# variant of criterion 1
# standardized with square root of product
# returns min(mean(dij^2/sqrt(dii^2*dij^2)))

  nbPoly=length(diag(matDistance))
  #set the initial value
  val=Inf

  #index of zone corresponding to minimum criterion value
  zoneVal=0

  for (i in 1:nbPoly)
  {
    nbNb=0
    tmp=0
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #for each pair of neighbor zones
        #add dij/sqrt(dii*djj)
        tmp=tmp+(matDistance[i,j]/sqrt(matDistance[j,j]*matDistance[i,i]))
        nbNb=nbNb+1
      }
    }
    #divide by number of neighbors to obtain mean
    if (nbNb !=0){
      tmp=tmp/nbNb
    }
    #keep smallest value
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
#' @details computes a quality criterion equal to min(min(dij^2/sqrt(dii^2*djj^2)))
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCrit4(resD$matDistanceCorr,K$zoneNModif)
#' # not run
calCrit4=function(matDistance,zoneNModif)
###############################################
{
#variant of criterion 2 standardized with square root of product of squares
#returns min(min(dij^2/sqrt(dii^2*djj^2)))

nbPoly=length(diag(matDistance))
  #set the initial value
  val=Inf

  #index of zone corresponding to minimum criterion value
  zoneVal=0

  #for each zone i
  for (i in 1:nbPoly)
  {
    tmpi=Inf

    #for each zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #for each pair of neighbor zones
        #compute dii/sqrt(dii*dij)
        tmpj=(matDistance[i,j]/sqrt(matDistance[j,j]*matDistance[i,i]))
        if(tmpj<tmpi)
        {
          #keep smallest value for zone i
          tmpi=tmpj
        }
      }
    }
    if(tmpi<val)
    {
      #keep smallest value for all zones
      zoneVal=i
      val=tmpi
    }
  }

  return(val)
}
###############################################
#' calCrit5
#'
#' @details computes a quality criterion equal to min(median(dij/sqrt(dii*dij))) 
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCrit5(resD$matDistanceCorr,K$zoneNModif)
#' # not run
calCrit5=function(matDistance,zoneNModif)
###############################################
{

## variant: use median instead of min or mean. Geometric standardization.

nbPoly=length(diag(matDistance))
  #set the initial value
  val=Inf

  #index of zone corresponding to minimum criterion value
  zoneVal=0

  mat= as.data.frame(matrix(0,nrow=nbPoly,ncol=nbPoly))
  v=list()
  #for each zone i
  for (i in 1:nbPoly)
  {
    v[[i]] = numeric()
    #for each zone j
    for (j in 1:nbPoly)
    {
      if(zoneNModif[i,j])
      {
        #for each pair of neighbor zones
        #compute dii/sqrt(dii*dij)
        v[[i]]=append(v[[i]],(matDistance[i,j]/sqrt(matDistance[j,j]*matDistance[i,i])))
      }
    }
  }
  a = numeric()
  for (i in 1:nbPoly )
  {
    a = append(a,median(v[[i]]))
  }
  return(min(a))
}

#################################################
#' calCritMinMean
#'
#' @details computes a quality criterion equal to min(mean(dij^2/sqrt(dii^2*djj^2)))
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCritMinMean(resD$matDistanceCorr,K$zoneNModif)
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
#' @details computes a quality criterion equal to mean(2*mean(dij/(dii+djj)))
#' @param matDistance zone distance matrix resulting from a call to calDistance
#' @param zoneNModif matrix of zone neigbors with FALSE on the diagonal}
#'
#' @return a numerical value 
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
#' calCrit7(resD$matDistanceCorr,K$zoneNModif)
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
