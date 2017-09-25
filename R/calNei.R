######################################################################
#' calNei
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param data xxxx
#' @param surfVoronoi xxxx
#' @param ptN xxxx
#' @param simplitol xxxx
#' @param remove xxxx
#' @param correct xxxx
#'
#' @return a ?
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' # not run
calNei=function(Z,data,surfVoronoi,ptN,simplitol,remove=TRUE,correct=FALSE)
######################################################################
{
  nbZ=length(Z)

  #list of pts in zones
  listZonePoint=zoneAssign(data,Z)

  #determine zone neighbors
  zoneN=matrix(logical(nbZ^2),nbZ,nbZ)
  vZ=calZoneN(ptN,zoneN,listZonePoint)
  zoneN = vZ$zoneN

  # remove zones with #pts <= n  from zoning Z or surf<minSizeNG
  if(remove)
	{
	res=removeFromZ(Z,zoneN,ptN,listZonePoint,data,simplitol,n=1)
	Z=res$Z
	zoneN=res$zoneN
	listZonePoint=res$listZonePoint
	}

  # correct zone neighbors
  # confusion when voronoi are close - so check again zone distance
  if(correct) voisinZone = correctN(Z,zoneN)
  #
  zoneNModif = zoneN
  diag(zoneNModif) = FALSE
  #
  meanTot=sum(data[[1]]*surfVoronoi)/sum(surfVoronoi)
  listSurf=sapply(Z,gArea)
  # compute (surface/perim^2)
  listPerim=sapply(Z,gLength)
  critSurf=listSurf/(listPerim^2)

  #zone mean values (each pt ponderated by its Voronoi surface)
  meanZone=wMean(2,listZonePoint,surfVoronoi,data)

  return(list(data=data,zoneN=zoneN,zoneNModif=zoneNModif,listZonePoint=listZonePoint,meanTot=meanTot,
              critSurf=critSurf,meanZone=meanZone,listSurf=listSurf,zonePolygone=Z))
 }


#########################################
#' labZone
#'
#' @details description, a paragraph
#' @param listK xxxx
#' @param qProb xxxx
#' @param matVal xxxx
#'
#' @return a ?
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' # not run
labZone=function(listK,qProb,matVal)
#########################################
  {
#create labels
    # input is zoning from calNei, quantile vector and data values
    # output has lab
  lab = rep(1,length(listK$zonePolygone))
  q1= quantile(matVal,na.rm=TRUE,prob=qProb)

  for (i in 1:length(listK$zonePolygone))
  {
    for (j in 1:length(q1))
    {
      if (listK$meanZone[i]>(q1[j]+0.1))
      {
        lab[i]=j+1
      }

    }
  }
listK$lab=lab
listK$qProb=qProb
return(listK)
}

######################################
#' correctN
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param zoneN xxxx
#'
#' @return a ?
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' # not run
correctN = function(Z,zoneN)
######################################
{
# confusion when voronoi are close -> false neighborhood
# so check again zone distance
# computationally costly - obsolete
# correct pt neighborhood instead (done in voronoiPolygons)

if (length(Z)>=2)
  {
    for (j in 1:(length(Z)-1))
    {
      for (k in (j+1):length(Z))
      {
      	a=gDistance(Z[[j]],Z[[k]])
  	if(a >0.001)
      	 zoneN[j,k]=zoneN[k,j]=FALSE
  	 } #end k loop
      }#end j loop
    } # end test more than one zone
return(zoneN)
}
