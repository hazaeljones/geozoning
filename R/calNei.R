######################################################################
#' calNei calculate zone neighborhood and assign data points to zones
#'
#' @details calNei first removes from zoning Z all zones with less than a minimum number of points. Then it calculates zone neighborhood, assigns each data point to a zone, computes zone mean values and areas. It does not assign zone labels (this is done by labZone function for the initial zoning, and by trLabZone function to transfer labels to corrected zonings).
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param spdata  SpatialPointsDataFrame containing the data pts and values
#' @param surfVoronoi Surfaces of the Voronoi polygons corresponding to data pts
#' @param ptN indices of data pts neighbours
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param remove if TRUE remove zones with less than nmin data points
#' @param correct if TRUE correct zone neighborhood
#' @param nmin number of points below wich a zone is removed from the zoning (default is 2)
#' @return a list with components
#' \describe{
#' \item{zoneN}{matrix of zone neigbors}
#' \item{zoneNModif}{modified matrix with FALSE on the diagonal}
#' \item{listZonePoint}{ indices of pts within each zone}
#' \item{meanTot}{zoning mean data value}
#' \item{meanZone}{vector of zone data mean values}
#' \item{listSurf}{vector of zone areas}
#' \item{critSurf}{vector of filiform zone characteristics}
#' \item{zonePolygone}{list of zones, each zone is a SpatialPolygons}
#' }
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' data(mapTest)
#' ptN=mapTest$krigN
#' spdata=mapTest$krigData
#' surfVoronoi=mapTest$surfVoronoi
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' K=calNei(Z,spdata,surfVoronoi,ptN)
#' names(K)
#' plotZ(K$zonePolygone)
#' K=calNei(Z,spdata,surfVoronoi,ptN,nmin=20) #keep only zones with a minimum of 20 data points
#' plotZ(K$zonePolygone)

calNei=function(Z,spdata,surfVoronoi,ptN,simplitol=1e-3,remove=TRUE,correct=FALSE,nmin=2)
#################################################################################
{
  nbZ=length(Z)

  #list of pts in zones
  listZonePoint=zoneAssign(spdata,Z)

  #determine zone neighbors
  zoneN=matrix(logical(nbZ^2),nbZ,nbZ)
  vZ=calZoneN(ptN,zoneN,listZonePoint)
  zoneN = vZ$zoneN

  # remove zones with #pts <= n  from zoning Z or surf<minSizeNG
  if(remove)
	{
	res=removeFromZ(Z,zoneN,ptN,listZonePoint,spdata,simplitol,n=max(0,nmin))
	Z=res$Z
	zoneN=res$zoneN
	listZonePoint=res$listZonePoint
	}

  # correct zone neighbors
  # confusion when voronoi are close - so check again zone distance
  if(correct) zoneN = correctN(Z,zoneN)
  #
  zoneNModif = zoneN
  diag(zoneNModif) = FALSE
  #
  meanTot=sum(spdata[[1]]*surfVoronoi)/sum(surfVoronoi)
  listSurf=sapply(Z,gArea)
  # compute (surface/perim^2)
  listPerim=sapply(Z,gLength)
  critSurf=listSurf/(listPerim^2)

  #zone mean values (each pt ponderated by its Voronoi surface)
  meanZone=wMean(2,listZonePoint,surfVoronoi,spdata)

  return(list(zoneN=zoneN,zoneNModif=zoneNModif,listZonePoint=listZonePoint,meanTot=meanTot,
              critSurf=critSurf,meanZone=meanZone,listSurf=listSurf,zonePolygone=Z))
 }


#########################################
#' labZone0
#'
#' @details assigns a class label (integer) to a zone depending on the zone mean value
#' and on the quantile values. Default label is 1, corresponding to mean value samller #' or equal to first quantile. For k ordered quantile values, if mean value is greater #' than quantile k plus 10% of the data range, zone label is k.
#' @param K zoning object, as returned by the calNei function
#' @param qProb probability vector used to generate quantile values for Z
#' @param dataF data used to generate labels and zoning
#'
#' @return a zoning object with labelled zones in lab component
#'
#' @keywords internal
#' @importFrom sp point.in.polygon
#'
#' @examples
#' data(mapTest)
#' dataF=mapTest$krigGrid
#' data(resZTest)
#' K=resZTest
#' p = K$qProb
#' geozoning:::labZone0(K,p,dataF)
#' # not run
labZone0=function(K,qProb,dataF)
#########################################
  {
    #create labels
    # input is zoning from calNei, quantile vector and data values
    # output has lab
  lab = rep(1,length(K$zonePolygone))
  # consider range of data values
  rate= max(dataF)-min(dataF)
  q1= quantile(dataF,na.rm=TRUE,prob=qProb)

  for (i in 1:length(K$zonePolygone))
  {
    for (j in 1:length(q1))
    {
      if (K$meanZone[i]>(q1[j]+0.01*rate))
      {
        lab[i]=j+1
      }

    }
  }
K$lab=lab
K$qProb=qProb
return(K)
  }

#################################################################################################
#' labZone
#'
#' @details assigns a class label (integer) to a zone depending on the zone mean value
#' and on the quantile values (as in PA paper). Default label is 1, corresponding to a mean value smaller or equal to first quantile. For p ordered quantile values, if mean value is greater than quantile k and smaller or equal to quantile k+1, zone label is k+1. if mean value is greater than quantile p, zone label is p+1.
#' @param K zoning object, as returned by the calNei function
#' @param qProb probability vector used to generate quantile values for Z
#' @param dataF data used to generate labels and zoning
#'
#' @return a zoning object with labelled zones in lab component
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' data(mapTest)
#' dataF=mapTest$krigGrid
#' data(resZTest)
#' K=resZTest
#' p = K$qProb
#' labZone(K,p,dataF)
labZone=function(K,qProb,dataF)
#########################################
  {
    #create labels
    #correct assignment according to PA paper, NEW VERSION OF labZone
    #input is zoning from calNei, quantile vector and data values
    #output has lab
  # consider range of data values
  qvec= quantile(dataF,na.rm=TRUE,prob=qProb)
  nL=length(qvec)+1
  dmin=min(dataF,na.rm=TRUE)
  dmax=max(dataF,na.rm=TRUE)
  lab=as.numeric(cut(K$meanZone,c(dmin,qvec,dmax),1:nL))
  K$lab=lab
  K$qProb=qProb
  return(K)
  }

######################################
#' correctN
#'
#' @details description, a paragraph
#' @param Z zoning geometry (list of SpatialPolygons)
#' @param zoneN zone neighborhood Logical matrix
#' @param dN maximum distance beyond which 2 zones cannot be considered as neighbors
#'
#' @return a new zone neighborhood Logical matrix
#'
#' @export
#' @importFrom sp point.in.polygon
#'
#' @examples
#' data(resZTest)
#' Z=resZTest$zonePolygone
#' H=correctN(Z,resZTest$zoneN,1e-8)
correctN = function(Z,zoneN,dN=1e-3)
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
  	if(a >dN)
      	 zoneN[j,k]=zoneN[k,j]=FALSE
  	 } #end k loop
      }#end j loop
    } # end test more than one zone
return(zoneN)
}
