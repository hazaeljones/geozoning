###########################################################################
#' randKmap:  Generate data for zoning or prepare real data
#'
#' @details Generates a map structure from simulated data or real data. Kriged data are normalized so that x-coordinates are between 0 and 1. y-coordinates are normalized with the same ratio used for x-coordinates. Kriging is either done with inverse distance interpolation, or with a variogram model.  It creates a structure that contains the data and parameters necessary to perform a zoning. The structure is identical wether the data are simulated or not.
#' @param DataObj =NULL: simulated data with seed or = a data frame with real data
#' @param seed numeric, seed used to randomly generate data points
#' @param nPoints numeric, number of raw data points, default 450
#' @param nPointsK numeric, number of kriged data points, default 2000
#' @param nSimuCond number of conditional simulations, reserved for future implementation.
#' @param typeMod type of variogram model (see vgm) "Gau", "Sph", "Exp"
#' @param Vpsill numeric, default 5
#' @param Vrange numeric, default 0.2
#' @param Vmean numeric, default 8
#' @param Vnugget numeric, default 0
#' @param boundary list contains x and y
#' @param manualBoundary logical, default FALSE
#' @param krig numeric
#' @param disp numeric
#' @param FULL logical, if TRUE the returned list is complete
#' @return a list
#' \describe{
#' \item{rawData}{simulated or real raw data within the boundary}
#' \item{step}{grid step}
#' \item{krigData}{kriged data as a SpatialPointsDataFrame}
#' \item{krigGrid}{kriged data in form of a grid-useful for image plots.}
#' \item{krigN}{list of neighbours of each kriged data point}
#' \item{krigSurfVoronoi}{list of areas of Voronoi polygons in the tesselation of kriged data}
#' \item{modelGen}{random fields model}
#' \item{VGMmodel}{vgm model}
#' \item{boundary}{(x,y) list of boundary points}
#' \item{ratio}{ratio used to normalize x data}
#' }
#' @export
#' @importFrom gstat krige
#' @importFrom sp coordinates
#' @importMethodsFrom sp coordinates
#' @examples
#' map = randKmap(NULL,nPointsK=200,Vmean=15,krig=2)
#' mean(map$krigGrid) # mean of generated kriged data
#' plotMap(map)
randKmap=function(DataObj,seed=NULL,nPoints=450,nPointsK=2000, nSimuCond=0,typeMod="Exp",Vpsill=5,Vrange=0.2,Vmean=8,Vnugget=0.2, boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),manualBoundary=FALSE,krig=2,disp=0,FULL=FALSE)
###########################################################################
{

    #simulation seed for random fields
    set.seed(seed)
    Vanis=0 # no anisotropy
    # genData reads real data in DataObj data frame and returns data and a vgm model
    # or simulates them if DataObj=NULL according to vgm model
    resGene=genData(DataObj,seed,nPoints,typeMod,Vpsill,Vrange,Vmean,Vnugget,Vanis,boundary,manualBoundary)
    if(is.null(resGene)) return(NULL)
    rawDataRaw=resGene$tabData
    boundary=resGene$boundary
    ratio=resGene$ratio

    x=as.numeric(rawDataRaw[,1])
    xsize=max(x)-min(x)
    y=as.numeric(rawDataRaw[,2])
    ysize=max(y)-min(y)
    z=as.numeric(rawDataRaw[,3])
    nPoints=length(z)

    #keep only pts within the boundary
    ptsInsRaw=point.in.polygon(x,y,boundary$x,boundary$y)#returns 0 if not within, 1 if strictly inside, 2 is inside an edge, 3 if vertex

    if(sum(ptsInsRaw!=0)==0)
    {
    print("no data points within boundary-please check boundary argument")
    return(NULL)
    }

    rawDataNa=cbind(rawDataRaw,ptsIns=ptsInsRaw) #all pts with index column
    rawDataNa$z[ptsInsRaw==0]=NA # replace pts outside boundary with NA
    rawData=rawDataNa[ptsInsRaw!=0,1:3] #exclude pts outside boundary
    z=as.numeric(rawData[,3])
    #make spatial object
    sp::coordinates(rawDataNa)=~x+y
    sp::coordinates(rawData)=~x+y
    ####define empty grid for kriging#####################
    # - compute step
    step=calStep(nPointsK,xsize,ysize)
    # - generate grid of kriged pt locations
    gridK=genEmptyGrid(step,xsize,ysize)
    tabK=gridK$tabEmpty
    xempty=as.numeric(tabK$x) #to save space in names
    yempty=as.numeric(tabK$y)
    vecZ=tabK$z
    nx=gridK$nx
    ny=gridK$ny
    xx=gridK$xx
    yy=gridK$yy
    # find future kriged pt locations within boundary
    ptsIns = point.in.polygon(xempty,yempty,boundary$x,boundary$y)
    maskIns=(ptsIns==1) #strictly interior pts

    # generate kriged pts at these locations
    # sometimes matrix is singular with krig=1
    if (krig == 2) #  inverse distance
    	 krigZ=krige(z~1,rawData,newdata=tabK[maskIns,])
    else if (krig == 1) # vgm model
    	 krigZ=krige(z~1,rawData,newdata=tabK[maskIns,],model=resGene$VGMmodel)
    else return(NULL)
    #transform into grid matrix
    vecZ[maskIns]=as.numeric(krigZ$var1.pred)
    krigData=data.frame(x=xempty,y=yempty,var1.pred=vecZ)
    sp::coordinates(krigData)=~x+y
    #NAs outside boundary
    #avoid side effects by using krigDataNa
    krigDataNa=cbind(data.frame(krigData),ptsIns)
    sp::coordinates(krigDataNa)=~x+y
    # effective number of kriged pt locations nx*ny
    krigGrid=matrix(vecZ,nx,ny)
    colnames(krigGrid)=round(yy,3)
    rownames(krigGrid)=round(xx,3)
    ##Voronoi on kriged pts
    #prepare matrix
    nbP= nx*ny
    neighBool=matrix(logical(nbP^2),nbP,nbP)

    resVoronoi=voronoiPolygons(krigDataNa,c(0,xsize,0,ysize),neighBool,FULL)
    neighNa=resVoronoi$neighBool
    surfVoronoiNa=resVoronoi$surfVoronoi
    surfVoronoi=surfVoronoiNa[maskIns] #remove pts outside boundary
    #pt neighborhood matrix
    neighBool=neighNa[maskIns,maskIns]
    # compute list of neighbor pts
    krigN = ptNei(neighBool)

    #if required (argument FULL) compute voronoi on raw pts also
    if(FULL)
    {
	    krigVoronoi=resVoronoi$voronoi
	    nbPB= nrow(rawDataNa)
      neighB=matrix(logical(nbPB^2),nbPB,nbPB)
      resVoronoiB=voronoiPolygons(rawDataNa,neighB)
      voronoiB=resVoronoiB$voronoi
      surfVoronoiNaB=resVoronoiB$surfVoronoi
      surfVoronoiB=surfVoronoiNaB[rawDataNa$ptsIns!=0]
      neighNaB=resVoronoiB$neighBool
      neighB=neighNaB[rawDataNa!=0,rawDataNa!=0]
      ptNB = ptNei(neighB)
    }

    # conditional simulation - to be added
    meanCondTab=NULL
    meanCondTabNa=NULL
    matMeanCond=NULL

    # reduced object size except if FULL argument was given
    if(FULL)
      return(list(rawData=rawData,ratio=ratio,step=step,xsize=xsize,ysize=ysize,krigData=krigData[maskIns,],
                  krigGrid=krigGrid,krigN=krigN,krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,
                  VGMmodel=resGene$VGMmodel,boundary=boundary,krigVoronoi=krigVoronoi,matMeanCond=matMeanCond,
                  boundary=boundary,meanCondTab=meanCondTab,meanCondTabNa=meanCondTabNa,surfVoronoiB=surfVoronoiB,
                  voronoiB=voronoiB,ptNB=ptNB))
    else
	     return(list(rawData=rawData,ratio=ratio,step=step,xsize=xsize,ysize=ysize,krigData=krigData[maskIns,],
	                 krigGrid=krigGrid,krigN=krigN,krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,
	                 VGMmodel=resGene$VGMmodel,boundary=boundary))
}
