##########################################################
#' randKmapGrid
#'
#' @details Prepare real data for zoning, data are already on a regular grid, hence no kriging is done.
#' @param DataObj =NULL: simulated data with seed or = a data frame with real data
#' @param seed numeric, seed
#' @param nPoints numeric, number of points, default 450
#' @param nPointsK numeric, default 2000
#' @param nSimuCond numeric
#' @param typeMod character, model type
#' @param Vpsill numeric, default 5
#' @param Vrange numeric, default 0.2
#' @param Vmean numeric, default 8
#' @param Vnugget numeric, default 0
#' @param boundary list contains x and y
#' @param manualBoundary logical, default FALSE
#' @param krig numeric
#' @param disp numeric
#' @param FULL logical, if TRUE the returned list is complete
#'
#' @importFrom sp coordinates
#'
#' @return a list
#' \describe{
#' \item{rawData}{rawData}
#' \item{step}{step}
#' \item{krigData}{krigData}
#' \item{krigGrid}{krigGrid}
#' \item{krigN}{krigN}
#' \item{krigSurfVoronoi}{krigSurfVoronoi}
#' \item{modelGen}{modelGen}
#' \item{VGMmodel}{VGMmodel}
#' \item{boundary}{boundary}
#' }
#' @export
#'
#' @examples
#' data(dataReg) #regular data on a square grid between 0 and 1
#' map = randKmapGrid(dataReg)
#' plotMap(map)
#' # not run
randKmapGrid=function(DataObj,nSimuCond=0,boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),manualBoundary=FALSE,disp=0,FULL=FALSE)
##########################################################
{
    # genData reads real data in DataObj data frame
    resGene=genData(DataObj=DataObj,Vnugget=0.5)
    if(is.null(resGene)) return(NULL)
    rawDataRaw=resGene$tabData
    boundary=resGene$boundary
 
    x=as.numeric(rawDataRaw[,1])
    xsize=max(x)-min(x)
    y=as.numeric(rawDataRaw[,2])
    ysize=max(y)-min(y)
    z=as.numeric(rawDataRaw[,3])

    #keep only pts within the boundary
    ptsInsRaw=point.in.polygon(x,y,boundary$x,boundary$y)#returns 0 if not within, 1 if strictly inside, 2 is inside an edge, 3 if vertex
    #
    rawDataNa=cbind(rawDataRaw,ptsIns=ptsInsRaw) #all pts with index column
    rawDataNa$z[ptsInsRaw==0]=NA # replace pts outside boundary with NA
    rawData=rawDataNa[ptsInsRaw!=0,1:3] #exclude pts outside boundary
    MatTest=rawData #data pts within boundary
    z=as.numeric(rawData[,3])
    #make spatial object
    sp::coordinates(rawDataNa)=~x+y
    ####define empty grid for kriging#####################
    # use empty grid of raw pts
    xempty= unique(MatTest$x)
    yempty= unique(MatTest$y)
    tempty = matrix(NA,nrow=length(xempty),ncol=length(yempty))
    colnames(tempty)= paste("Y",round(yempty,3),sep="")
    rownames(tempty)= paste("X",round(xempty,3),sep="")
    # fill matrix with data values
    # no kriging - pts are already on grid
    for (ii in 1:length(xempty))
    {
    	maskx=(MatTest$x==xempty[ii])
    	mati=MatTest[maskx,]
    	# complete mati to have a measure for each y value
    	masky=pmatch(mati$y,yempty)
    	tempty[ii,masky]=mati[,"z"]
    }
    # find pt locations within boundary
    ptsIns = point.in.polygon(xempty,yempty,boundary$x,boundary$y)
    maskIns=(ptsIns!=0)
    krigGrid=tempty[maskIns,] #  grid pt locations within boundary
    nx=nrow(krigGrid)
    ny=ncol(krigGrid)
    #Voronoi on data pts
    #prepare matrix
    nbP= nx*ny
    neighBool=matrix(logical(nbP^2),nbP,nbP)
    #
    resVoronoi=voronoiPolygons(rawDataNa,c(0,xsize,0,ysize),neighBool,FULL)
    neighNa=resVoronoi$neighBool
    surfVoronoiNa=resVoronoi$surfVoronoi
    surfVoronoi=surfVoronoiNa[maskIns] #remove pts outside boundary
    #pt neighborhood matrix
    neigh=neighNa[maskIns,maskIns]
    # compute list of neighbor pts
    krigN = ptNei(neigh) 
    #if required (argument FULL) compute voronoi on raw pts also
    if(FULL)
    {
	     krigVoronoi=resVoronoi$voronoi
	     nbPB= nrow(tabAleaNa)
      	     neighB=matrix(logical(nbPB^2),nbPB,nbPB)
     	     resVoronoiB=voronoiPolygons(tabAleaNa,neighB)
      	     voronoiB=resVoronoiB$voronoi
      	     surfVoronoiNaB=resVoronoiB$surfVoronoi
      	     surfVoronoiB=surfVoronoiNaB[tabAleaNa$ptsIns!=0]
      	     neighNaB=resVoronoiB$neighBool
      	     neighB=neighNaB[tabAleaNa$ptsIns!=0,tabAleaNa$ptsIns!=0]
     	      ptNB = vL(neighB)
    }

    # conditional simulation - to be added
    meanCondTab=NULL
    meanCondTabNa=NULL
    matMeanCond=NULL

    # reduced object size except if FULL argument was given
    if(FULL)
      return(list(rawData=rawData,step=step,krigData=rawData[ptsIns==1,],krigGrid=krigGrid,
                  krigN=krigN,krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,VGMmodel=resGene$VGMmodel,
                  boundary=boundary,krigVoronoi=krigVoronoi,matMeanCond=matMeanCond,boundary=boundary,
                  meanCondTab=meanCondTab,meanCondTabNa=meanCondTabNa,surfVoronoiB=surfVoronoiB,voronoiB=voronoiB,ptNB=ptNB))
    else
	    return(list(rawData=rawData,step=step,krigData=rawData[ptsIns==1,],krigGrid=krigGrid,krigN=krigN,
	                krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,VGMmodel=resGene$VGMmodel,
	                boundary=boundary))
}
