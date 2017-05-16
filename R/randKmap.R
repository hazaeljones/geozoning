###########################################################################
#' randKmap : todo title of the function in one line
#'
#' @details description, a paragraph
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
#' @return a list
#' \describe{
#' \item{rawData}{rawData}
#' \item{step}{step}
#' \item{krigData}{krigData}
#' \item{krigGrid}{krigGrid}
#' \item{krigN}{krigN}
#' \item{krigSurfVoronoi}{krigSurfVoronoi}
#' \item{modelGen}{modelGen}
#' \item{modelVGM}{modelVGM}
#' \item{boundary}{boundary}
#' }
#'
#' @export
#' @importFrom gstat krige
#'
#' @examples
#' # not run
#' # map<-randKmap(DataObj,seed=seed,Vpsill=5,Vrange=0.2,Vnugget=0,Vmean=8,krig=1,disp=0,FULL=FALSE)
#'
randKmap=function(DataObj,seed,nPoints=450,nPointsK=2000,nSimuCond=0,typeMod="Gau",Vpsill=5,Vrange=0.2,Vmean=8,Vnugget=0,
                  boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),manualBoundary=FALSE,krig=1,disp=0,FULL=FALSE)
###########################################################################
{

 #simulation seed for random fields
      set.seed(seed)

 # genData reads real data in DataObj data frame
 # or simulates them if DataObj=NULL
      resGene=genData(DataObj,seed,nPoints,typeMod,Vpsill,Vrange,Vmean,Vnugget,boundary,manualBoundary)
      if(is.null(resGene)) return(NULL)
      rawDataRaw=resGene$tabData
      boundary=resGene$boundary

      x=as.numeric(rawDataRaw[,1])
      y=as.numeric(rawDataRaw[,2])
      z=as.numeric(rawDataRaw[,3])
      nPoints=length(z)

      #keep only pts within the boundary
      ptsInsRaw=point.in.polygon(x,y,boundary$x,boundary$y)#returns 0 if not within, 1 if strictly inside, 2 is inside an edge, 3 if vertex
      #
      rawDataNa=cbind(rawDataRaw,ptsIns=ptsInsRaw) #all pts with index column
      rawDataNa$z[ptsInsRaw==0]=NA # replace pts outside boundary with NA
      rawData=rawDataNa[ptsInsRaw!=0,1:3] #exclude pts outside boundary
      #MatTest=rawData
      z=as.numeric(rawData[,3])
      #make spatial object
      coordinates(rawDataNa)=~x+y
      coordinates(rawData)=~x+y
      ####define empty grid for kriging#####################
      # - compute step
      step=calStep(nPointsK)
      # - Effective number of kriged pts according to step
      nKrigE=(1*step-1)*(1*step-1)
       # - generate square grid of kriged pt locations
      gridK=genEmptyGrid(step,nKrigE)
      xempty=as.numeric(gridK$x) #to save space in names
      yempty=as.numeric(gridK$y)
      vecTabAlea=gridK$z
      # find future kriged pt locations within boundary
      ptsIns = point.in.polygon(xempty,yempty,boundary$x,boundary$y)
      maskIns=(ptsIns!=0)
      # generate kriged pts at these locations
      # sometimes matrix is singular with krig=1
      if (krig == 2) #  inverse distance
      	 krigTabAleaPart=krige(z~1,rawData,newdata=gridK[maskIns,])
      else if (krig == 1) # vgm model
      	 krigTabAleaPart=krige(z~1,rawData,newdata=gridK[maskIns,],model=resGene$modelVGM)
	 else return(NULL)
      #transform into grid matrix
      vecTabAlea[maskIns]=as.numeric(krigTabAleaPart$var1.pred)
      krigData=data.frame(x=xempty,y=yempty,var1.pred=vecTabAlea)
      coordinates(krigData)=~x+y
      #NAs outside boundary
      krigGrid=matrix(vecTabAlea,1*step -1,1*step -1)
      colnames(krigGrid)=round(seq(1/step,1-1/step,by=1/step),3)
      rownames(krigGrid)=round(seq(1/step,1-1/step,by=1/step),3)

      #avoid side effects by using krigDataNa
      krigDataNa=cbind(data.frame(krigData),ptsIns)
      coordinates(krigDataNa)=~x+y

      ##Voronoi on kriged pts
      #prepare matrix
      nbP= nrow(gridK)
      neighBool=matrix(logical(nbP^2),nbP,nbP)
      #
      resVoronoi=voronoiPolygons(krigDataNa,neighBool,FULL)
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
            return(list(rawData=rawData,step=step,krigData=krigData[ptsIns==1,],krigGrid=krigGrid,krigN=krigN,
                        krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,modelVGM=resGene$modelVGM,boundary=boundary,
                        krigVoronoi=krigVoronoi,matMeanCond=matMeanCond,boundary=boundary,meanCondTab=meanCondTab,
                        meanCondTabNa=meanCondTabNa,surfVoronoiB=surfVoronoiB,voronoiB=voronoiB,ptNB=ptNB))
       else
	    return(list(rawData=rawData,step=step,krigData=krigData[ptsIns==1,],krigGrid=krigGrid,krigN=krigN,
	                krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,modelVGM=resGene$modelVGM,boundary=boundary))

}
