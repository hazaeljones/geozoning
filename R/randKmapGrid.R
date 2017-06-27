##########################################################
#' randKmapGrid
#'
#' @details description, a paragraph
#' @param DataObj xxxx
#' @param seed xxxx
#' @param nPointsK xxxx
#' @param nSimuCond xxxx
#' @param typeMod xxxx
#' @param Vpsill xxxx
#' @param Vrange xxxx
#' @param Vmean xxxx
#' @param boundary xxxx
#' @param manualBoundary xxxx
#' @param krig xxxx
#' @param disp xxxx
#' @param FULL xxxx
#'
#' @importFrom sp coordinates
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
randKmapGrid=function(DataObj,seed=NULL,nPointsK=200,nSimuCond=0,typeMod="Gau",Vpsill=5,Vrange=0.2,
                      Vmean=8,boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),
                      manualBoundary=FALSE,krig=1,disp=0,FULL=FALSE)
##########################################################
{
    # genData reads real data in DataObj data frame
    # or simulates them if DataObj=NULL
    # IS 19/05/2017 suppress Vang, Vanis and alphavario parameters not used!
    resGene=genData(DataObj=DataObj,seed=seed,nPoints=nPoints,typeMod=typeMod,Vpsill=Vpsill,
                    Vrange=Vrange,Vmean=Vmean,boundary=boundary,manualBoundary=manualBoundary)
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
    MatTest=rawData #data pts within boundary
    z=as.numeric(rawData[,3])
    #make spatial object
    sp::coordinates(rawDataNa)=~x+y
    sp::coordinates(rawData)=~x+y
    ####define empty grid for kriging#####################

    ################structure à kriger(definition de la grille)#################################
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
    krigGrid=tempty[maskIns,] #  pt locations within boundary

    krigDataNa=rawDataNa
    krigData=krigDataNa
    sp::coordinates(krigDataNa)=~x+y
    step=calStep(nrow(rawDataRaw))

    #Voronoi on kriged pts
    #prepare matrix
    nbP= nrow(gridK)
    neighBool=matrix(logical(nbP^2),nbP,nbP)
    #
    resVoronoi=voronoiPolygons(krigDataNa,neighBool,FULL)
    neighNa=resVoronoi$neighBool
    surfVoronoiNa=resVoronoi$surfVoronoi
    surfVoronoi=surfVoronoiNa[maskIns] #remove pts outside boundary
    #pt neighborhood matrix
    neigh=neighNa[maskIns,maskIns]
    # compute list of neighbor pts
    krigN = vL(neigh)
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
      return(list(rawData=rawData,step=step,krigData=krigData[ptsIns==1,],krigGrid=krigGrid,
                  krigN=krigN,krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,VGMmodel=resGene$VGMmodel,
                  boundary=boundary,krigVoronoi=krigVoronoi,matMeanCond=matMeanCond,boundary=boundary,
                  meanCondTab=meanCondTab,meanCondTabNa=meanCondTabNa,surfVoronoiB=surfVoronoiB,voronoiB=voronoiB,ptNB=ptNB))
    else
	    return(list(rawData=rawData,step=step,krigData=krigData[ptsIns==1,],krigGrid=krigGrid,krigN=krigN,
	                krigSurfVoronoi=surfVoronoi,modelGen=resGene$modelGen,VGMmodel=resGene$VGMmodel,
	                boundary=boundary))
}
