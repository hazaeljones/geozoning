##############################################################################
#' generate data
#'
#' @details description, a paragraph
#' @param DataObj =NULL: simulated data with seed or = a data frame with real data
#' @param seed numeric,
#' @param nPoints numeric,
#' @param typeMod character
#' @param Vpsill numeric,
#' @param Vrange numeric,
#' @param Vmean numeric,
#' @param Vnugget numeric,
#' @param boundary list, contains x and y
#' @param manualBoundary logical,
#'
#' @return a list
#' \describe{
#' \item{tabData}{tabData}
#' \item{boundary}{boundary}
#' \item{modelGen}{modelGen}
#' \item{modelVGM}{modelVGM}
#' }
#'
#' @export
#' @importFrom stats runif
#' @importFrom sp coordinates
#' @importFrom gstat vgm
#' @importFrom RandomFields RMtrend RMnugget RFsimulate
#'
#' @examples
#' # not run
#' # resGene=genData(DataObj,0,450,"Gau",5,0.2,8,0,list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),FALSE)
genData=function(DataObj=NULL,seed=0,nPoints=450,typeMod="Gau",Vpsill=5,Vrange=0.2,Vmean=8,Vnugget=0,
                 boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),manualBoundary=FALSE)
##############################################################################
{
  modelGen=NULL #variogram model
  modelVGM1=NULL

  # real or simulated data
  if(!is.null(DataObj))
    print(paste("reading DataObj,nrow(DataObj)=",nrow(DataObj),",ncol(DataObj)=",ncol(DataObj),collapse=","))
  else
    print(paste("DataObj=NULL, generating DataObj-seed=",seed))

  if(!is.null(DataObj)){
    #read data frame x y z
    #remove duplicated data pt locations
    tabData=DataObj[ ! duplicated(DataObj[,c(1,2)]),]
    names(tabData)=c("x","y","z")
    #draw boundary if required
    if(manualBoundary)
      {
        print("Draw boundary")
        plot(coordinates(tabData))
        boundary=locator(500,type="l")

        boundary$x[length(boundary$x)]=boundary$x[1]
        boundary$y[length(boundary$y)]=boundary$y[1]
      }

  #Normalize coordinates and boundary
    resN=datanorm(tabData,boundary)
   if(is.null(resN)) return(NULL)

  #Normalize boundary
    # eliminate pts outside field edges
    boundary$x[boundary$x<0]=0
    boundary$y[boundary$y<0]=0
    boundary$x[boundary$x>1]=1
    boundary$y[boundary$y>1]=1
  }
  # simulated data
  else{
     #Generate random (x,y) values within unit square
    set.seed(seed)
    x=runif(nPoints, min=0, max=1)
    y=runif(nPoints, min=0, max=1)

    #Generate z values according to (Gaussian) field
    #RMmodel starting from VGM model
    modelVGM=vgm(model=typeMod,range=Vrange,psill=Vpsill,mean=Vmean,ang1=Vang,anis1=Vanis)
    modelVGM1=vgm(model=typeMod,range=Vrange,psill=Vpsill,mean=Vmean,ang1=Vang,anis1=Vanis,nugget=Vnugget)
    modelGen=calRMmodel(modelVGM)
    modelGen=modelGen+RMtrend(mean=Vmean)
    if(Vnugget>1e-3)
	    modelGen=modelGen+RMnugget(var=Vnugget)

    testMap<-RFsimulate(modelGen,x,y)
     #store in dataframe
    tabData=data.frame(x=x,y=y,z=testMap$variable1)
    # normalize x,y coordinates
   tabData=datanormXY(tabData)
  }

  return(list(tabData=tabData,boundary=boundary,modelGen=modelGen,modelVGM=modelVGM1))
}



#####################################################
#' compute step for square grid
#'
#' @param nPointsK numeric
#'
#' @return a numeric
#' @export
#'
#' @examples
#' # not run
calStep=function(nPointsK)
#####################################################
{
  #compute step for square grid
  stepA=1
  stepB=-2/(1*1)
  stepC=-(nPointsK-1)/(1*1)
  stepApprox=(-stepB+sqrt(stepB^2 -4*stepA*stepC))/2
  step=trunc(stepApprox)
  return(step)
}

################################################################
#' generate grid from raw data
#'
#' @param step numeric
#' @param nKrigE numeric
#'
#' @return a dataframe that contains kriged positions based on original ones
#' @export
#'
#' @examples
#' # not run
genEmptyGrid=function(step,nKrigE)
################################################################
{
  # generate grid from raw data
  xempty=rep(seq(1/step, 1 -1/step, by=(1/step)),(1*step -1))
  yempty=as.vector( t( matrix( xempty,  1*step -1, (1*step -1)) ) )
  z=rep(NA,nKrigE)

  #dataframe des positions et valeurs des points à kriger à partir des originaux
  tabEmpty=data.frame(x=xempty,y=yempty,z)
  coordinates(tabEmpty)=~x+y
  return(tabEmpty)
}


####################################################################
#' returns list of pt neigbors for each pt
#'
#' @param neighBool numeric, boolean neighborhood matrix for pts
#'
#' @return a list of pt neigbors for each pt
#' @export
#'
#' @examples
#' # not run
ptNei=function(neighBool)
###################################################################
{
  # arg  = boolean neighborhood matrix for pts
  # returns list of pt neigbors for each pt
	vPt=apply(neighBool,1,function(x){return(grep(TRUE,x))})
	return(vPt)
}
