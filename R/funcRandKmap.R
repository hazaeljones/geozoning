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
#' \item{VGMmodel}{VGM model}
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
genData=function(DataObj=NULL,seed=0,nPoints=450,typeMod="Gau",Vpsill=5,Vrange=0.2,Vmean=8,Vnugget=0,Vanis=1,
                 boundary=list(x=c(0,0,1,1,0),y=c(0,1,1,0,0)),manualBoundary=FALSE)
##############################################################################
{
  modelGen=NULL #variogram model
  VGMmodel1=NULL

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
    resNorm=datanormX(tabData,boundary)
   if(is.null(resNorm)) return(NULL)
    tabData = resNorm$dataN 
    boundary = resNorm$boundaryN
    xmin=resNorm$xmin
    xmax=resNorm$xmax
    ymin=resNorm$ymin
    ymax=resNorm$ymax
    xyminmaxI=rbind(c(xmin,xmax),c(ymin,ymax))
    rownames(xyminmaxI)=c("InitialX","InitialY")
    colnames(xyminmaxI)=c("min","max")
    x=tabData$x
    xsize=max(x)-min(x)
    y=tabData$y
    ysize=max(y)-min(y)
    boundary$x[boundary$x<0]=0
    boundary$y[boundary$y<0]=0
    boundary$x[boundary$x>xsize]=xsize
    boundary$y[boundary$y>ysize]=ysize
    
# fit experimental variogram to model
  tabDataSp=tabData
  coordinates(tabDataSp)=~x+y
  expVario=variogram(z~1,data=tabDataSp)
  VGMmodel1=fit.variogram(expVario,vgm(c("Exp","Sph","Gau"))) # find best model to be fitted

  }
  # simulated data
  else{
     #Generate random (x,y) values within unit square
    set.seed(seed)
    x=runif(nPoints, min=0, max=1)
    y=runif(nPoints, min=0, max=1)

    #Generate z values according to (Gaussian) field
    #RMmodel starting from VGM model
    VGMmodel=vgm(model=typeMod,range=Vrange,psill=Vpsill,mean=Vmean,ang1=Vang,anis1=Vanis)
    VGMmodel1=vgm(model=typeMod,range=Vrange,psill=Vpsill,mean=Vmean,ang1=Vang,anis1=Vanis,nugget=Vnugget)

  modelGen=calRMmodel(VGMmodel)
    modelGen=modelGen+RMtrend(mean=Vmean)
    if(Vnugget>1e-3)
	    modelGen=modelGen+RMnugget(var=Vnugget)

    testMap<-RFsimulate(modelGen,x,y)
     #store in dataframe
    tabData=data.frame(x=x,y=y,z=testMap$variable1)
    # normalize x,y coordinates
   tabData=datanormXY(tabData)
    xyminmaxI=rbind(c(0,1),c(0,1))
    rownames(xyminmaxI)=c("InitialX","InitialY")
    colnames(xyminmaxI)=c("min","max")
  }

  return(list(tabData=tabData,boundary=boundary,xyminmaxI=xyminmaxI,modelGen=modelGen,VGMmodel=VGMmodel1))
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
calStep=function(nPointsK,xsize,ysize)
#####################################################
{
  #compute step to obtain square grid
  stepA=1
  stepB=-2/(xsize+ysize)
  stepC=-(nPointsK-1)/(xsize*ysize)
  stepApprox=(-stepB+sqrt(stepB^2 -4*stepA*stepC))/2
  step=1/trunc(stepApprox)
  return(step)
}

################################################################
#' generate grid from raw data
#'
#' @param step numeric
#' @param nKrigE numeric
#'
#' @importFrom sp coordinates
#'
#' @return a dataframe that contains kriged positions based on original ones
#' @export
#'
#' @examples
#' # not run
genEmptyGrid=function(step,xsize,ysize)
################################################################
{
 # generate grid from raw data

  xx=seq(from=step, to=xsize-step,by=step)
  yy=seq(from=step, to=ysize-step,by=step)
  nx=length(xx)
  ny=length(yy)
  xempty=rep(xx,times=ny)
  yempty=rep(yy,each=nx)
  z=rep(NA,nx*ny)
  
  # turn into dataframe 
  tabEmpty=data.frame(x=xempty,y=yempty,z=z)
  coordinates(tabEmpty)=~x+y
  return(list(tabEmpty=tabEmpty,xx=xx,yy=yy,nx=nx,ny=ny))
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
