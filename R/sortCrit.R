###############################################################################
#' sortCrit
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param crit xxxx
#' @param cost xxxx
#' @param costL xxxx
#' @param nz xxxx
#' @param mdist xxxx
#' @param listOfZ xxxx
#' @param map xxxx
#' @param disp xxxx
#' @param SAVE xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
sortCrit=function(qProb,crit,cost,costL,nz,mdist,listOfZ,map,disp,SAVE=FALSE)
###############################################################################
{

# for optim, single result required (crit)
# all results are saved  (if FULL=T)
le=length(crit)
if (disp>0)
  {
	for (ii in 1:le)
	{
		criti=crit[[ii]]

		print(paste("length(crit[[",ii,"]])=",length(criti),collapse=""))
  		print(unlist(criti))
	}

  } #end disp
  bestcrit=0
  bestK=listOfZ[[1]][[1]]
  bestZ=bestK$zonePolygone

  b=TRUE
  while (b && le>=2)
  {
  if (length(crit[[le]]) == 0 )
  {
	le=le-1
  }
  else
  {
    b=FALSE
    #number of quantiles
    nq=length(qProb)

  # last level is the right one to keep
    best=searchNODcrit(qProb,le,listOfZ,crit,cost,costL,nz)
    indList=best$ind
    critList=best$critList
    costList=best$costList
    costLList=best$costLList
    nzList=best$nzList
    #save all criterion results for final level
    #critList<<-critList
    #costList <<- costList
    #costLList <<- costLList
    #nzList <<- nzList
    #
  # best crit for non degenerated qs
  # if there is none, then best of best for all qs
    cr=sapply(critList,max)
    n=names(critList)
    labn=names(cr)
    #cri=which(cr==max(cr)) #best of best for all qs
    #
    labq=paste("q",nq,sep="")
    nodeg=match(labq,labn)
    if (!is.na(nodeg))
       cri=cr[nodeg]
    else
	cri=which(cr==max(cr))
    #
    cri=cri[1]
    labi=names(cri)
    ind=indList[[labi]]
    ind=ind[1]
    bestcrit=max(crit[[le]][[ind]])
    bestcrit=bestcrit[1]
    #
    bestK=listOfZ[[le]][[ind]]
    bestZ=bestK$zonePolygone
    bestmdist=mdist[[le]][[ind]]
    #
    if(SAVE) bestZ<<-bestZ
  }
  } #end while
# plot best corrected zoning
  if (disp==2)
  {
  x11()
  dispZ(map$step,map$krigGrid,zonePolygone=bestZ,boundary=map$boundary,nbLvl=0)
  }
  if (SAVE)
     return(list(bestcrit=round(bestcrit,3),critList=critList,costList=costList,costLList=costLList,nzList=nzList,zk=listOfZ,mdist=mdist,criterion=crit,cost=cost,costL=costL,nz=nz,qProb=qProb))
  else
	return(list(bestcrit=round(bestcrit,3),critList=critList,costList=costList,costLList=costLList,nzList=nzList,qProb=qProb))
  }
