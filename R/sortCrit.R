###############################################################################
#' sortCrit called by correctionTree
#'
#' @details sort last level criteria from list of zonings, return criteria and list of zonings if SAVE=TRUE, otherwise only return last level criteria
#' @param qProb probability vector used to generate quantile values
#' @param crit list of criteria
#' @param cost  list of costs
#' @param costL list of per label costs 
#' @param nz list of number of zones 
#' @param mdist list of distance matrices
#' @param listOfZ list of zoning objects
#' @param map object returned by function genMap or genMapR
#' @param disp 0: no info, 1: plot best corrected zoning
#' @param SAVE logical value, if TRUE function returns more elements
#'
#' @return a list with components
#'\describe{
#' \item{bestcrit}{best criterion value at last level}
#' \item{critList}{criterion values at last level}
#' \item{costList}{cost values at last level}
#' \item{costLList}{cost per label  values at last level}
#' \item{nzList}{vector of number of zones at last level}
#' \item{qProb}{vector of probabilities values used for quantiles}
#' \item{zk}{(SAVE=TRUE) list of zoning objects (such as returned by calNei function), first element corresponds to initial zoning, each other element is a list with each (last if ALL=FALSE) level zoning objects}
#' \item{mdist}{(SAVE=TRUE) list of initial distance matrix and all (last if ALL=FALSE) level distance matrices}
#' \item{crit}{(SAVE=TRUE) list of initial criterion and all (last if ALL=FALSE) level criteria }
#' \item{cost}{(SAVE=TRUE) list of initial cost and all (last if ALL=FALSE) level costs  }
#' \item{costL}{(SAVE=TRUE) list of initial cost per label and all (last if ALL=FALSE) level costs per label}
#' \item{nz}{(SAVE=TRUE) list of initial number of zones and all (last if ALL=FALSE) level number of zones}
#' }
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' qProb=c(0.4,0.7)
#' criti=correctionTree(qProb,mapTest)
# displays best criterion, corresponding costs and number of zones
#' sortCrit(qProb,criti$criterion,criti$cost,criti$costL,criti$nz,criti$mdist,criti$zk,mapTest) 
#' # not run
sortCrit=function(qProb,crit,cost,costL,nz,mdist,listOfZ,map,disp=0,SAVE=FALSE)
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
    }
  } #end while
  # plot best corrected zoning
  if (disp==2)
  {
    # IS 19/05/2017: add comment for x11
    #x11()
    dispZ(map$step,map$krigGrid,zonePolygone=bestZ,boundary=map$boundary,nbLvl=0)
  }
  if (SAVE)
     return(list(bestcrit=round(bestcrit,3),critList=critList,costList=costList,costLList=costLList,nzList=nzList,
                 zk=listOfZ,mdist=mdist,criterion=crit,cost=cost,costL=costL,nz=nz,qProb=qProb))
  else
	return(list(bestcrit=round(bestcrit,3),critList=critList,costList=costList,
	            costLList=costLList,nzList=nzList,qProb=qProb))
  }
