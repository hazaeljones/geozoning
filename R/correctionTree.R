####################################################################################
#' correctionTree
#'
#' @details description, a paragraph
#' @param qProb probability vector used to generate quantile values
#' @param map object returned by function genMap
#' @param pErr equality tolerance for distance calculations
#' @param optiCrit criterion choice
#' @param minSize zone area threshold under which a zone is too small to be manageable
#' @param minSizeNG zone area threshold under which a zone will be removed
#' @param distIsoZ threshold distance to next zone, above which a zone is considered to be isolated
#' @param simplitol tolerance for spatial polygons geometry simplification
#' @param LEQ length of quantile sequence used to grow isolated zone
#' @param MAXP quantile sequence maximum shift quantile sequence maximum shift
#' @param LASTPASS if TRUE, remove zones that are still too small at the last level of the correction tree
#' @param disp 0: no info, 1: some info, 2: detailed info
#' @param SAVE logical value, if TRUE function returns last level zonings, if FALSE function only returns best last level results
#' @param ONE logical value, if TRUE function returns only criterion value
#' @param ALL logical value, if TRUE function returns zonings at all levels
#'
#' @return a list with components
#'\describe{
#' \item{bestcrit}{best criterion value at last level (in all cases)}
#' \item{critList}{criterion values at last level (in all cases if ONE=FALSE)}
#' \item{costList}{cost values at last level (in all cases if ONE=FALSE)}
#' \item{costLList}{cost per label  values at last level (in all cases if ONE=FALSE)}
#' \item{nzList}{vector of number of zones at last level (in all cases if ONE=FALSE)}
#' \item{qProb}{vector of probabilities values used for quantiles (in all cases if ONE=FALSE)}
#' \item{zk}{list of zoning objects (such as returned by calNei function), first element corresponds to initial zoning, each other element is a list with each (last if ALL=FALSE) level zoning objects (only if SAVE=TRUE)}
#' \item{mdist}{list of initial distance matrix and all (last if ALL=FALSE) level distance matrices (only if SAVE=TRUE)}
#' \item{criterion}{list of initial criterion and all (last if ALL=FALSE) level criteria (only if SAVE=TRUE)}
#' \item{cost}{list of initial cost and all (last if ALL=FALSE) level costs   (only if SAVE=TRUE)}
#' \item{costL}{list of initial cost per label and all (last if ALL=FALSE) level costs per label (only if SAVE=TRUE)}
#' \item{nz}{list of initial number of zones and all (last if ALL=FALSE) level number of zones (only if SAVE=TRUE)}
#' }
#' @importFrom rgeos gArea
#'
#' @export
#'
#' @examples
#' data(mapTest)
# run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7
# saving initial zoning and last level zonings
#' criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE) 
#' plotZ(criti$zk[[1]][[1]]$zonePolygone)
#' plotZ(criti$zk[[2]][[1]]$zonePolygone) # zones 7 and 8 were handled
#'
correctionTree=function(qProb,map,pErr=0.9,optiCrit=2,minSize=0.012,minSizeNG=1e-3,distIsoZ=0.075,
                        simplitol=1e-3,LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=0,SAVE=TRUE,ONE=FALSE,ALL=FALSE)
####################################################################################
{
    # arguments
    # qProb=quantile probability vector
    # map=kriged data
    # choice of criterion = optiCrit
    # small zones handled by increasing size order
    # simplitol = tolerance for zone growing, polygone simplification
  

    # precaution !
  qProb=sort(unique(as.numeric(qProb)))
  #
  if(disp>0) print(paste("qProb=",qProb))
  #
  # results in global variable
  #######################################################
  #calcul du zonage initial avec le vecteur de quantiles donne en arg.
  resini = initialZoning(qProb,map,pErr,simplitol,optiCrit,disp)
  nbPoly = length(resini$resZ$zonePolygone)
  md0=resini$resDist$matDistanceCorr
  md0=normDistMat(md0,optiCrit)

  crit0=resini$resCrit
  cost0=resini$resDist$cost
  costL0=resini$cL
  nz0=nbPoly

  #######################################################
  # Detect smallZones
  ZI = detectSmallZones(resini$resZ$zonePolygone,minSize)

  listeZS = ZI$vectIndex
  if(disp>0)
	{
	print(paste(nbPoly,"zones,", length(listeZS)," small zones:"))
 	print(paste(listeZS,collapse=","))
	 }

  #######################################################
  #On recupere dans listOfZ les infos correspondant à ce decoupage
  #on recupere dans crit les criteres de tous les possibles zonages parcourus
  ###############################################################
     valRef= quantile(map$krigGrid,na.rm=TRUE,prob=qProb)
  ###############################################################
  listOfZ=list(list())
  crit=list(list())
  cost=list(list())
  costL=list(list())
  nz=list(list())
  mdist=list(list())
  listOfZ[[1]][[1]]=resini$resZ
  listOfZ[[1]][[1]]$qProb=qProb
  crit[[1]][[1]]=crit0
  cost[[1]][[1]]=cost0
  costL[[1]][[1]]=costL0
  nz[[1]][[1]]=nz0
  mdist[[1]][[1]]=md0

 if (disp>0) print(paste("level=1, initial crit=",round(crit0,3)))
 if(length(listeZS)==0) #no correction
  {

  if(ONE)
	return(round(crit0,3))
  else
  {
  	name = paste("q",length(qProb),sep="")
	li = list()
	li[name]=crit0
	critList=li
	li[name]=cost0
	costList=li
	li[name]=costL0
	costLList=li
	li[name]=nz0
	nzList=li

  	if (SAVE)
     	   return(list(bestcrit=round(crit0,3),critList=critList,costList=costList,costLList=costLList,nzList=nzList,zk=listOfZ,
                 mdist=mdist,criterion=crit,cost=cost,costL=costL,nz=nz))
  	else
		 return(list(bestcrit=round(crit0,3),critList=critList,costList=costList,costLList=costLList,nzList=nzList))
  }
  } # END NO CORRECTION CASE
  ###############################################################################################################

  # Tree traversal
  #Pour chaque zone trop petite soit on la supprime, soit on agrandit --> autant de resultats que de solutions possibles au probleme
  #on récupère le resultat
  #En toute rigueur, ordre de traitement des zones important.
  #pas important si les zones sont independantes.

  #etage arbre où lon se situe, cad quelle zone problematique on traite
  indCur=1
  
    #pour chaque zone a supprimer
    for (indZS in listeZS)
    {
     
      #Passage au prochain iter
      # Case of complete disparition
      curLen=length(listOfZ[[indCur]])
      
      if (curLen==0)
        {
	indCur=indCur-1
      	curLen=length(listOfZ[[indCur]])
	}
      else
      {
	#Add a stage
        listOfZ = append(listOfZ, list(list()))
      	crit=append(crit,list(list()))
	cost=append(cost,list(list()))
	costL=append(costL,list(list()))
	nz=append(nz,list(list()))
     	mdist=append(mdist,list(list()))
	
      }
      # make a copy for each branch

       if(disp>0)
	{
		cat("\n")
		print(paste("in loop level=",indCur+1,",zone to handle initial number (id)= ",resini$resZ$zonePolygone[[indZS]]@polygons[[1]]@ID,",",curLen, "branch(es) to examine "))
	}

      for (iter in (1:curLen))
      {

        checkSize=TRUE
        disparition = FALSE
        ##2 copies of current zoning, first for removal, 2nd for growing
	K=listOfZ[[indCur]][[iter]]
        zpCopy1 = K$zonePolygone
        zpCopy2 = zpCopy1
	
	#
        # iC=zone ne satisfaisant pas les contraintes
        iC = Identify(indZS,zpCopy1) # identique pour zpCopy2
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	# checking again because zone may now satisfy constraints
	#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      	if(disp>0) print(paste("iter=",iter," new zone number=",iC))
        #cas ou la zone a déja été integree dans une autre
        if (iC == 0 || gArea(zpCopy1[[iC]]) > minSize )
        {
          disparition = TRUE
	  if(disp>0)
		{
		print(paste("skipping zone current index",iC," initial index",indZS))
		if(iC !=0) print(paste("zone area=",gArea(zpCopy1[[iC]])))
		}
	  # keep only one copy for next level
	  zpCopy2 = NULL
        }
        else
        {
	       # 2 possibilities : include in zpCopy1 or grow in zpCopy2
	       # 1=merge zone indZS and zone near by
      	 Ns = getNs(K$zoneNModif,iC)
	       zpCopy1 = zoneFusion3(K,iC,Ns,map,minSize,simplitol,disp)
         if(disp>0) print(paste(length(zpCopy1)," polygons after zone merging"))
         # 2 = grow zone indZS

      	 zpCopy2 = zoneGrow(K,map,iC,optiCrit,minSizeNG,distIsoZ,LEQ,MAXP,simplitol,disp)
         if (disp>0) print(paste(length(zpCopy2)," zones after zone growing"))
         ###############################################################################################
        } # end else disparition

        #save infos for next iteration (indCur+1)
	Z=list(zpCopy1,zpCopy2)
	izk=0
	for (iz in 1:2)
	{
	izk=izk+1
	#only non NULL zonings are kept
	if (length(Z[[iz]])>0)
	   {
	# update crit[[indCur+1]], listOfZ, mdist
	   # keep only initial and current stages
       	   # except if ALL=TRUE, keep all stages
  	   
	   # saveZK appends a sublevel to listofZ[[indCur+1]]
	   resD=saveZK(map,K,Z[[iz]],qProb,listOfZ, indCur+1,crit,cost,costL,nz,mdist,pErr,optiCrit,simplitol)
      	   listOfZ=resD$listOfZ
	   mdist=resD$mdist
	   # save all criteria
	   crit=resD$crit
	   cost=resD$cost
	   costL=resD$costL
	   nz=resD$nz
	   }
	}#end for iz
	
     } # end for iter
     indCur=indCur+1
 # reuse allocated space for next level
  	   if((indCur>2)& !ALL)
	   {
	   listOfZ[indCur-1]=NULL
	   mdist[indCur-1]=NULL
	   crit[indCur-1]=NULL
	   cost[indCur-1]=NULL
	   costL[indCur-1]=NULL
	   nz[indCur-1]=NULL
	   indCur=indCur-1
	   }
    }# end for indZs

  # one more pass so that all zones are bigger than minSize in last level
  # simply remove zones of last level zonings that are too small and recalculate criteria
  if (LASTPASS)
	{
	  resPass=lastPass(map,qProb,listOfZ,crit,cost,costL,nz,mdist,pErr,optiCrit,minSize,simplitol,disp)
	  listOfZ = resPass$listOfZ
	  crit = resPass$crit
	  cost = resPass$cost
	  costL = resPass$costL
	  nz = resPass$nz
	  mdist = resPass$mdist
	}

  #  consider last step criteria
  #  sort last level criteria, return criteria and listOfZ if SAVE=TRUE, otherwise only return last level criteria
  #  and select the best one
  resC=sortCrit(qProb,crit,cost,costL,nz,mdist,listOfZ,map,disp,SAVE)
  #
  #garbage collection
  #gc()
 if (ONE)
        return(resC$bestcrit) #return single result (for optimization functions)
 else
	return(resC) #return full result
}


