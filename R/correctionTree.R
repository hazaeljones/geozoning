####################################################################################
#' correctionTree
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param map xxxx
#' @param pErr xxxx
#' @param optiCrit xxxx
#' @param minSize xxxx
#' @param minSizeNG xxxx
#' @param distIsoZ xxxx
#' @param simplitol xxxx
#' @param LEQ xxxx
#' @param MAXP xxxx
#' @param LASTPASS xxxx
#' @param disp xxxx
#' @param SAVE xxxx
#' @param ONE xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
correctionTree=function(qProb,map,pErr=0.9,optiCrit=2,minSize=0.012,minSizeNG=1e-3,distIsoZ=0.075,
                        simplitol=1e-3,LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=0,SAVE=TRUE,ONE=FALSE)
####################################################################################
{
    # arguments
    # qProb=quantile probability vector
    # map=kriged data
    # choice of criterion = optiCrit
    # small zones handled by increasing size order
    # simplitol = tolerance for zone growing, polygone simplification
    # disp =0 : no info
    # disp =1 : print info
    # disp =2 : detailed info
    # ONE=TRUE: returns only criterion value

    # precaution !
  qProb=sort(unique(qProb))
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
     return(list(bestcrit=round(crit0,3),critList=critList,costList=costList,costLList=costLList,nzList=nzList,zk=listOfZ,mdist=mdist,criterion=crit,cost=cost,costL=costL,nz=nz))
  else
	return(list(bestcrit=round(crit0,3),critList=critList,costList=costList,costLList=costLList,nzList=nzList))
	return(resC)
  }
  } # END NO CORRECTION CASE
  ###############################################################################################################

  # Tree traversal
  #Pour chaque zone trop petite soit on la supprime, soit on agrandit --> autant de resultats que de solutions possibles au probleme
  #on récupère le resultat
  #En toute rigueur, ordre de traitement des zones important.
  #pas important si les zones sont independantes.

  #etage arbre où lon se situe, cad quelle zone problematique on traite
  counter = 0

    #pour chaque zone a supprimer
    for (indZS in listeZS)
    {
      #Passage au prochain iter
      counter=counter +1
      # Cas de disparition complete
      curLen=length(listOfZ[[counter]])
      if (curLen==0)
        {
	counter=counter-1
      	curLen=length(listOfZ[[counter]])
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
		print(paste("in loop level=",counter+1,",zone to handle initial number= ",resini$resZ$zonePolygone[[indZS]]@polygons[[1]]@ID,",",curLen, "branch(es) to examine "))
	}

      for (iter in (1:curLen))
      {

        checkSize=TRUE
        disparition = FALSE
        ##On doit faire 2 copies du zonage, une pour la suppression, lautre pour lagrandissement
	K=listOfZ[[counter]][[iter]]
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
      	Ns = getNs(K,iC)
	zpCopy1 = zoneFusion3(zpCopy1,K,iC,Ns,map,minSize,simplitol,disp)
         if(disp>0) print(paste(length(zpCopy1)," polygons after zone merging"))
         # 2 = grow zone indZS

      	zpCopy2 = zoneGrow(zpCopy2,K,iC,Ns,map,optiCrit,valRef,qProb,minSizeNG,distIsoZ,LEQ,MAXP,simplitol,disp)
        if (disp>0) print(paste(length(zpCopy2)," polygons after zone growing"))
        ###############################################################################################
        } # end else disparition

        #save infos for next iteration (counter+1)
	Z=list(zpCopy1,zpCopy2)
	for (iz in 1:2)
	{
	#only non NULL zonings are kept
	if (length(Z[[iz]])>0)
		{
	   # update crit[[counter+1]], listOfZ, mdist
		resD=saveZK(map,K,Z[[iz]],qProb,listOfZ, counter,crit,cost,costL,nz,mdist,pErr,optiCrit,simplitol)
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

    }# end for indZs

# one more pass so that all zones are bigger than minSize in last level
# simply remove zones of last level zonings that are too small and recalculate criteria
if (LASTPASS)
	{
	resPass=lastPass(map,qProb,listOfZ,crit,cost,costL,nz,mdist,pErr,optiCrit,simplitol)
	listOfZ = resPass$listOfZ
	crit = resPass$crit
	cost = resPass$cost
	costL = resPass$costL
	nz = resPass$nz
	mdist = resPass$mdist
	}

#  consider last step criteria
#  sort criteria, assign zoning to global variables zf, zk, critere and critList (if SAVE=TRUE)
#  and select the best one
resC=sortCrit(qProb,crit,cost,costL,nz,mdist,listOfZ,map,disp,SAVE)
#

 if (ONE)
        return(resC$bestcrit) #return single result (for optimization functions)
 else
	return(resC) #return full result


}


