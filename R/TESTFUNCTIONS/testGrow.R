###############################################################################
testModifIso=function(qProb,map,Z,K,iC,optiCrit,minSize,minSizeNG,distIsoZ=0.075,LEQ=5,MAXP=0.1,simplitol=1e-12,disp=1)
###############################################################################
{
	# tester agrandissement zone isolee
	# env # 0 = dilate zone and use bbox as enveloppe
	print("zonage initial")
	printZid(Z)
	print(paste("growing isolated zone: ",iC))
	x11()
	dispZ(map$step,map$krigGrid,zonePolygone=Z,nbPoly=length(Z),boundary=map$boundary,nbLvl=0)
    	# zones proches
	Ns = getNs(K,iC)
	#
	zmi = NULL
	 ##############################################################
	 #fonction qui trouve le meilleur quantile pour agrandir la zone
         res = optiGrow(Z,K,iC,qProb,map,optiCrit,minSize,minSizeNG,distIsoZ,LEQ,MAXP,simplitol,disp)
	 ##############################################################
	 if (!is.null(res))
	 {
	 x11()
	 zmi=res$Zopti
   	 dispZ(map$step,map$krigGrid,zonePolygone=zmi,nbPoly=length(zmi),boundary=map$boundary,nbLvl=0)
	 }
	 return(zmi)
	 }
##############################################################
testModifNonIso=function(qProb,map,Z,K,iC,simplitol)
##############################################################
{
	# test growing non isolated zone
	print("initial zoning IDs")
	for (ii in 1:length(Z)){print(paste("ii=",ii," ID=", Z[[ii]]@polygons[[1]]@ID))}
	x11()
	dispZ(map$step,map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
	# agrandissement zone iC
	resP = detZoneClose(iC,Z,K) # renvoie FALSE si zone trop proche dune autre, TRUE sinon
        ##############################################################
        InterZoneSpace = resP$InterZoneSpace
        zoneClose = resP$zoneClose
	
	zmi = zoneModifnonIso(Z,K,qProb,map,zoneClose,iC,simplitol,disp=TRUE)
	 
   	 if (!is.null(zmi))
	    {
	    x11()	
	    dispZ(map$step,map$krigGrid,zonePolygone=zmi,boundary=map$boundary,nbLvl=0)
	    }
	 return(zmi)
	
}	
##############################################################
testOptiRegroup=function(Z,K,i1,i2,simplitol,map)
##############################################################
{
	zmi = optiRG(Z,K,map,i1,i2,simplitol,disp=1)	
	 x11()
   	 dispZ(map$step,map$krigGrid,zonePolygone=zmi,boundary=map$boundary,nbLvl=0)
	 return(zmi)
}

##############################################################
testOptiRegroupDiff=function(Z,K,i1,i2,map,qProb)
##############################################################
# !!!!! pb !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##############################################################
{
	valRef= quantile(map$krigGrid,na.rm=TRUE,prob=qProb)
	refPoint = gCentroid(Z[[i1]])
	 x11()
	 dispZ(map$step,map$krigGrid,zonePolygone=Z,nbPoly=length(Z),boundary=map$boundary,nbLvl=0)
	 delim=  calFrame(i1,K, 0.15, 0.075)
	plotCad(delim)
	new = funOptiRegroupDiff(Z,K,i1,i2,qProb,map,delim,valRef,refPoint)
	if (is.null(new))  return(NULL)
	
	vRef = new$optimum
        polyUni = new$polyUni
	indiceZE = new$indiceZE
				                    	
        if(vRef !=0)
	{
                 Zoptizg=crZones(Z,i1,i2,indiceZE,map,delim,vRef,polyUni)
                          	
	 	x11()
   	 	dispZ(map$step,map$krigGrid,zonePolygone=zg,nbPoly=length(zg),boundary=map$boundary,nbLvl=0)
	 	return(zg)
	}
	else return(NULL)
}

##############################################################
testCad = function(qProb,map,Z,K,iC,wcad=0.15)
##############################################################
{
	Ns = getNs(K,iC)
	#
	refPoint = gCentroid(Z[[iC]])
	#surface de la petite zone
        refSurf = getSurf(Z,iC)
        #On recupere les points qui sont proches de la zone pour pouvoir les utiliser pour generer de nouveaux contours!!!!!!!!!
	# calculCadre2 reduit le cadre pour ne pas intersecter avec
	# les zones non voisines
	
        cad = calculCadre2(iC,Ns,K, wcad, map,Z,intersec=TRUE)
	if (is.null(cad)) return (NULL)
	#
	plotCad(cad)
	return(cad)
	
}

##############################################################
testGrowM = function(iC,Z,K,map,qProb,LEQ=5,MAXP=0.1,distIsoZ=0.075,simplitol=1e-3,optiCrit=2)
##############################################################
{
#source("prep_simu.R") #to get map
#qProb=c(0.1,0.9)
#load("testGr") #resu optim
#Z=li$Z
#K=li$K
plotZ(Z)

iE=detZoneEng(iC,Z,K)

if (iE == 0)
{
print("no englobing zone")
return(NULL)
}
else
   {
   envel=calFrame(iC,Z,K,distIsoZ)
   if (is.null(envel))
   {
   print("no envelope")
   return(NULL)
   }
   linesSp(envel)
   #ptsN1 = ptsInSp(envel,map$krigTabAlea,hole=FALSE) 
   #points(ptsN1@coords,col="red")
   #valRef= quantile(map$krigGrid,na.rm=TRUE,prob=qProb)

   Qseq = genQseq(qProb,K,map,iC,iE,LEQ,MAXP,disp=1)

   critG=rep(0,length(Qseq))
   area=critG
   Zopt=list()
   refPoint = gCentroid(Z[[iC]])
   
   for (i in 1:length(Qseq))
       {
       resi = findCinZ(iC,Z,K,map,Qseq[i],envel)
       # returns NULL if no grow
       Zopti=NULL
       if(!is.null(resi))
       {
	resp = checkContour(resi$contourSp,map$step,refPoint,minSizeNG)
	# if  condition not met try next contour
	if (is.null(resp)) next
	Zopti=zoneQ(resi$contourSp,iC,iE,Z,K,map,simplitol) # current zone is now last zone
	 if (!is.null(Zopti))
      	 {
	 # create comments for holes
	 Zopti = crComment(Zopti)
	 #
	 criti = calcDCrit(Zopti,map,optiCrit)
	 if (!is.null(criti))
	   {
	   critG[i]=criti$resCrit
	   Zopt[[i]]=Zopti
	   area[i]=gArea(Zopti[[findNumZ(Zopti,iC)]])
	   } # end criti not null
       	 } # end Zopti not null
	} #end contour found for ith quantile 
       
       } #end loop on Qseq
       
   # find best criterion for zone area > minSize
       n = rev(order(critG))
       mask = area[n]>minSize
       if (any(mask))
       {
       nm = n[mask]
       nm = nm[rev(order(area[nm]))] # sort by area
       iM = nm[1]
       } else # area condition not satisfied - take biggest area
	{
	 n = rev(order(area))
	 iM = n[1]
	}

   } #end englobing zone exists
   
if (length(Zopt)>0)
{
	numiC = findNumZ(Zopt[[iM]],iC)
	linesSp(Zopt[[iM]][[numiC]],col="blue")
	}
else print(paste("impossible to grow zone",iC))

return(list(crit=critG,Zopt=Zopt, area=area, iM=iM))
}