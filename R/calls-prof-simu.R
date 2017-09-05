rm(list=ls())
source("srcZ.R") # source libraries and functions - params are in initParam.R
#

################################################################
# step 1 - generate 2D map from simulation - default = 450 pts
################################################################
# kriging default 2000 pt grid
seed=80
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2,Vnugget=Vnugget)

# display 2D map
plotMap(map)
#

################################################################
# step 2 - generates zoning Z from map, quantile vector (0.4 here)-no display
################################################################
ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
# initialZoning calls zoneGeneration, calNei, calDistance, calCrit
# zoneGeneration calls contourAuto
# contourAuto calls contourLines and extensionLine
#
# plot zoning (8 zones in this case)
Z=ZK$resZ$zonePolygone
K=ZK$resZ
DC0=ZK$resD$matDistanceCorr
DC0N=normDistMat(DC0,2)

plotZ(Z) #only zones
# more detailed plot
valRef=quantile(map$krigGrid,na.rm=TRUE,prob=c(0.4,0.7))
dispZ(map$step,map$krigGrid,zonePolygone=Z,K=K,boundary=map$boundary,nbLvl=0,id=FALSE,mu=2)
 title(paste(" q=[",toString(round(valRef,2)),"]   crit=",round(ZK$resCrit,2),sep=""))
# print zoning labels
printLabZ(list(ZK$resZ))
# print zoning surfaces
printZsurf(ZK$resZ$zonePolygone)
# print zoning ids
printZid(ZK$resZ$zonePolygone)

p=profvis({K=calNei(Z,map$krigData,map$krigSurfVoronoi,map$krigN,simplitol)})
p

################################################################
# step 3 - generate tree of possible corrections for small zones (remove, i.e. merge into englobing zone, r grow - if grow, isolated zone -> grows bigger but remains isolated from others (distIsoZ parameter) - non isolated zone -> joins closest zone with same labels
################################################################
p=profvis({criti=correctionTree(c(0.4,0.7),map,pErr=0.9,optiCrit=2,minSize=0.012,minSizeNG=1e-3,distIsoZ=0.075,simplitol=1e-3,LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=0,SAVE=TRUE,ONE=FALSE,ALL=TRUE)
})# argument SAVE and ALL will keep results at all branches in zf,zk,critere,mdist
p # 

zk=criti$zk

# here 1 small zone (#3), hence 2 levels (level 1 is initial zoning, level 2 has 2 branches, one for each corrected zoning-first branch=zone removal, second branch=zone junction, because zone 3 is not isolated and close to zone 5 (BOTH ZONES 3 and 5 have same lab)
plotZ(zk[[2]][[1]]$zonePolygone) # result of removal of zones 8 and 7
# or
K=zk[[2]][[1]]
Z=K$zonePolygone
dispZ(map$step,map$krigGrid,zonePolygone=Z,K=K,boundary=map$boundary,nbLvl=0,id=FALSE)
plotZ(zk[[2]][[2]]$zonePolygone) # result of junction of zones 3 and 5
# or
K=zk[[2]][[2]]
Z=K$zonePolygone
dispZ(map$step,map$krigGrid,zonePolygone=Z,K=K,boundary=map$boundary,nbLvl=0,id=FALSE)

# test distance calculation
crit=correctionTree(qProb=c(0.4,0.7),map,SAVE=T)
zk=crit$zk
K=zk[[1]][[1]]
Z=K$zonePolygone
mdist=crit$mdist
md=mdist[[1]][[1]]
sum(abs(md-DC0N)) #check
profvis({
for (i in 1:100){
resD=calDistance(typedist=1,map$krigData,K$listZonePoint,K$zoneN,map$krigSurfVoronoi,K$meanZone,pErr)}
})

sum(abs(resD$matDistanceCorr-DC0))

crit2 = calCrit(resD$matDistanceCorr,K$zoneNModif,optiCrit)
#critC = calCritC(crit2,resD$cost)
normMat=normDistMat(resD$matDistanceCorr,optiCrit)#dist(Z4,Z6)=0.59 


profvis({
qProb=c(0.1,0.2);criti=correctionTree(qProb,map,pErr=0.9,optiCrit=2,minSize=0.012,minSizeNG=1e-3,distIsoZ=0.075,simplitol=1e-3,LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=0,SAVE=TRUE,ONE=FALSE)
})

zk=criti$zk

sig2=rep(0,length(Z))
SI=sig2
for(k in 1:length(Z))
{
res=Sigmai2(k,K$listZonePoint,map$krigData,map$krigSurfVoronoi,K$meanZone)
sig2[k]=res$sigmai2
SI[k]=res$SI
}

qProb=c(0.1,0.2);criti=correctionTree(qProb,map)
res=searchNODcrit1(qProb,criti)

ZK=initialZoning(qProb=c(0.4,0.7),map)
Z=ZK$resZ$zonePolygone

iSmall=detectSmallZones(Z,minSize) # 2 small zones 7 and 8

################################################################
# step 4 - optimizing loops - step=0.075 and diffQ=0.2
################################################################
p=profvis({
ro1=loopQ1(map,disp=0,step=0.2) # 1 quantile vector
})
ro1
p

ro2=loopQ2(map,disp=0,step=0.2) # 2 quantile vector
ro2

ro3=loopQ3(map,disp=0,step=0.2) # 3 quantile vector
ro3


ro4=loopQ4(map,disp=0,step=0.2) # 4 quantile vector
ro4


# test optiRG
###############
qProb=c(0.1,0.2);ZK=initialZoning(qProb,map)
#criti=correctionTree(qProb,map) # pb vers la fin
Z=ZK$resZ$zonePolygone
K=ZK$resZ

plotZ(Z) #only zones
printZsurf(Z) # zones 4 3 9
# merge zones 3 and 2
detZoneClose(3,Z,K$zoneNModif,0.01)
kmi=optiRG(K,map,3,2,disp=1) #geometry intersection pb solved -> kmi=NULL
