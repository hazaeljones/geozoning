rm(list=ls())
library(profvis)

source("srcZ.R") # source libraries and functions - params are in initParam.R

# simulate map - write data in data directory with day date
# parametres 
krig=2
Vpsill=c(10)
Vrange=c(0.25)
Vmean=25
Vnugget=0
#
map=genMap(DataObj=NULL,seed=30,krig=krig,Vpsill=Vpsill,Vrange=Vrange,Vnugget=Vnugget,Vmean=Vmean,disp=0)
os(map) # object size

# generate zoning
p=profvis({
ZK=initialZoning(qProb=c(0.4,0.7), map,disp=FALSE)
})
p # profiling results in browser

os(ZK) # object size
#
# get results
K=ZK$resZ
Z=K$zonePolygone
qProb=K$qProb
valRef=quantile(map$krigGrid,na.rm=TRUE,prob=qProb) #data values for quantiles
#
# test Fusion (zone removal = merging into englobing zone) - does not test zone size, so can be applied to any zone
testFusion(qProb,map,Z,K,3)
testFusion(qProb,map,Z,K,4)
# get englobing zone
detZoneEng(4,Z,K) # returns 2
detZoneEng(3,Z,K) # returns 2
detZoneEng(2,Z,K) # returns 0 = no englobing zone - fusion will return NULL
#
# test isolated zone Growing: 1- generate envelope into which future zone must fit
# 2- try quantiles corresponding to quantile sequence (LEQ and MAXP parameters in initParam file
valRef=quantile(map$krigGrid,na.rm=TRUE,prob=qProb) # just to see

testGrowM(3,Z,K,map,qProb) # impossible (no envelope) as zone 3 is not isolated
testGrowM(7,Z,K,map,qProb) # impossible - quantile sequence corresponding zone is not included into envelope
LEQ=50
testGrowM(7,Z,K,map,qProb,LEQ) # possible - opt quantile = 23.18728
# note: after growing, isolated zone is still isolated, and may be still too small - last pass in funcCritereCNO

#
#
# generate other zoning
ZK=initialZoning(qProb=c(0.1,0.4,0.7), map,disp=FALSE)
# get results
K=ZK$resZ
Z=K$zonePolygone
qProb=K$qProb
plotZ(Z,id=FALSE)
# test grow zone 9
LEQ=50
testGrowM(9,Z,K,map,qProb,LEQ)

# test zone junction (#4 and #6 - same label)
testModifNonIso(qProb,map,Z,K,4,simplitol)

# test isolated zone grow (zone8 gets just a little bigger)
LEQ=50
testGrowM(8,Z,K,map,qProb)
#
# tree of corrected zonings
minSize=0.023 # 2 small zones #4, #8
LASTPASS=FALSE
resC=correctionTree(qProb=c(0.1,0.4,0.7),map,disp=1,SAVE=T)
critList=resC$critList
# examine critList (sorted by number of effective quantiles)
# critList written in sortCrit
nq=getNq(critList)
#
zk=resC$resZ
crit=resC$criterion
md=resC$mdist
# no change with lastPass (no small zone in last tree level)
resPass=lastPass(map,qProb,resZ,crit,md,disp=F)
Z=resPass$listOfZonage
K=resPass$listOfKrig
crit=resPass$crit
md=resPass$mdist

# other tree of corrected zonings
minSize=0.025 # 3 small zones #4, #8, #9
LASTPASS=FALSE

p=profvis({
resC=correctionTree(qProb=c(0.1,0.4,0.7),map,disp=1,SAVE=T)
})
p #view profiling in browser
critList=resC$critList
# examine critList (sorted by number of effective quantiles)
# critList written in sortCrit
nq=getNq(critList)
#
ZK=initialZoning(qProb=c(0.1,0.9), map,disp=FALSE)
K=ZK$resZ
Z=K$zonePolygone
qProb=K$qProb
plotZ(Z)
LEQ=5
res=testGrowM(6,Z,K,map,qProb) # grows  a lot
#
res=testGrowM(3,Z,K,map,qProb)# grows a little
#
LEQ=50
res=testGrowM(1,Z,K,map,qProb)

map=genMap(DataObj=NULL,seed=80,disp=FALSE)
ZK=initialZoning(qProb=c(0.1,0.4,0.7), map,disp=FALSE) #lots of small zones
K=ZK$resZ
Z=K$zonePolygone
qProb=K$qProb
LEQ=5
res=testGrowM(7,Z,K,map,qProb)

qProb=c(0.1,0.3,0.9)
ZK=initialZoning(qProb=qProb, map,disp=FALSE) #seed=80
K=ZK$resZ
Z=K$zonePolygone
qProb=K$qProb
plotZ(Z)
res=testGrowM(5,Z,K,map,qProb)

