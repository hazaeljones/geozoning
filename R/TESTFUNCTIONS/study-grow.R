source("prep_simu.R") #seed 89

# small zones (isolated or not)
qProb=c(0.1,0.5,0.9)
crit2=correctionTree(qProb,map,pErr=0.9,optiCrit=2,minSize=0.012,minSizeNG=1e-3,distIsoZ=0.075,simplitol=1e-3,LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=1,SAVE=TRUE,ONE=FALSE)
#[1] "12 zones, 5  small zones:"
#[1] "12,10,3,8,11"


# last level
calczf(qProb,map,optiCrit,minSize,minSizeNG,disp=0,pdf1=NULL,FULL=FALSE,data=NULL)

zk=crit2$zk
K=zk[[1]][[1]]
Z=K$zonePolygone
# 2 isolated zones 8 10 

testGrowM(8,Z,K,map,qProb) # grow isolated zone 8
testGrowM(2,Z,K,map,qProb) # no englobing zone
testGrowM(10,Z,K,map,qProb) # grow isolated zone 10
testModifNonIso(qProb,map,Z,K,11,simplitol) # join zone 11 and zone 12
testModifNonIso(qProb,map,Z,K,12,simplitol) # idem


