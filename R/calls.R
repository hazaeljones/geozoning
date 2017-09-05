rm(list=ls())
source("srcZ.R") # source libraries and functions - params are in initParam.R
#
################################################################
# step 1 - generate 2D map from simulation - default = 450 pts
################################################################
# kriging default 2000 pt grid
#seed=30
seed=80
map=genMap(DataObj=NULL,seed=seed,disp=FALSE)

# display 2D map
plotMap(map)
#
################################################################
# or build 2D map from real data (+kriging)
source("prep_real_data-vig.R") # example of call for Arnel vigor data (and frame) - generates map, normalizes boundary and minSize
################################################################
#
#
#UnzonageFromSimu(seed=80,qProb=c(0.4,0.7),display=TRUE) # example of call
# calls genMap to generate 2D map, then calls func1Zonage 
# generates zoning Z and saves it in R object data/Z.Rdata: polys,nb of zones, seed, raw data, result of genMap and calNei called by func1Zonage
# objects can be reloaded later on (to calculate criterion for instance)
# function returns void
#
#
################################################################
# step 2 - generates zoning Z from map, quantile vector (0.4 here)-no display
################################################################
ZK=initialZoning(qProb=c(0.4,0.7), map,disp=FALSE) # names(Z) "zonePolygone" "resCritereA"  "resDistanceA" "resKrig"
# initialZoning calls zoneGeneration, calNei, calDistance, calCrit
# zoneGeneration calls contourAuto
# contourAuto calls contourLines and extensionLigne
#
# plot zoning (7 zones in this case)
plotZ(ZK$resZ$zonePolygone) #only zones
# more detailed plot
valRef=quantile(map$krigGrid,na.rm=TRUE,prob=c(0.4,0.7))
dispZ(map$step,map$krigGrid,ZK$resZ$zonePolygone,boundary=map$boundary,nbLvl=0,id=FALSE)
 title(stepte(" valRef=[",toString(round(valRef,2)),"]   critere=",round(ZK$resCritereA,2),sep=""))
# print zoning labels
printLabZ(list(ZK$resZ))
# print zoning surfaces
printZsurf(ZK$resZ$zonePolygone)
# print zoning ids
printZid(ZK$resZ$zonePolygone)
#
################################################################
# step 3 - generate tree of possible corrections for small zones (remove, i.e. merge into englobing zone, r grow - if grow, isolated zone -> grows bigger but remains isolated from others (distIsoZ parameter) - non isolated zone -> joins closest zone with same labels
################################################################
resC=correctionTree(qProb=c(0.4,0.7),map,disp=1,SAVE=T) # argument SAVE will keep results at all branches
zk=resC$zk
#
# here 1 small zone (#3), hence 2 levels (level 1 is initial zoning, level 2 has 2 branches, one for each corrected zoning-first branch=zone removal, second branch=zone junction, because zone 3 is not isolated and close to zone 5 (BOTH ZONES 3 and 5 have same lab)
plotZ(zk[[2]][[1]]$zonePolygone) # result of removal of zone 3
# or
dispZ(map$step,map$krigMatTest,zonePolygone=zk[[2]][[1]]$zonePolygone,boundary=map$boundary,nbLvl=0,id=FALSE)
plotZ(zk[[2]][[2]]) # result of junction of zones 3 and 5
# or
dispZ(map$step,map$krigMatTest,zonePolygone=zk[[2]][[2]]$zonePolygone,boundary=map$boundary,nbLvl=0,id=FALSE)

################################################################
# step 4 - optimizing loops - step=0.1 and diffQ=0.2
################################################################
# results for seed=30
ro1=loopQ1(map,disp=0) 1 quantile vector
ro1
#  crit  iq nq
#9 5.76 0.9  1
#1 5.20 0.1  1
#2 4.75 0.2  1
#3 4.35 0.3  1
#5 4.33 0.5  1
#6 3.86 0.6  1
#7 3.63 0.7  1
#8 3.52 0.8  1
#4 3.07 0.4  1

ro2=loopQ2(map,disp=0)# 2 quantile vector
ro2
#       crit  iq  kq nq
#7  6.269950 0.1 0.9  2
#13 5.746957 0.2 0.9  2
#10 5.261905 0.2 0.6  2
#18 5.059242 0.3 0.9  2
#3  4.992160 0.1 0.5  2
#4  4.823662 0.1 0.6  2
#11 4.609160 0.2 0.7  2
#9  4.589413 0.2 0.5  2
#15 4.518982 0.3 0.6  2
#16 4.435907 0.3 0.7  2
#5  4.348327 0.1 0.7  2
#25 4.328800 0.5 0.9  2
#12 4.180170 0.2 0.8  2
#1  4.148492 0.1 0.3  2
#6  4.045923 0.1 0.8  2
#17 4.011524 0.3 0.8  2
#27 3.862833 0.6 0.9  2
#28 3.541946 0.7 0.9  2
#14 3.465107 0.3 0.5  2
#21 3.458674 0.4 0.8  2
#22 3.364550 0.4 0.9  2
#8  3.066124 0.2 0.4  2
#2  3.066124 0.1 0.4  2
#24 2.807293 0.5 0.8  2
#20 2.630430 0.4 0.7  2
#23 2.135092 0.5 0.7  2
#26 1.720284 0.6 0.8  2
#19 1.518963 0.4 0.6  2
#

ro3=loopQ3(map,disp=0) # 3 quantile vector
ro3
#      crit  iq  jq  kq nq
#24 5.261905 0.2 0.6 0.9  3
#12 5.173814 0.1 0.5 0.9  3
#5  5.024095 0.1 0.3 0.9  3
#22 4.906687 0.2 0.5 0.9  3
#14 4.823662 0.1 0.6 0.9  3
#30 4.518982 0.3 0.6 0.9  3
#2  4.518982 0.1 0.3 0.6  3
#3  4.435907 0.1 0.3 0.7  3
#25 4.397772 0.2 0.7 0.9  3
#31 4.201150 0.3 0.7 0.9  3
#15 4.189143 0.1 0.7 0.9  3
#4  4.011524 0.1 0.3 0.8  3
#1  3.539323 0.1 0.3 0.5  3
#28 3.465107 0.3 0.5 0.9  3
#18 3.458674 0.2 0.4 0.8  3
#8  3.458674 0.1 0.4 0.8  3
#19 3.364550 0.2 0.4 0.9  3
#9  3.364550 0.1 0.4 0.9  3
#27 2.807293 0.3 0.5 0.8  3
#21 2.807293 0.2 0.5 0.8  3
#11 2.807293 0.1 0.5 0.8  3
#34 2.630430 0.4 0.7 0.9  3
#17 2.630430 0.2 0.4 0.7  3
#7  2.630430 0.1 0.4 0.7  3
#35 2.135092 0.5 0.7 0.9  3
#26 2.135092 0.3 0.5 0.7  3
#20 2.135092 0.2 0.5 0.7  3
#10 2.135092 0.1 0.5 0.7  3
#29 1.720284 0.3 0.6 0.8  3
#23 1.720284 0.2 0.6 0.8  3
#13 1.720284 0.1 0.6 0.8  3
#33 1.518963 0.4 0.6 0.9  3
#32 1.518963 0.4 0.6 0.8  3
#16 1.518963 0.2 0.4 0.6  3
#6  1.518963 0.1 0.4 0.6  3
#


ro4=loopQ4(map,disp=0)# 4 quantile vector
ro4
#      crit  iq  jq  kq  pq nq
#5  4.518982 0.1 0.3 0.6 0.9  4
#6  4.201150 0.1 0.3 0.7 0.9  4
#3  3.539323 0.1 0.3 0.5 0.9  4
#2  2.807293 0.1 0.3 0.5 0.8  4
#13 2.630430 0.2 0.4 0.7 0.9  4
#9  2.630430 0.1 0.4 0.7 0.9  4
#15 2.135092 0.3 0.5 0.7 0.9  4
#14 2.135092 0.2 0.5 0.7 0.9  4
#10 2.135092 0.1 0.5 0.7 0.9  4
#1  2.135092 0.1 0.3 0.5 0.7  4
#4  1.720284 0.1 0.3 0.6 0.8  4
#12 1.518963 0.2 0.4 0.6 0.9  4
#11 1.518963 0.2 0.4 0.6 0.8  4
#8  1.518963 0.1 0.4 0.6 0.9  4
#7  1.518963 0.1 0.4 0.6 0.8  4

# results for seed=80
ro1=loopQ1(map,disp=0) 1 quantile vector
ro1
#  crit  iq nq
#8 4.31 0.8  1
#4 4.29 0.4  1
#9 3.95 0.9  1
#3 3.88 0.3  1
#1 3.37 0.1  1
#7 3.30 0.7  1
#2 2.90 0.2  1
#6 2.72 0.6  1
#5 2.41 0.5  1
#
ro2=loopQ2(map,disp=0)# 2 quantile vector
ro2
#        crit  iq  kq nq
#17 3.7860396 0.3 0.8  2
#7  3.5796879 0.1 0.9  2
#18 3.5462219 0.3 0.9  2
#6  3.4483083 0.1 0.8  2
#22 3.2562323 0.4 0.9  2
#5  3.2211852 0.1 0.7  2
#16 3.0683691 0.3 0.7  2
#21 3.0579565 0.4 0.8  2
#13 2.8957611 0.2 0.9  2
#12 2.8957611 0.2 0.8  2
#11 2.7633934 0.2 0.7  2
#20 2.7092086 0.4 0.7  2
#4  2.2623312 0.1 0.6  2
#25 2.0245842 0.5 0.9  2
#3  1.8979156 0.1 0.5  2
#14 1.6274908 0.3 0.5  2
#27 1.5442379 0.6 0.9  2
#10 1.4767690 0.2 0.6  2
#9  1.4228514 0.2 0.5  2
#2  1.2957988 0.1 0.4  2
#26 1.1384033 0.6 0.8  2
#24 1.1223975 0.5 0.8  2
#28 1.0301268 0.7 0.9  2
#8  1.0143685 0.2 0.4  2
#15 0.9508294 0.3 0.6  2
#1  0.8148798 0.1 0.3  2
#23 0.6983560 0.5 0.7  2
#19 0.5729509 0.4 0.6  2
#
ro3=loopQ3(map,disp=0) # 3 quantile vector
ro3
#        crit  iq  jq  kq nq
#12 1.8979156 0.1 0.5 0.9  3
#28 1.6274908 0.3 0.5 0.9  3
#14 1.5442379 0.1 0.6 0.9  3
#24 1.4767690 0.2 0.6 0.9  3
#22 1.4228514 0.2 0.5 0.9  3
#9  1.2957988 0.1 0.4 0.9  3
#8  1.2957988 0.1 0.4 0.8  3
#7  1.2957988 0.1 0.4 0.7  3
#23 1.1384033 0.2 0.6 0.8  3
#13 1.1384033 0.1 0.6 0.8  3
#27 1.1223975 0.3 0.5 0.8  3
#21 1.1223975 0.2 0.5 0.8  3
#11 1.1223975 0.1 0.5 0.8  3
#34 1.0301268 0.4 0.7 0.9  3
#31 1.0301268 0.3 0.7 0.9  3
#25 1.0301268 0.2 0.7 0.9  3
#15 1.0301268 0.1 0.7 0.9  3
#19 1.0143685 0.2 0.4 0.9  3
#18 1.0143685 0.2 0.4 0.8  3
#17 1.0143685 0.2 0.4 0.7  3
#30 0.9508294 0.3 0.6 0.9  3
#29 0.9508294 0.3 0.6 0.8  3
#5  0.8148798 0.1 0.3 0.9  3
#4  0.8148798 0.1 0.3 0.8  3
#3  0.8148798 0.1 0.3 0.7  3
#2  0.8148798 0.1 0.3 0.6  3
#1  0.8148798 0.1 0.3 0.5  3
#35 0.6983560 0.5 0.7 0.9  3
#26 0.6983560 0.3 0.5 0.7  3
#20 0.6983560 0.2 0.5 0.7  3
#10 0.6983560 0.1 0.5 0.7  3
#33 0.5729509 0.4 0.6 0.9  3
#32 0.5729509 0.4 0.6 0.8  3
#16 0.5729509 0.2 0.4 0.6  3
#6  0.5729509 0.1 0.4 0.6  3
#

ro4=loopQ4(map,disp=0) # 4 quantile vector
ro4
#        crit  iq  jq  kq  pq nq
#9  1.0301268 0.1 0.4 0.7 0.9  4
#13 1.0143685 0.2 0.4 0.7 0.9  4
#6  0.8148798 0.1 0.3 0.7 0.9  4
#5  0.8148798 0.1 0.3 0.6 0.9  4
#4  0.8148798 0.1 0.3 0.6 0.8  4
#3  0.8148798 0.1 0.3 0.5 0.9  4
#2  0.8148798 0.1 0.3 0.5 0.8  4
#15 0.6983560 0.3 0.5 0.7 0.9  4
#14 0.6983560 0.2 0.5 0.7 0.9  4
#10 0.6983560 0.1 0.5 0.7 0.9  4
#1  0.6983560 0.1 0.3 0.5 0.7  4
#12 0.5729509 0.2 0.4 0.6 0.9  4
#11 0.5729509 0.2 0.4 0.6 0.8  4
#8  0.5729509 0.1 0.4 0.6 0.9  4
#7  0.5729509 0.1 0.4 0.6 0.8  4
#
seed=2
map=genMap(DataObj=NULL,seed=seed,disp=FALSE)
qProb= c(0.05,0.2,0.35,0.5)
resC=correctionTree(qProb,map,disp=1,SAVE=T)
zk=resC$zk
zk0=zk[[1]][[1]]
Z0=zk0$zonePolygone
K0=zk0
printZsurf(Z0) #16,5,6,14=small zones
testGrowM(16,Z0,K0,map,qProb) #grow isolated zone 16 - no space
testGrowM(16,Z0,K0,map,qProb,distIsoZ=1e-3) # works
testGrowM(5,Z0,K0,map,qProb) #grow isolated zone 5 - works
testGrowM(14,Z0,K0,map,qProb) #grow isolated zone 14 - no space

zg=optiRG(Z0,K0,map,6,3,simplitol,disp=0)#grow non isolated zone 6 close to zone 3-works

Z=zk[[3]][[1]]$zonePolygone
K=zk[[3]][[1]]
zg=optiRG(Z,K,map,5,3,simplitol,disp=0)