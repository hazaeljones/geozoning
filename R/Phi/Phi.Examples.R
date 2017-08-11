
##################################################################################################################
################################### TEST for EXAMPLE documentation ###############################################
##################################################################################################################


load("mapTest.rda")
map = mapTest
plotM(map)


res2 = loopQ2(map,disp=0,step=0.075,QUIET=FALSE)
res3 = loopQ3(map,disp=0,step=0.075,QUIET=FALSE)
res4 = loopQ4(map,disp=0,step=0.075,QUIET=FALSE)


epsilon = 1
len = 11 # length(list_nested_quantile) , initialized to 6 to enter the loop
while(len >= 7){
  list_nested_quantile = list()
  res2_best = res2[which( (res2[1,1] - res2[,1])<=epsilon ), ]
  res3_best = res3[which( (res3[1,1] - res3[,1])<=epsilon ), ]
  res4_best = res4[which( (res4[1,1] - res4[,1])<=epsilon ), ]
  index = 0
  for(i in 1: nrow(res2_best)){
    q2 = as.matrix(res2_best[i,5:6])
    for(j in 1:nrow(res3_best)){
      q3 = as.matrix(res3_best[j, 5:7])
      if(sum(q2%in%q3) == 2){
        for(k in 1:nrow(res4_best)){
          q4 = as.matrix(res4_best[k, 5:8])
          if(sum(q3%in%q4) == 3){
            index = index+1
            list_nested_quantile[[index]] = list(q2,q3,q4)
          }
        }
      }
    }
  }
  epsilon = epsilon - 0.05
  len = length(list_nested_quantile)
}


criti = correctionTree(qProb = c(0.575,0.725),map = map,disp = 1)




# Points_Near_Boundary() ########################################################################################################
load("mapTest.rda")
map = mapTest
Points_Near_Boundary(map = map)



# Transition_Zone_Near_Boundary() ########################################################################################################
seed=2
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
ZK=initialZoning(qProb=c(0.45,0.85),map)
Z=ZK$resZ$zonePolygone # list of zones
lab = ZK$resZ$lab # label of zones
plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
# zone 6 is a transition zone that has commun boundary with the map
numZ = 6
Estimation = Transition_Zone_Near_Boundary(map = map, Z = Z, numZ = numZ)
# compute the cost
cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
print(cL$cost_Laplace)
print(cM$cost_Mean)
# zone 6 is a zone with gradient



# Transition_Zone_Far_Boundary() ########################################################################################################
seed=9
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
ZK=initialZoning(qProb=c(0.65,0.8),map)
Z=ZK$resZ$zonePolygone # list of zones
lab = ZK$resZ$lab # label of zones
plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
# zone 7 is a transition zone that is far from map boundary
numZ = 7
Estimation = Transition_Zone_Far_Boundary(map = map, Z = Z, numZ = numZ)
# compute the cost
cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
print(cL$cost_Laplace)
print(cM$cost_Mean)
# zone 7 is a zone with gradient


# Extreme_Zone() ########################################################################################################

seed=6
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
ZK=initialZoning(qProb=c(0.8),map)
Z=ZK$resZ$zonePolygone # list of zones
lab = ZK$resZ$lab # label of zones
plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
# zone 2 is a zone with maximum label
numZ = 2
Estimation = Extreme_Zone(map = map, Z = Z, numZ = numZ, label.is.min = FALSE)
# compute the cost
cL = Cost_By_Laplace(map = map, Z = Z, numZ = numZ, Estimation = Estimation)
cM = Cost_By_Mean(map = map, Z = Z, numZ = numZ)
print(cL$cost_Laplace)
print(cM$cost_Mean)
# zone 2 is homogeneous


# list_Zone_2_Neighbour()

seed=6
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
ZK=initialZoning(qProb=c(0.67,0.8),map)
Z=ZK$resZ$zonePolygone # list of zones
lab = ZK$resZ$lab # label of zones
plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
# zone 4 and 6 are transition zones and have exactly 2 neighbours with different labels.
list_Zone_2_Neighbours(Z = Z, lab = lab)





# new_krigGrid_for_visualisation() ########################################################################################################
seed=2
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
ZK=initialZoning(qProb=c(0.55,0.85),map)
Z=ZK$resZ$zonePolygone # list of zones
lab = ZK$resZ$lab # label of zones
plotM(map = map,Z = Z,lab = lab, byLab = FALSE)
# zone 6 is a transition zone that has commun boundary with the map
numZ = 6
Estimation = Transition_Zone_Near_Boundary(map = map, Z = Z, numZ = numZ)

result = new_krigGrid_for_visualisation(map = map, Z = Z, numZ = numZ, solution = Estimation)
new_krigGrid = result$new_krigGrid
new_data = result$new_data
quant1 = quantile(map$krigData@data$var1.pred,probs = 0.55)
quant2 = quantile(map$krigData@data$var1.pred,probs = 0.85)

# plot initial isocontours
plotM(map = map,Z = Z,lab = lab, byLab = TRUE)
listContours = contourBetween(map = map, krigGrid = map$krigGrid, q1 = quant1, q2 = quant2)
for (i in 1:length(listContours)){
  plot(listContours[[i]]$contour,add=TRUE,col = "red")
}
# plot modified isocontours
plotM(map = map,Z = Z,lab = lab, byLab = TRUE)
listContours = contourBetween(map = map, krigGrid = new_krigGrid, q1 = quant1, q2 = quant2)
for (i in 1:length(listContours)){
  plot(listContours[[i]]$contour,add=TRUE,col = "red")
}


#




# correctBoundaryMap ########################################################################################################
seed=1
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
criti = correctionTree(qProb = c(0.5), map = map)
Z = criti$zk[[1]][[1]]$zonePolygone
lab = criti$zk[[1]][[1]]$lab
plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
class(gIntersection(Z[[1]],Z[[2]])) [1]
class(gIntersection(Z[[1]],Z[[5]])) [1]
class(gIntersection(Z[[2]],Z[[3]])) [1]
class(gIntersection(Z[[2]],Z[[4]])) [1]
res = correctBoundaryMap(Zi = Z, map = map)
Z = res$Z
class(gIntersection(Z[[1]],Z[[2]])) [1]
class(gIntersection(Z[[1]],Z[[5]])) [1]
class(gIntersection(Z[[2]],Z[[3]])) [1]
class(gIntersection(Z[[2]],Z[[4]])) [1]


# smoothingZone ########################################################################################################

seed=1
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
criti = correctionTree(qProb = c(0.5), map = map)
Z = criti$zk[[1]][[1]]$zonePolygone
lab = criti$zk[[1]][[1]]$lab
# zones' correction
res = correctBoundaryMap(Zi = Z, map = map)
Z = res$Z
# map boundary after correction
boundary = Z[[1]]
for(i in 2:length(Z)){
  boundary = gUnion(boundary, Z[[i]])
}
# plot map
plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
# smoothing
zone = Z[[2]]
newZone = smoothingZone(z = zone, width = 0.05, boundary = boundary)
plot(zone)
X11()
plot(newZone)


# touch.border ########################################################################################################
seed=1
map = genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
criti = correctionTree(qProb = c(0.5), map = map)
Z = criti$zk[[1]][[1]]$zonePolygone
lab = criti$zk[[1]][[1]]$lab
# zones' correction
res = correctBoundaryMap(Zi = Z, map = map)
Z = res$Z
# map boundary after correction
boundary = Z[[1]]
for(i in 2:length(Z)){
  boundary = gUnion(boundary, Z[[i]])
}
# plot map
plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
# verification
for(i in 1:length(Z)){
  print(touch.border(z = Z[[i]], boundary = boundary))
}


# zone.extended ########################################################################################################
seed=1
map = genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
criti = correctionTree(qProb = c(0.5), map = map)
Z = criti$zk[[1]][[1]]$zonePolygone
lab = criti$zk[[1]][[1]]$lab
# zones' correction
res = correctBoundaryMap(Zi = Z, map = map)
Z = res$Z
# map boundary after correction
boundary = Z[[1]]
for(i in 2:length(Z)){
  boundary = gUnion(boundary, Z[[i]])
}
# plot map
plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
# extend zone
z = zone.extended(z = Z[[1]], boundary = boundary)
plot(z)
plot(Z[[1]],add=TRUE)


# cal.max.width.Zone ########################################################################################################
seed=1
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
criti = correctionTree(qProb = c(0.4,0.6), map = map)
Z = criti$zk[[2]][[1]]$zonePolygone
lab = criti$zk[[2]][[1]]$lab
# zones' correction
res = correctBoundaryMap(Zi = Z, map = map)
Z = res$Z
# map boundary after correction
boundary = Z[[1]]
for(i in 2:length(Z)){
  boundary = gUnion(boundary, Z[[i]])
}
# plot map
plotM(map = map, Z = Z, lab = lab, byLab = FALSE)
widthMax = cal.max.width.Zone2(z = Z[[3]], step = 0.001, widthMax = 0.05, boundary = boundary, erosion = TRUE)
zone = zone.extended(z = Z[[3]], boundary = boundary)
erosion1 = gBuffer(zone ,width = - (widthMax + 0.002) ,joinStyle="ROUND",capStyle = "ROUND")
erosion2 = gBuffer(zone ,width = - (widthMax - 0.002) ,joinStyle="ROUND",capStyle = "ROUND")
plot(erosion1)
X11()
plot(erosion2)



# smoothingMap ########################################################################################################
seed=1
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
criti = correctionTree(qProb = c(0.4,0.6), map = map)
Z = criti$zk[[2]][[1]]$zonePolygone
newZ = smoothingMap2(Z = Z, width = 0.05, map = map, disp = TRUE)
plotM(map = map, Z = Z, lab = 1:length(Z))
plotM(map = map, Z = newZ, lab = 1:length(Z), new_Window = TRUE)



load("yieldMap+Z.Rdata")
width = 0.01
Z1.1 = smoothingMap2(Z = Z1, width = width, map = map, disp = TRUE)
Z2.1 = smoothingMap2(Z = Z2, width = width, map = map, disp = TRUE)
Z3.1 = smoothingMap2(Z = Z3, width = width, map = map, disp = TRUE)
Z4.1 = smoothingMap2(Z = Z4, width = width, map = map, disp = TRUE)
Z5.1 = smoothingMap2(Z = Z5, width = width, map = map, disp = TRUE)

plotM(map = map, Z = Z1.1, lab = 1:length(Z1.1))
plotM(map = map, Z = Z2.1, lab = 1:length(Z2.1), new_Window = TRUE)
plotM(map = map, Z = Z3.1, lab = 1:length(Z3.1), new_Window = TRUE)
plotM(map = map, Z = Z4.1, lab = 1:length(Z4.1), new_Window = TRUE)
plotM(map = map, Z = Z5.1, lab = 1:length(Z5.1), new_Window = TRUE)

plotM(map = map, Z = Z5, lab = 1:length(Z5), new_Window = TRUE)



res = correctBoundaryMap(Zi = Z3, map = map)
Z = res$Z
z = Z[[7]]
boundary = Z[[1]]
for(i in 2:length(Z)){
  boundary = gUnion(boundary, Z[[i]])
}
zone = zone.extended(z = z, boundary = boundary)

z= zone
step = 0.001
widthMax = 0.02
erosion = FALSE
width = 0.02
cal.max.width.Zone2(z = z, step = 0.001, widthMax = width, boundary = boundary, erosion = FALSE)



dilatation1 = gBuffer(zone,width = widthExt,joinStyle="ROUND",capStyle = "ROUND")
erosion1 = gBuffer(dilatation1,width = -widthExt,joinStyle="ROUND",capStyle = "ROUND")




























