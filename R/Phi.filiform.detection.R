# this script is about detecting filiform zone

##############################################################################################################
############### DÉTECTION DES ZONES FILIFORMES PAR LE RAPPORT PÉRIMÈTRE/AIRE #################################
##############################################################################################################

load("yieldMap+Z.Rdata")
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

X = correctBoundaryMap(Z3,boundary)

Z = X[[1]]
plotZ(Z)



z = Z[[5]]

plot(z)
width = 0.02

z2 = gBuffer(z,width = -width)
#z2 = gBuffer(z2,width = width)
#z2 = gBuffer(z2,width = -width)
#z2 = gBuffer(z2,width = width)


plot(gDifference(z,z2),add=TRUE,col="yellow")
plot(z2,add=TRUE,col="green")


diff = gDifference(z,z2)

listPolygon = diff@polygons[[1]]@Polygons
listArea = c()

for (i in 1:length(listPolygon)){
  listArea = c(listArea, listPolygon[[i]]@area)
}


# classement des aires
table = matrix(rep(0,2*length(listPolygon)), ncol = 2)
newlistArea = sort(listArea,decreasing = TRUE)
for(i in 1:length(newlistArea)){
  table[i,1] = newlistArea[i]
  table[i,2] = which(listArea == newlistArea[i])
}

# classement des (aires / périmètres)
listRapport = c()
for(i in 1:length(listPolygon)){
  p = listPolygon[[i]]
  ps = Polygons(list(p), ID = "p")
  sp = SpatialPolygons(list(ps))
  perimeter = gLength(sp)
  rapport = perimeter/listArea[i]
  listRapport = c(listRapport, rapport)
}

Z = list()

table = matrix(rep(0,3*length(listPolygon)), ncol = 3)
newListRapport = sort(listRapport,decreasing = TRUE)
for(i in 1:length(newListRapport)){
  table[i,1] = which(listRapport == newListRapport[i])
  table[i,2] = newListRapport[i]
  table[i,3] = listArea[which(listRapport == newListRapport[i])]

  if (i<=10){
    j = which(listRapport == newListRapport[i])
    p = listPolygon[[j]]
    ps = Polygons(list(p), ID = "p")
    Z = c(Z, SpatialPolygons(list(ps)))

  }
}
plotZ(Z)



##############################################################################################################
######################## DÉTECTION ZONES FILIFORME  ##########################################################
##############################################################################################################

seed = 2
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
ZK=initialZoning(qProb=c(0.5, 0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
plotZ(Z)

x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

X = correctBoundaryMap(Z,boundary)
Z = X[[1]]
plotZ(Z)

# les zones sur les quelles on va travailler: Union(Z[[2]],Z[[4]],Z[[5]],Z[[9]]) et Z[[3]], Z[[10]]
z3 = Z[[3]]
z10 = Z[[10]]


z2 = gUnion(gUnion(gUnion(Z[[2]],Z[[4]]), Z[[5]]),Z[[9]])

plot(z2)
width = 0.03

z22 = gBuffer(z2,width = -width)
z22 = gBuffer(z22,width = width)
#z2 = gBuffer(z2,width = -width)
#z2 = gBuffer(z2,width = width)


plot(gDifference(z2,z22),add=TRUE,col="yellow")





########## LES OPTIONS DE gBuffer #####################################################################

width = 0.030
z22 = gBuffer(z2,width = -width, capStyle = "ROUND",joinStyle = "BEVEL",mitreLimit = 4)
z22 = gBuffer(z22,width = width, capStyle = "ROUND",joinStyle = "BEVEL",mitreLimit = 4)
plot(z2)
plot(gDifference(z2,z22),add=TRUE,col="yellow")





##############################################################################################################
######################## APPROCHE QUANTILES EMBOITEES  #######################################################
##############################################################################################################


seed = 2

map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)

ZK=initialZoning(qProb=c(0.4,0.55, 0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))

X = correctBoundaryMap(Z,boundary)
Z = X[[1]]
plotZ(Z)
################################### VISUALISATION QUANTILE EMBOITÉE ##########################################

# fonction pour le plot Z ave couleur et label
plotZ.color = function(Z,lab,boundary) {
  plot(boundary)
  nblab = length(levels(as.factor(lab)))
  for (i in 1:length(Z)){
    plot(Z[[i]],add= TRUE, col = terrain.colors(nblab, alpha = 1)[lab[i]])
  }
  for (i in 1:length(Z)){
    text(coordinates(Z[[i]]), labels = i)
  }
}

seed = 2

map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)


ZK=initialZoning(qProb=c(0.5, 0.7),map,pErr,simplitol,optiCrit,disp=0)
Z=ZK$resZ$zonePolygone # zone


ZK=initialZoning(qProb=c(0.3,0.5, 0.7),map,pErr,simplitol,optiCrit,disp=0)
Z=ZK$resZ$zonePolygone # zone


lab = ZK$resZ$lab
x = map$boundary$x
y = map$boundary$y
z = cbind(x,y)
p = Polygon(z)
ps = Polygons(list(p),ID = "p")
boundary = SpatialPolygons(list(ps))



plotZ.color(Z,lab = lab,boundary = boundary)
















