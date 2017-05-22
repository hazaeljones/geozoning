#install.packages("rgeos")
#library(rgeos)
#install.packages("sp")
#library(sp)

##############################################################
# Exemple gBuffer
##############################################################


?gBuffer
?sp
?readWKT

# ?????????????????????? QUESTIONS PRELIMINAIRES ?????????????

# Si label(situé au centroid de la zone) se situe dans l'autre zone, on fait quoi?
# Morphologie : comment on calcule la moyenne de la dilatation et erosion
# Quel est le critère pour choisir le paramètre "width" dans "gBuffer" pour faire la morphologie
# Comment détecter les régions d'un zone où on veut effectuer la morphologie locale

rm(list=ls())
source("srcZ.R") # source libraries and functions - params are in initParam.R



##############################################################
# TEST TEST TEST
# jusqu'à la ligne 165, quelque test des fonctions de "rgeos"
##############################################################

seed=81
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
plotMap(map)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone

# plot of dilatation and erosion
plotZ(Z) # ????????????? Si label(situé au centroid de la zone) se situe dans l'autre zone, on fait quoi?
# dilatation
buffer1 = gBuffer(Z[[6]],width = 0.02,joinStyle="ROUND",capStyle = "ROUND")
plot(buffer1, add=TRUE)
# erosion
buffer2 = gBuffer(Z[[6]],width = -0.02, joinStyle = "ROUND",capStyle = "ROUND")
plot(buffer2,add=TRUE)


##############################################################
# search the centroid of polygon
##############################################################
center = gCentroid(Z[[6]])
center1 = gCentroid(buffer1)
center2 = gCentroid(buffer2)

plot(center, col = "blue",add= TRUE)
plot(center1, col = "red", add = TRUE)
plot(center2, col = "black",add = TRUE)


##############################################################
# example on a modified map
##############################################################
#  generate tree of possible corrections for small zones
criti=correctionTree(c(0.4,0.7),map,pErr=0.9,optiCrit=2,minSize=0.012,minSizeNG=1e-3,
                     distIsoZ=0.075,simplitol=1e-3,LEQ=5,MAXP=0.1,LASTPASS=TRUE,disp=0,
                     SAVE=TRUE,ONE=FALSE)

zk=criti$zk

# here 1 small zone (#3), hence 2 levels (level 1 is initial zoning, level 2 has 2 branches, one for each corrected zoning-first branch=zone removal, second branch=zone junction, because zone 3 is not isolated and close to zone 5 (BOTH ZONES 3 and 5 have same lab)
plotZ(zk[[2]][[1]]$zonePolygone) # result of removal of zone 3
# or
plotZ(zk[[2]][[2]]$zonePolygone) # result of junction of zones 3 and 5



##############################################################
# transform an "sp" object to "spDataFrame"
##############################################################
plotZ(Z)
plot(buffer1,add=TRUE)
plot(buffer2,add=TRUE)

# utilisation de "SpatialPolygonsDataFrame" semble non convenable

# on essaie le package "broom" : Convert Statistical Analysis Objects into Tidy Data Frames
# ou la fonction "fortify" du package ggplot2
zone6 = zk[[2]][[1]]$zonePolygone [[6]]

zone6.Data.Frame = fortify(zone6) [,1:2]
buffer1.Data.Frame = fortify(buffer1) [,1:2]
buffer2.Data.Frame = fortify(buffer2) [,1:2]

zone6.Data.Frame = tidy(zone6)[,1:2]
buffer1.Data.Frame = tidy(buffer1) [,1:2]
buffer2.Data.Frame = tidy(buffer2) [,1:2]

zone6.Data.Frame = geom(zone6)[,5:6]
buffer1.Data.Frame = geom(buffer1) [,5:6]
buffer2.Data.Frame = geom(buffer2) [,5:6]

center.Data.Frame = as.data.frame(center) # tidy doesn't work with spatial point object, use as.data.frame instead

plot(zone6.Data.Frame,type = "l", xlim= c(0.5,1), ylim = c(0,0.4))
lines(buffer1.Data.Frame)
lines(buffer2.Data.Frame)
points(center.Data.Frame)







##############################################################
# use gBuffer : dilatation -> erosion -> dilatation -> erosion
# in order to see the effect on the smoothing
##############################################################
plot(zone6)
plot(buffer1,add=TRUE)
plot(buffer2,add=TRUE)

# initialize the param "width" in "gBuffer"
width = 0.022
# dilatation -> erosion -> dilatation -> erosion
dilatation1 = gBuffer(Z[[6]],width = width,joinStyle="ROUND",capStyle = "ROUND")
erosion1 = gBuffer(dilatation1,width = -width,joinStyle="ROUND",capStyle = "ROUND")
dilatation2 = gBuffer(erosion1,width = width,joinStyle="ROUND",capStyle = "ROUND")
erosion2 = gBuffer(dilatation2,width = -width,joinStyle="ROUND",capStyle = "ROUND")

plot(zone6, col = "yellow")
plot(erosion1,add=TRUE)
plot(erosion2,add=TRUE)

# erosion -> dilatation -> erosion -> dilatation
erosion1 = gBuffer(Z[[6]],width = -width,joinStyle="ROUND",capStyle = "ROUND")
dilatation1 = gBuffer(erosion1,width = width,joinStyle="ROUND",capStyle = "ROUND")
erosion2 = gBuffer(dilatation1,width = -width,joinStyle="ROUND",capStyle = "ROUND")
dilatation2 = gBuffer(erosion2,width = width,joinStyle="ROUND",capStyle = "ROUND")

plot(zone6, col = "yellow")
plot(dilatation1,add=TRUE)
plot(dilatation2,add=TRUE)



##############################################################
# Calculer la distance entre un point et un polygon
##############################################################
# test
point1 = readWKT("POINT(0 0)")
polygon1 = readWKT("POLYGON((1 1,2 2,3 1,1 1))")
gDistance(point1, polygon1)
##############################################################
# chercher l'intersection entre une droite et un polygone
##############################################################
# test
droite1 = readWKT("LINESTRING(0 0,4 4)")
polygon1 = readWKT("POLYGON((2 2,2 3,3 3,3 2,2 2))")
intersection1 = gIntersection(droite1, polygon1)
plot(droite1)
plot(polygon1,add=TRUE)
intersection.Data.Frame = geom(intersection1)[,4:5]



###################################################################################################
############################ C'EST PARTI ##########################################################
# toutes les approches ci-dessus ont pour but de trouver le milieu entre la dilatation et l'érosion
###################################################################################################
seed=81
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)


ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone



##############################################################
### PREMIERE APPROCHE ########################################
# manipuler l'intersection et la différence entre les polygones
##############################################################
plot(Z[[6]])
width = 0.03

# lissage extérieur
dilatation = gBuffer(Z[[6]],width = width,joinStyle="ROUND",capStyle = "ROUND")
smooth.ext = gBuffer(dilatation,width = -width,joinStyle="ROUND",capStyle = "ROUND")
# lissage intérieur
erosion = gBuffer(Z[[6]],width = -width,joinStyle="ROUND",capStyle = "ROUND")
smooth.int = gBuffer(erosion,width = width,joinStyle="ROUND",capStyle = "ROUND")

plot(smooth.ext)
plot(smooth.int,add=TRUE)

# difference
difference = gDifference(smooth.ext,smooth.int)
plot(difference,add=TRUE,col="yellow")

# intersection
intersection = gIntersection(smooth.ext,difference)
plot(intersection, add=TRUE,col="red")

# DIFFICULTE DE SOUSTRAIRE ET AJOUTER DES REGION QU'ON VEUT MODIFIER
# ????????????????????????? PISTE A CONTINUER ??????????????????????????


########################################################################################
### DEUXIEME APPROCHE ##################################################################
# Calculer le polygone moyenne entre la dilatation et la erosion par rapport au centroid
########################################################################################

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
zone6 = Z[[6]]
width = 0.02

buffer1 = gBuffer(zone6,width = width,joinStyle="ROUND",capStyle = "ROUND")
buffer2 = gBuffer(zone6,width = -width,joinStyle="ROUND",capStyle = "ROUND")

buffer1.Data.Frame = geom(buffer1)[,5:6]
buffer2.Data.Frame = geom(buffer2)[,5:6]

zone6.modif = matrix(data = rep(0,2*nrow(buffer1.Data.Frame)),ncol = 2)

centroid = gCentroid(zone6)
centroid.Data.Frame = geom(centroid)[,2:3]

for (i in 1:nrow(buffer1.Data.Frame)){
  droit = readWKT(paste("LINESTRING(", toString(centroid.Data.Frame[1]), " ", toString(centroid.Data.Frame[2]), ",",
                        toString(buffer1.Data.Frame[i,1]), " ",toString(buffer1.Data.Frame[i,2]), ")"))
  intersection = gIntersection(droit, gBoundary(buffer2))
  middle.point = 1/2 * ( geom(intersection)[nrow(geom(intersection)),2:3] + buffer1.Data.Frame[i,])
  zone6.modif[i,] = middle.point
}
plot(buffer1.Data.Frame,type="l",col="blue")
lines(buffer2.Data.Frame,col="red")
lines(zone6.modif,col="green")
points(centroid.Data.Frame)

# DIFFICULTE QUAND LA DROITE COUPE 2 FOIS LA FRONTIERE DU POLYGONE
# ????????????????????????? PISTE A CONTINUER ??????????????????????????


########################################################################################
### TROISIEME APPROCHE #################################################################
# Calculer le polygone moyenne entre la (dilatation-erosion) et la (erosion-dialatation
# ce calcul est basé sur la projection des points de l'éroson sur le polygone  dilatation
########################################################################################
zone6 = Z[[6]]
width = 0.02

# lissage extérieur
dilatation = gBuffer(zone6,width = width,joinStyle="ROUND",capStyle = "ROUND")
smooth.ext = gBuffer(dilatation,width = -width,joinStyle="ROUND",capStyle = "ROUND")
# lissage intérieur
erosion = gBuffer(zone6,width = -width,joinStyle="ROUND",capStyle = "ROUND")
smooth.int = gBuffer(erosion,width = width,joinStyle="ROUND",capStyle = "ROUND")

plot(zone6)
plot(smooth.ext,add=TRUE)
plot(smooth.int,add=TRUE)

smooth.int.D.F = geom(smooth.int)[,5:6]
smooth.ext.D.F = geom(smooth.ext)[,5:6]
boundary.ext = gBoundary(smooth.ext)

zone6.modif = matrix(data = rep(0,2*nrow(smooth.int.D.F)),ncol = 2)

for (i in 1:nrow(smooth.int.D.F)){
  point = readWKT(paste("POINT(",toString(smooth.int.D.F[i,1])," ", toString(smooth.int.D.F[i,2]),")" ))
  nearest.Point = geom(gNearestPoints(point, boundary.ext) [2,] ) [2:3]
  middle.point = 1/2 * (geom(point)[2:3]+nearest.Point)
  zone6.modif[i,] = middle.point
}

plot(smooth.ext.D.F,col="blue",type="l")
lines(zone6.modif,type = "l")
lines(smooth.int.D.F,col = "red")


# DIFFICULTE QUAND ON A DES RÉGIONS TRÈS CONCAVES
# ????????????????????????? PISTE A CONTINUER ??????????????????????????



########################################################################################
### QUATRIEME APPROCHE #################################################################
# Calculer le polygone moyenne entre la dilatation et la erosion
# ce calcul est basé sur la projection des points de l'éroson sur le polygone  dilatation
########################################################################################

zone6 = Z[[6]]
width = 0.02


dilatation = gBuffer(zone6,width = width,joinStyle="ROUND",capStyle = "ROUND")
erosion = gBuffer(zone6,width = -width,joinStyle="ROUND",capStyle = "ROUND")


plot(dilatation)
plot(erosion,add=TRUE)

zone6.D.F = geom(zone6)[,5:6]
dilatation.D.F = geom(dilatation)[,5:6]
erosion.D.F = geom(erosion)[,5:6]
boundary.ext = gBoundary(dilatation)

zone6.modif = matrix(data = rep(0,2*nrow(erosion.D.F)),ncol = 2)

for (i in 1:nrow(erosion.D.F)){
  point = readWKT(paste("POINT(",toString(erosion.D.F[i,1])," ", toString(erosion.D.F[i,2]),")" ))
  nearest.Point = geom(gNearestPoints(point, boundary.ext) [2,] ) [2:3]
  middle.point = 1/2 * (geom(point)[2:3]+nearest.Point)
  zone6.modif[i,] = middle.point
}

plot(dilatation.D.F,col="blue",type="l")
lines(zone6.modif,type = "l")
lines(erosion.D.F,col = "red")
lines(zone6.D.F,col="green")


########################################################################################
### CINQUIEME APPROCHE #################################################################
# dilatation-erosion-erosion-dilatation OU erosion-dilatation-dilatation-erosion
########################################################################################

zone = Z[[1]]
width = 0.03
capStyle = "ROUND"
joinStyle = "ROUND"
mitreLimit = 1
nb.Iter = 5

dilatation = gBuffer(zone,width = width,joinStyle=joinStyle, capStyle = capStyle,mitreLimit = mitreLimit)
erosion = gBuffer(dilatation,width = -width,joinStyle=joinStyle, capStyle = capStyle, mitreLimit = mitreLimit)
erosion = gBuffer(erosion,width = -width,joinStyle=joinStyle, capStyle = capStyle,  mitreLimit = mitreLimit)
dilatation = gBuffer(erosion,width = width,joinStyle=joinStyle, capStyle = capStyle, mitreLimit = mitreLimit)

for (i in 2:nb.Iter){
  dilatation = gBuffer(dilatation,width = width,joinStyle=joinStyle, capStyle = capStyle,mitreLimit = mitreLimit)
  erosion = gBuffer(dilatation,width = -width,joinStyle=joinStyle, capStyle = capStyle, mitreLimit = mitreLimit)
  erosion = gBuffer(erosion,width = -width,joinStyle=joinStyle, capStyle = capStyle,  mitreLimit = mitreLimit)
  dilatation = gBuffer(erosion,width = width,joinStyle=joinStyle, capStyle = capStyle, mitreLimit = mitreLimit)
}
plot(zone)
plot(dilatation,add=TRUE)
# La qualité de lissage est superbe
# deplus, on n'a qu'à regler la parammètre "width" pour calculer la critère
# méthode simple et efficace.

# ???????? PROBLEME : LISSAGE est fait sur les frontières de la carte, or on veut que ce soit inchangé


zone.df = geom(zone)[,5:6]
dilatation.df = geom(dilatation)[,5:6]
plot(zone.df,type = "l")
lines(dilatation.df,col = "red")

plot(zone)
plot(dilatation,add= TRUE)

difference = gDifference(zone,dilatation)
plot(difference,col="yellow")

difference.df = geom(difference)[,5:6]


plot(difference.df[1:20,],type="l")

which(difference.df[,1] == 0)
which(difference.df[,1] == 1)
which(difference.df[,2] == 0)
which(difference.df[,2] == 1)




newZ = smoothingZone(zone, 0.02)
plot(newZ)


########################################################################################
### TESTER LES MODIFICATION SUR 2 ZONES VOISINNES ######################################
########################################################################################

zone1 = Z[[2]]
zone2 = Z[[3]]
intersection=gIntersection(zone1,zone2)
plot(zone1)
plot(zone2,add=TRUE)
plot(intersection,col="red",add=TRUE)

width = 0.03

newZ1 = smoothingZone(zone1, width,order =1)
newZ2 = smoothingZone(zone2, width,order =0)
newIntersection = gIntersection(newZ1,newZ2)

plot(newZ2)
plot(newZ1,add=TRUE)
plot(newIntersection,add=TRUE,col="yellow")
########################################################################################
### TESTER LES MODIFICATION SUR 1 ZONE AVEC L'ORDRE DIFFERENT ##########################
########################################################################################
width = 0.03
newZ1.1 = smoothingZone(zone1, width, order=1)
newZ1.2 = smoothingZone(zone1, width, order=0)
difference = gDifference(newZ1.1,newZ1.2)


plot(newZ1.1)
plot(newZ1.2,add=TRUE)
plot(difference,add= TRUE,col="yellow")


########################################################################################
### GERER LES FRONTIERES COMMUNES AVEC LA CARTE  #######################################
########################################################################################

zone = Z[[2]]
plot(zone)

zone.df = geom(zone)
which(zone.df[,5] == 0)
which(zone.df[,5] == 1)
which(zone.df[,6] == 0)
which(zone.df[,6] == 1)



# certains points très proches de la frontière ne sont pas pris en compte : ici le point à la ligne 138
zone.df[138,]
# ont donne un seuil pour les détecters

epsilon = 0.001
for (i in 1:nrow(zone.df)){
  if (zone.df[i,5]<=epsilon)
    zone.df[i,5] = -0.2
  if (zone.df[i,5]>=1-epsilon)
    zone.df[i,5] = 1.2
  if (zone.df[i,6]<= epsilon)
    zone.df[i,6] = -0.2
  if (zone.df[i,6]>=1-epsilon)
    zone.df[i,6] = 1.2
}
plot(zone.df[,5:6],type="l")

level = zone.df[,2]
level = levels(as.factor(level))


# transform data frame to Spatial polygon

for (i in 1:length(level)){
  if (i==1){
    P = paste("polygon",i,sep = "")
    assign(P, Polygon(zone.df[which(zone.df[,2]==i), 5:6],hole = FALSE))
  }
  else{
    P = paste("polygon",i,sep = "")
    assign(P, Polygon(zone.df[which(zone.df[,2]==i), 5:6], hole = TRUE))
  }
}

listPolygons = list()
for (i in 1:length(level)){
  listPolygons = c(listPolygons, get(paste("polygon",i,sep = "")))
}

polygons = Polygons(listPolygons,"p")
zone.extended = SpatialPolygons(list(polygons))

width = 0.01
newZ = smoothingZone(zone.extended, width = width,order = 1)

# searche the intersection between the new smoothed zone and the map
map = g3=readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
newzone = gIntersection(newZ,map)


########################################################################################
### TESTER LES FONCTIONS smoothingZone, zone.extended et touch.border ##################
########################################################################################

rm(list=ls())
source("srcZ.R") # source libraries and functions - params are in initParam.R


seed=81
map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=2)
plotMap(map)

ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone

zone = Z[[1]]
plot(zone)
#touch.border
touch.border(zone)

#zone.extended
zone.etendue = zone.extended(zone)
plot(zone.etendue)

#smoothingZone
zone.smoothed = smoothingZone(zone,width = 0.04)
plot(zone.smoothed)

plot(zone,add=TRUE)





