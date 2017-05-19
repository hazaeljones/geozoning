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



##############################################################
# CAS PRATIQUE
##############################################################
# real case
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


##############################################################
# Calculer le polygone moyenne entre la dilatation et l'erosion
##############################################################
ZK=initialZoning(qProb=c(0.4,0.7),map,pErr,simplitol,optiCrit,disp=0) # names(Z)  "resCrit"  "resDist" "resZ"
Z=ZK$resZ$zonePolygone # zone
zone6 = Z[[6]]
width = 0.02

# Cas 1 : moyenne de la dilatation et de l'erosion sur la zone
buffer1 = gBuffer(zone6,width = width,joinStyle="ROUND",capStyle = "ROUND")
buffer2 = gBuffer(zone6,width = -width,joinStyle="ROUND",capStyle = "ROUND")

buffer1.Data.Frame = geom(buffer1)[,5:6]
buffer2.Data.Frame = geom(buffer2)[,5:6]

zone6.modif = matrix(data = rep(0,2*nrow(buffer1.Data.Frame)),ncol = 2)

for (i in 1:nrow(buffer1.Data.Frame)){
  point = readWKT()
}




# ???????????????????????????????????????????????????????????
# Si label(situé au centroid de la zone) se situe dans l'autre zone, on fait quoi?
# Morphologie : comment on calcule la moyenne de la dilatation et erosion
# Quel est le critère pour choisir le paramètre "width" dans "gBuffer" pour faire la morphologie
# Comment détecter les régions d'un zone où on veut effectuer la morphologie locale























