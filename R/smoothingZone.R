##################################################################
smoothingZone = function (Z,width,order=1)
##################################################################
# description : fucntion return a new smoothed zone

# input:
# Z : zone to be modified
# width : smoothing parameter in gBuffer
# order = 1 if dilatation-erosion-erosion-dilatation
#       = 0 if erosion-dilatation-dilatation-erosion

# output:
# newZ : smoothed zone



{
source("touch.border.R")
source("zone.extended.R")

if (touch.border(Z)){
  Z = zone.extended(Z)
}

if (order == 1) {
  dilatation1 = gBuffer(Z,width = width,joinStyle="ROUND",capStyle = "ROUND")
  erosion1 = gBuffer(dilatation1,width = -width,joinStyle="ROUND",capStyle = "ROUND")
  erosion2 = gBuffer(erosion1,width = -width,joinStyle="ROUND",capStyle = "ROUND")
  newZ = gBuffer(erosion2,width = width,joinStyle="ROUND",capStyle = "ROUND")
}
else{
  erosion1 = gBuffer(Z,width = -width,joinStyle="ROUND",capStyle = "ROUND")
  dilatation1 = gBuffer(erosion1,width = width,joinStyle="ROUND",capStyle = "ROUND")
  dilatation2 = gBuffer(dilatation1,width = width,joinStyle="ROUND",capStyle = "ROUND")
  newZ = gBuffer(dilatation2,width = -width,joinStyle="ROUND",capStyle = "ROUND")
}
  # searche the intersection between the new smoothed zone and the map
  map = g3=readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")
  newZ = gIntersection(newZ,map)
  return(newZ)
}

