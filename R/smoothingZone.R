##################################################################
smoothingZone = function (z,width,boundary, order=0)
##################################################################
# description : fucntion return a new smoothed zone

# input:
# z : zone to be modified
# width : smoothing parameter in gBuffer
# order = 1 if dilatation-erosion-erosion-dilatation
#       = 0 if erosion-dilatation-dilatation-erosion

# output:
# newZ : smoothed zone



{
  

  z = zone.extended(z,boundary)
  
  widthMax = cal.max.width.Zone(z,errMax = 10^-6)
  width2 = width
  if(width>widthMax){
    width2 = widthMax
  }

  if (order !=0) {
    dilatation1 = gBuffer(z,width = width,joinStyle="ROUND",capStyle = "ROUND")
    erosion1 = gBuffer(dilatation1,width = -width,joinStyle="ROUND",capStyle = "ROUND")
    erosion2 = gBuffer(erosion1,width = -width2,joinStyle="ROUND",capStyle = "ROUND")
    newZ = gBuffer(erosion2,width = width2,joinStyle="ROUND",capStyle = "ROUND")
  }else{
    erosion1 = gBuffer(z,width = -width2,joinStyle="ROUND",capStyle = "ROUND")
    dilatation1 = gBuffer(erosion1,width = width2,joinStyle="ROUND",capStyle = "ROUND")
    dilatation2 = gBuffer(dilatation1,width = width,joinStyle="ROUND",capStyle = "ROUND")
    newZ = gBuffer(dilatation2,width = -width,joinStyle="ROUND",capStyle = "ROUND")
  }
  # search the intersection between the new smoothed zone and the map
  # 2 ways to do: (we have to check if the geometry is valid)
  # 1st way : we search directly the intersection
  # 2nd way : diff = boundary - newZ , then , newZ = boundary-diff
  
  # 1st way
  newZ1 = gIntersection(newZ,boundary)
  newZ1 = gBuffer(newZ1,width = 0)
  
  # 2nd way
  diff = gDifference(boundary, newZ)
  diff = gBuffer(diff,width = 0)
  newZ2 = gDifference(boundary,diff)
  
  if(is.null(newZ1) == FALSE){
    newZ = newZ1
  }else{
    newZ = newZ2
  }
  
  return(newZ)
}


