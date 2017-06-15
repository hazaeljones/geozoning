
##################################################################
zone.extended = function (z ,boundary)
  ##################################################################
# description : fonction that returns a zone with extended border if the border is in common with the map

# input: 
# z : zone to be extended si touch.border(Z) = TRUE

# output: 
# z : new zone extended

{
  
  boundaryLineExtend = gBoundary(gConvexHull(gBuffer(boundary,width = 0.2)))
  
  if(touch.border(z, boundary)){
    lineInter = gIntersection(gBoundary(boundary),z) # intersection of zone and the boundary of the map
    li = geom(lineInter) # transform lineInter into a dataframe
    level = length(levels(as.factor(li[,2]))) # compute the number of pieces of lines in the spatial line
    for(i in 1:level){ 
      numPoint = which(li[,2]==i) # index of the point that belongs to the piece of the line i
      text = "POLYGON(("
      for (j in numPoint){
        text = paste(text,li[j,4],li[j,5],",")
        point = readWKT(paste("POINT(",li[j,4],li[j,5],")"))
        ptProject = gNearestPoints(point,boundaryLineExtend)@coords[2,]
        text = paste(text,ptProject[1],ptProject[2],",")
      }
      text = paste(text,li[numPoint[1],4],li[numPoint[1],5],"))")
      smallZ.extend = gConvexHull(readWKT(text))
      z = gUnion(z,smallZ.extend)
      z = gBuffer(z,width = 0)
    }
  }
  
  return(z)
}






