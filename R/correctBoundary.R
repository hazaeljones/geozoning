correctBoundary = function(z1,z2)

# description : fix the problem linked to the border between two neighbour zones

# input :
  # z1, z2 : neighbour zones to be modified


# output :
  # newZ1 : z1 modified version
  # newZ2 : z2 modified version

{
  border =readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")

  newZ = gDifference(border,gUnion(z1,z2))

  newZ.df = geom(newZ) # data frame corresponding to the spatial polygon Z
  level = newZ.df[,2]
  level = levels(as.factor(level))

  area = rep(0,length(level))
  for (i in 1:length(level)){

    P = paste("P",i,sep = "")
    assign(P, Polygon(newZ.df[which(newZ.df[,2]==i), 5:6], hole = FALSE))

    Ps = paste("Ps",i,sep="")
    assign(Ps,Polygons(list(get(paste("P",i,sep=""))), ID = "a"))

    sPs = paste("sPs",i,sep="")
    assign(sPs,SpatialPolygons(list(get(paste("Ps",i,sep="")))))

    ind  = which(newZ.df[,2]==i)
    if (newZ.df[ind[1],4] == 1){
      assign(paste("sPs",i,sep=""), gDifference(get(paste("sPs",i-1,sep="")),get(paste("sPs",i,sep=""))))
    }

    area[i] = gArea(get(paste("sPs",i,sep="")))
  }

  index = which(area<10^-3)

  if(length(index)>0){
    if (is.null(gIntersection(z1,z2))){
      intersection = get(paste("sPs",index[1],sep=""))
    }else{
      intersection = gIntersection(z1,z2)
    }
    for(i in index){
      intersection = gUnion(intersection,get(paste("sPs",i,sep="")))
    }
  }

  union = gUnion(gUnion(z1,z2),intersection)
  newZ2 = z2
  newZ1 = gDifference(union,z2)
  gIntersection(newZ1,newZ2)

  return(list(newZ1,newZ2))
}










