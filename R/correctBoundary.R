# This script was an old version that worked on the square simulated map

# Please use the new version written in correctBoundaryMap.R





correctBoundary1 = function(z1,z2,boundary)

# description : fix the problem linked to the border between two neighbour zones

# input :
  # z1, z2 : neighbour zones to be modified


# output :
  # newZ1 : z1 modified version
  # newZ2 : z2 modified version

{
  #border =readWKT("POLYGON((0 0,1 0,1 1,0 1,0 0))")

  newZ = gDifference(boundary,gUnion(z1,z2))

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


#########################################################################################################################

correctBoundaryMap1 = function(Z,zN)

  # description : fix the problem linked to the border between neighbour zones of the map

  # input :
  # Z : a map with zones that have problem
  # zN : matrix of neighbourhood

  # output :
  # Z

{
  nbZ = length(Z)

  for (i in 1:(nbZ-1)){
    for (j in (i+1):nbZ){
      if (zN[i,j]==TRUE){
        listZ = correctBoundary1(Z[[i]],Z[[j]])
        Z[[i]] = listZ[[1]]
        Z[[j]] = listZ[[2]]
        print(paste(i,j))
      }
    }
  }
  return(Z)
}



#########################################################################################################################
######################### VERSION n°2 of CORRECTBOUNDARY : NOT WORKING ##################################################
# THIS VERSION TAKES IN PARAMETER THE BOUNDARY OF THE MAP
# PLEASE USE THE VERSION WRITTEN IN correctBoundaryMap.R
#########################################################################################################################


correctBoundary = function(z1,z2,boundary)
  
  # description : fix the problem linked to the border between two neighbour zones
  
  # input :
  # z1, z2 : neighbour zones to be modified
  # boundary : new Boundary of the map
  
  
  # output :
  # newZ1 : z1 modified version
  # newZ2 : z2 modified version

{
  
  newZ = gDifference(boundary,gUnion(z1,z2))
  
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



############################################################################################################################################

correctBoundaryMap = function(Z,boundary)
  
  # description : function that correct the problem linked to the boundary between Zones, and between Zone and the map boundary
  
  # INPUT :
  # Z : list of Zones that have problem
  # boundary : initial boundary of the map (spatial polygon)
  
  # OUPUT :
  # Z : new list of Zone with proper boundary
  # newBoundary : new Boundary of the map
  
{
  
  # first : correct the boundary between 1 zone and the map
  newBoundary = gBuffer(boundary,width=-0.0001)
  for (i in 1:length(Z)){
    Z[[i]] = gIntersection(newBoundary, Z[[i]])
  }
  
  # second : correct the border between 2 zones
  nbZ = length(Z)
  # matrix of neighbourhood
  zN = matrix(rep(0,nbZ*nbZ),ncol=nbZ)
  for(i in 1:nbZ){
    for(j in 1:nbZ){
      if (gDistance(Z[[i]],Z[[j]])<10^-3){
        zN[i,j] = TRUE
      }else{
        zN[i,j] = FALSE
      }
    }
  }
  
  for (i in 1:(nbZ-1)){
    for (j in (i+1):nbZ){
      if (zN[i,j]==TRUE){
        listZ = correctBoundary(Z[[i]],Z[[j]], newBoundary)
        Z[[i]] = listZ[[1]]
        Z[[j]] = listZ[[2]]
      }
    }
  }
  return(list(Z,newBoundary))
}



















