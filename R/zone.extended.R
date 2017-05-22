
##################################################################
zone.extended = function (Z)
  ##################################################################
# description : fonction that returns a zone with extended border if the border is in common with the map

# input:
# Z : zone to be extended si touch.border(Z) = TRUE

# output:
# Z.extended : new zone extended

{
  Z.df = geom(Z)
  epsilon = 0.001
  for (i in 1:nrow(Z.df)){
    if (Z.df[i,5]<=epsilon)
      Z.df[i,5] = -0.2
    if (Z.df[i,5]>=1-epsilon)
      Z.df[i,5] = 1.2
    if (Z.df[i,6]<= epsilon)
      Z.df[i,6] = -0.2
    if (Z.df[i,6]>=1-epsilon)
      Z.df[i,6] = 1.2
  }


  level = Z.df[,2]
  level = levels(as.factor(level))


  # transform data frame to Spatial polygon

  for (i in 1:length(level)){
    if (i==1){ # if i=1 then polygon1 contains all others polygons : hole = FALSE
      P = paste("polygon",i,sep = "")
      assign(P, Polygon(Z.df[which(Z.df[,2]==i), 5:6],hole = FALSE))
    }
    else{ # if i!=1 then polygoni is contained in polygon1 : hole = TRUE
      P = paste("polygon",i,sep = "")
      assign(P, Polygon(Z.df[which(Z.df[,2]==i), 5:6], hole = TRUE))
    }
  }

  listPolygons = list()
  for (i in 1:length(level)){
    listPolygons = c(listPolygons, get(paste("polygon",i,sep = "")))
  }

  polygons = Polygons(listPolygons,ID = "p")
  #comment(polygons) = createPolygonsComment(polygons)

  Z.extended = SpatialPolygons(list(polygons))

  return(Z.extended)
}
