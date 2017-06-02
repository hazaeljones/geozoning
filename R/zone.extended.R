
##################################################################
zone.extended = function (z)
  ##################################################################
# description : fonction that returns a zone with extended border if the border is in common with the map

# input:
# z : zone to be extended si touch.border(z) = TRUE

# output:
# z.extended : new zone extended

{

  z.df = geom(z)
  level = z.df[,2]
  level = levels(as.factor(level))


  z.df = geom(z)
  epsilon = 0.0001
  for (i in 1:nrow(z.df)){
    if (z.df[i,5]<=epsilon)
      z.df[i,5] = -0.5
    if (z.df[i,5]>=1-epsilon)
      z.df[i,5] = 1.5
    if (z.df[i,6]<= epsilon)
      z.df[i,6] = -0.5
    if (z.df[i,6]>=1-epsilon)
      z.df[i,6] = 1.5
  }



  # transform data frame to Spatial polygon

  for (i in 1:length(level)){
    if (i==1){ # if i=1 then polygon1 contains all others polygons : hole = FALSE
      P = paste("polygon",i,sep = "")
      assign(P, Polygon(z.df[which(z.df[,2]==i), 5:6],hole = FALSE))
    }
    else{ # if i!=1 then polygoni is contained in polygon1 : hole = TRUE
      P = paste("polygon",i,sep = "")
      assign(P, Polygon(z.df[which(z.df[,2]==i), 5:6], hole = TRUE))
    }
  }

  listPolygons = list()
  for (i in 1:length(level)){
    listPolygons = c(listPolygons, get(paste("polygon",i,sep = "")))
  }

  polygons = Polygons(listPolygons,ID = "p")
  comment(polygons) = createPolygonsComment(polygons)

  z.extended = SpatialPolygons(list(polygons))

  return(z.extended)
}
