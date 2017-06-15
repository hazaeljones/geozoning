cal.max.width.Zone = function(zone, errMax)

# description : function that return the maximal value of the parameter "width" in function gBuffer
#               in order not to make zone disappear (problem that occurs when we realise an interior erosion)

# method used : dichrotomy (between 0 and 1)

# input :
  # zone : zone to be smoothed
  # errMax : maximal error between the best value of width and its approximation

# output :
  # width : approxition of the biggest value of "width"

{
  min = 0
  max = 0.05
  width = 1/2*(min+max)
  while((max-min)>errMax){
    buff = gBuffer(zone,width = -width)
    
    if (is.null(buff)){
      max = width
    }else{
      if (length(buff@polygons[[1]]@Polygons)>1){
        max = width
      }else{
        min = width
      }
    }
    width = 1/2*(max+min)
  }
  width = min
  return(width)
}


##############################################################################################################################


cal.max.width.Map = function(Z, errMax)

  # description : function that return the maximal value of the parameter "width"
  #               in order not to make zones disappear from a map (problem that occurs when we realise an interior erosion)

  # input :
    # Z : map containing zones to be smoothed
    # errMax : maximal error between the best value of width and its approximation

  # output :
    # width : approxition of the biggest value of "width" that could be used in smoothingMap()

{
  listWidth = rep(0,length(Z))
  for (i in 1:length(Z)){
    listWidth[i] = cal.max.width.Zone(Z[[i]], errMax = errMax,boundary = boundary)
  }
  width = min(listWidth)
  return(width)
}
