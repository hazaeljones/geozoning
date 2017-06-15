##################################################################
touch.border = function (Z, boundary)
  ##################################################################
# description : test if the zone Z has a common border with the map

# input: 
# Z : zone to be tested
# boundary : a polygon which represent the map 

# output: 
# TRUE : if Z has a common border with the map
# FALSE : otherwise

{
  lineBoundary = gBoundary(boundary) # transform the polygon to a line
  
  if (gDistance(Z,lineBoundary) > 10^-3) {
    res = FALSE
  }
  else{
    res = TRUE
  }
  return(res)
}

