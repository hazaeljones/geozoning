##################################################################
touch.border = function (Z)
##################################################################
# description : test if the zone Z has a common border with the map

# input: 
# Z : zone to be tested

# output: 
# TRUE : if Z has a common border with the map
# FALSE : otherwise

{
  Z.df = geom(Z)
  if (length(which(Z.df[,5]==0))!=0 & length(which(Z.df[,5]==1))!=0 &
      length(which(Z.df[,6]==1))!=0 & length(which(Z.df[,6]==1))!=0) {
    res = FALSE
  }
  else{
    res = TRUE
  }
  return(res)
}

