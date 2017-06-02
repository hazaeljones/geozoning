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
  epsilon = 10^-4
  if (length(which(Z.df[,5]<=epsilon))==0 & length(which(Z.df[,5]>=1-epsilon))==0 &
      length(which(Z.df[,6]>=1-epsilon))==0 & length(which(Z.df[,6]<=epsilon))==0) {
    res = FALSE
  }
  else{
    res = TRUE
  }
  return(res)
}

