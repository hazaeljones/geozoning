correctBoundaryMap = function(Z,zN)

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
        listZ = correctBoundary(Z[[i]],Z[[j]])
        Z[[i]] = listZ[[1]]
        Z[[j]] = listZ[[2]]
        print(paste(i,j))
      }
    }
  }
  return(Z)
}
