#####################################################
testFusion=function(qProb,map,Z,K,iC)
#####################################################
{
# 
print("Initial zoning")
printZid(Z)
x11()
dispZ(map$step,map$krigGrid,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
# merge small zone into neighbour zone
Ns = getNs(K,iC)
zg=zoneFusion3(Z,K,iC,Ns,map,minSize,simplitol,disp=1)
K2=calNei(zg,map$krigData,map$krigSurfVoronoi,map$krigN)	
Z2=K2$zonePolygone
# transfer zone labels in K2$lab from K$lab 
# use zone id correspondence
# do not reassign labels
K2=trLabZone(K,K2,Z,Z2,map,qProb,disp=0)

# afficher ids
print("Final zoning")

if (!is.null(zg))
{
for (ii in 1:length(zg)){print(paste("ii=",ii," ID=", zg[[ii]]@polygons[[1]]@ID))}
x11()
  dispZ(map$step,map$krigGrid,zonePolygone=zg,boundary=map$boundary,nbLvl=0)
  }
  else
  print("NULL")
#
return(list(zg=zg,kg=K2))
}