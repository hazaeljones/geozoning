##################################################################
#' saveZoningFromSimu
#'
#' @details computes initial zoning obtained from simulation and saves it with map data in R object
#' @param seed simulation seed
#' @param qProb vector of probabilities values used for quantiles
#' @param disp 0: no info, 1: plot zoning
#' @param file: filename for saved R object
#'
#' @return an empty value
#'
#' @export
#'
#' @examples
#' # not run
saveZoningFromSimu=function(seed=80,qProb=0.4,disp=disp,file="Z.Rdata")
##################################################################
{

  # seed=NULL pour renouveler le germe
  #seed=80
  #simulation, DataObj=NULL, de nbPoints=450 (defini dans initParam), coordonnees entre 0 et 1
  # sinon DataObj=data frame des donnees x y z
  # si spatial data point, le transformer en data frame avant l'appel
  map=genMap(DataObj=NULL,seed=seed,disp=disp)
  if(disp) plotMap(map) # affiche les donnees brutes de la map qui doit etre zonee
  # IS 19/05/2017: add comment for x11
  #if(disp) x11()

  ZK=initialZoning(qProb, map,disp=disp)
  #save zoning Z in R objects for reusing later

  resZ=ZK$resZ
  ZPS=resZ$zonePolygone
  NBZ=length(ZPS)

  save(ZPS,NBZ, seed,map,resZ,file=file)
  print("results saved in Z.Rdata")
  print(paste("Crit=",round(ZK$resCrit,3)))#valeur du critere

  print(names(resZ))
  print(paste(length(ZPS),"zones")) # nb of zones

}
