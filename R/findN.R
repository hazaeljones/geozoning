###############################################
#' findN
#'
#' @details Find the neighbor zone into which to merge the current zone.
#' It must be a neighbor in the sense of Voronoi polygons. In case of ties, choose the smallest zone for merging into
#' @param K zoning object, as returned by the calNei function
#' @param listN list of neighbor zones
#' @param iZ index of current zone in zoning
#' @param minSize minimum admissible zone size
#'
#' @return the index of the zone into which to merge the current zone
#'
#' @export
#'
#' @examples
#' data(mapTest)
#' data(resZTest)
#' K=resZTest
#' Ns=getNs(K$zoneNModif,4) # neighbors of zone 4
#' listN =  grep( TRUE , Ns) # zones 2 and 5
#' findN(K,listN,4) # zone 4 will be merged into zone 5
#' # not run
findN=function(K,listN,iZ,minSize=0.012)
###############################################
  {
       #Find the neighbor zone with which to merge the current zone
       #It must be a neighbor in the sense of Voronoi polygons
       
  if (length(listN) == 0) return(0)

  if (length(listN) == 1) return(listN[1])

   Z=K$zonePolygone
   iZN=0
   #Check all neighbor zones
   potN = numeric()
   for (i in listN)
            {
              iZNp = i
              tempo = gDistance(Z[[iZNp]],Z[[iZ]])

              #potN useful when current zone is between 2 zones
              if (tempo < 0.001)
              {
                potN = append(potN,i)
              }
             } #end for i

    mina=1
    for (v in potN)
        {
	 area=getSurf(Z,v)
  	 if (area<= mina) # choose the close zone with the smallest area
                  {
                    mina = area
                    iZN= v
                  }
         } #end for v
  
  return(iZN)
}
