###############################################
#' findN
#'
#' @details description, a paragraph
#' @param K zoning object, as returned by the calNei function
#' @param listN list of neighbor zones
#' @param iC index of current zone in zoning
#' @param minSize minimum admissible zone size
#'
#' @return the index of the zone with which to merge the current zone
#'
#' @export
#'
#' @examples

#' # not run
findN=function(K,listN,iC,minSize=0.012)
###############################################
  {
       #Find the neighbor zone with which to merge the current zone
       #It must be a neighbor in the sense of Voronoi polygons
       
          if (length(listN) == 0) return(0)

          if (length(listN) == 1)
          {
            indZV = listN[1]
          }
          else
          {
	    Z=K$zonePolygone
	    indZV=0
            #Check all neighbor zones
	    potN = numeric()
            for (i in listN)
            {
              indZVp = i
              tempo = gDistance(Z[[indZVp]],Z[[iC]])

	      # modif bch october 2016
              #potN = numeric()
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
                  if (area<= mina)
                  {
                    mina = area
                    indZV= v
                  }
                } #end for v
		if (mina >=minSize)
		   {
		   diffmu=abs(K$meanZone[potN]-K$meanZone[iC])
		   iVP=which(diffmu==min(diffmu))
		   indZV=potN[iVP]
		   }
    } #end else

  return(indZV)
}
