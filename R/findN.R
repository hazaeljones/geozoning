###############################################
#' findN
#'
#' @details description, a paragraph
#' @param Z xxxx
#' @param K xxxx
#' @param listN xxxx
#' @param iC xxxx
#' @param minSize xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
findN=function(Z,K,listN,iC,minSize=0.012)
###############################################
  {
       #On détermine la zone qui va absorber notre petite zone
       #elles doivent se toucher (voisines par Voronoi)
       # si plusieurs choix, fusionner avec la plus petite si sa surface < surface minimum
       # sinon fusionner avec la zone qui a la moyenne la plus proche

          if (length(listN) == 0) return(0)

          if (length(listN) == 1)
          {
            indZV = listN[1]
          }
          else
          {
	    indZV=0
            #On parcourt tous les voisins de notre zone
	    # modif bch october 2016
	    potN = numeric()
            for (i in listN)
            {
              indZVp = i
              tempo = gDistance(Z[[indZVp]],Z[[iC]])

	      # modif bch october 2016
              #potN = numeric()
              #potN est utile lorsque notre petite zone est entre 2 zones
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
