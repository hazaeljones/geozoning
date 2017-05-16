###############################################################
#' readS returns coords, ranges for x and y of a shapefile
#'
#' @details description, a paragraph
#' @param file file name
#' @param rep directory
#'
#' @return a ?
#' @importFrom raster shapefile
#'
#' @export
#'
#' @examples
#' # not run
###############################################################
readS = function(file, rep)
{
  #read shp file
  nom = paste(rep,file, sep = "")

  if(file.exists(nom))
  {
    #shape1 <- readOGR(dsn = rep, layer = file)
    #shape1 <-readOGR(dsn=path.expand(nom), layer=file)

    # to avoid bugs with readOGR...
    shape1 <- shapefile(nom)

    #obtention of coords
    p = shape1@polygons[[1]]
    p1 = p@Polygons[[1]]
    p2 = p1@coords

    xlim=range(p2[,1])
    ylim=range(p2[,2])
    # return list with coord, ranges for x and y
    return(list(p=p2, xlim = xlim, ylim = ylim,sp=shape1))
  }
  else
  {
    print("Read error: File or folder does not exist!")
  }
}
