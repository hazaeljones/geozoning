###############################################################
#' readS returns coords, ranges for x and y of a shapefile
#'
#' @details reads a polygon shp file in a directory and extracts coordinates and x and y ranges.
#' @param file file name
#' @param dir directory
#'
#' @return a list with coords, ranges for x and y
#' @importFrom raster shapefile
#'
#' @export
#'
#' @examples
#' # not run
###############################################################
readS = function(file, dir)
{
  #read shp file
  name = paste(dir,file, sep = "")

  if(file.exists(name))
  {
    #shape1 <- readOGR(dsn = dir, layer = file)
    #shape1 <-readOGR(dsn=path.expand(name), layer=file)

    # to avoid bugs with readOGR...
    shape1 <- shapefile(name)

    #obtention of coords
    p = shape1@polygons[[1]]
    p1 = p@Polygons[[1]]
    p2 = p1@coords

    xlim=range(p2[,1])
    ylim=range(p2[,2])
    # return list with coords, ranges for x and y
    return(list(p=p2, xlim = xlim, ylim = ylim,sp=shape1))
  }
  else
  {
    print("Read error: File or folder does not exist!")
  }
}
