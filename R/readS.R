###############################################################
#' readS returns coords, ranges for x and y of a shapefile
#'
#' @details reads a polygon shp file in a directory and extracts coordinates and x and y ranges.
#' @param file file name
#' @param dir directory
#'
#' @return a list with components
#' SpatialPolygonsDataFrame, ranges for x and y
#' @importFrom raster shapefile
#'
#' @export
#' @examples
#' \donttest{
#' #  readS was used to create the shape1 object in geozoning package
#' z=readS("Field_8_zones.shp",dir="../data/")
#' plot(z$sp)
#' }
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
    xlim=shape1@bbox[1,]
    ylim=shape1@bbox[2,]
    # return list with coords, ranges for x and y
    return(list(sp=shape1,xlim = xlim, ylim = ylim))
  }
  else
  {
    print("Read error: File or folder does not exist!")
  }
}
