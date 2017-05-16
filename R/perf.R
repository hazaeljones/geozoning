#########################################################
#' os
#'
#' @details description, a paragraph
#' @param li xxxx
#'
#' @return a ?
#' @importFrom stats na.omit var
#' @importFrom utils object.size write.table
#' @importFrom sp over
#'
#' @export
#'
#' @examples
#' # not run
os=function(li)
#########################################################
{
s=0.0
if(is.null(names(li))) names(li)=1:length(li)
for (n in names(li))
{
o=object.size(li[n])
s=s+o
print(paste(n,":",o))
}

return(s)
}

#########################################################
#' evans
#'
#' @details description, a paragraph
#' @param pts xxxx
#' @param ply xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
evans <- function(pts,ply){
  #########################################################
  prid <- over(pts,ply)
  ptid <- na.omit(prid)
  pt.poly <- pts[as.numeric(as.character(row.names(ptid))),]
  return(pt.poly)
}

#########################################################
#' rowlings
#'
#' @details description, a paragraph
#' @param pts xxxx
#' @param ply xxxx
#'
#' @return a ?
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' # not run
rowlings <- function(pts,ply){
#########################################################
  return(pts[!is.na(over(pts,as(ply,"SpatialPolygons"))),])
}

#########################################################
#' rowlings2
#'
#' @details description, a paragraph
#' @param pts xxxx
#' @param ply xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
rowlings2 <- function(pts,ply){
#########################################################
  class(ply) <- "SpatialPolygons"
  return(pts[!is.na(over(pts,ply)),])
}

#########################################################
#' obrien
#'
#' @details description, a paragraph
#' @param pts xxxx
#' @param ply xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
obrien <- function(pts,ply){
#########################################################
  return(pts[apply(gIntersects(ply,pts,byid=TRUE),1,sum)==1,])
}

#########################################################
#' obrien2
#'
#' @details description, a paragraph
#' @param pts xxxx
#' @param ply xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
obrien2 <- function(pts,ply){
#########################################################
  m=gIntersects(ply,pts,byid=TRUE)
  npts=1:nrow(m)
  return(npts[m])
}
