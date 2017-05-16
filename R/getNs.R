##################################################################
#' getNs
#'
#' @details description, a paragraph
#' @param resZ xxxx
#' @param iC xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
##################################################################
getNs=function(resZ,iC)
{
  Ns=resZ$zoneNModif[iC,]
  return(Ns)
}
