#################################################################################
#' generate empty grid
#'
#' @details description, a paragraph
#' @param mat xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
#################################################################################
# generate empty grid
gridXY=function(mat)
{
  # x en col 1, y en col 2, z en col3
  #
      xvide= unique(mat[,1])
      yvide= unique(mat[,2])
      tvide = matrix(NA,nrow=length(xvide),ncol=length(yvide))
      colnames(tvide)= round(yvide,3)
      rownames(tvide)= round(xvide,3)
      # fill matrix with data values
      for (ii in 1:length(xvide))
      {
	      maskx=(mat[,1]==xvide[ii])
      	mati=mat[maskx,]
      	# complete mati to have a measure for each y value
      	masky=pmatch(mati[,2],yvide)
      	tvide[ii,masky]=mati[,3]
      }
  return(tvide)
}
