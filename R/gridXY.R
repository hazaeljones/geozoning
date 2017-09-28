#################################################################################
#' generate empty grid
#'
#' @details generate rectangular empty grid corresponding to x and y values in matrix
#' @param mat matrix with x and y coordinates in the first two columns, data in third column
#'
#' @return a grid 
#'
#' @keywords internal
#'
#' @examples
#' data(dataReg)
#' geozoning:::gridXY(dataReg)
#' # not run
#################################################################################
# generate empty grid
gridXY=function(mat)
{
  # x in col 1, y in col 2, z in col3
  #
      xempty= unique(mat[,1])
      yempty= unique(mat[,2])
      tempty = matrix(NA,nrow=length(xempty),ncol=length(yempty))
      colnames(tempty)= round(yempty,3)
      rownames(tempty)= round(xempty,3)
      # fill matrix with data values
      for (ii in 1:length(xempty))
      {
	      maskx=(mat[,1]==xempty[ii])
      	      mati=mat[maskx,]
      	      # complete mati to have a measure for each y value
      	      masky=pmatch(mati[,2],yempty)
      	      tempty[ii,masky]=mati[,3]
      }
  return(tempty)
}
