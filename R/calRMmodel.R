###############################
#' transform VGM model into model usable by RandomFields
#'
#' @param vgmodel à remplir
#'
#' @return anisotropy matrix
#'
#' @export
#' @seealso http://www.techmat.vgtu.lt/~art/proc/file/BudrLi.pdf
#' @importFrom RandomFields RMgauss RMspheric
#'
#' @examples
#' # not run
calRMmodel=function(vgmodel)
###############################
{
  matAniso=NULL
  modelVar=NULL
  #if anisotropic
  if(vgmodel$anis1!=1)
  {
    #rotation
    angle=(vgmodel$ang1)*pi/180
    matR=matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),ncol=2,nrow=2)
    #transformation
    amax=vgmodel$range
    amin=vgmodel$anis1*amax
    matT=matrix(c(1/amax,0,0,1/amin),2,2)
    #on obtient la matrice d'anisotropie
    matAniso=matT%*%matR
  }
  else
  {
    matAniso=matrix(c(1,0,0,1),ncol=2,nrow=2)
  }

  if(vgmodel$model=="Gau")
  {
    if(vgmodel$anis1!=1)
	modelVar=RMgauss(var=vgmodel$psill,scale=vgmodel$range,Aniso=matAniso)
    else
	modelVar=RMgauss(var=vgmodel$psill,scale=vgmodel$range)
  }
  if(vgmodel$model=="Sph")
  {
    modelVar=RMspheric(var=vgmodel$psill,scale=vgmodel$range,Aniso=matAniso)

  }
  return(modelVar)

}
