#########################################################
#' figCritN
#'
#' @details description, a paragraph
#' @param seed xxxx
#' @param m1 xxxx
#' @param m2 xxxx
#' @param m3 xxxx
#' @param m4 xxxx
#' @param m5 xxxx
#' @param NEW xxxx
#' @param ONE xxxx
#' @param title xxxx
#' @param pdf xxxx
#'
#' @return a plot
#'
#' @export
#'
#' @examples
#' # not run
figCritN=function(seed=89,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,NEW=FALSE,ONE=FALSE,title="Gaussian field simulation",pdf=NULL)
#########################################################
{
  if(is.null(m1)) m1=read.table(paste("res-simuseed",seed,"-1q-pE",pErr,".csv",sep=""))
  if(is.null(m2)) m2=read.table(paste("res-simuseed",seed,"-2q-pE",pErr,".csv",sep=""))
  if(is.null(m3)) m3=read.table(paste("res-simuseed",seed,"-3q-pE",pErr,".csv",sep=""))
  if(is.null(m4)) m4=read.table(paste("res-simuseed",seed,"-4q-pE",pErr,".csv",sep=""))
  if(is.null(m5)) m5=read.table(paste("res-simuseed",seed,"-5q-pE",pErr,".csv",sep=""))
  # simus- recup matrices 1q 2q 3q 4q - plot avec
  # en x nq
  # en y pour chaque q, crit (noir) et cost(rouge)
  # tels que (max(crit)-crit) <=1
  if (NEW)
  {
    if(!is.null(pdf)) pdf(pdf)
    if (!ONE) par(mfrow=c(3,2))
  }
  q=1:5
  nl=length(q)+1
  maxy=0
  best=c()

  for (k in q)
  {
    mk=get(paste("m",k,sep=""))
    mbkname=paste("mb",k,sep="")
    mask=(-mk[,"crit"]+ mk[1,"crit"])<0.5 & mk[,"nq"]==k
    mbk=mk[mask,]
    assign(mbkname,mbk)
    maxy=max(maxy,max(mbk[1,"crit"]))
  }

  plot(2:nl,rep(mb1[1,"crit"],nl-1),type="n",xlim=c(1.5,6),ylim=c(0,ceiling(maxy)),ylab="Criteria",xlab="Number of labels",xaxt="n",main=title)
  axis(1,at=2:6)
  step=0.05

  for (k in q)
  {
    mbkname=paste("mb",k,sep="")
    mbk=get(mbkname)
    n=nrow(mbk)
    kx=k+1+seq(-n/2,n/2,length=n)*step
    points(kx,mbk[,"crit"],pch=20,cex=1.5)
    mdb=mbk[,"cost"]
    points(kx,mdb,col="blue")
    crit=mbk[,"crit"]-1

    numic=(ncol(mbk)-1-k+1):(ncol(mbk)-1)
    best=c(best,mbk[1,numic])
  }
  legend(x=5.5,y=3,leg=c("Crit","Cost"),col=c("black","blue"),pch=c(20))
  if(ONE & !is.null(pdf)) dev.off()

  return(unlist(best))
}

