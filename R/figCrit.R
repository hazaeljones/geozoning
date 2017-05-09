#########################################################
#' figCrit
#'
#' @details description, a paragraph
#' @param seed xxxx
#' @param gr xxxx
#' @param m1 xxxx
#' @param m2 xxxx
#' @param m3 xxxx
#' @param m4 xxxx
#' @param NEW xxxx
#' @param ONE xxxx
#' @param title xxxx
#'
#' @return a plot
#' @importFrom graphics axis legend
#'
#' @export
#'
#' @examples
#' # not run
figCrit=function(seed=89,gr=1,m1=NULL,m2=NULL,m3=NULL,m4=NULL,NEW=FALSE,ONE=FALSE,title=NULL)
#########################################################
{
  if(is.null(m1)) m1=read.table(paste("res-simuseed",seed,"-1q.csv",sep=""))
  if(is.null(m2)) m2=read.table(paste("res-simuseed",seed,"-2q.csv",sep=""))
  if(is.null(m3)) m3=read.table(paste("res-simuseed",seed,"-3q.csv",sep=""))
  if(is.null(m4)) m4=read.table(paste("res-simuseed",seed,"-4q.csv",sep=""))

  # simus- recup matrices 1q 2q 3q 4q - plot avec
  # en x nq
  # en y pour chaque q, crit (noir) et critM(rouge)
  # tels que (max(crit)-crit) <=1

  mask1=(-m1[,"crit"]+ m1[1,"crit"])<1
  mb1=m1[mask1,]
  mask2=((-m2[,"crit"]+ m2[1,"crit"])<1) & (m2[,"nq"]==2)
  mb2=m2[mask2,]
  mask3=((-m3[,"crit"]+ m3[1,"crit"])<1) & (m3[,"nq"]==3)
  mb3=m3[mask3,]
  mask4=((-m4[,"crit"]+ m4[1,"crit"])<1) & (m4[,"nq"]==4)
  mb4=m4[mask4,]

  if (is.null(title)) title="Gaussian field simulation"
  nq=2:5
  best=rep(0,10)
  if (NEW)
  {
    pdf(paste("figCrit2-seed",seed,".pdf",sep=""))
    if (!ONE) par(mfrow=c(3,2))
  }
  maxy=max(mb1[1,"crit"],mb2[1,"crit"],mb3[1,"crit"],mb4[1,"crit"])
  if(gr==1)
  plot(nq,rep(mb1[1,"crit"],4),type="n",xlim=c(1.5,6),ylim=c(0,ceiling(maxy)),ylab="Criterion",
       xlab="Number of labels",xaxt="n",main=title)
  else
  plot(nq,rep(mb1[1,"crit"],4),type="n",xlim=c(1.5,6),ylim=c(0,ceiling(maxy))*0.5,ylab="criteria",
       xlab="number of labels",xaxt="n",main=paste("seed=",seed,sep=""))

  axis(1,at=2:5)
  n=nrow(mb1)
  step=0.05
  if(gr==1){
  points(2+seq(-n/2,n/2,length=n)*step,mb1[,"crit"],pch=20,cex=1.5)
  #points(2+c(0:(n-1))*step,mb1[,"critM"],col="red")
  }
  else{
    mdb=mb1[1,"crit"]- mb1[1,"critM"]
    points(2,mdb,col="blue")
  }
  best[1]=mb1[1,"iq"]
  n=nrow(mb2)
  if(gr==1){
    points(3+seq(-n/2,n/2,length=n)*step,mb2[,"crit"],pch=20,cex=1.5)
    #points(3+c(0:(n-1))*step,mb2[,"critM"],col="red")
  } else{
    mdb=mb2[1,"crit"]- mb2[1,"critM"]
    points(3,mdb,col="blue")
  }
  best[2:3]=mb2[1,c("iq","kq")]
  n=nrow(mb3)
  if(gr==1){
    points(4+seq(-n/2,n/2,length=n)*step,mb3[,"crit"],pch=20,cex=1.5)
    #points(4+c(0:(n-1))*step,mb3[,"critM"],col="red")
  } else{
    mdb=mb3[1,"crit"]- mb3[1,"critM"]
    points(4,mdb,col="blue")
  }
  best[4:6]=mb3[1,c("iq","jq","kq")]
  n=nrow(mb4)
  if(gr==1){
    points(5+seq(-n/2,n/2,length=n)*step,mb4[,"crit"],pch=20,cex=1.5)
    #points(5+c(0:(n-1))*step,mb4[,"critM"],col="red")
  } else{
    mdb=mb4[1,"crit"]- mb4[1,"critM"]
    points(5,mdb,col="blue")
  }
  best[7:10]=mb4[1,c("iq","jq","kq","pq")]

  #legend(x=5,y=5,leg=c("Crit"),col=c("black"),pch=c(20))
  if(ONE) dev.off()

  return(unlist(best))
}

