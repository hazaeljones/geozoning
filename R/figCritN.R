#########################################################
#' figCritN
#'
#' @details  reads loopQ1-5 results, filters results by keeping th best criteria ) and plots them together with corresponding costs.
#' @param seed 
#' @param m1 dataset with loopQ1 results
#' @param m2 dataset with loopQ2 results
#' @param m3 dataset with loopQ3 results
#' @param m4 dataset with loopQ4 results
#' @param m5 dataset with loopQ5 results
#' @param NEW new plot
#' @param ONE single plot
#' @param title plot title
#' @param pdf pdf file name
#'
#' @return a vector of probabilities corresponding to best results
#'
#' @export
#'
#' @examples
#' assuming that loopQ results for seed 33 were saved in RESD directory
#' figCritN(seed=33,basefile="RESD/res-simuseed")
#' # not run
figCritN=function(seed=89,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,NEW=FALSE,ONE=FALSE,title="Gaussian field simulation",pdf=NULL,basefile="res-simuseed")
#########################################################
{
 if(is.null(m1)) m1=read.table(paste(basefile,seed,"-1q-pE",pErr,".csv",sep=""))
if(is.null(m2)) m2=read.table(paste(basefile,seed,"-2q-pE",pErr,".csv",sep=""))
if(is.null(m3)) m3=read.table(paste(basefile,seed,"-3q-pE",pErr,".csv",sep=""))
if(is.null(m4)) m4=read.table(paste(basefile,seed,"-4q-pE",pErr,".csv",sep=""))
if(is.null(m5)) m5=read.table(paste(basefile,seed,"-5q-pE",pErr,".csv",sep=""))
# en x nq
# en y pour chaque q, crit (noir) et cost(rouge)
# tels que (max(crit)-crit) <=maxd

if (NEW)
{
if(!is.null(pdf)) pdf(pdf)
if (!ONE) par(mfrow=c(3,2))
}
maxd=0.8
q=1:4
nl=length(q)+1
maxy=0
miny=Inf
best=c()

for (k in q)
{
mk=get(paste("m",k,sep=""))
mbkname=paste("mb",k,sep="")
mask=(-mk[,"crit"]+ mk[1,"crit"])<maxd & mk[,"nq"]==k

mbk=mk[mask,]
if(nrow(mbk)>5) mbk=mbk[1:6,]
assign(mbkname,mbk)
maxy=max(maxy,max(mbk[1,"crit"]))
miny=min(miny,min(mbk[,"crit"]))
}
miny=min(miny,4)

plot(2:nl,rep(mb1[1,"crit"],nl-1),type="n",xlim=c(1.5,6),ylim=c(miny,ceiling(maxy)),ylab="Criteria",xlab="Number of labels",xaxt="n",cex=2,main=title)
axis(1,at=2:nl,cex=2)
step=0.05

for (k in q)
{
mbkname=paste("mb",k,sep="")
mbk=get(mbkname)
n=nrow(mbk)
kx=k+1+seq(-n/2,n/2,length=n)*step
points(kx,mbk[,"crit"],pch=20,cex=1.5)
if (k!=q[length(q)]) lines(c(k+1.5,k+1.5),c(miny*0.8,ceiling(maxy)*1.2),lty=2)
mdb=mbk[,"cost"]
#points(kx,mdb,col="blue")
crit=mbk[,"crit"]-1

numic=(ncol(mbk)-1-k+1):(ncol(mbk)-1)
best=c(best,mbk[1,numic])
}

#legend(x=5.5,y=3,leg=c("Crit","Cost"),col=c("black","blue"),pch=c(20))
if(ONE & !is.null(pdf)) dev.off()

return(unlist(best))
}

