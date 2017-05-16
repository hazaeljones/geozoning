#########################################################
#' optCost
#'
#' @details description, a paragraph
#' @param seed xxxx
#' @param map xxxx
#' @param Cmin xxxx
#' @param f xxxx
#' @param disp xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
optCost=function(seed,map,Cmin=0,f=1,disp=1)
#########################################################
{
# make seed data

print(paste("seed=",seed))
# prepare simu
#map=genMap(DataObj=NULL,seed,disp=0)
d=map$krigData
sigma2G=as.numeric(var(d@data))
m1=read.table(paste("res-simuseed",seed,"-1q-pE",pErr,".csv",sep=""))
m2=read.table(paste("res-simuseed",seed,"-2q-pE",pErr,".csv",sep=""))
m3=read.table(paste("res-simuseed",seed,"-3q-pE",pErr,".csv",sep=""))
m4=read.table(paste("res-simuseed",seed,"-4q-pE",pErr,".csv",sep=""))
m5=read.table(paste("res-simuseed",seed,"-5q-pE",pErr,".csv",sep=""))
costi=c(m1[,"cost"],m2[,"cost"],m3[,"cost"],m4[,"cost"],m5[,"cost"])
nzi=c(m1[,"nz"],m2[,"nz"],m3[,"nz"],m4[,"nz"],m5[,"nz"])
nqi=c(m1[,"nq"],m2[,"nq"],m3[,"nq"],m4[,"nq"],m5[,"nq"])
criti=c(m1[,"crit"],m2[,"crit"],m3[,"crit"],m4[,"crit"],m5[,"crit"])

numq1=c()
numqf=c()
for (kk in 1:5)
{
numq1[kk]=5
numqf[kk]=numq1[kk]+kk-1
}
vq1=m1[,c(numq1[1]:numqf[1])]
mq1=cbind(vq1,matrix(NA,nrow=length(vq1),ncol=4))
vq2=m2[,c(numq1[2]:numqf[2])]
mq2=cbind(vq2,matrix(NA,nrow=nrow(vq2),ncol=3))
vq3=m3[,c(numq1[3]:numqf[3])]
mq3=cbind(vq3,matrix(NA,nrow=nrow(vq3),ncol=2))
vq4=m4[,c(numq1[4]:numqf[4])]
mq4=cbind(vq4,matrix(NA,nrow=nrow(vq4),ncol=1))
mq5=m5[,c(numq1[5]:numqf[5])]
colnames(mq1)=paste("q",1:5,sep="")
colnames(mq2)=paste("q",1:5,sep="")
colnames(mq3)=paste("q",1:5,sep="")
colnames(mq4)=paste("q",1:5,sep="")
colnames(mq5)=paste("q",1:5,sep="")
mqi=rbind(mq1,mq2,mq3,mq4,mq5)
colnames(mqi)=paste("q",1:5,sep="")

mask=order(nzi) # sort by number of zones

nzk=nzi[mask]
costk=costi[mask]
critk=criti[mask]
nqk=nqi[mask]
mqk=mqi[mask,]

Cseq=f*seq(Cmin,10/sigma2G,length=1000)

p=list()
nq=list()
crit=list()
cost=list()
lab=list()
mq=list()
C=list()
costAdj=list()
j=0
for (Cval in Cseq)
        {
	vec=costk+Cval*nzk
	ind=which(vec==min(vec))
	for (jk in ind)
	{
	j=j+1
	pmin=as.numeric(nzk[jk])
	lab[[j]]=round(Cval,5)
	p[[j]]=pmin
	nq[[j]]=nqk[jk]
	crit[[j]]=critk[jk]
	cost[[j]]=costk[jk]
	mq[[j]]=as.numeric(mqk[jk,])
	C[[j]]=Cval
	costAdj[[j]]=vec[jk]
	}
	}
pu=unlist(p)
n=table(pu)
nb=table(nzi)

mqm=c()
for (k in 1:length(mq))
mqm=rbind(mqm,mq[[k]])

nam=unlist(lab)
rownames(mqm)=nam
crit=as.matrix(unlist(crit))
cost=as.matrix(unlist(cost))
costAdj=as.matrix(unlist(costAdj))
nq=as.matrix(unlist(nq))
pu=as.matrix(pu)
C=as.matrix(unlist(C))
mat=cbind(C,pu,costAdj,crit,cost,mqm,nq)
mode(mat)="numeric"
rownames(mat)=NULL
colnames(mat)=c("C","p","costAdj","crit","cost",paste("q",1:5,sep=""),"nq")
if(disp) plot(lab,pu,ylim=c(2,14),main=paste("seed=",seed,"-Minimize adj. cost",sep=""),xlab="")

return(list(opt=mat,n=n,nb=nb))
}
