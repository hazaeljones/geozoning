#########################################################
#' optCostL
#'
#' @details description, a paragraph
#' @param seed xxxx
#' @param map xxxx
#' @param Cmin xxxx
#' @param f xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
optCostL=function(seed,map,Cmin=0,f=2)
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
costi=c(m1[,"costL"],m2[,"costL"],m3[,"costL"],m4[,"costL"],m5[,"costL"])

nzi=c(m1[,"nz"],m2[,"nz"],m3[,"nz"],m4[,"nz"],m5[,"nz"])
nqi=c(m1[,"nq"],m2[,"nq"],m3[,"nq"],m4[,"nq"],m5[,"nq"])
nli=nqi+1
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
nlk=nli[mask]
costk=costi[mask]
critk=criti[mask]
nqk=nqi[mask]
mqk=mqi[mask,]

Cseq=f*seq(Cmin,1/sigma2G,length=100)

p=list()
nq=list()
crit=list()
cost=list()
lab=list()
mq=list()
C=list()
j=0
for (Cval in Cseq)
        {
	vec=costk+Cval*nlk
	ind=which(vec==min(vec))
	for (jk in ind)
	{
	j=j+1
	pmin=as.numeric(nlk[jk])
	lab[[j]]=round(Cval,3)
	label=paste("C",round(Cval,3),"-q",paste(mqk[jk,],collapse="-"),sep="")
	p[[label]]=pmin
	nq[[label]]=nqk[jk]
	crit[[label]]=critk[jk]
	cost[[label]]=costk[jk]
	mq[[label]]=as.numeric(mqk[jk,])
	C[[label]]=Cval
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
nq=as.matrix(unlist(nq))
pu=as.matrix(pu)
C=as.matrix(unlist(C))
mat=cbind(C,pu,crit,cost,mqm,nq)
mode(mat)="numeric"
rownames(mat)=NULL
colnames(mat)=c("C","p","crit","costL",paste("q",1:5,sep=""),"nq")
plot(lab,pu)

return(list(opt=mat,n=n,nb=nb))
}
