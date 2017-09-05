###################################################################
loopQ2=function(map,disp=1,step=0.075,QUIET=FALSE)
###################################################################
{
#loop with correction-2 quantiles
iseq=seq(0.05,0.95,step)
kseq=seq(0.125,0.95,step)
#diffQ=0.1
diffQ=0.15

r=data.frame(NULL)

for (i in iseq)
{
  for(k in kseq)
  {
    if ((k-i) < (diffQ-1e-3)) next
    resC=correctionTree(c(i,k),map,disp=disp,SAVE=F)
    critList=resC$critList
    costList=resC$costList
    costLList=resC$costLList
    nzList=resC$nzList
    criti=resC$bestcrit
  # examine critList (sorted by number of effective quantiles)
    # critList written in sortCrit
    n=names(critList)
    for (qq in 1:length(critList))
    	{
	nq=sapply(strsplit(n[qq],"q"),function(x){return(x[2])})
	nq=as.numeric(nq)
	crit=critList[[qq]][1]
	co=costList[[qq]][1]
	coL=costLList[[qq]][1]
	nz=nzList[[qq]][1]
	r=rbind(r,c(round(crit,3),round(co,3),round(coL,3),nz,i,k,nq))
	if(!QUIET) print(paste(i,k,"criterion=",round(crit,3),"cost=",round(co,3),"costL=",round(coL,3),"nz=",nz,"nq=",nq))
	}
    
  }
}

colnames(r)=c("crit","cost","costL","nz","iq","kq","nq")
ro=r[rev(order(r[,"nq"],r[,"crit"])),]


return(ro)
}
