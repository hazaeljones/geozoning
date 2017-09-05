#####################################################
loopQ5=function(map,disp=1,step=0.075,SAVE=F,QUIET=F)
#####################################################
{
#loop with correction-5 quantiles
#sink(paste("res4Q-sT",minSize,"-STNG",minSizeNG,sep=""))
iseq=seq(0.05,0.95,step)
jseq=seq(0.125,0.95,step)
kseq=seq(0.2,0.95,step)
pseq=seq(0.275,0.95,step)
qseq=seq(0.325,0.95,step)

diffQ=0.15
r=data.frame(NULL)
ro=NULL

for (i in iseq)
{
	for(j in jseq)
	{
		if((j-i) < (diffQ-1e-3)) next
		for(k in kseq)
    	  	{
	  	if ((k-j) < (diffQ-1e-3)) next
		   for(p in pseq)
		   {
		   if((p-k) < (diffQ-1e-3)) next
		     for (q in qseq)
		     {
		      if((q-p) < (diffQ-1e-3)) next
    	  	   resC=correctionTree(c(i,j,k,p,q),map,disp=disp,SAVE=F)
		   critList=resC$critList
    		   costList=resC$costList
    		   costLList=resC$costLList
    		   nzList=resC$nzList
    		   critijkpq=resC$bestcrit
		   # examine critList (sorted by number of effective quantiles)
    		   n=names(critList)
    		   for (qq in 1:length(critList))
    		   {
		   nq=sapply(strsplit(n[qq],"q"),function(x){return(x[2])})
		   nq=as.numeric(nq)
		   crit = critList[[qq]][1]
		   co = costList[[qq]][1]
		   coL = costLList[[qq]][1]
		   nz=nzList[[qq]][1]
		   r=rbind(r,c(round(crit,3),round(co,3),round(coL,3),nz,i,j,k,p,q,nq))
		   }
    		   if(!QUIET) print(paste(i,j,k,p,q,"criterion=",critijkpq,"cost=",round(costList[[qq]][1],3),"costL=",round(costLList[[qq]][1],3),"nz=",nz,"nq=",nq))
		   }
		   }
		}
 	  }
 }
if(nrow(r)>0)
	{
	colnames(r)=c("crit","cost","costL","nz","iq","jq","kq","pq","qq","nq")
	# 
	ro=r[rev(order(r[,"nq"],r[,"crit"])),]
	}

return(ro)
}
