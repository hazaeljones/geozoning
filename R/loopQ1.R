####################################################################
loopQ1=function(map,disp=1,step=0.075,QUIET=FALSE)
####################################################################
{
#loop with correction - 1 quantile

iseq=seq(0.05,0.95,step)
r=data.frame(NULL)
for (i in iseq)
{
  resC=correctionTree(i,map,disp=disp,SAVE=F,ONE=F)
  critList=resC$critList
  costList=resC$costList
  costLList=resC$costLList
  nzList=resC$nzList
  criti=resC$bestcrit
  #
  co=costList[[1]][1]
  coL=costLList[[1]][1]
  crit=critList[[1]][1]
  nz=nzList[[1]][[1]]
  r=rbind(r,c(round(crit,3),round(co,3),round(coL,3),nz,i,1))
  
  if(!QUIET) print(paste(i,"criterion=",round(criti,3),"cost=",round(co,3),"costL=",round(coL,3),"nz=",nz))
}
colnames(r)=c("crit","cost","costL","nz","iq","nq")
# best ones first
ro=r[rev(order(r[,1])),]

#garbage collection
gc()
return(ro)
}
