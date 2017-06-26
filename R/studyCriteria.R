###############################################################################
studyCriteria=function(vseed=89,krig=1,Vpsill=5,Vrange=0.2,Vnugget=0,Vmean=8)
###############################################################################
{
if (is.null(vseed)) vseed=floor(runif(1,100,10000))
for (seed in vseed)
    {
    print(paste("seed=",seed))
    # prepare simu
    map=genMap(DataObj=NULL,seed=seed,krig=krig,Vpsill=Vpsill,Vrange=Vrange,Vnugget=Vnugget,Vmean=Vmean,disp=0)
    # run optim 1q
    print(paste("loopQ1 starts"))
    m1=loopQ1(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ1 ended"))
    write.table(m1,paste("res-simuseed",seed,"-1q-pE",pErr,".csv",sep=""))
    # run optim 2q
    print(paste("loopQ2 starts"))
    m2=loopQ2(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ2 ended"))
    write.table(m2,paste("res-simuseed",seed,"-2q-pE",pErr,".csv",sep=""))
    # run optim 3q
    print(paste("loopQ3 starts"))
    m3=loopQ3(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ3 ended"))
    write.table(m3,paste("res-simuseed",seed,"-3q-pE",pErr,".csv",sep=""))
    # run optim 4q
    print(paste("loopQ4 starts"))
    m4=loopQ4(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ4 ended"))
    write.table(m4,paste("res-simuseed",seed,"-4q-pE",pErr,".csv",sep=""))
    print(paste("loopQ5 starts"))
    m5=loopQ5(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ5 ended"))
    write.table(m5,paste("res-simuseed",seed,"-5q-pE",pErr,".csv",sep=""))
    pdf=paste("figCrit2-seed",seed,".pdf",sep="")
    figCritN(seed=seed,m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,NEW=TRUE,ONE=TRUE,pdf=pdf)
     }
    return()
}
