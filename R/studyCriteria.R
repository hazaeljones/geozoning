###############################
#' studyCriteria
#'
#' @details description, a paragraph
#' @param vseed xxxx
#'
#' @return a ?
#' @importFrom graphics boxplot persp
#'
#' @export
#'
#' @examples
#' # not run
studyCriteria=function(vseed=89)
###############################
{
  if (is.null(vseed)) vseed=floor(runif(1,100,10000))
  for (seed in vseed)
  {
    print(paste("seed=",seed))
    # prepare simu
    map=genMap(DataObj=NULL,seed,disp=0)
    # run optim 1q
    m1=loopQ1(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ1 ended"))
    write.table(m1,paste("res-simuseed",seed,"-1q-pE",pErr,".csv",sep=""))
    # run optim 2q
    m2=loopQ2(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ2 ended"))
    write.table(m2,paste("res-simuseed",seed,"-2q-pE",pErr,".csv",sep=""))
    # run optim 3q
    m3=loopQ3(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ3 ended"))
    write.table(m3,paste("res-simuseed",seed,"-3q-pE",pErr,".csv",sep=""))
    # run optim 4q
    m4=loopQ4(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ4 ended"))
    write.table(m4,paste("res-simuseed",seed,"-4q-pE",pErr,".csv",sep=""))
    m5=loopQ5(map,disp=0,step=0.075,QUIET=TRUE)
    print(paste("loopQ5 ended"))
    write.table(m5,paste("res-simuseed",seed,"-5q-pE",pErr,".csv",sep=""))
    pdf=paste("figCrit2-seed",seed,".pdf",sep="")
    figCritN(seed=seed,m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,NEW=TRUE,ONE=TRUE,pdf=pdf)
  }
  return()
}

##############################
#' studyCriteria2
#'
#' @details description, a paragraph
#' @param vseed xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
studyCriteria2=function(vseed)
###############################
{


  for (seed in vseed)
  {
    pdf=paste("figCritZ-seed",seed,".pdf",sep="")
    best=figCritN(seed,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,NEW=TRUE,ONE=FALSE,pdf=pdf,title=paste("seed=",seed,sep=""))
    best=figCritN(seed,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,NEW=FALSE,ONE=FALSE)
    # best 1q
    map=genMap(DataObj=NULL,seed=seed,disp=0)
    visuZ(best[1],map,nq=1)
    visuZ(best[2:3],map,nq=2)
    visuZ(best[4:6],map,nq=3)
    visuZ(best[7:10],map,nq=4)
    dev.off()
  }

  return()
}

#################################################################
#' studyCriteria3
#'
#' @details description, a paragraph
#' @param vseed xxxx
#' @param thr xxxx
#' @param med xxxx
#' @param medM xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
studyCriteria3=function(vseed=NULL,thr=0.75,med=NULL,medM=NULL)
#################################################################
{

  if (is.null(vseed))
  {
    vseed=listSeeds()
  }
  keep=rep(FALSE,length(vseed))
  cr=rep(0,length(vseed))
  mcrit=matrix(0,nrow=length(vseed),ncol=4)
  mcritM=mcrit

  for (k in 1:length(vseed))
  {
    # IS 19/05/2017 medM doesn't exist in selMaps definition. Change by medC...
    # note that studyCriteria3() not used...
    res=selMaps(seed=vseed[k],thr=thr,med=med,medC=medM,m1=NULL,m2=NULL,m3=NULL,m4=NULL)
    keep[k]=res$keep
    cr[k]=res$cr
    mcrit[k,]=res$crit
    mcritM[k,]=res$critM
  }
  #print(keep)
  #print(mcrit)
  keeps=vseed[which(!keep)]
  crmat=matrix(0,nrow=length(vseed),ncol=2)
  crmat[,2]=cr
  crmat[,1]=vseed
  #print(crmat[order(crmat[,2]),])

  boxplot(mcrit,names=paste("nL=",2:5,sep=""))
  mcrit=cbind(vseed,mcrit)

  mcritM=cbind(vseed,mcritM)

  return(list(keeps=keeps,mcrit=mcrit,mcritM=mcritM,med=apply(mcrit[,2:5],2,median),medM=apply(mcritM[,2:5],2,median)))
}



#################################################################
#' studyCriteria4
#'
#' @details description, a paragraph
#' @param vseed xxxx
#' @param thr xxxx
#' @param med xxxx
#' @param medC xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
studyCriteria4=function(vseed=NULL,thr=0.75,med=NULL,medC=NULL)
#################################################################
{
  if (is.null(vseed))
  {
    vseed=listSeeds()
  }
  keep=rep(FALSE,length(vseed))
  mcrit=matrix(0,nrow=length(vseed),ncol=5)

  for (k in 1:length(vseed))
  {
    res=selMaps(seed=vseed[k],thr=thr,med=med,medC=medC,m1=NULL,m2=NULL,m3=NULL,m4=NULL)
    keep[k]=res$keep
    mcrit[k,]=res$crit-res$wcrit
  }

  keeps=vseed[which(!keep)]

  # diff best-worst
  boxplot(mcrit,names=paste("nL=",2:6,sep=""))
  mcrit=cbind(vseed,mcrit)

  return(list(keeps=keeps,mcrit=mcrit,med=apply(mcrit[,2:5],2,median)))
}


####################################################
#' studyCriteria5
#'
#' @details description, a paragraph
#' @param map xxxx
#' @param seed xxxx
#' @param m1 xxxx
#' @param m2 xxxx
#' @param m3 xxxx
#' @param m4 xxxx
#' @param m5 xxxx
#' @param title xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
studyCriteria5=function(map,seed,m1,m2,m3,m4,m5,title="")
####################################################
{

  pdf=paste("figCritZ-yield.pdf",sep="")
  best=figCritN(seed=seed,m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,NEW=TRUE,ONE=FALSE,pdf=pdf,title=title)
  # best 1q 2q 3q 4q 5q
  visuZ(best[1],map,nq=1)
  visuZ(best[2:3],map,nq=2)
  visuZ(best[4:6],map,nq=3)
  visuZ(best[7:10],map,nq=4)
  visuZ(best[11:15],map,nq=4)
  dev.off()

  return()
}
