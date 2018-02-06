####################################################################
#' loopQ1
#' @details exploratory loop on probability values associated to quantiles. Performs map zonings for each value of the 1 quantile loop (yielding a 2-label zoning).
#' see also \code{\link{loopQ2}},\code{\link{loopQ3}}, \code{\link{loopQ4}}, \code{\link{loopQ5}} for loops with a different number of labels
#' @param map object returned by function genMa
#' @param disp 0: no info, 1: some info, 2: detailed info
#' @param step loop increment
#' @param minSize zone area threshold under which a zone is too small to be manageable
#' @param minSizeNG zone area threshold under which a zone will be removed
#' @param QUIET run in silence-no display
#'
#' @return a matrix with 6 columns and as many rows as loop elements. Columns contain the following values calculated for each quantile vector:  criterion, cost, cost per label, number of zones, quantile associated probability values and number of non degenerated quantiles.
#'
#' @export
#'
#' @examples
#' # not run
loopQ1=function(map,disp=1,step=0.075,minSize=0.012,minSizeNG=1e-3,QUIET=FALSE)
####################################################################
{
#loop with correction - 1 quantile

iseq=seq(0.05,0.95,step)
r=data.frame(NULL)
for (i in iseq)
{
  resC=correctionTree(i,map,minSize=minSize,minSizeNG=minSizeNG,disp=disp,SAVE=F,ONE=F)
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

return(ro)
}


###################################################################
#' loopQ2
#'
#' @details exploratory loop on probability values associated to quantiles. Performs map zonings for each value of the 1 quantile loop (yielding a 42-label zoning).
#' see also \code{\link{loopQ1}},\code{\link{loopQ3}}, \code{\link{loopQ4}}, \code{\link{loopQ5}} for loops with a different number of labels
#' @param map object returned by function genMa
#' @param disp 0: no info, 1: some info, 2: detailed info
#' @param step loop increment
#' @param minSize zone area threshold under which a zone is too small to be manageable
#' @param minSizeNG zone area threshold under which a zone will be removed
#' @param QUIET run in silence-no display
#'
#' @return a matrix with 7 columns and as many rows as loop elements. Columns contain the following values calculated for each quantile vector:  criterion, cost, cost per label, number of zones, quantile associated probability values and number of non degenerated quantiles.

#' @export
#'
#' @examples
#' # not run
loopQ2=function(map,disp=1,step=0.075,minSize=0.012,minSizeNG=1e-3,QUIET=FALSE)
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
      resC=correctionTree(c(i,k),map,minSize=minSize,minSizeNG=minSizeNG,disp=disp,SAVE=F)
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

#################################################################
#' loopQ3
#'
#' @details exploratory loop on probability values associated to quantiles. Performs map zonings for each value of the 3 quantile loop (yielding a 4-label zoning).
#' see also \code{\link{loopQ1}}, \code{\link{loopQ2}}, \code{\link{loopQ4}}, \code{\link{loopQ5}} for loops with a different number of labels
#' @param map object returned by function genMa
#' @param disp 0: no info, 1: some info, 2: detailed info
#' @param step loop increment
#' @param minSize zone area threshold under which a zone is too small to be manageable
#' @param minSizeNG zone area threshold under which a zone will be removed
#' @param QUIET run in silence-no display
#'
#' @return a matrix with 8 columns and as many rows as loop elements. Columns contain the following values calculated for each quantile vector:  criterion, cost, cost per label, number of zones, quantile associated probability values and number of non degenerated quantiles.
#' @export
#'
#' @examples
#' \donttest{
#' # not run, take a while - >5s CPU
#' seed=10
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=1)
#' loopQ3(map,step=0.1,disp=0,QUIET=TRUE)
#' }
loopQ3=function(map,disp=1,step=0.075,minSize=0.012,minSizeNG=1e-3,QUIET=F)
  #################################################################################
{
  #loop with correction-3 quantiles
  iseq=seq(0.05,0.95,step)
  jseq=seq(0.125,0.95,step)
  kseq=seq(0.200,0.95,step)
  diffQ=0.15
  r=data.frame(NULL)

  for (i in iseq)
  {
    for(j in jseq)
    {
      if((j-i) < (diffQ-1e-3)) next
      for(k in kseq)
      {
        if ((k-j) < (diffQ-1e-3)) next
        resC=correctionTree(c(i,j,k),map,minSize=minSize,minSizeNG=minSizeNG,disp=disp,SAVE=F)
        critList=resC$critList
        costList=resC$costList
        costLList=resC$costLList
        nzList=resC$nzList
        criti=resC$bestcrit
        # examine critList (sorted by number of effective quantiles)
        n=names(critList)
        for (qq in 1:length(critList))
        {
          nq=sapply(strsplit(n[qq],"q"),function(x){return(x[2])})
          nq=as.numeric(nq)
          crit=critList[[qq]][1]
          co=costList[[qq]][1]
          coL=costLList[[qq]][1]
          nz=nzList[[qq]][1]
          r=rbind(r,c(round(crit,3),round(co,3),round(coL,3),nz,i,j,k,nq))
          if(!QUIET) print(paste(i,j,k,"criterion=",round(crit,3),"cost=",round(co,3),"costL=",round(coL,3),"nz=",nz,"nq=",nq))

        }

      }
    }
  }
  colnames(r)=c("crit","cost","costL","nz","iq","jq","kq","nq")
  # best ones
  ro=r[rev(order(r[,"nq"],r[,"crit"])),]

  return(ro)
}

########################################################
#' loopQ4
#'
#' @details exploratory loop on probability values associated to quantiles. Performs map zonings for each value of the 1 quantile loop (yielding a 42-label zoning).
#' see also \code{\link{loopQ1}},\code{\link{loopQ2}}, \code{\link{loopQ3}}, \code{\link{loopQ5}} for loops with a different number of labels
#' @param map object returned by function genMa
#' @param disp 0: no info, 1: some info, 2: detailed info
#' @param step loop increment
#' @param minSize zone area threshold under which a zone is too small to be manageable
#' @param minSizeNG zone area threshold under which a zone will be removed
#' @param QUIET run in silence-no display
#'
#' @return a matrix with 9 columns and as many rows as loop elements. Columns contain the following values calculated for each quantile vector:  criterion, cost, cost per label, number of zones, quantile associated probability values and number of non degenerated quantiles.
#' @export
#'
#' @examples
#' # not run
loopQ4=function(map,disp=1,step=0.075,minSize=0.012,minSizeNG=1e-3,QUIET=F)
########################################################
{
  #loop with correction-4 quantiles

  #sink(paste("res4Q-sT",minSize,"-STNG",minSizeNG,sep=""))
  iseq=seq(0.05,0.95,step)
  jseq=seq(0.125,0.95,step)
  kseq=seq(0.200,0.95,step)
  pseq=seq(0.275,0.95,step)
  diffQ=0.15
  r=data.frame(NULL)

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
          resC=correctionTree(c(i,j,k,p),map,minSize=minSize,minSizeNG=minSizeNG,disp=disp,SAVE=F)
          critList=resC$critList
          costList=resC$costList
          costLList=resC$costLList
          nzList=resC$nzList
          critijkp=resC$bestcrit
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
            r=rbind(r,c(round(crit,3),round(co,3),round(coL,3),nz,i,j,k,p,nq))
          }
          if(!QUIET) print(paste(i,j,k,p,"criterion=",critijkp,"cost=",round(costList[[qq]][1],3),"costL=",round(costLList[[qq]][1],3),"nz=",nz,"nq=",nq))
        }
      }
    }
  }
  colnames(r)=c("crit","cost","costL","nz","iq","jq","kq","pq","nq")
  #
  ro=r[rev(order(r[,"nq"],r[,"crit"])),]


  return(ro)
}

#####################################################
#' loopQ5
#'
#' @details exploratory loop on probability values associated to quantiles. Performs map zonings for each value of the 1 quantile loop (yielding a 42-label zoning).
#' see also \code{\link{loopQ1}},\code{\link{loopQ2}}, \code{\link{loopQ3}}, \code{\link{loopQ4}} for loops with a different number of labels
#' @param map object returned by function genMa
#' @param disp 0: no info, 1: some info, 2: detailed info
#' @param step loop increment
#' @param minSize zone area threshold under which a zone is too small to be manageable
#' @param minSizeNG zone area threshold under which a zone will be removed
#' @param QUIET run in silence-no display
#'
#' @return a matrix with 9 columns and as many rows as loop elements. Columns contain the following values calculated for each quantile vector:  criterion, cost, cost per label, number of zones, quantile associated probability values and number of non degenerated quantiles.
#' @export
#'
#' @examples
#' # not run
loopQ5=function(map,disp=1,step=0.075,minSize=0.012,minSizeNG=1e-3,QUIET=F)
####################################################################################
{
  #loop with correction-5 quantiles
  #sink(paste("res4Q-sT",minSize,"-STNG",minSizeNG,sep=""))
  iseq=seq(0.05,0.95,step)
  jseq=seq(0.125,0.95,step)
  kseq=seq(0.2,0.95,step)
  pseq=seq(0.275,0.95,step)
  qseq=seq(0.35,0.95,step)

  diffQ=0.15
  r=data.frame(NULL)

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
            resC=correctionTree(c(i,j,k,p,q),map,minSize=minSize,minSizeNG=minSizeNG,disp=disp,SAVE=F)
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
  colnames(r)=c("crit","cost","costL","nz","iq","jq","kq","pq","qq","nq")
  #
  ro=r[rev(order(r[,"nq"],r[,"crit"])),]

  return(ro)
}
