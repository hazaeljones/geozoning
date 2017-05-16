#########################################################
#' optCrit
#'
#' @details description, a paragraph
#' @param vseed xxxx
#' @param k xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
optRank=function(vseed=NULL,k=5)
#########################################################
{
  if(is.null(vseed)) vseed=listSeeds()

  for (seed in vseed)
  {
    file=paste("opt-seed",seed,"crit-cost.pdf",sep="")
    pdf(file,width=7,height=3.5)
    par(mfrow=c(1,2))
    map=genMap(DataObj=NULL,seed,disp=0)
    resCrit=optCrit(seed,map,Cmin=0,f=10,disp=1)
    matCrit=resCrit$opt
    matCrit=matCrit[,c("p","crit","cost","q1","q2","q3","q4","q5","nq")]
    matCrit=matCrit[,c("p","crit","cost","q1","q2","q3","q4","q5","nq")]
    matCrit=matCrit[!duplicated(matCrit),,drop=F]
 #
    resCost=optCost(seed,map,Cmin=0,f=1,disp=1)
    dev.off()

    matCost=resCost$opt
    matCost=matCost[,c("p","crit","cost","q1","q2","q3","q4","q5","nq")]
    matCost=matCost[!duplicated(matCost),,drop=F]
 #
    file=paste("critAdj-seed",seed,".csv",sep="")
    write.table(matCrit,file)
#
    file=paste("costAdj-seed",seed,".csv",sep="")
    write.table(matCost,file)
 #
    file=paste("opt-seed",seed,"-critA.pdf",sep="")
    plotmat(matCrit,pdf=file)
    file=paste("opt-seed",seed,"-costA.pdf",sep="")
    plotmat(matCost,pdf=file)
 #
    # k best criteria for nq=1 to 5
    m=getKMloop(seed=seed,k=k)
    file=paste("opt-seed",seed,"-m.pdf",sep="")
    plotmat(m,pdf=file)

    file=paste("Allopt-seed",seed,".pdf",sep="")
    plotOpt(seed,file,k=6)
 }
 return(vseed)
}
