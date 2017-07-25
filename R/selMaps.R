#########################################################
#' selMaps
#'
#' @details reads loopQ1-5 results, removes degenerate results and returns criteria and costs for each quantile probability vector. 
#' @param seed simulation seed
#' @param med criterion median vector, if not NULL, these values will be subtracted from criterion values
#' @param medC cost median vector, if not NULL, these values will be subtracted from cost values
#' @param m1 dataset with loopQ1 results
#' @param m2 dataset with loopQ2 results
#' @param m3 dataset with loopQ3 results
#' @param m4 dataset with loopQ4 results
#' @param m5 dataset with loopQ5 results
#'
#' @return a list with components
#'\describe{
#' \item{crit}{vector of best criterion for each loopQ}
#' \item{crit}{vector of cost corresponding to best criterion for each loopQ}
#' \item{wcrit}{vector of worst criterion for each loopQ}
#' }
#'
#' @export
#'
#' @examples
#' assuming that loopQ results for seed 33 were saved in RESD directory
#' selMaps(seed=33,basefile="RESD/res-simuseed")

#' # not run
selMaps=function(seed=89,thr=0.5,med=NULL,medC=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,pErr=0.9,basefile="res-simuseed")
#########################################################
{
if(is.null(m1)) m1=read.table(paste(basefile,seed,"-1q-pE",pErr,".csv",sep=""))
if(is.null(m2)) m2=read.table(paste(basefile,seed,"-2q-pE",pErr,".csv",sep=""))
if(is.null(m3)) m3=read.table(paste(basefile,seed,"-3q-pE",pErr,".csv",sep=""))
if(is.null(m4)) m4=read.table(paste(basefile,seed,"-4q-pE",pErr,".csv",sep=""))
if(is.null(m5)) m5=read.table(paste(basefile,seed,"-5q-pE",pErr,".csv",sep=""))
# remove degenerate quantiles
mask2=m2[,"nq"]==2
mb2=m2[mask2,]
mask3=m3[,"nq"]==3
mb3=m3[mask3,]
mask4=m4[,"nq"]==4
mb4=m4[mask4,]
mask5=m5[,"nq"]==5
mb5=m5[mask5,]

if(is.null(med)) med=rep(0,5)
if(is.null(medC)) medC=rep(0,5)

crit1=m1[1,"crit"] - med[1]
crit2=mb2[1,"crit"] - med[2]
crit3=mb3[1,"crit"] - med[3]
crit4=mb4[1,"crit"] - med[4]
crit5=mb5[1,"crit"] - med[5]

cost1=m1[1,"cost"] - medC[1]
cost2=mb2[1,"cost"] - medC[2]
cost3=mb3[1,"cost"] - medC[3]
cost4=mb4[1,"cost"] - medC[4]
cost5=mb5[1,"cost"] - medC[5]

wcrit1=m1[nrow(m1),"crit"] - med[1]
wcrit2=mb2[nrow(mb2),"crit"] - med[2]
wcrit3=mb3[nrow(mb3),"crit"] - med[3]
wcrit4=mb4[nrow(mb4),"crit"] - med[4]
wcrit5=mb5[nrow(mb5),"crit"] - med[5]

keep=TRUE

return(list(keep=keep,crit=c(crit1,crit2,crit3,crit4,crit5),
            cost=c(cost1,cost2,cost3,cost4,cost5),
            wcrit=c(wcrit1,wcrit2,wcrit3,wcrit4,wcrit5)))
}
