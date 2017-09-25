#################################################################
#' loopQ3
#'
#' @details exploratory loop on probability values associated to quantiles. Performs map zonings for each value of the 3 quantile loop (yielding a 4-label zoning).
#' see also \code{\link{loopQ1, loopQ2, loopQ4, loopQ5}} for loopw with a different number of labels
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
#' seed=10
#' map=genMap(DataObj=NULL,seed=seed,disp=FALSE,krig=1)
#' # not run
#' loopQ3(map,step=0.1,disp=0,QUIET=TRUE)
#' 
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
