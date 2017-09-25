#####################################################
#' nextPerm
#'
#' @details description, a paragraph
#' @param a xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
nextPerm=function(a)
#####################################################
{
# gives all permutations in lexicographic order
# call with loop as follows:
# a0=a;a=1:length(a0);res=nextPerm(a);a=res$a;print(a0[a])
# init
N=length(a)
k=1
if(N<=1) return(list(a=a,MORE=FALSE))
#
for (i in 2:(N))
    {
    if (a[i] > a[i-1]) k = i #find index of last element just next to the one to swap, keep in k
    }
# reverse partial order from k index until the end of a
revpa=rev(a[k:length(a)])
if (k>1) a=c(a[1:(k-1)],revpa)
#
if (k!=1) # more perms to do
   {
   s=k # find index of element to swap with a(k-1)
   while(s<N)
	{
	if (a[s] > a[k-1]) break;
	s=s+1
	}

   # swap a(k-1) et a(s)
   tmp=a[k-1]
   a[k-1]=a[s]
   a[s]=tmp
   #

   return(list(a=a,MORE=TRUE))
   }
   #if k=1 reached last permutation
    return(list(a=NULL,MORE=FALSE))
}
