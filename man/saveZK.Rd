% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveZK.R
\name{saveZK}
\alias{saveZK}
\title{saveZK function called by correctionTree}
\usage{
saveZK(map, K1, Z2, qProb, listOfZ, indCur, crit, cost, costL, nz, mdist,
  pErr = 0.9, optiCrit = 2, simplitol = 0.001)
}
\arguments{
\item{map}{object returned by function genMap}

\item{K1}{previous zoning}

\item{Z2}{current zoning geometry (list of SpatialPolygons)}

\item{qProb}{probability vector used to generate quantile values}

\item{listOfZ}{list of zoning objects}

\item{indCur}{index of new list element}

\item{crit}{list of criteria}

\item{cost}{list of costs}

\item{costL}{list of per label costs}

\item{nz}{list of number of zones}

\item{mdist}{list of distance matrices}

\item{pErr}{equality tolerance for distance calculations}

\item{optiCrit}{criterion choice}

\item{simplitol}{tolerance for spatial polygons geometry simplification}
}
\value{
a  list with components
\describe{
\item{listOfZ}{updated list of zoning objects, first element corresponds to initial zoning, each other element is a list with each (last if ALL=FALSE) level zoning objects}
\item{mdist}{list of initial distance matrix and all (last if ALL=FALSE) level distance matrices}
\item{crit}{list of initial criterion and all (last if ALL=FALSE) level criteria }
\item{cost}{list of initial cost and all (last if ALL=FALSE) level costs }
\item{costL}{list of initial cost per label and all (last if ALL=FALSE) level costs per label}
\item{nz}{list of initial number of zones and all (last if ALL=FALSE) level number of zones}
}
}
\description{
saveZK function called by correctionTree
}
\details{
Given a map object, a list of zonings, a current and a previous zoning, adds the current zoning to the list of zonings if it has at least 2 zones,after recalculating zone neighborhood and transferring zone labels.
}
\examples{
data(mapTest)
# run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7
criti=correctionTree(c(0.4,0.7),mapTest,LASTPASS=FALSE,SAVE=TRUE)
K1=criti$zk[[1]][[1]]#initial zoning
Z1=K1$zonePolygone
printZsurf(Z1) # 8 zones with 2 small zones (7 and 8)
Z2 = geozoning:::remove1FromZ(Z1,7,K1$zoneN)
printZsurf(Z2) #7 zones
indCur=2
newRes=geozoning:::saveZK(mapTest,K1,Z2,c(0.4,0.7),criti$zk,indCur,
        criti$criterion,criti$cost,criti$costL,criti$nz,criti$mdist)
newZ=newRes$listOfZ[[2]][[1]]$zonePolygone
printZsurf(newZ) #6 zones
}
\keyword{internal}
