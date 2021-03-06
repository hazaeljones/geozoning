% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{getPolySp}
\alias{getPolySp}
\title{getPolySp}
\usage{
getPolySp(sp, k = 1)
}
\arguments{
\item{sp}{SpatialPolygons object}

\item{k}{polygon number}
}
\value{
a polygon (object of class Polygon)
}
\description{
getPolySp
}
\details{
get the kth polygon of the current SpatialPolygons
}
\examples{
ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
Z=ZK$resZ$zonePolygone
sp=Z[[5]]
P1=geozoning:::getPolySp(sp,1)
P2=geozoning:::getPolySp(sp,2) # second polygon is a hole
plot(P1@coords,type="l")
lines(P2@coords,type="l",col="blue")
}
\keyword{internal}
