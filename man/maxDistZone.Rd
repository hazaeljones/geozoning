% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilZ.R
\name{maxDistZone}
\alias{maxDistZone}
\title{maxDistZone}
\usage{
maxDistZone(Z, iZ, k)
}
\arguments{
\item{Z}{zoning geometry (list of SpatialPolygons)}

\item{iZ}{current zone index}

\item{k}{polygon number within current zone}
}
\value{
the maximum distance within kth polygon of the current zone
}
\description{
maxDistZone
}
\details{
maximum distance within kth polygon of current zone
}
\examples{
\donttest{
ZK=initialZoning(qProb=c(0.4,0.2,0.7),mapTest)
Z=ZK$resZ$zonePolygone
geozoning:::maxDistZone(Z,5,1)
}
}
\keyword{internal}
