% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcCalCrit.R
\name{calCritMinMean}
\alias{calCritMinMean}
\title{calCritMinMean}
\usage{
calCritMinMean(matDistance, zoneNModif)
}
\arguments{
\item{matDistance}{zone distance matrix resulting from a call to calDistance}

\item{zoneNModif}{matrix of zone neigbors with FALSE on the diagonal}
}
\value{
a numerical value
}
\description{
calCritMinMean
}
\details{
computes a quality criterion equal to min(mean(dij^2/sqrt(dii^2*djj^2)))
see also  \code{\link{calCrit1}}, \code{\link{calCrit2}}, \code{\link{calCrit3}}, \code{\link{calCrit4}}, \code{\link{calCrit5}} for other criteria
}
\examples{
data(mapTest)
data(resZTest)
K=resZTest
resD = calDistance(typedist=1,mapTest$krigData,K$listZonePoint,K$zoneN,
       mapTest$krigSurfVoronoi,K$meanZone,pErr=0.9)
calCritMinMean(resD$matDistanceCorr,K$zoneNModif)
}
