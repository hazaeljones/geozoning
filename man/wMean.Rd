% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcCalNei.R
\name{wMean}
\alias{wMean}
\title{wMean}
\usage{
wMean(type, listZonePoint, surfVoronoi, data)
}
\arguments{
\item{type}{1-squared mean, 2-mean}

\item{listZonePoint}{list of data points belonging to zone}

\item{surfVoronoi}{areas of Voronoi polygon corresponding to data points}

\item{data}{SpatialPointsDataFrame}
}
\value{
a vector of mean zone values
}
\description{
wMean
}
\details{
computes weighted mean or squared mean of zone data
}
\examples{
data(mapTest)
data(resZTest)
K=resZTest
wMean(1,K$listZonePoint,mapTest$krigSurfVoronoi,mapTest$krigData)
}
