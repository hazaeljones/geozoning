% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilZ.R
\name{setIds}
\alias{setIds}
\title{setIds}
\usage{
setIds(Z)
}
\arguments{
\item{Z}{zoning geometry (list of SpatialPolygons)}
}
\value{
a zoning geometry
}
\description{
setIds
}
\details{
set all zone identifiers in a zoning by assigning zone number to each identifier.
}
\examples{
\donttest{
data(mapTest)
criti=correctionTree(c(0.4,0.5),mapTest,SAVE=TRUE)
Z=criti$zk[[2]][[1]]$zonePolygon
Z1=setIds(Z)
}
}
