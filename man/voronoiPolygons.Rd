% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcCalNei.R
\name{voronoiPolygons}
\alias{voronoiPolygons}
\title{voronoiPolygons}
\usage{
voronoiPolygons(spdata, gridLim = c(0, 1, 0, 1), neighBool,
  PTJUNCTION = FALSE, FULL = FALSE)
}
\arguments{
\item{spdata}{SpatialPointsDataFrame}

\item{gridLim}{list of boundary coordinates}

\item{neighBool}{empty point neighborhood Logical matrix}

\item{PTJUNCTION}{logical value, if FALSE (default): pts are not neighbors if their Voronoi polygons only have a vertex in common}

\item{FULL}{logical value, if FALSE (default): do not return Vornoi polygons}
}
\value{
a list with components
\describe{
\item{surfVoronoi}{Voronoi polygons areas}
\item{neighBool}{Voronoi point neighborhood Logical matrix}
if FULL=TRUE (warning: uses a lot of memory space), also:
\item{voronoi}{Voronoi polygons}
}
}
\description{
voronoiPolygons
}
\details{
determines the Voronoi neighborhood of data points
}
\examples{
\donttest{
 data(mapTest)
 rx=range(mapTest$krigData$x)
 ry=range(mapTest$krigData$y)
 nx=nrow(mapTest$krigGrid)
 ny=ncol(mapTest$krigGrid)
 nB=matrix(logical((nx*ny)^2),nx*ny,nx*ny) # big matrix
 vP=voronoiPolygons(mapTest$krigData,c(rx,ry),nB)
 length(vP$surfVoronoi) #as many as kriged data points
}
}
\seealso{
http://www.carsonfarmer.com/2009/09/voronoi-polygons-with-r/
}
