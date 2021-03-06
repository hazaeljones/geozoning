% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{costLab}
\alias{costLab}
\title{costLab}
\usage{
costLab(K, map)
}
\arguments{
\item{K}{zoning object, as returned by the calNei function}

\item{map}{object returned by genMap function}
}
\value{
the sum of per label costs
}
\description{
costLab
}
\details{
description, a paragraph
}
\examples{
data(mapTest)
# run zoning with 2 quantiles corresponding to probability values 0.4 and 0.7,
# saving initial zoning and last level zonings
criti=correctionTree(c(0.4,0.7),mapTest,SAVE=TRUE)
K=criti$zk[[1]][[1]] # initial zoning
costLab(K,mapTest) #identical to criti$costL[[1]][[1]]
}
