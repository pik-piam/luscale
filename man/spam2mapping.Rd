% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spam2mapping.R
\name{spam2mapping}
\alias{spam2mapping}
\title{Spam2Mapping converter}
\usage{
spam2mapping(spam, cellregions = NULL, clusterregions = NULL)
}
\arguments{
\item{spam}{spam object}

\item{cellregions}{a region vector with names cell of the spam object}

\item{clusterregions}{a region vector with names cluster of the spam object}
}
\value{
data.frame containing the mapping and weights
}
\description{
Converts a spam object into a mapping and in addition returns statistics about cells
and clusters per region
}
\examples{
spam <- spam::as.spam(matrix(c(0.2,0,0.8,0,0,1),2,3))
spam2mapping(spam,paste0("CE",1:3), paste0("CL",1:2))
}
\seealso{
\code{\link{clusterspam}}
}
\author{
Jan Philipp Dietrich
}
