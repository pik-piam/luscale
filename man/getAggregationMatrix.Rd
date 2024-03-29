% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAggregationMatrix.R
\name{getAggregationMatrix}
\alias{getAggregationMatrix}
\title{getAggregationMatrix}
\usage{
getAggregationMatrix(rel, from=NULL, to=NULL)
}
\arguments{
\item{rel}{file name of a region mapping (".csv" as file ending, ";" used as
separator and with columns containing the region names) or a mapping matrix
(matrix with region names in columns)}

\item{from}{Name of the first column to be used (if not set the first or
second column will be used)}

\item{to}{Name of the second column to be used (if not set the second or
third column will be used)}
}
\value{
A matrix nregions1 x nregions2 with 1s and 0s showing the mapping of
countries to regions.
}
\description{
Function which converts the supplied regionmapping file to a transformation
matrix which then can be used for aggregation with
\code{\link{speed_aggregate}}.
}
\examples{
\dontrun{
  x <- cbind(reg=c("REG1","REG2","REG2","REG3"),country=c("C1","C2","C3","C4")) 

  getAggregationMatrix(x)
  
  getAggregationMatrix(x,from="reg",to="reg")
}
}
\author{
Jan Philipp Dietrich
}
