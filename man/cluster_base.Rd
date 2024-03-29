% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_base.R
\name{cluster_base}
\alias{cluster_base}
\title{Cluster Base}
\usage{
cluster_base(
  ifolder = ".",
  cfiles = c("lpj_yields_rf", "lpj_yields_ir", "lpj_airrig", "transport_distance"),
  years2use = "y1995",
  spatial_header = NULL,
  use_cache = TRUE
)
}
\arguments{
\item{ifolder}{input folder where the MAgPIE input files are located}

\item{cfiles}{a vector containin the names of the MAgPIE input files
(beginning of the name is enough)}

\item{years2use}{A vector with years with should be taken into account for
the clustering}

\item{spatial_header}{A vector of the form c("REG.1","REG.2") (region name,
cell number) with entries for each spatial entity of the MAgPIE input files
which should be used to replace the names given in the inputs (required for
flexible region aggregation as here region names might change.). If set to
NULL the original information is used.}

\item{use_cache}{Read data from cache file if available (dangerous as changes
in settings will not be considered if an existing cache file is found).}
}
\value{
A matrix containing the data
}
\description{
Reads a series of MAgPIE files and combines them to a matrix which is then
used by mag_kmeans or mag_hierarchical for calculating a clustering
}
\details{
As this procedure is typically quite time consuming the function is a the
begining looking for a file cluster.Rdata in the input folder and is loading
the data from this file instead of cfiles if it exists. If it does not exist
data is read from cfiles but written to cluster.Rdata so that it is already
available in a second attempt.
}
\seealso{
\code{\link{mag_kmeans}}, \code{\link{mag_hierarchical}},
\code{\link{clusterspam}}
}
\author{
Jan Philipp Dietrich
}
