% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_data.R
\name{normalize.serp_data}
\alias{normalize.serp_data}
\title{Perform library-size normalization}
\usage{
\method{normalize}{serp_data}(data)
}
\arguments{
\item{data}{A \code{serp_data} object.}
}
\value{
A \code{serp_data} object.
}
\description{
Normalize raw SeRP counts to library size. Each value is divided by the total number of reads
in the respective sample and multiplied by \eqn{10^6}.
}
