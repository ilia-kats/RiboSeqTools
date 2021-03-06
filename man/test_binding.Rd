% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binding_scores.R
\name{test_binding}
\alias{test_binding}
\title{Statistical test for position-wise enrichment in selective ribosome profiling data}
\usage{
test_binding(data, window_size, bpparam = BiocParallel::bpparam())
}
\arguments{
\item{data}{A \code{serp_data} object. \code{\link{fit_background_model}} must have been run on the data.}

\item{window_size}{Neighborhood size for the confidence interval calculation in nucleotides. If missing, the default
window size of the data set will be used.}

\item{bpparam}{A \code{\link[BiocParallel]{BiocParallelParam-class}} object.}
}
\value{
A \code{serp_data} object.
}
\description{
Using the background model estimated by \code{\link{fit_background_model}}, performs significance testing for enrichment of reads in
\code{sample1} compared to \code{sample2} for all genes and positions.
}
\details{
Enrichment spikes at individual positions are likely not indicative of real chaperone binding. Nevertheless, they frequently occur,
probably due to ligation and/or sequencing bias. To alleviate this problem somewhat, read counts within a \code{window_size}-wide
neighborhood are summed up individually for \code{sample1} and \code{sample2} and used for significance testing. This results in
strongly correlatd p-values, therefore the Benjamini-Yekutieli-procedure is used for FDR control. \code{sample1} and \code{sample2}
from \code{\link{fit_background_model}} are used. Note that the binning mode (bynuc or byaa) previously given to \code{\link{fit_background_model}}
will be used.
}
\seealso{
\code{\link{get_binding_pvalues}}, \code{\link{fit_background_model}}, \code{\link{get_binding_positions}}, \code{\link{plot_binding_positions}}
}
