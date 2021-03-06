% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_enrichment.R
\name{plot.serp_data}
\alias{plot.serp_data}
\title{Plot a profile along a gene}
\usage{
\method{plot}{serp_data}(
  data,
  gene,
  type = c("enrichment", "rpm"),
  samples,
  sample1,
  sample2,
  exp,
  rep,
  bin,
  window_size,
  conf.level = 0.95,
  color = exp,
  highlightregion = list(),
  highlightargs = list()
)
}
\arguments{
\item{data}{A \code{serp_data} object. Must contain raw (unnormalized) read counts.}

\item{gene}{Name of the gene/ORF to plot.}

\item{type}{Plot type. One of \code{enrichment}, \code{rpm}.}

\item{samples}{Samples to plot. If missing, all samples will be plotted.}

\item{sample1}{Name of the first sample (the numerator). If missing, the default sample1 of the data set
will be used.}

\item{sample2}{Name of the second sample (the denominator). If missing, the default sample2 of the data set
will be used.}

\item{exp}{Character vector of experiments to plot. If missing, all experiments are plotted.}

\item{rep}{Character vector of replicates to plot. If missing, all replicates will be plotted.}

\item{bin}{Bin level (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
will be used.}

\item{window_size}{Neighborhood size for the confidence interval calculation in nucleotides. If missing, the default
window size of the data set will be used.}

\item{conf.level}{Confidence level.}

\item{highlightargs}{Named list of additional parameters to pass to \code{\link[ggplot2]{annotate}}.
Each list element will be recycled to match the number of annotations.}

\item{colaes}{Variable to use for the color scale.}

\item{higlightregion}{Regions in the plot to be highlighted. Can be either a list containing two-element
numeric vectors or a two-column matrix. Rectangles spanning the entire height of the plot will be
drawn between the given x positions.}
}
\value{
A \code{\link[ggplot2]{ggplot}} object.
}
\description{
For each position, a confidence interval is plotted. Transparency reflects the total number of reads
contributing to the confidence interval.
}
\details{
If \code{type} is \code{enrichment}, a binomial confidence interval for enrichment of \code{sample1}
compared to \code{sample2} is calculated using \code{\link{binom_ci_profile}}. The \code{samples}
argument is ignored.

If \code{type} is \code{samples}, a Poisson confidence interval for RPM in samples \code{samples}
is calculated using \code{\link{pois_ci_profile}}. The upper and lower bounds are then divided by
\code{window_size} to indicate confidence in the local smoothed read density. The \code{sample1}
and \code{sample2} arguments are ignored.
}
\seealso{
\code{\link{defaults}}
}
