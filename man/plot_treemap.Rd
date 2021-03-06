% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treemaps.R
\name{plot_treemap}
\alias{plot_treemap}
\title{Treemap plot of read counts per gene with an additional grouping variable}
\usage{
plot_treemap(data, exp, rep, sample, geneclass, title = "", palette = "Set2")
}
\arguments{
\item{data}{A \code{serp_data} object.}

\item{exp}{Experiment name.}

\item{rep}{Replicate name. If missing or NULL, read counts will be averaged over all replicates.}

\item{sample}{Sample type.}

\item{geneclass}{Data frame with columns \code{gene} and \code{class}. \code{class} must be an
ordered factor.}

\item{title}{Plot title.}

\item{palette}{RColorBrewer palette to color-code the class variable}
}
\value{
An invisible list from \code{\link[treemap]{treemap}}
}
\description{
This function tries really hard to make sure that the ordering of class groups in the plot
corresponds to the ordering of the factor levels. The treemap is plotted using \code{\link[treemap]{treemap}}.
}
\seealso{
\code{\link[treemap]{treemap}}
}
