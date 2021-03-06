% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_data.R
\name{get_elbow_threshold}
\alias{get_elbow_threshold}
\title{Threshold data using the elbow method}
\usage{
get_elbow_threshold(xvals, yvals, upper_plateau = FALSE)
}
\arguments{
\item{xvals}{X values}

\item{yvals}{Y values}

\item{upper_plateau}{If \code{TRUE}, instead of connecting the first and last data points, the
left support is chosen such that no points are above the line. This is useful if there are
a few outliers with very low X values.}
}
\value{
The threshold, i.e. the value of the X axis where the distance between the Y value and
     the line connecting the first and last points is maximal.
}
\description{
Threshold data using the elbow method
}
\examples{
     xvals <- 1:20
     yvals <- c(20:11, -0.2 * 1:10 + 11)
     get_elbow_threshold(xvals, yvals)
}
