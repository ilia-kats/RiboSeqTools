#' @param data A \code{serp_data} object.
#' @param sample1 Name of the first sample (the numerator). If missing, the default sample1 of the data set
#'      will be used.
#' @param sample2 Name of the second sample (the denominator). If missing, the default sample2 of the data set
#'      will be used.
#' @param bin Binning mode (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
#'      will be used.
#' @param window_size Neighborhood size for the confidence interval calculation in nucleotides. If missing, the default
#'      window size of the data set will be used.
#' @param skip_5prime How many nucleotides to skip at the 5' end of the ORF. Useful if you know that the 5' end
#'      contains artifacts.
#' @param skip_3prime How many nucleotides to skip at the 3' end of the ORF. useful if you know that the 3' end
#'      contains artifacts.
#' @param conf.level Confidence level.
#' @param bpparam A \code{\link[BiocParallel]{BiocParallelParam-class}} object.
