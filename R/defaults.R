get_default_param <- function(serp_data, param, error=TRUE) {
    pname <- as.character(rlang::ensym(param))
    if (missing(param)) {
        param <- serp_data$defaults[[pname]]
    }
    if (is.null(param) && error)
        stop(sprintf("invalid %s argument", pname))
    param
}

.defaults <- list(bin='byaa',
                  window_size=45,
                  plot_ylim=2^c(-2.8,6.8),
                  plot_ybreaks=2^(-2:6))

#' Default parameters for a data set
#'
#' Named list of default values used in various functions. Usually, the default can be overriden by an argument
#' to the respective function. Currently, the following defaults are recognized:
#' \describe{
#'      \item{bin}{Bin level, one of \code{byaa} for codon-level binning or \code{bynuc} for nucleotide-level
#'          binning (corresponds to no binning).}
#'      \item{sample1}{Sample to use in the nominator of enrichment calculations.}
#'      \item{sample2}{Sample to use in the denumerator of enrichment calculations.}
#'      \item{window_size}{Neighborhood size to use in enrichment CI calculations.}
#'      \item{plot_ylim}{Y axis limit for \link[=plot.serp_data]{enrichment plots}.}
#'      \item{plot_ybreaks}{Y axis breaks for \link[=plot.serp-data]{enrichment plots}.}
#' }
#'
#' @name defaults
NULL

