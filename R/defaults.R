get_default_param <- function(serp_data, param, error=TRUE) {
    pname <- as.character(rlang::ensym(param))
    if (missing(param) || !rlang::env_has(rlang::caller_env(), nms=pname, inherit=TRUE)) {
        dflts <- get_defaults(serp_data)
        if (error && !(pname %in% names(dflts)))
            rlang::abort(sprintf("no value for %s", pname))
        else
            param <- dflts[[pname]]
    }
    param
}

.defaults <- list(bin='byaa',
                  window_size=45,
                  plot_ylim=c(2e-2, NA),
                  plot_ybreaks=scales::log_breaks(base=2))

#' Default parameters for a data set
#'
#' Named list of default values used in various functions. Usually, the default can be overriden by an argument
#' to the respective function. Currently, the following defaults are recognized:
#' \describe{
#'      \item{bin}{Bin level, one of \code{byaa} for codon-level binning or \code{bynuc} for nucleotide-level
#'          binning (corresponds to no binning).}
#'      \item{sample1}{Sample to use in the numerator of enrichment calculations.}
#'      \item{sample2}{Sample to use in the denominator of enrichment calculations.}
#'      \item{window_size}{Neighborhood size to use in enrichment CI calculations.}
#'      \item{plot_ylim}{Y axis limit for \link[=plot.serp_data]{enrichment plots}.}
#'      \item{plot_ybreaks}{Y axis breaks for \link[=plot.serp_data]{enrichment plots}.}
#'      \item{plot_fill_scale}{A \link[ggplot2:scale_fill_discrete]{ggplot2 fill scale} to use for
#'          \link[=plot.serp_data]{enrichment plots}. If this is a \link[ggplot2:scale_fill_manual]{manual scale},
#'          value names must correspond to experiment names.}
#' }
#'
#' @name defaults
NULL

