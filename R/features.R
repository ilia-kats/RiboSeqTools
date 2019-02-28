#' Create a \code{serp_features} object
#'
#' This class holds position-wise features that are independent of experiment and replicate
#'
#' @param ... Name-value pairs of matrices. The name of each argument denotes the feature type,
#'      the value is a matrix with each row corresponding to a gene/ORF.
#' @param ref ref Reference data frame containing at least the following columns:
#'      \describe{
#'          \item{gene}{Gene/ORF name. Must match the names given in the read count tables.}
#'          \item{length}{ORF length in nucleotides.}
#'      }
#' @param bin Whether the features are per codon (\code{byaa}) or nucleotide (\code{bynuc}).
#' @param defaults Default parameters for the feature set.
#' @return An object of class \code{serp_features}
#' @seealso \link{serp_feature_accessors}, \link{defaults}
#' @export
serp_features <- function(..., ref, bin=c('bynuc', 'byaa'), defaults=list()) {
    features <- rlang::list2(...)
    stopifnot('gene' %in% colnames(ref) && 'length' %in% colnames(ref))
    bin <- match.arg(bin)

    data <- sapply(features, function(f) {
        if (!is.list(f)) {
            ret <- list()
            ret[[bin]] <- f
            ret
        } else {
            stopifnot(all(names(f) %in% c('bynuc', 'byaa')))
            f
        }
    })
    ref$cds_length <- ref$length %/% 3
    ret <- list(ref=ref, data=data)
    ret <- structure(ret, class='serp_features')

    # TODO: extra defaults for features
    defaults <- purrr::list_modify(.defaults, !!!defaults)
    if (defaults$bin == 'byaa' && !('byaa' %in% bin))
        defaults$bin <- 'bynuc'
    ret$defaults <- defaults

    ret
}
