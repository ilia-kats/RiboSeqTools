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
    }, simplify=FALSE)
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

#' @export
print.serp_features <- function(data) {
    cat(sprintf('A ribosome profiling feature set with %d features\n', length(get_data(data))))
    print_list_name("", " .", "", get_data(data))
    cat('\n')
}

#' @export
`[.serp_features` <- function(data, i) {
    set_data(data, get_data(data)[i])
}

#' @export
c.serp_features <- function(...) {
    dat <- list(...)
    sapply(dat, check_serp_features_class)
    outdata <- get_data(dat[[1]])
    for (d in dat[-1]) {
        cdata <- get_data(d)
        cdnames <- names(cdata)

        for (i in cdnames) {
            if (i %in% names(outdata)) {
                for (j in names(d[[i]])) {
                    if (j %in% names(outdata[[i]])) {
                        if (d[[i]][[j]] != outdata[[i]][[j]])
                            rlang::abort(sprintf("Different features for %s %s", i, j))
                    } else {
                        outdata[[i]][[j]] <- d[[i]][[j]]
                    }
                }
            } else {
                x[[i]] <- y[[i]]
            }
        }
    }

    outref <- purrr::reduce(purrr::map(dat, get_reference), dplyr::union)
    outdefaults <- purrr::reduce(purrr::map(dat, get_defaults), combine_defaults)
    ret <- structure(list(ref=outref, data=outdat, defaults=list()), class='serp_features')
    set_defaults(ret, !!!outdefaults)
}
