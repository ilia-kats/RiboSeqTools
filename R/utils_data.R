#' Perform library-size normalization
#'
#' Normalize raw SeRP counts to library size. Each value is divided by the total number of reads
#' in the respective sample and multiplied by \eqn{10^6}.
#'
#' @param data A \code{serp_data} object.
#' @param exclude ORFs to exclude. Note that excluded ORFs still contribute to the total read number.
#' @return A \code{serp_data} object.
#' @rdname normalize
#' @export
normalize.serp_data <- function(data, exclude=c()) {
    stopifnot(!is_normalized(data))
    set_normalized(set_data(data, mapply(function(exp, texp) {
        mapply(function(rep, trep) {
            mapply(function(sample, tsample) {
                sapply(sample, function(bin) {
                    toinclude <- rownames(bin)
                    toinclude <- toinclude[!(toinclude %in% exclude)]
                    bin[toinclude,] / tsample * 1e6
                })
            }, rep, trep, SIMPLIFY=FALSE)
        }, exp, texp, SIMPLIFY=FALSE)
    }, get_data(data), get_total(data), SIMPLIFY=FALSE)), TRUE)
}

#' @export
normalize <- function(data, exclude=c()) {
    UseMethod("normalize")
}

#' Threshold data using the elbow method
#'
#' @param xvals X values
#' @param yvals Y values
#' @return The threshold, i.e. the value of the X axis where the distance between the Y value and
#'      the diagonal connecting the first and last points is maximal
#' @examples
#'      xvals <- 1:20
#'      yvals <- c(20:11, -0.2 * 1:10 + 11)
#'      get_elbow_threshold(xvals, yvals)
#' @export
get_elbow_threshold <- function(xvals, yvals) {
    if (missing(yvals) && is.matrix(xvals) && ncol(xvals) == 2) {
        yvals <- xvals[,2]
        xvals <- xvals[,1]
    } else if (missing(yvals)) {
        stop("y coordinates required")
    }
    v <- c(diff(range(xvals)), -diff(range(yvals)))
    w <- mapply(function(x, y, xmin, ymax)c(x - xmin, y - ymax), xvals, yvals, min(xvals), max(yvals))
    x <- apply(w, 2, function(w)sqrt(sum((w - as.vector(v %*% w) * v / sum(v^2))^2)))
    min(xvals) + w[1, which.max(x)]
}

#' Accessors for \code{serp_data} objects
#'
#' @param data A \code{serp_data} object.
#' @return \describe{
#'      \item{get_data}{Nested named list, with first level representing the experiment, second level
#'          the replicate, third level the sample type, and fourth level the binning. Count tables are
#'          sparse matrices with each row corresponding to an ORF.}
#'      \item{get_reference}{Reference data frame. Guaranteed to contain at least the following columns:
#'          \describe{
#'              \item{gene}{Gene/ORF name. Must match the names given in the read count tables.}
#'              \item{length}{ORF length in nucleotides.}
#'              \item{cds_length}{ORF length in codons.}}
#'          }
#'      \item{get_total_counts}{Nested named list containing toal read counts for each sample.}
#'      \item{get_defaults}{Named list of default parameters for this data set.}
#'      \item{is_normalized}{A logical value indicating whether the data object contains raw or normalized
#'          read counts.}
#'      }
#' @rdname serp_accessors
#' @export
get_data <- function(data) {
    check_serp_class(data)
    data$data
}

#' @rdname serp_accessors
#' @export
get_reference <- function(data) {
    check_serp_class(data)
    data$ref
}

#' @rdname serp_accessors
#' @export
get_total_counts <- function(data) {
    check_serp_class(data)
    data$total
}

#' @rdname serp_accessors
#' @export
get_defaults <- function(data) {
    check_serp_class(data)
    data$defaults
}

#' @rdname serp_accessors
#' @export
is_normalized <- function(data) {
    check_serp_class(data)
    data$normalized
}

#' Set default parameters of a \code{serp_data} object
#'
#' @param data A \code{serp_data} object.
#' @param ... Name-value-pairs of parameters
#' @export
set_defaults <- function(data, ...) {
    check_serp_class(data)
    data$defaults <- purrr::list_modify(data$defaults, ...)
    data
}

set_data <- function(data, newdata) {
    check_serp_class(data)
    data$data <- newdata
    data
}

set_normalized <- function(data, normalized) {
    check_serp_class(data)
    data$normalized <- normalized
    data
}
