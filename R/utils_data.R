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
    }, get_data(data), get_total_counts(data), SIMPLIFY=FALSE)), TRUE)
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
#' @rdname serp_data_accessors
#' @name serp_data_accessors
NULL

#' Accessors for \code{serp_features} objects
#'
#' @param data A \code{serp_features} object.
#' @return \describe{
#'      \item{get_data}{Nested named list, with first level representing the feature type and second level
#'          the binning. Data tables are matrices with each row corresponding to an ORF.}
#'      \item{get_reference}{Reference data frame. Guaranteed to contain at least the following columns:
#'          \describe{
#'              \item{gene}{Gene/ORF name. Must match the names given in the read count tables.}
#'              \item{length}{ORF length in nucleotides.}
#'              \item{cds_length}{ORF length in codons.}}
#'          }
#'      \item{get_defaults}{Named list of default parameters for this feature set.}
#'      }
#' @rdname serp_feature_accessors
#' @name serp_feature_accessors
NULL

#' @export
print.serp_data <- function(data) {
    cat(sprintf('A ribosome profiling data set with %d experiments\n', length(get_data(data))))
    cat(sprintf('Counts are normalized to library size: %s\n', ifelse(is_normalized(data), 'Yes', 'No')))
    print_list_name("", " .", "", get_data(data))
    cat('\n')
    cat('Defaults:\n')
    print(get_defaults(data))
}

#' @export
get_data <- function(data) {
    UseMethod("get_data")
}

#' @rdname serp_data_accessors
#' @export
get_data.serp_data <- function(data) {
    data$data
}

#' @rdname serp_feature_accessors
#' @export
get_data.serp_features <- function(data) {
    data$data
}

#' @export
get_reference <- function(data) {
    UseMethod("get_reference")
}

#' @rdname serp_data_accessors
#' @export
get_reference.serp_data <- function(data) {
    data$ref
}

#' @rdname serp_feature_accessors
#' @export
get_reference.serp_features <- function(data) {
    data$ref
}

#' @export
get_total_counts <- function(data) {
    UseMethod("get_total_counts")
}

#' @rdname serp_data_accessors
#' @export
get_total_counts.serp_data <- function(data) {
    data$total
}

#' @export
get_defaults <- function(data) {
    UseMethod("get_defaults")
}

#' @rdname serp_data_accessors
#' @export
get_defaults.serp_data <- function(data) {
    data$defaults
}

#' @rdname serp_feature_accessors
#' @export
get_defaults.serp_features <- function(data) {
    data$defaults
}

#' @export
is_normalized <- function(data) {
    UseMethod("is_normalized")
}

#' @rdname serp_data_accessors
#' @export
is_normalized.serp_data <- function(data) {
    data$normalized
}

#' Set default parameters
#'
#' @param data The data
#' @param ... Name-value-pairs of parameters
#' @export
set_defaults <- function(data, ...) {
    UseMethod("set_defaults")
}

#' @rdname set_defaults
#' @export
set_defaults.serp_data <- function(data, ...) {
    data$defaults <- purrr::list_modify(data$defaults, ...)
    data
}

#' @rdname set_defaults
#' @export
set_defaults.serp_features <- function(data, ...) {
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
