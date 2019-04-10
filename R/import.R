load_experiment <- function(..., bin=c('bynuc', 'byaa'), exclude=NULL) {
    bin <- match.arg(bin, several.ok=TRUE)
    ret <- mapply(function(...) {
        paths <- rlang::list2(...)
        ret <- sapply(paths, function(path) {
            ret <- list()
            m = as.matrix(read.table(path , fill=TRUE, header=TRUE, sep=",", row.names=1))
            acols <- ncol(m) %% 3
            if (acols > 0) {
                acols <- ncol(m):(ncol(m) - acols + 1)
                if (all(is.na(m[,acols])))
                m <- m[,-acols]
            }
            m[is.na(m)] <- 0
            libsize <- sum(m, na.rm=TRUE)
            if ('bynuc' %in% bin)
                ret$bynuc <- m
            if ('byaa' %in% bin) {
                ret$byaa <- m[,seq(1, ncol(m), by=3)] + m[,seq(2, ncol(m), by=3)] + m[,seq(3, ncol(m), by=3)]
            }
            if (!is.null(exclude) && length(exclude) > 0) {
                ret <- lapply(ret, function(x)Matrix::Matrix(x[!(rownames(x) %in% exclude),], sparse=TRUE))
            }
            else {
                ret <- lapply(ret, function(x)Matrix::Matrix(x, sparse=TRUE))
            }
            ret
        }, simplify=FALSE)
        ret
    },  ..., SIMPLIFY=FALSE)
    names(ret) <- 1:unique(purrr::map_int(list(...), length))
    ret
}

#' Import ribosome profiling data
#'
#' Import read count tables. Expects each count table to be a csv file with each row representing a different
#' ORF and the first column containing the ORF names. Other columns represent positions from the 5' end of of
#' the ORF in nucleotdes and must contain integer-valued read counts. A header must be present.
#'
#' @param ... Name-value pairs of lists. The name of each argument will be the name of an experiment.
#'      The name of each element will be the sample type (e.g. TT for total translatome), the value
#'      of each element must be a character vector of file paths, where each file is a read count
#'      table of one replicate experiment. Replicate order must match between sample types.
#' @param ref Reference data frame containing at least the following columns:
#'      \describe{
#'          \item{gene}{Gene/ORF name. Must match the names given in the read count tables.}
#'          \item{length}{ORF length in nucleotides.}
#'      }
#' @param normalize Normalize the read counts to library size? Output will then be in RPM.
#' @param bin Bin the data. \code{bynuc}: No binning (i.e. counts per nucleotide). \code{byaa}: Bin by residue.
#' @param exclude Genes to exclude in all future analyses. This genes will also be excluded from total read count
#'      calculation. Note that the raw count tables will not be modified. Named list with names corresponding to
#'      experiments.
#' @param defaults Default parameters of the data set.
#' @return An object of class \code{serp_data}
#' @seealso \code{\link{serp_data_accessors}}, \code{\link{defaults}}
#' @examples \dontrun{
#'      data <- load_serp(DnaK=list(ip=c('data/dnak1_ip.csv', 'data/dnak2_ip.csv'),
#'                                  tt=c('data/dnak1_tt.csv', 'data/dnak2_tt.csv')),
#'                        TF=list(ip=c('datatf1_ip.csv', 'data/tf2_ip.csv'),
#'                                tt=c('data/tf1_tt.csv', 'data/tf2_tt.csv')),
#'                        ref=reference_df,
#'                        bin='byaa')
#'      }
#' @export
load_serp <- function(..., ref, normalize=FALSE, bin=c('bynuc', 'byaa'), exclude=list(), defaults=list()) {
    experiments <- rlang::list2(...)
    stopifnot('gene' %in% colnames(ref) && 'length' %in% colnames(ref))
    bin <- match.arg(bin)
    data <- sapply(experiments, purrr::lift_dl(load_experiment), bin=bin, simplify=FALSE)
    what <- ifelse('bynuc' %in% bin, 'bynuc', 'byaa')

    ref$cds_length <- ref$length %/% 3
    ret <- list(ref=ref, data=data, normalized=FALSE)
    ret <- structure(ret, class="serp_data")
    ret <- set_excluded(ret, exclude)
    ret <- set_total_counts(ret, calc_total_counts(ret, what))
    if (normalize)
        ret <- normalize(ret)
    defaults <- purrr::list_modify(.defaults, !!!defaults)

    if (defaults$bin == 'byaa' && !('byaa' %in% bin))
        defaults$bin <- 'bynuc'
    ret$defaults <- defaults

    set_downsampled(ret, FALSE)
}
