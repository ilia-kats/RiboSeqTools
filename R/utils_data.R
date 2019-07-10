#' Perform library-size normalization
#'
#' Normalize raw SeRP counts to library size. Each value is divided by the total number of reads
#' in the respective sample and multiplied by \eqn{10^6}.
#'
#' @param data A \code{serp_data} object.
#' @return A \code{serp_data} object.
#' @rdname normalize
#' @export
normalize.serp_data <- function(data) {
    stopifnot(!is_normalized(data))
    exclude <- excluded(data)
    data <- set_data(data, purrr::pmap(list(exp=get_data(data), texp=get_total_counts(data)[names(get_data(data))]), function(exp, texp) {
        purrr::pmap(list(rep=exp, trep=texp[names(exp)]), function(rep, trep) {
            purrr::pmap(list(sample=rep, tsample=trep[names(rep)]), function(sample, tsample) {
                sapply(sample, function(bin) {
                    bin / tsample * 1e6
                })
            })
        })
    }))
    what <- get_default_param(data, bin, error=FALSE)
    set_total_counts(data, calc_total_counts(data, what)) %>%
        set_normalized(TRUE)
}

#' @export
normalize <- function(data) {
    UseMethod("normalize")
}

#' Calculate total read counts per gene
#'
#' @param data A \code{serp_data} object.
#' @return A \link[tibble]{tibble} with columns \code{exp}, \code{rep}, \code{sample}, \code{counts}, and \code{RPM}
#' @export
get_genecounts <- function(data) {
    check_serp_class(data)
    if (is_normalized(data))
        rlang::abort("normalized data given")

    purrr::map2_dfr(get_data(data), get_total_counts(data)[names(get_data(data))], function(exp, texp) {
        purrr::map2_dfr(exp, texp[names(exp)], function(rep, trep) {
            purrr::map2_dfr(rep, trep[names(rep)], function(sample, tsample) {
                touse <- names(sample)[1]
                Matrix::rowSums(sample[[touse]]) %>%
                    tibble::enframe(name="gene", value="counts") %>%
                    dplyr::mutate(RPM=counts / tsample * 1e6)
            }, .id="sample")
        }, .id="rep")
    }, .id="exp")
}

#' Threshold data using the elbow method
#'
#' @param xvals X values
#' @param yvals Y values
#' @param upper_plateau If \code{TRUE}, instead of connecting the first and last data points, the
#'      left support is chosen such that no points are above the line. This is useful if there are
#'      a few outliers with very low X values.
#' @return The threshold, i.e. the value of the X axis where the distance between the Y value and
#'      the line connecting the first and last points is maximal.
#' @examples
#'      xvals <- 1:20
#'      yvals <- c(20:11, -0.2 * 1:10 + 11)
#'      get_elbow_threshold(xvals, yvals)
#' @export
get_elbow_threshold <- function(xvals, yvals, upper_plateau=FALSE) {
    if (missing(yvals) && is.matrix(xvals) && ncol(xvals) == 2) {
        yvals <- xvals[,2]
        xvals <- xvals[,1]
    } else if (missing(yvals)) {
        stop("y coordinates required")
    }
    if (upper_plateau) {
        lastx <- max(xvals)
        lasty <- min(yvals)
        idx <- which(!purrr::map2_lgl(xvals, yvals, function(x, y) {
            any(yvals > (lasty - y) / (lastx - x) * (xvals - x) + y)
        }))
        minx <- min(xvals[idx])
        idx <- which(xvals >= minx)
        xvals <- xvals[idx]
        yvals <- yvals[idx]
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
#'      \item{get_background_model}{Background model estimated by \code{\link{fit_background_model}}.}
#'      \item{get_binding_pvalues}{Binding p-values calculated by \code{\link{test_binding}}.}
#'      \item{is_normalized}{A logical value indicating whether the data object contains raw or normalized
#'          read counts.}
#'      }
#'      \item{is_downsampled}{A logical value indicating whether the data object has been downsampled to
#'          the lowest total read counts by \code{\link{downsample}}.}
#'      \item{excluded}{Character vector of genes excluded from analyses.}
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

    cat('Excluded genes:')
    if (length(excluded(data)) == 0) {
        cat(' None\n')
    } else {
        cat('\n  ', paste(purrr::imap(excluded(data), function(exc, exp)sprintf("%s: %s", exp, paste(exc, collapse=', '))), collapse='\n  '), '\n', sep='')
    }
    print_list_name("", " .", "", get_data(data))
    cat('\n')
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

#' @rdname serp_data_accessors
#' @export
is_downsampled <- function(data) {
    check_serp_class(data)
    data$downsampled
}

#' @export
excluded <- function(data) {
    UseMethod("excluded")
}

#' @rdname serp_data_accessors
#' @export
excluded.serp_data <- function(data) {
    data$exclude
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
    ndefs <- rlang::list2(...)
    if (is.null(data$defaults)) {
        data$defaults <- ndefs
    } else {
        data$defaults[names(ndefs)] <- ndefs
    }
    data
}

#' @rdname set_defaults
#' @export
set_defaults.serp_features <- function(data, ...) {
    ndefs <- rlang::list2(...)
    if (is.null(data$defaults)) {
        data$defaults <- ndefs
    } else {
        data$defaults[names(ndefs)] <- ndefs
    }
    data
}

#' @rdname serp_data_accessors
#' @export
get_background_model <- function(data) {
    check_serp_class(data)
    data$background_model
}

#' @rdname serp_data_accessors
#' @export
get_binding_pvalues <- function(data) {
    check_serp_class(data)
    data$binding_pvals
}

set_reference <- function(data, ref) {
    UseMethod("set_reference")
}

set_reference.serp_data <- function(data, ref) {
    data$ref <- ref
    data
}

set_background_model <- function(data, bgmodel) {
    check_serp_class(data)
    data$background_model <- bgmodel
    data
}

set_binding_pvalues <- function(data, pvals) {
    check_serp_class(data)
    data$binding_pvals <- pvals
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

set_downsampled <- function(data, downsampled) {
    check_serp_class(data)
    data$downsampled <- downsampled
    data
}

set_excluded <- function(data, exclude) {
    check_serp_class(data)
    data$excluded <- exclude
    data
}

set_total_counts <- function(data, newdata) {
    check_serp_class(data)
    data$total <- newdata
    data
}

calc_total_counts <- function(data, what=NULL) {
    check_serp_class(data)
    purrr::map2(get_data(data), excluded(data)[names(get_data(data))], function(exp, exc) {
        purrr::map(exp, function(rep) {
            purrr::map(rep, function(sample) {
                if (!is.null(what)) {
                    genes <- rownames(sample[[what]])
                    genes <- genes[!(genes %in% exc)]
                    sum(sample[[what]][genes,])
                } else {
                    median(unlist(purrr::map(sample, function(bintype) {
                        genes <- rownames(bintype)
                        genes <- genes[!(genes %in% exc)]
                        sum(bintype[genes,])
                    })))
                }
            })
        })
    })
}

#' @export
`[.serp_data` <- function(data, i) {
    data <- set_data(data, get_data(data)[i])
    data <- set_total_counts(data, get_total_counts(data)[i])
    data <- set_excluded(data, excluded(data)[i])

    bgmodel <- get_background_model(data)
    if (!is.null(bgmodel)) {
        bgmodel$model <- bgmodel$model[i]
        data <- set_background_model(data, bgmodel)
    }

    pvals <- get_binding_pvalues(data)
    if (!is.null(pvals))
        data <- set_binding_pvalues(data, dplyr::filter(pvals, exp %in% i))
    data
}

#' @export
c.serp_data <- function(...) {
    dat <- rlang::list2(...)
    sapply(dat, check_serp_class)
    nrml <- unique(sapply(dat, is_normalized))
    if (length(nrml) > 1)
        rlang::abort("Mixture of normalized and unnormalized data given.")
    dsmpld <- unique(sapply(dat, is_downsampled))
    if (length(dsmpld) > 1)
        rlang::abort("Mixture of downsampled and not downsampled data given.")
    outdat <- get_data(dat[[1]])
    outtotal <- get_total_counts(dat[[1]])
    outbgmodel <- NULL
    if (any(!sapply(dat, function(x)is.null(get_background_model(x)))))
        outbgmodel <- list()
    outpvals <- NULL
    if (any(!sapply(dat, function(x)is.null(get_binding_pvalues(x)))))
        outpvals <- tibble::tibble()
    outexclude <- excluded(dat[[1]])
    for (d in dat[-1]) {
        cdata <- get_data(d)
        ctotal <- get_total_counts(d)
        cbgmodel <- get_background_model(d)
        cpvals <- get_binding_pvalues(d)
        cdnames <- names(cdata)
        for (i in cdnames) {
            if (i %in% names(outdat)) {
                cdrepnames <- names(cdata[[i]])
                for (j in cdrepnames) {
                    if (j %in% names(outdat[[i]])) {
                        nrep <- list(cdata[[i]][[j]])
                        nreptotal <- list(ctotal[[i]][[j]])
                        newrepname <- max(names(outdat[[i]])) + 1
                        names(nrep) <- names(nreptotal) <- newrepname
                        outdat[[i]] <- c(outdat[[i]], nrep)
                        outtotal[[i]] <- c(outtotal[[i]], nreptotal)

                        if (!is.null(cbgmodel)) {
                            nrepbg <- list(cbgmodel[[i]][[j]])
                            names(nrepbg) <- newrepname
                            outbgmodel[[i]] <- c(outbgmodel[[i]], nrepbg)
                        }
                        if (!is.null(cpvals)) {
                            outpvals <- dplyr::filter(cpvals, exp == i, rep == j) %>%
                                mutate(rep=nrepname) %>%
                                bind_rows(outpvals)
                        }
                    } else {
                        outdat[[i]] <- c(outdat[[i]], cdata[[i]][j])
                        outtotal[[i]] <- c(outtotal[[i]], ctotal[[i]][[j]])
                        if (!is.null(cbgmodel)) {
                            outbgmodel[[i]] <- c(outbgmodel[[i]], cbgmodel[[i]][[j]])
                        }
                        if (!is.null(cpvals)) {
                            outpvals <- dplyr::bind_rows(outpvals, dplyr::filter(cpvals, exp == i, rep == j))
                        }
                    }
                }
            } else {
                outdat <- c(outdat, cdata[i])
                outtotal <- c(outtotal, ctotal[i])
                if (!is.null(cbgmodel))
                    outbgmodel <- c(outbgmodel, cbgmodel[i])
                if (!is.null(cpvals))
                    outpvals <- dplyr::bind_rows(outpvals, dplyr::filter(cpvals, exp == i))
            }

            if (i %in% names(outexclude)) {
                outexclude[[i]] <- union(outexclude[[i]], excluded(d)[[i]])
            } else if (!is.null(excluded(d)[[i]])) {
                outexclude[[i]] <- excluded(d)[[i]]
            }
        }
    }

    outref <- purrr::reduce(purrr::map(dat, get_reference), dplyr::union)
    outdefaults <- purrr::reduce(purrr::map(dat, get_defaults), combine_defaults)

    structure(list(defaults=list()), class='serp_data') %>%
        set_data(outdat) %>%
        set_reference(outref) %>%
        set_total_counts(outtotal) %>%
        set_normalized(nrml) %>%
        set_downsampled(dsmpld) %>%
        set_background_model(outbgmodel) %>%
        set_defaults(!!!outdefaults) %>%
        set_excluded(outexclude)
}

combine_defaults <- function(x, y) {
    out <- list()
    if (isTRUE(x$bin != y$bin))
        out$bin <- 'byaa'
    else
        out$bin <- x$bin

    if (isTRUE(x$sample1 == y$sample1))
        out$sample1 <- x$sample1
    if (isTRUE(x$sample2 == y$sample2))
        out$sample2 <- x$sample2

    if (isTRUE(x$window_size == y$window_size))
        out$window_size <- x$window_size

    if (isTRUE(x$plot_ylim == y$plot_ylim))
        out$plot_ylim <- x$plot_ylim
    else if (!is.null(x$plot_ylim) || !is.null(y$plot_ylim))
        out$plot_ylim <- range(x$plot_ylim, y$plot_ylim)

    out$plot_ybreaks <- dplyr::union(x$plot_ybreaks, y$plot_ybreaks)

    out
}

#' Downsample all samples to the same total counts
#'
#' The minimum total count of every sample and sample type is used as downsampling target.
#'
#' @param data A \code{serp_data} object
#' @return A \code{serp_data} object
#' @export
downsample <- function(data) {
    check_serp_class(data)
    if (is_normalized(data))
        rlang::abort("normalized data given")
    msizes <- purrr::map(get_total_counts(data), purrr::transpose) %>%
              purrr::transpose() %>%
              purrr::map(unlist) %>%
              purrr::map(min)
    ndata <- get_data(data)
    ndata <- purrr::map2(ndata, get_total_counts(data)[names(ndata)], function(exp, nexp) {
        purrr::map2(exp, nexp[names(exp)], function(rep, nrep) {
            purrr::pmap(list(sample=rep, nsample=nrep[names(rep)], samplename=names(rep)), function(sample, nsample, samplename) {
                msize <- msizes[[samplename]]
                if (nsample == msize)
                    sample
                else {
                    purrr::map(sample, function(bintype) {
                        Matrix(rbinom(length(bintype), as.vector(bintype), msize / nsample), nrow=nrow(bintype), ncol=ncol(bintype), dimnames=list(rownames(bintype), NULL))
                    })
                }
            })
        })
    })
    data <- set_data(data, ndata)
    set_total_counts(data, calc_total_counts(data)) %>%
        set_downsampled(TRUE)
}
