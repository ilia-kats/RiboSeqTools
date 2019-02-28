align_stop <- function(d, lengths) {
    Matrix::Matrix(t(mapply(function(r, l) {
        x <- d[r,]
        if (length(x) >= l) {
            x[(length(x) - l + 1):length(x)] <- x[1:l]
            x[1:(length(x) - l)] <- 0
        }
        x
    }, rownames(d), lengths[rownames(d)])), sparse=TRUE)
}

do_boot <- function(n, profilefun, mats, bpparam=BiocParallel::bpparam(), ...) {
    pars <- rlang::list2(...)
    res <- BiocParallel::bpmapply(function(...) {
        nsample <- unique(sapply(mats, nrow))
        stopifnot(length(nsample) == 1)
        s <- sample(nsample, nsample, replace=TRUE)
        m <- lapply(mats, function(x)x[s,,drop=FALSE])
        rlang::exec(profilefun, !!!m, !!!pars)
    }, 1:n, SIMPLIFY=FALSE, BPPARAM=bpparam)
    if (all(sapply(res, is.list))) {
        res <- transpose(res)
        res <- lapply(res, function(x)rlang::exec(rbind, !!!x, deparse.level=0))
    } else {
        res <- rlang::exec(rbind, !!!res, deparse.level=0)
    }
    res
}

mat_to_df <- function(mat, boot) {
    if (is.vector(mat)) {
        tibble::tibble(id=0, pos=1:length(mat), summary=mat, boot=boot)
    } else if (is.matrix(mat)) {
        r <- tibble::as_tibble(reshape2::melt(mat, varnames=c('id', 'pos'), value.name='summary'))
        r$boot <- boot
        r
    }
}

#' Align data matrices
#'
#' Centers each gene at a given position within the gene.
#'
#' @param data Names list of data matrices.
#' @param pos Vector of positions.
#' @param lengths Named vector of gene/ORF lengths.
#' @param pwidth Width of the final matrices.
#' @param filter Genes to include. Defaults to all genes.
#' @param binwidth Bin width, if binning is desired.
#' @param binmethod Binning method. \code{sum}: Read counts within each bin are summed up; \code{mean}:
#'      Read counts are averaged.
#' @return Named list of matrices. Matrices are of width \eqn{2\cdot \textrm{pwidth} + 1}, with \code{pos} in
#'      column \eqn{\textrm{pwidth} + 1}.
#' @export
make_aligned_mats <- function(data, pos, lengths, pwidth, filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), bpparam=BiocParallel::bpparam()) {
    binmethod <- match.arg(binmethod)
    binmethod <- switch(binmethod, sum=sum, mean=mean)
    BiocParallel::bplapply(data, function(d) {
        .filter <- filter
        stopifnot(is.null(.filter) || length(.filter) == length(pos))
        if (is.null(.filter))
            .filter <- rownames(d)
        .lengths <- lengths[.filter];
        t(mapply(function(g, p, l) {
            tlen <- 2 * pwidth + 1
            efflen <- ceiling(tlen / binwidth)
            out <- rep(NA_real_, efflen)
            nstart <- max(1, pwidth - p + 2)
            nstop <- 2 * pwidth + 1 - (p + pwidth - min(l, p + pwidth))
            if (binwidth > 1) {
                bins <- rle(1:tlen / tlen * efflen)
                starts <- c(1, cumsum(bins$lengths)[1:(length(bins$lengths))-1]+1)
                ends <- cumsum(bins$lengths)
                firstbin <- 1
                lastbin <- length(bins$values)
                if (pwidth > p)
                    firstbin <- floor((pwidth - p) / binwidth)
                if (p + pwidth > l)
                    lastbin <- lastbin - floor((p + pwidth - l) / binwidth)
                start <- max(1, p - pwidth)
                out[firstbin:lastbin] <- mapply(function(s, e)binmethod(d[g, s:e]), start + starts[firstbin:lastbin] - 1, start + ends[firstbin:lastbin] - 1)
            }
            else {
                out[nstart:nstop] <- d[g, max(1, p - pwidth):min(p + pwidth, l)]
            }
            out
        }, .filter, pos, .lengths))
    }, BPPARAM=bpparam)
}

metagene_profile <- function(d, profilefun, len, bin, refs, extrapars=list(), filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), normalizefun=NULL, align=c('start', 'stop'), nboot=100, bpparam=BiocParallel::bpparam()) {
    mats <- lapply(d, function(m) {
        m <- m[[bin]]
        cfilter <- if(is.null(filter)) TRUE else filter[filter %in% rownames(m)]
        m[cfilter,, drop=FALSE]
    })
    pars <- list(binwidth=binwidth, binmethod=binmethod, align=align)

    if (!is.null(normalizefun)) {
        mats <- rlang::exec(normalizefun, !!!mats, !!!extrapars, !!!pars)
    }

    if (length(align) == 1 && align %in% c('start', 'stop')) {
        if (align == 'stop')
            mats <- lapply(mats, align_stop, refs)

        len <- sapply(mats, function(m)min(ncol(m), len))

        if (binwidth > 1) {
            rval <- ifelse(binmethod == 'mean', 1 / binwidth, 1)
            if (length(unique(sapply(mats, ncol))) == 1) {
                len <- ceiling(len[1] / binwidth)

                r <- rep(c(rep(rval, binwidth), rep(0, ncol(mats[[1]]))), length.out=ncol(mats[[1]]) * len)
                if (align == 'stop')
                    r <- rev(r)
                r <- matrix(r, ncol=len, byrow=FALSE)
                mats <- lapply(mats, function(m) tidyr::replace_na(m, 0) %*% r)
            } else {
                mats <- mapply(function(m, l) {
                    len <- ceiling(l / binwidth)

                    r <- rep(c(rep(rval, binwidth), rep(0, ncol(m))), length.out=ncol(m) * len)
                    if (align == 'stop')
                        r <- rev(r)
                    tidyr::replace_na(m, 0) %*% matrix(r, ncol=len, byrow=FALSE)
                }, mats, len, SIMPLIFY=FALSE)
            }
        }

        mats <- mapply(function(m, len) {
            coords <- if (align == 'stop') (ncol(m) - len + 1):ncol(m) else 1:len
            as.matrix(m[,coords, drop=FALSE])
        }, mats, len, SIMPLIFY=FALSE)
    } else {
        mats <- make_aligned_mats(mats, align, refs, len, binwidth=1, binmethod='sum')
    }
    all <- rlang::exec(profilefun, !!!mats, !!!extrapars, !!!pars)
    boot <- rlang::exec(do_boot, n=nboot, profilefun=profilefun, mats=mats, bpparam=bpparam, !!!extrapars, !!!pars)

    if (!is.list(all) && !is.list(boot)) {
        all <- mat_to_df(all, FALSE)
        boot <- mat_to_df(boot, TRUE)
    } else {
        all <- lapply(all, mat_to_df, FALSE)
        boot <- lapply(boot, mat_to_df, TRUE)

        all <- bind_rows(all, .id='type')
        boot <- bind_rows(boot, .id='type')
    }
    d <- bind_rows(all, boot)
    d$pos <- d$pos - 1
    if (length(align) == 1 && align %in% c('start', 'stop')) {
        if (align == 'stop')
            d$pos <- d$pos - (len - 1)
        d$pos <- d$pos * binwidth
    } else {
        d$pos <- (d$pos - len) %/% binwidth * binwidth

        # this is faster than binning in make_aligned_mats
        bm <- switch(binmethod, sum=sum, mean=mean)
        d <- group_by_at(d, vars(-summary)) %>%
            summarize(summary=bm(summary)) %>%
            ungroup()
    }
    d
}

#' Calculate metagene profiles
#'
#' Calculates a metagene profile from the full data set as well as bootstrapping samples (sampling genes).
#' Profiles are calculated separately for each experiment and replicate.
#'
#' Data matrices are first filtered to contain only genes in \code{filter}. They are then passed to
#' \code{normalizefun} as named arguments, with names corresponding to sample type. If \code{align}
#' is one of \code{start} or \code{stop}, normalized data matrices are first aligned, then binned
#' and trimmed to \code{len} columns. If \code{align} is a vector of positions, data matrices are
#' centered using \code{\link{make_aligned_mats}} without binning for performance reasons, binning
#' is performed on the final profiles. The centered and trimmed data matrices are passed as named
#' arguments to \code{profilefun}, which calculates the final profiles. \code{profilefun} can return
#' either a single numeric vector, representing a single profile calculated from all data matrices,
#' or a named list of numeric vectors.
#'
#' @param data The data.
#' @param profilefun Function that calculates a profile. Must accept named arguments for all sample types
#'      present in the data set as well as \code{exp} (experiment name), \code{rep} (replicate name),
#'      \code{binwidth} (bin width), \code{binmethod} (binning method), \code{align} (alignment). Must return
#'      either a single numeric vector or a named list of numeric vectors.
#' @param len Length of the profile.
#' @param bin Bin level (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
#'      will be used
#' @param filter List of genes to include. Defaults to all genes.
#' @param binwidth Bin width.
#' @param binmethod How to bin the data. \code{sum}: Sums all read counts, \code{mean}: Averages read counts
#' @param normalizefun Function that normalizes the data pefore binning and profile calculation. Must accept
#'      the same arguments as \code{profilefun}. Must return a named list of matrices.
#' @param align Alignment of the metagene profile, one of \code{start} (5' end) or \code{stop} (3' end).
#'      Alternatively, a numeric vector of alignment positions for each gene in code{filter} can be given. In
#'      this case, genes will be aligned to the given positions and the profile will be centered at 0, spanning
#'      \code{len} positions in either direction.
#' @param nboot Number of bootstrap samples.
#' @param bpparam A \code{\link[BiocParallel]{BiocParallelParam-class}} object.
#' @return A \link[tibble]{tibble} with the following columns: \describe{
#'      \item{id}{ID of the bootstrap sample}
#'      \item{pos}{Distance from the \code{align} position. If \code{bin == 'byaa'} this is measured in codons,
#'          otherwise in nucleotides.}
#'      \item{boot}{Logical, indicates whether this profile was generated from a bootstrap sample or from the
#'          complete data set.}
#'      \item{type}{Only if \code{profilefun} returns a list. Corresponds to the name of the list element
#'          containing the profile.}
#'      \item{summary}{The value returned by \code{profilefun}.}
#'}
#' @seealso \link{defaults}
#' @export
metagene_profiles <- function(data, ...) {
    UseMethod("metagene_profiles")
}

#' @rdname metagene_profiles
#' @export
metagene_profiles.serp_data <- function(data, profilefun, len, bin, filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), normalizefun=NULL, align=c('start', 'stop'), nboot=100, bpparam=BiocParallel::bpparam()) {
    bin <- get_default_param(data, bin)
    binmethod <- match.arg(binmethod)
    refs <- setNames(get_reference(data)$length, get_reference(data)$gene)
    if (bin == 'byaa')
        refs <- refs / 3
    d <- mapply(function(d, exp) {
        ret <- mapply(function(d, rep) {
                        .filter <- filter
                        if (is.list(.filter)) {
                            .filter <- .filter[[exp]]
                            if (is.list(.filter))
                            .filter <- .filter[[rep]]
                        }
                        metagene_profile(d,
                                         profilefun,
                                         len,
                                         bin,
                                         refs,
                                         extrapars=list(exp=exp, rep=rep),
                                         filter=.filter,
                                         binwidth=binwidth,
                                         binmethod=binmethod,
                                         normalizefun=normalizefun,
                                         align=align,
                                         nboot=nboot,
                                         bpparam=bpparam)
                    },
                    d,
                    if(!is.null(names(d))) names(d) else 1:length(d),
                    SIMPLIFY=FALSE)
        bind_rows(ret, .id='rep')
    }, get_data(data), names(get_data(data)), SIMPLIFY=FALSE)
    d <- mutate(bind_rows(d, .id='exp'), id=as.factor(id))
    d
}

#' @rdname metagene_profiles
#' @export
metagene_profiles.serp_features <- function(data, profilefun, len, bin, filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), normalizefun=NULL, align=c('start', 'stop'), nboot=100, bpparam=BiocParallel::bpparam()) {
    bin <- get_default_param(data, bin)
    binmethod <- match.arg(binmethod)
    refs <- setNames(get_reference(data)$length, get_reference(data)$gene)
    if (bin == 'byaa')
        refs <- refs / 3
    metagene_profile(d,
                     profilefun,
                     len,
                     bin,
                     refs,
                     filter=.filter,
                     binwidth=binwidth,
                     binmethod=binmethod,
                     normalizefun=normalizefun,
                     align=align,
                     nboot=nboot,
                     bpparam=bpparam)
}

#' Plot metagene profiles
#'
#' Plots a metagene profile previously calculated with \code{\link{metagene_profiles}}. Lines show
#' the real profile using all data, shading indicates a bootstrapping-based confdience interval.
#'
#' @param df A data frame created by \code{\link{metagene_profiles}}.
#' @param ylab Y axis label
#' @param exp Experiments to plot. Defaults to all experiments.
#' @param colaes Variable to use for the color scale.
#' @param align Alignment to use in X axis label.
#' @template plot_annotations
#' @param conf.level Confidence level to plot.
#' @param ci.alpha Transparency level for the CI shading.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
plot_metagene_profiles <- function(df, ylab, exp=NULL, colaes=exp, align='start', highlightregion=list(), highlightargs=list(), conf.level=0.95, ci.alpha=0.3) {
    colaes <- rlang::enexpr(colaes)
    if (!is.null(exp))
        df <- dplyr::filter(df, exp %in% !!exp)
    p <- ggplot2::ggplot(df, ggplot2::aes(pos, summary, fill=!!colaes, color=!!colaes, group=interaction(rep, !!colaes))) +
        annotate_profile(highlightregion, !!!highlightargs) +
        ggplot2::stat_summary(ggplot2::aes(color=NULL), data=function(x)dplyr::filter(x, boot), geom='ribbon', fun.ymin=function(x)quantile(x, 0.5 * (1 - conf.level)), fun.ymax=function(x)quantile(x, 1 - 0.5 * (1 - conf.level)), alpha=ci.alpha) +
        ggplot2::geom_line(ggplot2::aes(y=summary), data=function(x)dplyr::filter(x, !boot)) +
        ggplot2::scale_y_continuous(trans='log2') +
        ggplot2::labs(x=sprintf('distance from %s / codons', align), y=ylab) +
        ggplot2::scale_x_continuous(expand=ggplot2::expand_scale()) +
        ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(alpha=1)), color=FALSE)
    p
}
