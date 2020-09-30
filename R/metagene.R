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
        res <- purrr::transpose(res)
        res <- lapply(res, function(x)rlang::exec(rbind, !!!x, deparse.level=0))
    } else {
        res <- rlang::exec(rbind, !!!res, deparse.level=0)
    }
    res
}

mat_to_df <- function(mat, boot) {
    if (is.vector(mat) && length(mat) > 0) {
        tibble::tibble(id=0, pos=1:length(mat), summary=mat, boot=boot)
    } else if (is.matrix(mat) && all(dim(mat) > 0)) {
        r <- tibble::as_tibble(reshape2::melt(unname(mat), varnames=c('id', 'pos'), value.name='summary'))
        r$boot <- boot
        r
    } else {
        tibble::tibble()
    }
}

#' Generate the default function for metagene enrichment profiles
#'
#' Position-wise meta-enrichment is calculated as\eqn{E_i = \frac{\sum_g \textrm{sample1}_{gi}}{\sum_g \textrm{sample2}_{gi}}}{E_i = sum(sample[,i]) / sum(sample[,i])}
#' where \eqn{i} is the position and \eqn{g} the gene.
#'
#' @param data A \code{serp_data} object
#' @param sample1 Name of the first sample (the numerator). If missing, the default sample1 of the data set
#'      will be used.
#' @param sample2 Name of the second sample (the denominator). If missing, the default sample2 of the data set
#'      will be used.
#' @return A function that can be passed as \code{profilefun} to \code{\link{metagene_profiles}}
#' @seealso \link{defaults}, \code{\link{metagene_profiles}}, \code{\link{make_average_profilefun}},
#'          \code{\link{make_expression_normalizefun}}
#' @export
make_enrichment_profilefun <- function(data, sample1, sample2) {
    check_serp_class(data)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)

    fbody <- substitute(
        if (all(dim(sample1) > 0) && all(dim(sample2) > 0))
            colSums(sample1, na.rm=TRUE) / colSums(sample2, na.rm=TRUE)
        else
            numeric()
    , list(sample1=rlang::sym(sample1), sample2=rlang::sym(sample2)))

    profilefun <- function(){}
    args <- alist(sample1=, sample2=, ...=)
    names(args)[1] <- sample1
    names(args)[2] <- sample2
    formals(profilefun) <- args
    body(profilefun) <- fbody
    environment(profilefun) <- parent.env(environment())
    profilefun
}

make_function <- function(argnames, body) {
    args <- alist()
    for (arg in argnames) {
        tmp <- alist(tmp=)
        names(tmp) <- arg
        args <- c(args, tmp)
    }
    args <- c(args, alist(...=))
    fun <- function(){}
    formals(fun) <- args
    body(fun) <- as.call(c(as.name("{"), body))
    environment(fun) <- parent.env(environment())
    fun
}

make_average_profilefun_impl <- function(featurenames) {
    fbody <- sapply(featurenames, function(x)as.expression(substitute(ret[[x]] <- if(!missing(xx) && all(dim(xx) > 0)) colMeans(xx, na.rm=TRUE) else numeric(), env=list(x=x, xx=as.name(x)))))
    fbody <- c(expression(ret <- list()), fbody, expression(ret))
    make_function(featurenames, fbody)
}

#' Generate a default function for metagene profiles that computes the position-wise arithmetic mean.
#'
#' @param data The data
#' @return A function that can be passed as \code{profilefun} to \code{\link{metagene_profiles}}
#' @seealso \code{\link{metagene_profiles}}, \link{make_enrichment_profilefun},
#'          \code{\link{make_expression_normalizefun}}
#' @export
make_average_profilefun <- function(data) {
    UseMethod("make_average_profilefun")
}

#' @rdname make_average_profilefun
#' @export
make_average_profilefun.serp_data <- function(data) {
    make_average_profilefun_impl(samples(data))
}

#' @rdname make_average_profilefun
#' @export
make_average_profilefun.serp_features <- function(data) {
    make_average_profilefun_impl(features(data))
}

#' Generate the default function expression-normalized for metagene profiles
#'
#' Each gene is normalized to its total exression divided by its length.
#'
#' @param data A \code{serp_data} object
#' @param sample2 Name of the sample which all other samples will be normalized against. If missing,
#'      the default sample2 of the data set (usually the total translatome) will be used. If missing
#'      and no default is set, each sample will be normalized to its own expression.
#' @return A function that can be passed as \code{normalizefun} to \code{\link{metagene_profiles}}
#' @seealso \link{defaults}, \code{\link{metagene_profiles}}, \code{\link{make_average_profilefun}},
#'          \code{\link{make_enrichment_profilefun}}
#' @export
make_expression_normalizefun <- function(data, sample2) {
    check_serp_class(data)
    sample2 <- get_default_param(data, sample2, error=FALSE)
    samplenames <- samples(data)
    if (is.null(sample2)) {
        rlang::warn("Could not determine a total translatome sample. Each sample will be normalized to its own expression. Note: Since only genes with at least one read are included, the result may contain different genes in each sample, make certain that this is handled by profilefun.")
        fbody <- sapply(samplenames, function(x)as.expression(substitute(
            if(!missing(xx) && all(dim(xx) > 0)) {
                norm <- lengths / Matrix::rowSums(xx, na.rm=TRUE)
                idx <- which(is.finite(norm))
                ret[[x]] <- xx[idx]  * norm[idx]
            }
            , env=list(x=x, xx=as.name(x)))))
        fbody <- c(fbody, as.expression(substitute()))
    } else {
        rlang::inform(sprintf("Normalizing all samples to %s.", sample2))
        fbody <- sapply(samplenames, function(x)as.expression(substitute(
            if(!missing(xx) && all(dim(xx) > 0))
                ret[[x]] <- xx[idx,]  * norm
            , env=list(x=x, xx=as.name(x)))))
        fbody <- c(as.expression(substitute(norm <- lengths / Matrix::rowSums(tt, na.rm=TRUE), env=list(tt=as.name(sample2)))),
                   expression(idx <- which(is.finite(norm))),
                   expression(norm <- norm[idx]),
                   fbody)
    }
    fbody <- c(expression(ret <- list()), fbody, expression(ret))
    make_function(c(samplenames, "lengths"), fbody)
}

#' Align data matrices
#'
#' Centers each gene at a given position within the gene.
#'
#' @param data Names list of data matrices.
#' @param align Named vector of positions.
#' @param lengths Named vector of gene/ORF lengths. Names must correspond to names of \code{align}.
#' @param pwidth Width of the final matrices.
#' @param filter Genes to include. Defaults to all genes.
#' @param binwidth Bin width, if binning is desired.
#' @param binmethod Binning method. \code{sum}: Read counts within each bin are summed up; \code{mean}:
#'      Read counts are averaged.
#' @return Named list of matrices. Matrices are of width \eqn{2\cdot \textrm{pwidth} + 1}, with \code{align} in
#'      column \eqn{\textrm{pwidth} + 1}.
#' @export
make_aligned_mats <- function(data, align, lengths, pwidth, filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), bpparam=BiocParallel::bpparam()) {
    binmethod <- match.arg(binmethod)
    binmethod <- switch(binmethod, sum=sum, mean=mean)
    if (is.null(names(align))) {
        rlang::abort("align must be a named vector")
    }
    BiocParallel::bplapply(data, function(d) {
        .filter <- filter
        if (is.null(.filter))
            .filter <- rownames(d)
        .filter <- intersect(.filter, names(align))
        alidx <- names(align) %in% .filter # this allows to have the same gene multiple times in align
        .filter <- names(align)[alidx]
        align <- align[alidx]
        .lengths <- lengths[.filter];
        t(mapply(function(g, p, l) {
            tlen <- 2 * pwidth + 1
            efflen <- ceiling(tlen / binwidth)
            out <- rep(NA_real_, efflen)
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
                nstart <- max(1, pwidth - p + 2)
                nstop <- 2 * pwidth + 1 - (p + pwidth - min(l, p + pwidth))
                out[nstart:nstop] <- d[g, max(1, p - pwidth):min(p + pwidth, l)]
            }
            out
        }, .filter, align, .lengths))
    }, BPPARAM=bpparam)
}

metagene_profile <- function(d, profilefun, len, bin, refs, extrapars=list(), exclude=c(), filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), normalizefun=NULL, align='start', nboot=100, bpparam=BiocParallel::bpparam()) {
    genesintersect <- intersect(names(refs), purrr::reduce(lapply(d, function(x)rownames(x[[bin]])), intersect))
    cfilter <- if(is.null(filter)) genesintersect else filter[filter %in% genesintersect]
    if (!is.null(exclude))
        cfilter <- cfilter[!(cfilter %in% exclude)]
    refs <- refs[cfilter]
    maxlen <- max(refs)
    if (maxlen < 1)
        return(tibble())
    mats <- lapply(d, function(m) m[[bin]][cfilter, 1:maxlen, drop=FALSE])

    pars <- list(binwidth=binwidth, binmethod=binmethod, align=align, lengths=refs)

    if (!is.null(normalizefun)) {
        mats <- rlang::exec(normalizefun, !!!mats, !!!extrapars, !!!pars)
        if (!is.list(mats) || is.null(names(mats)))
            rlang::abort("normalizefun must return a named list of matrices")
    }

    if (length(align) == 1 && align %in% c('start', 'stop')) {
        if (align == 'stop')
            mats <- lapply(mats, align_stop, refs)

        len <- min(sapply(mats, function(m)min(ncol(m), len)))

        if (binwidth > 1) {
            rval <- ifelse(binmethod == 'mean', 1 / binwidth, 1)

            binmat <- function(m, r, len) {
                na_ind <- which(as.matrix(is.na(m)), arr.ind=TRUE)
                m[na_ind] <- 0
                m <- as.matrix(m %*% r)
                na_ind[,2] <- na_ind[,2] %/% binwidth
                m[na_ind[na_ind[,2] <= len %/% binwidth], ] <- NA_real_
                m
            }
            len <- ceiling(len[1] / binwidth)

            r <- rep(c(rep(rval, binwidth), rep(0, ncol(mats[[1]]))), length.out=ncol(mats[[1]]) * len)
            if (align == 'stop')
                r <- rev(r)
            r <- matrix(r, ncol=len, byrow=FALSE)
            mats <- lapply(mats, binmat, r, len)
        }

        mats <- lapply(mats, function(m, len) {
            coords <- if (align == 'stop') (ncol(m) - len + 1):ncol(m) else 1:len
            m <- m[,coords, drop=FALSE]
            if (align == 'start') {
                t(mapply(function(r, l) {
                    x <- m[r,]
                    if (length(x) > l %/% binwidth)
                        x[(l %/% binwidth + 1):length(x)] <- NA_real_
                    x
                }, rownames(m), refs[rownames(m)]))
            } else {
                t(mapply(function(r, l) {
                    x <- m[r,]
                    if (length(x) > l %/% binwidth)
                        x[1:(length(x) - l %/% binwidth)] <- NA_real_
                    x
                }, rownames(m), refs[rownames(m)]))
            }
        }, len)
    } else if (is.numeric(align)) {
        mats <- make_aligned_mats(mats, align, refs, len, binwidth=1, binmethod='sum', bpparam=bpparam)
    } else {
        rlang::abort("unrecognized align parameter")
    }
    all <- rlang::exec(profilefun, !!!mats, !!!extrapars, !!!pars)
    if (nboot > 0)
        boot <- rlang::exec(do_boot, n=nboot, profilefun=profilefun, mats=mats, bpparam=bpparam, !!!extrapars, !!!pars)

    if (!is.list(all) && !is.list(boot)) {
        all <- mat_to_df(all, FALSE)
        if (nboot > 0) boot <- mat_to_df(boot, TRUE)
    } else {
        all <- lapply(all, mat_to_df, FALSE)
        if (nboot > 0) boot <- lapply(boot, mat_to_df, TRUE)

        all <- dplyr::bind_rows(all, .id='type')
        if (nboot > 0) boot <- dplyr::bind_rows(boot, .id='type')
    }
    if (nboot > 0)
        d <- dplyr::bind_rows(all, boot)
    else
        d <- all
    d$pos <- d$pos - 1
    if (length(align) == 1 && align %in% c('start', 'stop')) {
        if (align == 'stop')
            d$pos <- d$pos - (len - 1)
        d$pos <- d$pos * binwidth
    } else {
        d$pos <- (d$pos - len) %/% binwidth * binwidth

        # this is faster than binning in make_aligned_mats
        bm <- switch(binmethod, sum=sum, mean=mean)
        d <- dplyr::group_by_at(d, vars(-summary)) %>%
            dplyr::summarize(summary=bm(summary)) %>%
            dplyr::ungroup()
    }
    d
}

#' Calculate metagene profiles
#'
#' Calculates a metagene profile from the full data set as well as bootstrapping samples (sampling genes).
#' Profiles are calculated separately for each experiment and replicate.
#'
#' Data matrices are first filtered to contain the same set of genes (i.e. only the intersection of
#' genes contained in all matrices is retained). For \code{serp_data} objects, this filtering is performed
#' separately for each experiment and replicate. Data matrices are the filtered to contain only genes in
#' \code{filter}, if \code{filter is given}. Matrices are then passed to \code{normalizefun}
#' as named arguments, with names corresponding to sample type. If \code{align} is one of \code{start}
#' or \code{stop}, normalized data matrices are first aligned, then binned and trimmed to \code{len}
#' columns. If \code{align} is a vector of positions, data matrices are centered using
#' \code{\link{make_aligned_mats}} without binning for performance reasons, binning
#' is performed on the final profiles. The centered and trimmed data matrices are passed as named
#' arguments to \code{profilefun}, which calculates the final profiles. \code{profilefun} can return
#' either a single numeric vector, representing a single profile calculated from all data matrices,
#' or a named list of numeric vectors.
#'
#' @param data The data.
#' @param profilefun Function that calculates a profile. Must accept named arguments for all sample types
#'      present in the data set (matrices with genes in rows cropped to the profile region.
#'      Missing values, e.g. if a gene ends within the profile region, are encoded as \code{NA}.)
#'      as well as \code{exp} (experiment name), \code{rep} (replicate name),
#'      \code{binwidth} (bin width), \code{binmethod} (binning method), \code{align} (alignment),
#'      \code{lengths} (named vector of gene lengths). Must return either a single numeric vector or
#'      a named list of numeric vectors.
#' @param len Length of the profile.
#' @param bin Bin level (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
#'      will be used
#' @param filter List of genes to include. Defaults to all genes.
#' @param binwidth Bin width.
#' @param binmethod How to bin the data. \code{sum}: Sums all read counts, \code{mean}: Averages read counts
#' @param normalizefun Function that normalizes the data pefore binning and profile calculation. Must accept
#'      the same arguments as \code{profilefun}. Data matrices will contain data for entire genes.
#'      The \code{lengths} argument is guaranteed to contain only genes present in the data matrices
#'      in the same order.
#'      Must return a named list of matrices.
#' @param align Alignment of the metagene profile, one of \code{start} (5' end) or \code{stop} (3' end).
#'      Alternatively, a named numeric vector of alignment positions for each gene in \code{filter} can be given. In
#'      this case, names must correspond to gene names and genes will be aligned to the given positions
#'      and the profile will be centered at 0, spanning \code{len} positions in either direction.
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
#' @seealso \link{defaults}, \code{\link{make_average_profilefun}},
#'          \code{\link{make_enrichment_profilefun}}, \code{\link{make_expression_normalizefun}}
#' @export
metagene_profiles <- function(data, ...) {
    UseMethod("metagene_profiles")
}

#' @rdname metagene_profiles
#' @export
metagene_profiles.serp_data <- function(data, profilefun, len, bin, filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), normalizefun=NULL, align='start', nboot=100, bpparam=BiocParallel::bpparam()) {
    bin <- get_default_param(data, bin)
    binmethod <- match.arg(binmethod)
    ref <- get_reference(data)
    refs <- setNames(ref$length, ref$gene)
    exclude <- excluded(data)
    if (bin == 'byaa')
        refs <- refs / 3
    if (is.numeric(align))
        names(align) <- map_genenames(data, names(align))
    d <- purrr::imap_dfr(get_data(data), function(d, exp) {
        purrr::imap_dfr(d, function(d, rep) {
            .filter <- filter
            if (is.list(.filter)) {
                .filter <- .filter[[exp]]
                if (is.list(.filter))
                    .filter <- .filter[[rep]]

            }
            .filter <- map_genenames(data, .filter)

            metagene_profile(d,
                             profilefun,
                             len,
                             bin,
                             refs,
                             extrapars=list(exp=exp, rep=rep),
                             exclude=exclude[[exp]],
                             filter=.filter,
                             binwidth=binwidth,
                             binmethod=binmethod,
                             normalizefun=normalizefun,
                             align=align,
                             nboot=nboot,
                             bpparam=bpparam)
        }, .id='rep')
    }, .id='exp')
    if (nrow(d) > 0)
        d <- dplyr::mutate(d, id=as.factor(id))
    d
}

#' @rdname metagene_profiles
#' @export
metagene_profiles.serp_features <- function(data, profilefun, len, bin, filter=NULL, binwidth=1, binmethod=c('sum', 'mean'), normalizefun=NULL, align='start', nboot=100, bpparam=BiocParallel::bpparam()) {
    bin <- get_default_param(data, bin)
    binmethod <- match.arg(binmethod)
    refs <- setNames(get_reference(data)$length, get_reference(data)$gene)
    if (bin == 'byaa')
        refs <- refs %/% 3
    metagene_profile(get_data(data),
                     profilefun,
                     len,
                     bin,
                     refs,
                     filter=filter,
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
#' @param ylab Y axis label.
#' @param ytrans Y axis transformation.
#' @param exp Experiments to plot. Defaults to all experiments.
#' @param color Variable to use for the color scale. Defaults to \code{exp}, if a column \code{exp}
#'      exists in the data frame.
#' @param group Additional grouping variable. Geoms will be grouped by \code{color} and \code{group}.
#'      Defaults to \code{rep}, if a column \code{rep} exists in the data frame.
#' @param align Alignment to use in X axis label.
#' @template plot_annotations
#' @param conf.level Confidence level to plot.
#' @param ci.alpha Transparency level for the CI shading.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @examples
#' \dontrun{
#'      plot_metagene_profiles(df, 'enrichment', exp, rep)
#' }
#' @export
plot_metagene_profiles <- function(df, ylab, color, group, align='start', ytrans='identity', highlightregion=list(), highlightargs=list(), conf.level=0.95, ci.alpha=0.3) {
    pars <- list()
    interactions <- list()
    col <- rlang::enexpr(color)
    grp <- rlang::enexpr(group)
    if (!rlang::is_missing(col)) {
        pars$fill <- pars$color <- col
        interactions <- c(interactions, col)
    } else if ('exp' %in% colnames(df)) {
        col <- rlang::expr(exp)
        pars$fill <- pars$color <- col
        interactions <- c(interactions, col)
    }
    if (!rlang::is_missing(grp)) {
        interactions <- c(interactions, grp)
    } else if ('rep' %in% colnames(df)) {
        grp <- rlang::expr(rep)
        interactions <- c(interactions, grp)
    }
    if (length(interactions) > 0)
        pars$group <- rlang::expr(interaction(!!!interactions))

    p <- ggplot2::ggplot(df, ggplot2::aes(pos, summary, !!!pars)) +
        annotate_profile(highlightregion, !!!highlightargs) +
        ggplot2::stat_summary(ggplot2::aes(color=NULL), data=function(x)dplyr::filter(x, boot), geom='ribbon', fun.ymin=function(x)quantile(x, 0.5 * (1 - conf.level)), fun.ymax=function(x)quantile(x, 1 - 0.5 * (1 - conf.level)), alpha=ci.alpha) +
        ggplot2::geom_line(ggplot2::aes(y=summary), data=function(x)dplyr::filter(x, !boot)) +
        ggplot2::scale_y_continuous(trans=ytrans) +
        ggplot2::labs(x=sprintf('distance from %s / codons', align), y=ylab) +
        ggplot2::scale_x_continuous(expand=ggplot2::expand_scale()) +
        ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(alpha=1)), color=FALSE)
    p
}
