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

do_boot <- function(n, profilefun, mats, ...) {
    pars <- list(...)
    res <- BiocParallel::bpmapply(function(...) {
        nsample <- unique(sapply(mats, nrow))
        stopifnot(length(nsample) == 1)
        s <- sample(nsample, nsample, replace=TRUE)
        m <- lapply(mats, function(x)x[s,,drop=FALSE])
        rlang::exec(profilefun, !!!m, !!!pars)
    }, 1:n, SIMPLIFY=FALSE)
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
        tibble::tibble(id=1, pos=1:length(mat), counts=mat, boot=boot)
    } else if (is.matrix(mat)) {
        r <- tibble::as_tibble(reshape2::melt(mat, varnames=c('id', 'pos'), value.name='counts'))
        r$boot <- boot
        r
    }
}

#' @export
make_aligned_mats <- function(data, pos, lengths, pwidth, filter=NULL, bin=1, binmethod=c('sum', 'mean')) {
    binmethod <- match.arg(binmethod)
    binmethod <- switch(binmethod, sum=sum, mean=mean)
    BiocParallel::bplapply(data, function(d) {
        .filter <- filter
        stopifnot(is.null(.filter) || length(.filter) == length(pos))
        if (is.null(.filter))
            .filter <- rownames(d)
        .lengths <- lengths[filter];
        t(mapply(function(g, p, l) {
            tlen <- 2 * pwidth + 1
            efflen <- ceiling(tlen / .bin)
            out <- rep(NA_real_, efflen)
            nstart <- max(1, pwidth - p + 2)
            nstop <- 2 * pwidth + 1 - (p + pwidth - min(l, p + pwidth))
            if (bin > 1) {
                bins <- rle(1:tlen / tlen * efflen)
                starts <- c(1, cumsum(bins$lengths)[1:(length(bins$lengths))-1]+1)
                ends <- cumsum(bins$lengths)
                firstbin <- 1
                lastbin <- length(bins$values)
                if (pwidth > p)
                    firstbin <- floor((pwidth - p) / bin)
                if (p + pwidth > l)
                    lastbin <- lastbin - floor((p + pwidth - l) / bin)
                start <- max(1, p - pwidth)
                out[firstbin:lastbin] <- mapply(function(s, e)binmethod(d[g, s:e]), start + starts[firstbin:lastbin] - 1, start + ends[firstbin:lastbin] - 1)
            }
            else {
                out[nstart:nstop] <- d[g, max(1, p - pwidth):min(p + pwidth, l)]
            }
            out
        }, .filter, pos, .lengths))
    })
}

#' @export
metagene_profiles <- function(data, profilefun, len, filter=NULL, bin=1, binmethod=c('sum', 'mean'), normalizefun=NULL, align=c('start', 'stop'), nboot=100, what='byaa') {
    check_serp_class(data)
    binmethod <- match.arg(binmethod)
    refs <- setNames(get_reference(data)$length, get_reference(data)$gene)
    if (what == 'byaa')
        refs <- refs / 3
    d <- mapply(function(d, exp) {
        ret <- mapply(function(d, rep) {
            .filter <- filter
            if (is.list(.filter)) {
                .filter <- .filter[[exp]]
                if (is.list(.filter))
                    .filter <- .filter[[rep]]
            }

            mats <- lapply(d, function(m) {
                m <- m[[what]]
                cfilter <- if(is.null(.filter)) TRUE else .filter[.filter %in% rownames(m)]
                m[cfilter,, drop=FALSE]
            })
            pars <- list(exp=exp, rep=rep, bin=bin, binmethod=binmethod, align=align)

            if (!is.null(normalizefun)) {
                mats <- rlang::exec(normalizefun, !!!mats, !!!pars)
            }

            if (length(align) == 1 && align %in% c('start', 'stop')) {
                if (align == 'stop')
                    mats <- lapply(mats, align_stop, refs)

                len <- sapply(mats, function(m)min(ncol(m), len))

                if (bin > 1) {
                    rval <- ifelse(binmethod == 'mean', 1 / bin, 1)
                    if (length(unique(sapply(mats, ncol))) == 1) {
                        len <- ceiling(len[1] / .bin)

                        r <- rep(c(rep(rval, bin), rep(0, ncol(mats[[1]]))), length.out=ncol(mats[[1]]) * len)
                        if (align == 'stop')
                            r <- rev(r)
                        r <- matrix(r, ncol=len, byrow=FALSE)
                        mats <- lapply(mats, function(m) tidyr::replace_na(m, 0) %*% r)
                    } else {
                        mats <- mapply(function(m, l) {
                            len <- ceiling(l / bin)

                            r <- rep(c(rep(rval, bin), rep(0, ncol(m))), length.out=ncol(m) * len)
                            if (align == 'stop')
                                r <- rev(r)
                            tidyr::replace_na(m, 0) %*% matrix(r, ncol=len, byrow=FALSE)
                        }, mats, len, SIMPLIFY=FALSE)
                    }
                }

                mats <- mapply(function(m, len) {
                    coords <- if (align == 'stop') (ncol(m) - len + 1):ncol(m) else 1:len
                    m[,coords, drop=FALSE]
                }, mats, len, SIMPLIFY=FALSE)
            } else {
                mats <- make_aligned_mats(mats, align, refs, len, bin=1, binmethod='sum')
            }
            all <- rlang::exec(profilefun, !!!mats, !!!pars)
            boot <- rlang::exec(do_boot, n=nboot, profilefun=profilefun, mats=mats, !!!pars)

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
                d$pos <- d$pos * bin
            } else {
                d$pos <- (d$pos - len) %/% bin * bin

                # this is faster than binning in make_aligned_mats
                bm <- switch(binmethod, sum=sum, mean=mean)
                d <- group_by_at(d, vars(-counts)) %>%
                    summarize(counts=bm(counts)) %>%
                    ungroup()
            }
            d
        }, d, if(!is.null(names(d))) names(d) else 1:length(d), SIMPLIFY=FALSE)
        bind_rows(ret, .id='rep')
    }, get_data(data), names(get_data(data)), SIMPLIFY=FALSE)
    d <- mutate(bind_rows(d, .id='exp'), id=as.factor(id))
    d
}

#' @export
plot_metagene_profiles <- function(df, exp=NULL, colaes=exp, align='start', highlightregion=list(), conf.level=0.95) {
    colaes <- enexpr(colaes)
    if (!is.null(exp))
        df <- filter(df, exp %in% !!exp)
    p <- ggplot(df, aes(pos, counts, fill=!!colaes, color=!!colaes, group=interaction(rep, !!colaes))) +
        annotate_profile(highlightregion) +
        stat_summary(aes(color=NULL), data=function(x)filter(x, boot), geom='ribbon', fun.ymin=function(x)quantile(x, 0.5 * (1 - conf.level)), fun.ymax=function(x)quantile(x, 1 - 0.5 * (1 - conf.level)), alpha=0.2) +
        geom_line(aes(y=counts), data=function(x)filter(x, !boot)) +
        scale_y_continuous(trans='log2') +
        labs(x=sprintf('distance from %s / codons', align), y='enrichment') +
        scale_x_continuous(expand=expand_scale()) +
        guides(fill=guide_legend(override.aes=list(alpha=1)), color=FALSE)
    p
}
