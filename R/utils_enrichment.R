prob2odds <- function(p) p / ( 1 - p )

#' @export
binom_ci <- function(sample1, sample2, sample1_total, sample2_total, conf.level=0.95) {
    probnorm <- sample1_total / sample2_total
    cidf <- binom::binom.agresti.coull(sample1, sample1 + sample2, conf.level = 0.95)
    cidf$lower <- (prob2odds(pmax(0, cidf$lower)) %:% 0) / probnorm
    cidf$upper <- (prob2odds(pmin(1, cidf$upper)) %:% Inf) / probnorm
    cidf$mean <- (prob2odds(cidf$mean) %:% 0) / probnorm
    tibble::tibble(cidf)
}

#' @export
binom_ci_profile <- function(data, gene, sample1, sample2, exp, rep, bin, window_size, conf.level=0.95) {
    check_serp_class(data)
    idx <- which(get_reference(data)$gene == gene)
    genelen <- get_reference(data)$length[idx]
    cdslen <- get_reference(data)$cds_length[idx]

    df <- mapply(function(exp, nexp) {
        df <- mapply(function(rep, nrep) {
            if (!missing(bin) && !(is.null(rep[[sample1]][[bin]]) || !is.null(rep[[sample2]][[bin]]))) {
                stop(sprintf("requested binning level not found for experiment %s replicate %s samples %s and %s", nexp, nrep, sample1, sample2))
            } else if (!is.null(rep[[sample1]]$byaa) && !is.null(rep[[sample2]]$byaa)) {
                winsize <- window_size %/% 3
                bin <- 'byaa'
                len <- genelen
            } else if (!is.null(rep[[sample1]]$bynuc) && !is.null(rep[[sample2]]$bynuc)) {
                winsize <- window_size
                bin <- 'bynuc'
                len <- cdslen
            } else {
                stop(sprintf("binning does not match for experiment %s replicate %s samples %s and %s", nexp, nrep, sample1, sample2))
            }

            s1_data <- rep[[sample1]][[bin]][gene, 1:len]
            s2_data <- rep[[sample2]][[bin]][gene, 1:len]
            win_s1 <- as.integer(round(convolve(s1_data, rep(1, winsize), type='open')[ceiling(winsize/2):(len+floor(winsize/2))]))
            win_s2 <- s.integer(round(convolve(s2_data, rep(1, winsize), type='open')[ceiling(winsize/2):(len+floor(winsize/2))]))

            cidf <- binom_ci(win_s1, win_s2, get_total_counts(data)[[nexp]][[nrep]][[sample1]], get_total_counts(data)[[nexp]][[nrep]][[sample2]], conf.level=conf.level, bin=bin)
            tibble::tibble(winmid=1:len, !!sample1 := s1_data, !!sample2 := s2_data, !!paste0('win_', sample1) := win_s1, !!paste0('win_', sample2) := win_s2, ratio_mean = cidf$mean, lo_CI = cidf$lower, hi_CI = cidf$upper)
        }, exp[rep], names(exp[rep]), SIMPLIFY=FALSE)
        dplyr::bind_rows(df, .id='rep')
    }, get_data(data)[exp], names(get_data(data)[exp]), SIMPLIFY=FALSE)
    dplyr::bind_rows(df, .id='exp')
}

#' @export
make_label_fun <- function(...) {
    dots <- rlang::list2(...)
    refids <- tolower(names(dots))
    names <- unlist(dots)

    function(ids) {
        unname(sapply(tolower(ids), function(i) {
            n <- names[refids == i]
            if (length(n) == 0)
                i
            else
                n
        }))
    }
}

annotate_profile <- function(highlightregion, ...) {
    l <- rlang::list2(...)
    l$geom <- 'rect'
    l$ymin <- 0
    l$ymax <- Inf
    if (is.null(l$fill))
        l$fill <- 'black'
    if (is.null(l$alpha))
        l$alpha <- 0.1
    ann <- function(x)rlang::exec(annotate, xmin=x[1], xmax=x[2], !!!l)
    if (is.list(highlightregion)) {
        lapply(highlightregion, ann)
    } else if (is.matrix(highlightregion) && ncol(highlightregion) == 2) {
        apply(highlightregion, 1, ann)
    }
}

#' @export
plot.serp_data <- function(data, gene, sample1, sample2, exp, rep, bin, window_size, conf.level=0.95, colaes=exp, highlightregion=list(), highlightargs=list()) {
    bin <- get_default_param(bin)
    sample1 <- get_default_param(sample1)
    sample2 <- get_default_param(sample2)
    window_size <- get_default_param(window_size)

    colaes <- rlang::enexpr(colaes)

    df <- binom_ci_profile(data, gene, sample1, sample2, exp, rep, bin, window_size=window_size, conf.level=conf.level)

    xunit <- ifelse(bin == 'byaa', 'codons', 'nucleotides')

    dplyr::group_by(df, rep, exp) %>%
        dplyr::mutate(xmin= (winmid - 0.5) ,
               xmax=(winmid + 0.5),
               alpha=1/((1/!!rlang::sym(paste0('win_', sample1)) + 1/!!rlang::sym(paste0('win_', sample2))) * window_size),
               next.lo=lead(lo_CI, default=0),
               next.hi=lead(hi_CI, default=0),
               overlap=if_else(lo_CI - next.hi > 0, 1L, if_else(next.lo - hi_CI > 0, 2L, 0L)),
               overlap.xmin=xmax - (xmax - xmin) * 0.5,
               overlap.xmax=lead(xmin) + (lead(xmax) - lead(xmin)) * 0.5,
               overlap.ymin=recode(overlap, `0`=NA_real_, `1`=next.hi, `2`=hi),
               overlap.ymax=recode(overlap, `0`=NA_real_, `1`=lo, `2`=next.lo),
               mean_alpha=mean(c(alpha, lead(alpha)), na.rm=TRUE)) %>%
        dplyr::ungroup() %>%
        ggplot2::ggplot(aes(fill=!!colaes)) +
            ggplot2::scale_y_continuous( trans="log2", limits = 2^c(-2.8,6.8), oob=scales::squish, expand=expand_scale(), breaks=2^(-2:6)) +
            ggplot2::scale_x_continuous(expand=expand_scale()) +
            ggplot2::scale_alpha_continuous(trans="sqrt", breaks=c(1,3,10,30), limits=c(0,30), name="prec", range=c(0.02,1)) +
            #expand_limits(x=0) +
            rlang::exec(annotate_profile, highlightregion=highlightregion, !!!highlightargs) +
            ggplot2::geom_rect(aes(xmin= xmin, xmax=xmax, ymin=lo, ymax=hi, alpha=alpha)) +
            ggplot2::geom_rect(aes(xmin=overlap.xmin, xmax=overlap.xmax, ymin=overlap.ymin, ymax=overlap.ymax, alpha=mean_alpha), data=function(y)filter(y, overlap > 0, !is.na(overlap.xmax))) +
            ggplot2::labs(title=gene, x=sprintf("position / %s", xunit), y="enrichment")
}
