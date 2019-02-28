prob2odds <- function(p) p / (1 - p)

#' Calculate enrichment confidence interval.
#'
#' Given total read counts in two samples, calculates library size-normalized confidence intervals
#' for the ratio \eqn{\frac{\textrm{sample1}}{\textrm{sample2}}}{sample1/sample2}.
#'
#' Uses the Agresti-Coull approximation to calculate binomial confidence intervals of
#' \eqn{\frac{\textrm{sample1}}{\textrm{sample1} + \textrm{sample2}}}{sample1/(sample1 + sample2)}.
#' These values are then converted to CIs of \eqn{\frac{\textrm{sample1}}{\textrm{sample2}}}{sample1/sample2}
#' and normalized for library size.
#'
#' @param sample1 Read counts in sample1.
#' @param sample2 Read counts in sample2.
#' @param sample1_total Total number of reads in sample1 (for library size normalization).
#' @param sample2_total Total number of reads in sample2 (for library size normalization).
#' @param conf.level Confidence level.
#' @return A \link[tibble]{tibble} containing the observed proportions and the lower
#' and upper bounds of the confidence interval
#' @export
binom_ci <- function(sample1, sample2, sample1_total, sample2_total, conf.level=0.95) {
    probnorm <- sample1_total / sample2_total
    cidf <- binom::binom.agresti.coull(sample1, sample1 + sample2, conf.level = 0.95)
    cidf$lower <- (prob2odds(pmax(0, cidf$lower)) %:% 0) / probnorm
    cidf$upper <- (prob2odds(pmin(1, cidf$upper)) %:% Inf) / probnorm
    cidf$mean <- (prob2odds(cidf$mean) %:% 0) / probnorm
    tibble::as_tibble(cidf)
}

#' Calculate position-wise enrichment confidence intervals
#'
#' For a given gene, confidence intervals of enrichment (read count ratio) are calculated for each position.
#'
#' At each position within the gene, read counts within a \code{window_size}-wide neighborhood are summed up
#' and used for CI calculation. A confidence interval for the ratio \eqn{\frac{\textrm{sample1}}{\textrm{sample2}}}{sample1/sample2}
#' is calculated using \link{binom_ci}.
#'
#' @param data A \code{serp_data} object.
#' @param gene Name of the gene/ORF.
#' @param sample1 Sample name to use in the numerator.
#' @param sample2 Sample name to use in the denominator.
#' @param exp Experiment name(s). If missing, all experiments will be used.
#' @param rep Replicate name(s). If missing, all replicates will be used.
#' @param bin Binning to use. If missing, will be automatically determined.
#' @param window_size Neighborhood size in nucleotides. Will be automatically converted to codons if binning
#'      by codoons
#' @param conf.level Confidence level.
#' @return A \link[tibble]{tibble} with the following columns: \describe{
#'      \item{winmid}{Center of the neighborhood for which the CI was calculated.}
#'      \item{sample1}{Raw counts of sample1 at position \code{winmid}. Note that the actual column name
#'              is the value of \code{sample1}.}
#'      \item{sample2}{Raw counts of sample2 at position \code{winmid}. Note that the actual column name
#'              is the value of \code{sample2}.}
#'      \item{win_sample1}{Total read counts of sample1 within the neighborhood centered at \code{winmid}.
#'              Note that the actual column name is the value of \code{sample1}.}
#'      \item{win_sample2}{Total read counts of sample2 within the neighborhood centered at \code{winmid}.
#'              Note that the actual column name is the value of \code{sample2}.}
#'      \item{ratio_mean}{Estimated ratio \eqn{\frac{sample1}{sample2}}.}
#'      \item{lo_CI}{Lower confidence bound.}
#'      \item{hi_CI}{Upper confidence bound.}
#'}
#' @export
binom_ci_profile <- function(data, gene, sample1, sample2, exp, rep, bin, window_size, conf.level=0.95) {
    check_serp_class(data)
    stopifnot(!is_normalized(data))

    if (missing(exp))
        exp <- TRUE
    if (missing(rep))
        rep <- TRUE
    binmissing <- missing(bin)

    idx <- which(get_reference(data)$gene == gene)
    genelen <- get_reference(data)$length[idx]
    cdslen <- get_reference(data)$cds_length[idx]

    df <- mapply(function(exp, nexp) {
        df <- mapply(function(rep, nrep) {
            if (!binmissing && !(is.null(rep[[sample1]][[bin]]) || !is.null(rep[[sample2]][[bin]]))) {
                stop(sprintf("requested binning level not found for experiment %s replicate %s samples %s and %s", nexp, nrep, sample1, sample2))
            } else if (!is.null(rep[[sample1]]$byaa) && !is.null(rep[[sample2]]$byaa)) {
                bin <- 'byaa'
            } else if (!is.null(rep[[sample1]]$bynuc) && !is.null(rep[[sample2]]$bynuc)) {
                bin <- 'bynuc'
            } else {
                stop(sprintf("binning does not match for experiment %s replicate %s samples %s and %s", nexp, nrep, sample1, sample2))
            }
            len <- genelen
            winsize <- window_size
            if (bin == 'byaa') {
                len <- cdslen
                winsize <- window_size %/% 3
            }

            s1_data <- rep[[sample1]][[bin]][gene, 1:len]
            s2_data <- rep[[sample2]][[bin]][gene, 1:len]
            win_s1 <- as.integer(round(convolve(s1_data, rep(1, winsize), type='open')[ceiling(winsize/2):(len+floor(winsize/2))]))
            win_s2 <- as.integer(round(convolve(s2_data, rep(1, winsize), type='open')[ceiling(winsize/2):(len+floor(winsize/2))]))

            cidf <- binom_ci(win_s1, win_s2, get_total_counts(data)[[nexp]][[nrep]][[sample1]], get_total_counts(data)[[nexp]][[nrep]][[sample2]], conf.level=conf.level)
            tibble::tibble(winmid=1:len, !!sample1 := s1_data, !!sample2 := s2_data, !!paste0('win_', sample1) := win_s1, !!paste0('win_', sample2) := win_s2, ratio_mean = cidf$mean, lo_CI = cidf$lower, hi_CI = cidf$upper)
        }, exp[rep], names(exp[rep]), SIMPLIFY=FALSE)
        dplyr::bind_rows(df, .id='rep')
    }, get_data(data)[exp], names(get_data(data)[exp]), SIMPLIFY=FALSE)
    dplyr::bind_rows(df, .id='exp')
}

#' Create a mapping between experiment identifiers and names
#'
#' @param ... Name-value pairs of strings. Argument names correspond to experiment IDs and values to
#'      experiment names to be used e.g. in plot legends.
#' @return A function that takes a character vector of experiment IDs and returns a character vector
#'      of experiment names. This function can be used e.g. as labeller in \code{\link[ggplot2]{facet_wrap}}
#'      and\code{\link[ggplot2]{facet_grid}}
#' @examples
#'      labelfun <- make_label_fun(id1='DnaK', id2='DnaK -tig')
#'      labelfun('id1')
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
    ann <- function(x)rlang::exec(ggplot2::annotate, xmin=x[1], xmax=x[2], !!!l)
    if (is.list(highlightregion)) {
        lapply(highlightregion, ann)
    } else if (is.matrix(highlightregion) && ncol(highlightregion) == 2) {
        apply(highlightregion, 1, ann)
    }
}

#' Plot an enrichment profile
#'
#' For each position, a binomial confidence interval calculated using \code{\link{binom_ci_profile}} is
#' plotted. Transparency reflects the total number of reads contributing to the confidence interval.
#'
#' @param data A \code{serp_data} object. Must contain raw (unnormalized) read counts.
#' @param gene Name of the gene/ORF to plot.
#' @param sample1 Name of the first sample (the numerator). If missing, the default sample1 of the data set
#'      will be used
#' @param sample2 Name of the second sample (the denominator). If missing, the default sample2 of the data set
#'      will be used
#' @param exp Character vector of experiments to plot. If missing, all experiments are plotted.
#' @param rep Character vector of replicates to plot. If missing, all replicates will be plotted.
#' @param bin Bin level (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
#'      will be used
#' @param window_size Window size for the confidence interval calculation. If missing, the default window size
#'      of the data set will be used
#' @param conf.level Confidence level.
#' @param colaes Variable to use for the color scale.
#' @template plot_annotations
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \link{defaults}
#' @export
plot.serp_data <- function(data, gene, sample1, sample2, exp, rep, bin, window_size, conf.level=0.95, colaes=exp, highlightregion=list(), highlightargs=list()) {
    stopifnot(!is_normalized(data))

    bin <- get_default_param(data, bin)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)
    window_size <- get_default_param(data, window_size)

    ylim <- get_default_param(data, plot_ylim)
    ybreaks <- get_default_param(data, plot_ybreaks)

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
               overlap.ymin=recode(overlap, `0`=NA_real_, `1`=next.hi, `2`=hi_CI),
               overlap.ymax=recode(overlap, `0`=NA_real_, `1`=lo_CI, `2`=next.lo),
               mean_alpha=mean(c(alpha, lead(alpha)), na.rm=TRUE)) %>%
        dplyr::ungroup() %>%
        ggplot2::ggplot(ggplot2::aes(fill=!!colaes)) +
            ggplot2::scale_y_continuous(trans="log2", limits=ylim, oob=scales::squish, expand=ggplot2::expand_scale(), breaks=ybreaks) +
            ggplot2::scale_x_continuous(expand=ggplot2::expand_scale()) +
            ggplot2::scale_alpha_continuous(trans="sqrt", breaks=c(1,3,10,30), limits=c(0,30), name="prec", range=c(0.02,1)) +
            rlang::exec(annotate_profile, highlightregion=highlightregion, !!!highlightargs) +
            ggplot2::geom_rect(ggplot2::aes(xmin=xmin, xmax=xmax, ymin=lo_CI, ymax=hi_CI, alpha=alpha)) +
            ggplot2::geom_rect(ggplot2::aes(xmin=overlap.xmin, xmax=overlap.xmax, ymin=overlap.ymin, ymax=overlap.ymax, alpha=mean_alpha), data=function(y)filter(y, overlap > 0, !is.na(overlap.xmax))) +
            ggplot2::labs(title=gene, x=sprintf("position / %s", xunit), y="enrichment")
}
