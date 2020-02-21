convolve_selection <- function(winsize, len) {
    selectstart <- floor(0.5 * winsize + 1)
    selectstop <- winsize - ceiling(0.5 * winsize)
    c(selectstart, selectstop + len)
}

binding_scores_per_position <- function(data, sample1, sample2, bin, window_size, skip_5prime=0, skip_3prime=0, conf.level=0.95, bpparam=BiocParallel::bpparam()) {
    check_serp_class(data)
    stopifnot(!is_normalized(data))

    bin <- get_default_param(data, bin)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)
    window_size <- get_default_param(data, window_size)

    exclude <- excluded(data)

    fref <- dplyr::filter(get_reference(data), length > (skip_5prime + window_size + skip_3prime))
    lencol <- 'length'
    if (bin == 'byaa') {
       skip_5prime <- skip_5prime %/% 3
       skip_3prime <- skip_3prime %/% 3
       lencol <- 'cds_length'
    }

    BiocParallel::bplapply(rlang::set_names(fref$gene),
                                    function(gene, ...)binom_ci_profile(data, gene, ...),
                                    sample1=sample1,
                                    sample2=sample2,
                                    bin=bin,
                                    window_size=window_size,
                                    conf.level=conf.level,
                                    ignore_genecol=TRUE,
                                    BPPARAM=bpparam) %>%
        dplyr::bind_rows(.id='gene') %>%
        {suppressWarnings(dplyr::inner_join(., get_reference(data), by='gene'))} %>%
        dplyr::group_by(exp) %>%
        dplyr::filter(!(gene %in% exclude[[exp[1]]]))
}

#' Calculate binding scores for each gene in a SeRP experiment
#'
#' A binding score for a gene is defined as the highest value of the position-wise confidence interval for the
#' ratio \eqn{\frac{\textrm{sample1}}{\textrm{sample2}}}{sample1/sample2}. Enrichment confidence intervals are
#' calculated with \link{binom_ci_profile} separately for each experiment and replicate. Calculated values for
#' all replicates are averaged to create a new replicate \code{avg} which is returned along with values for the
#' individual replicates.
#'
#' @param data A \code{serp_data} object.
#' @param sample1 Name of the first sample (the numerator). If missing, the default sample1 of the data set
#'      will be used.
#' @param sample2 Name of the second sample (the denominator). If missing, the default sample2 of the data set
#'      will be used.
#' @param bin Binning mode (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
#'      will be used.
#' @param window_size Neighborhood size for the confidence interval calculation in nucleotides. If missing, the default
#'      window size of the data set will be used.
#' @param skip_5prime How many nucleotides to skip at the 5' end of the ORF. Useful if you know that the 5' end
#'      contains artifacts.
#' @param skip_3prime How many nucleotides to skip at the 3' end of the ORF. useful if you know that the 3' end
#'      contains artifacts.
#' @param conf.level Confidence level.
#' @param bpparam A \code{\link[BiocParallel]{BiocParallelParam-class}} object.
#' @return A \link[tibble]{tibble} with the following columns: \describe{
#'      \item{gene}{The gene/ORF name.}
#'      \item{exp}{The experiment name.}
#'      \item{rep}{The replicate name.}
#'      \item{score_position}{Position within the gene for which the confidence interval is returned. Corresponds
#'          to the position where the highest value of the lower CI bound was observed. If \code{bin == 'byaa'}
#'          this is measured in codons, otherwise in nucleotides.}
#'      \item{lo_CI}{Lower confidence bound at the position \code{score_position}.}
#'      \item{hi_CI}{Upper confidence bound at the position \code{score_position}.}
#'      \item{mean}{Estimated enrichment at the position \code{score_position}.}
#'      \item{total_enrichment}{Total enrichment of the gene, calculated using all reads mapped anywhere within
#'          the ORF}
#'      \item{sample1_total_counts}{Total number of reads mapped to this ORF. Note that the actual column name
#'          is the value of \code{sample1}.}
#'      \item{sample2_total_counts}{Total number of reads mapped to this ORF. Note that the actual column name
#'          is the value of \code{sample2}.}
#'      \item{sample1_total_RPM}{\code{sample1_total_counts} normalized to the number of mapped reads. Note that
#'          the actual column name is the value of \code{sample1}.}
#'      \item{sample2_total_RPM}{\code{sample2_total_counts} normalized to the number of mapped reads. Note that
#'          the actual column name is the value of \code{sample2}.}
#'      \item{sample1_avg_read_density}{Average read density in sample 1, calculated as
#'          \eqn{\frac{\sum_{i=1}^n \textrm{sample1}_i}{n}}{sum(sample1)/lengh}. If \code{bin} is \code{byaa},
#'          density is in reads/nucleotide, otherwise in reads/codon. Note that the actual column
#'          name is the value of \code{sample1}.}
#'      \item{sample2_avg_read_density}{Average read density in sample 2, calculated as
#'          \eqn{\frac{\sum_{i=1}^n \textrm{sample2}_i}{n}}{sum(sample2)/lengh}. If \code{bin} is \code{byaa},
#'          density is in reads/nucleotide, otherwise in reads/codon. Note that the actual column
#'          name is the value of \code{sample2}.}
#'      \item{rank}{Ranking of the gene within the experiment and replicate. Genes are ranked by \code{lo_CI}
#'          in descending order.}
#'}
#' @seealso \code{\link{defaults}}
#' @export
binding_scores <- function(data, sample1, sample2, bin, window_size, skip_5prime=0, skip_3prime=0, conf.level=0.95, bpparam=BiocParallel::bpparam()) {
    check_serp_class(data)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)

    binding_scores_per_position(data, sample1, sample2, bin, window_size, skip_5prime, skip_3prime, conf.level, bpparam) %>%
        dplyr::group_by(exp, rep, gene, !!rlang::sym(lencol)) %>%
        dplyr::group_modify(function(.x, .y) {
            pos <- which.max(.x$lo_CI[skip_5prime:(nrow(.x) - skip_3prime)]) + skip_5prime - 1
            gene <- as.character(.y$gene)
            exp <- .y$exp
            rep <- .y$rep
            len <- .y[[lencol]]

            s1 <- get_data(data)[[exp]][[rep]][[sample1]][[bin]][gene, 1:len]
            s2 <- get_data(data)[[exp]][[rep]][[sample2]][[bin]][gene, 1:len]

            total <- get_total_counts(data)[[exp]][[rep]]

            tibble::tibble(lo_CI=.x$lo_CI[pos],
                           hi_CI=.x$hi_CI[pos],
                           mean=.x$ratio_mean[pos],
                           score_position=pos,
                           total_enrichment=sum(s1) / sum(s2) * total[[sample2]] / total[[sample1]],
                           !!paste0(sample1, '_total_counts') := sum(s1),
                           !!paste0(sample2, '_total_counts') := sum(s2),
                           !!paste0(sample1, '_total_RPM') := sum(s1) / total[[sample1]] * 1e6,
                           !!paste0(sample2, '_total_RPM') := sum(s2) / total[[sample2]] * 1e6,
                           !!paste0(sample1, '_avg_read_density') := sum(s1) / len,
                           !!paste0(sample2, '_avg_read_density') := sum(s2) / len)
        }) %>%
        dplyr::ungroup()
    avgscores <- dplyr::group_by(scores, exp, gene) %>%
        dplyr::filter_if(is.numeric, dplyr::all_vars(is.finite(.))) %>%
        dplyr::filter(n() > 1) %>%
        dplyr::group_trim() %>%
        dplyr::summarize_if(is.numeric, mean) %>%
        dplyr::mutate(rep='avg')
    dplyr::bind_rows(scores, avgscores) %>%
        dplyr::group_by(exp, rep) %>%
        dplyr::mutate(rank=dplyr::dense_rank(dplyr::desc(lo_CI))) %>%
        dplyr::ungroup() %>%
        map_df_genenames(data) %>%
        dplyr::mutate(exp=as.factor(exp), rep=as.factor(rep))
}

#' @importFrom rmutil dbetabinom
betabinom_ll <- function(pars, x, n) {
    sum(dbetabinom(x, n, m=pars['m'], s=pars['s'], log=TRUE))
}

#' @importFrom rmutil dbetabinom
# with workaround to handle choose() returning Inf for large values
betabinom_gradient <- function(pars, x, n) {
    f <- x + pars['s'] * pars['m']
    g <- n - x + pars['s']*(1 - pars['m'])

    i <- rep(pars['s'] * pars['m'], length(x))
    j <- rep(pars['s'] * (1 - pars['m']), length(x))

    J_h <- matrix(c(pars['m'], pars['s'], 1 - pars['m'], -pars['s']), nrow=2, byrow=TRUE)
    J_k <- matrix(c(pars['m'], pars['s'], 1 - pars['m'], -pars['s']), nrow=2, byrow=TRUE)

    J_b <- function(a,b)matrix(c(beta(a,b) * (digamma(a) - digamma(a + b)), beta(a,b) * (digamma(b) - digamma(a+b))), ncol=2, byrow=FALSE)

    J_B1 <- J_b(f,g) %*% J_h
    J_B2 <- J_b(i,j) %*% J_k

    J_p <- (beta(i,j) * J_B1 - beta(f,g) * J_B2) / beta(i,j)^2
    sgn <- sign(J_p)
    J_p <- exp(log(abs(J_p)) + lchoose(n,x)) * sgn

    grad <- as.vector(apply(J_p, 2, function(y)sum(y / dbetabinom(x, n, m=pars['m'], s=pars['s']))))
    grad
}

#' Fit a background model to a selective ribosome profiling data set
#'
#' Reads in \code{sample1} at ribosome positions for which the nascent chain is completely obscured
#' by the ribosome exit tunnel (typically the first 30 codons of a CDS) are assumed to represent background
#' binding. A background distribution is fitted to read counts within the tunnel.
#'
#' For each gene, read counts within \code{tunnelcoords} are summed up (separately for \code{sample1} and
#' \code{sample2}). Counts in \code{sample1} conditioned on \code{sample1} + \code{sample2} are assumed
#' to follow a beta-binomial distribution. Parameter estimation is performed by maximum likelihood separately
#' for each sample and replicate.
#'
#' @param data A \code{serp_data} object.
#' @param sample1 Name of the first sample (the numerator). If missing, the default sample1 of the data set
#'      will be used.
#' @param sample2 Name of the second sample (the denominator). If missing, the default sample2 of the data set
#'      will be used.
#' @param bin Binning mode (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
#'      will be used.
#' @param tunnelcoords Ribosome positions in nucleotites for which it is assumed that the nascent chain is
#'      fully obscured by the exit tunnel.
#' @param quantile Only genes whith \code{sample2 >= quantile(sample2, use_quantile)} will be used for
#'      parameter estimation. This reduces the effect of random drop-outs on the final model.
#' @return A \code{serp_data} object.
#' @seealso \code{\link{get_background_model}}, \code{\link{test_binding}}, \code{\link{get_binding_positions}}, \code{\link{plot_binding_positions}}
#' @export
fit_background_model <- function(data, sample1, sample2, bin, tunnelcoords=18:90, use_quantile=0.5) {
    check_serp_class(data)
    if(is_normalized(data))
        rlang::abort("Need count data")
    bin <- get_default_param(data, bin)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)

    coords <- tunnelcoords
    exclude <- excluded(data)

    if (bin == 'byaa') {
        coords <- coords %/% 3
    }

    ret <- purrr::imap(get_data(data), function(exp, nexp) {
        purrr::map(exp, function(rep) {
            s1 <- Matrix::rowSums(rep[[sample1]][[bin]][,coords], na.rm=TRUE)
            s2 <- Matrix::rowSums(rep[[sample2]][[bin]][,coords], na.rm=TRUE)
            genes <- intersect(names(s1), names(s2)[s2 >= quantile(s2, use_quantile)])
            genes <- genes[!(genes %in% exclude[[nexp]])]
            opt <- optim(c('s'=2,'m'=0.5), betabinom_ll, gr=betabinom_gradient, x=s1[genes], n=s1[genes] + s2[genes], method='L-BFGS-B', control=list(fnscale=-1), lower=rep(.Machine$double.eps, 2), upper=c(Inf, 1 - .Machine$double.eps))
            ret <- list(samples=list())
            ret$samples[[sample1]] <- s1
            ret$samples[[sample2]] <- s2
            ret$quantile <- use_quantile
            ret$used_genes <- genes
            ret$fit <- opt
            ret
        })
    })
    set_background_model(data, list(model=ret, bin=bin, sample1=sample1, sample2=sample2, tunnelcoords=tunnelcoords)) %>%
        set_binding_pvalues(NULL)
}

#' Quality control plots for background model estimation
#'
#' Creates one plot per sample and replicate in a faceted display.
#'
#' In \code{density} mode, plots a histogram of \eqn{\frac{\textrm{sample1}}{\textrm{sample1} + \textrm{sample2}}}{sample1/(sample1 + sample2)}
#' overlaid with the probability density function of a \link[=dbeta]{Beta distribution} using the
#' parameters estimated by \code{\link{fit_background_model}}.
#'
#' In \code{ratio} mode, plots a scatterplot of \eqn{\frac{\textrm{sample1}}{\textrm{sample1} + \textrm{sample2}}}{sample1/(sample1 + sample2)}
#' over \code{sample2} overlaid with points randomly drawn from a beta-binomial distribution using
#' the parameters estimated by \code{\link{fit_background_model}}.
#'
#' @param data A \code{serp_data} object.
#' @param type Type of plot, one of \code{density} or \code{ratio}.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{fit_background_model}}
#' @importFrom rmutil rbetabinom
#' @export
plot_background_model <- function(data, type=c('density', 'ratio')) {
    check_serp_class(data)
    bgmodel <- get_background_model(data)
    if(is.null(bgmodel))
        rlang::abort("No background model present. Run fit_background_model first.")
    type <- match.arg(type)
    s1 <- bgmodel$sample1
    s2 <- bgmodel$sample2
    counts <- purrr::map_dfr(bgmodel$model, function(exp) {
        purrr::map_dfr(exp, function(rep) {
            dplyr::inner_join(tibble::enframe(rep$samples[[s1]], name='gene', value=s1), tibble::enframe(rep$samples[[s2]], name='gene', value=s2), by='gene')
        }, .id='rep')
    }, .id='exp') %>%
        mutate(type='observed')

    if (type == 'density') {
        densities <- purrr::map_dfr(bgmodel$model, function(exp) {
            purrr::map_dfr(exp, function(rep) {
                x <- seq(from=0, to=1, by=0.01)
                y <- dbeta(x, rep$fit$par['m'] * rep$fit$par['s'], rep$fit$par['s'] * (1 - rep$fit$par['m']))
                tibble::tibble(x=x, density=y)
            }, .id='rep')
        }, .id='exp')
        ggplot2::ggplot(counts) +
            ggplot2::geom_histogram(ggplot2::aes(!!rlang::sym(s1) / (!!rlang::sym(s1) + !!rlang::sym(s2)), stat(density)), binwidth=0.01) +
            ggplot2::geom_line(ggplot2::aes(x, density), data=densities, color='red') +
            ggplot2::facet_grid(rep~exp)
    } else if (type == 'ratio') {
        dplyr::mutate(counts, denom=!!rlang::sym(s1)+!!rlang::sym(s2)) %>%
            dplyr::group_by(exp, rep) %>%
            dplyr::group_modify(function(.x, .y) {
                pars <- bgmodel$model[[.y$exp]][[.y$rep]]$fit$par
                ret <- .x
                ret[[s1]] <- rbetabinom(nrow(ret), ret$denom, m=pars['m'], s=pars['s'])
                ret$type='random'
                dplyr::bind_rows(.x, ret)
            }) %>%
            ggplot2::ggplot(ggplot2::aes(!!rlang::sym(s2), !!rlang::sym(s1)/denom, shape=type, color=type)) +
                ggplot2::geom_point(alpha=0.5) +
                ggplot2::scale_x_log10() +
                ggplot2::scale_shape_manual(values=c(observed=1, random=20)) +
                ggplot2::scale_color_manual(values=c(observed='black', random='red')) +
                ggplot2::facet_grid(rep~exp)
    }
}

windowed_readcounts <- function(mat, windowlen) {
    t(apply(mat, 1, function(x) {
        len <- sum(!is.na(x))
        total <- rep(NA_integer_, length(x))
        cwindow <- convolve_selection(windowlen, len)
        total[1:len] <- as.integer(round(convolve(x[1:len], rep(1, windowlen), type='open')[cwindow[1]:cwindow[2]]))
        total
    }))
}

#' Statistical test for position-wise enrichment in selective ribosome profiling data
#'
#' Using the background model estimated by \code{\link{fit_background_model}}, performs significance testing for enrichment of reads in
#' \code{sample1} compared to \code{sample2} for all genes and positions.
#'
#' Enrichment spikes at individual positions are likely not indicative of real chaperone binding. Nevertheless, they frequently occur,
#' probably due to ligation and/or sequencing bias. To alleviate this problem somewhat, read counts within a \code{window_size}-wide
#' neighborhood are summed up individually for \code{sample1} and \code{sample2} and used for significance testing. This results in
#' strongly correlatd p-values, therefore the Benjamini-Yekutieli-procedure is used for FDR control. \code{sample1} and \code{sample2}
#' from \code{\link{fit_background_model}} are used. Note that the binning mode (bynuc or byaa) previously given to \code{\link{fit_background_model}}
#' will be used.
#'
#' @param data A \code{serp_data} object. \code{\link{fit_background_model}} must have been run on the data.
#' @param window_size Neighborhood size for the confidence interval calculation in nucleotides. If missing, the default
#'      window size of the data set will be used.
#' @param bpparam A \code{\link[BiocParallel]{BiocParallelParam-class}} object.
#' @return A \code{serp_data} object.
#' @seealso \code{\link{get_binding_pvalues}}, \code{\link{fit_background_model}}, \code{\link{get_binding_positions}}, \code{\link{plot_binding_positions}}
#' @export
#' @importFrom rmutil pbetabinom dbetabinom
test_binding <- function(data, window_size, bpparam=BiocParallel::bpparam()) {
    check_serp_class(data)
    window_size <- get_default_param(data, window_size)

    bgmodel <- get_background_model(data)
    if (is.null(bgmodel))
        rlang::abort("No background model present. Run fit_background_model first.")

    ref <- get_reference(data)
    exclude <- excluded(data)
    rawdata <- get_data(data)

    winsize <- window_size
    maxtun <- max(bgmodel$tunnelcoords)
    if (bgmodel$bin == 'byaa') {
        winsize <- winsize %/% 3
        maxtun <- maxtun %/% 3
        ref$length <- ref$length %/% 3
    }

    bgdata <- bgmodel$model
    if (!all(names(rawdata) == names(bgdata))) {
        rlang::warn("Sample mismatch between data and model. Restricting p-value calculation to intersection of samples.")
        smpls <- intersect(names(data), names(bgdata))
        rawdata <- rawdata[smpls]
        bgdata <- bgdata[smpls]
    }

    pvals <- BiocParallel::bpmapply(function(dexp, mexp, nexp) {
        purrr::map2_dfr(dexp, mexp, function(drep, mrep) {
            s1 <- windowed_readcounts(drep[[bgmodel$sample1]][[bgmodel$bin]], winsize)[, -(1:maxtun)]
            s2 <- windowed_readcounts(drep[[bgmodel$sample2]][[bgmodel$bin]], winsize)[, -(1:maxtun)]

            genes <- intersect(rownames(s1), rownames(s2))
            fref <- dplyr::filter(ref, gene %in% genes, !(gene %in% exclude[[nexp]]))

            purrr::map2_dfr(rlang::set_names(as.character(fref$gene)), fref$length, function(g, l) {
                l <- l - maxtun
                if (l > 0) {
                    pvals <- 1 - pbetabinom(s1[g,1:l], s1[g,1:l] + s2[g,1:l], m=mrep$fit$par['m'], s=mrep$fit$par['s']) + dbetabinom(s1[g,1:l], s1[g,1:l] + s2[g,1:l], m=mrep$fit$par['m'], s=mrep$fit$par['s'])
                    tibble::tibble(pos=1:l + maxtun, pval=pvals)
                } else {
                    tibble::tibble()
                }
            }, .id='gene')
        }, .id='rep')
    }, rawdata, bgdata, names(rawdata), SIMPLIFY=FALSE, BPPARAM=bpparam) %>%
        dplyr::bind_rows(.id='exp') %>%
        dplyr::mutate(exp=as.factor(exp), rep=as.factor(rep)) %>%
        dplyr::group_by(exp, rep) %>%
        dplyr::mutate(p.adj=p.adjust(pval, method='BY')) %>%
        dplyr::ungroup()
    set_binding_pvalues(data, pvals)
}


pos_to_peak_id <- function(pos) {
    diffs <- rle(c(1, diff(pos)))
    cumsum(rep(diffs$values > 1, times=diffs$lengths)) + 1
}

pos_to_start_end <- function(pos) {
    pidx <- rle(pos_to_peak_id(pos))
    start <- 1
    if (length(pidx$lengths) > 1)
        start <- cumsum(c(start, pidx$lengths[1:(length(pidx$lengths) - 1)]))
    end <- cumsum(pidx$lengths)
    tibble::tibble(start=pos[start], end=pos[end], width=end - start + 1)
}

#' Extract statistically significant binding regions
#'
#' Calculates continuous binding regions at given FDR.
#'
#' @param data A \code{serp_data} object. \code{\link{test_binding}} must have been run on the data.
#' @param fdr False discovery rate.
#' @templateVar sample_colname_suffix  given to \code{\link{fit_background_model}}
#' @template binding_positions
#' @seealso \code{\link{fit_background_model}}, \code{\link{test_binding}}, \code{\link{plot_binding_positions}}, \code{\link{get_binding_positions_by_threshold}}
#' @importFrom rlang %@%<-
#' @export
get_binding_positions <- function(data, fdr=0.01) {
    check_serp_class(data)
    bgmodel <- get_background_model(data)
    pvals <- get_binding_pvalues(data)
    if (is.null(pvals))
        rlang::abort("No p-values present. Run test_binding first.")
    df <- dplyr::filter(pvals, p.adj < fdr) %>%
        dplyr::group_by(exp, rep, gene) %>%
        dplyr::group_modify(function(.x, .y) {
            s1 <- as.integer(get_data(data)[[.y$exp]][[.y$rep]][[bgmodel$sample1]][[bgmodel$bin]][.y$gene,])
            s2 <- as.integer(get_data(data)[[.y$exp]][[.y$rep]][[bgmodel$sample2]][[bgmodel$bin]][.y$gene,])
            pos_to_start_end(.x$pos) %>%
                dplyr::mutate(!!bgmodel$sample1 := purrr::map2_int(start, end, function(s,e)sum(s1[s:e])), !!bgmodel$sample2 := purrr::map2_int(start, end, function(s,e)sum(s2[s:e])))
        }) %>%
        map_df_genenames(data) %>%
        dplyr::ungroup()
    df %@% ref <- get_reference(data)
    df
}

#' Extract binding regions based on a threshold
#'
#' Calculates binding regions given a threshold. For each gene, enrichment confidence intervals over
#' smoothed read counts are calculated by \code{\link{binom_ci_profile}}. If the lower confidence bound
#' is above the threshold at a position, this position is considered bound with high confidence. If the
#' threshold lies within the confidence interval for a position, this position is considered bound with
#' low confidence. Positions with the upper bound of the confidence interval below the threshold are
#' considered unbound.
#'
#' @template scores
#' @param threshold The enrichment threshold.
#' @templateVar additional_columns \item{binding_class}{Binding class of the binding region. 2 for high-confidence regions, 1 for low-confidence regions.}
#' @template binding_positions
#' @seealso \code{\link{get_binding_positions}}, \code{\link{plot_binding_positions}}
#' @importFrom rlang %@%<-
#' @export
get_binding_positions_by_threshold <- function(data, sample1, sample2, bin, window_size, skip_5prime=0, skip_3prime=0, conf.level=0.95, bpparam=BiocParallel::bpparam(), threshold=1.0) {
    check_serp_class(data)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)

    df <- binding_scores_per_position(data, sample1, sample2, bin, window_size, skip_5prime, skip_3prime, conf.level, bpparam) %>%
        dplyr::rename(pos=winmid) %>%
        dplyr::mutate(binding_class=dplyr::if_else(lo_CI > threshold, 2L, dplyr::if_else(hi_CI < threshold, 0L, 1L))) %>%
        dplyr::select(exp, rep, gene, pos, binding_class, !!sample1, !!sample2) %>%
        dplyr::filter(binding_class > 0) %>%
        dplyr::group_by(binding_class, exp, rep, gene) %>%
        dplyr::mutate(pid=pos_to_peak_id(pos)) %>%
        dplyr::group_by(pid, add=TRUE) %>%
        dplyr::summarize(!!sample1 := sum(!!rlang::sym(sample1)), !!sample2 := sum(!!rlang::sym(sample2)), start=pos[1], end=pos[dplyr::n()]) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(width=end - start + 1, binding_class=as.factor(binding_class)) %>%
        dplyr::select(exp, rep, gene, start, end, width, !!sample1, !!sample2, binding_class)

    df %@% ref <- get_reference(data)
    df
}

#' Plot binding regions
#'
#' Plots chaperone binding regions in a heat map. Color intensity is proportional to the fraction of replicates
#' in which binding at the respective position was observed.
#'
#' @param df A data frame created by \code{\link{get_binding_positions}} or \code{\link{get_binding_positions_by_threshold}}.
#' @param bgcolor Background color to visualize unbound regions.
#' @param ylabels Whether to show Y axis labels (gene names). Defaults to suppressing Y axis labels if more than
#'      10 genes are plotted.
#' @param lowconf If \code{df} was created by \code{\link{get_binding_positions_by_threshold}}, this
#'      determines the color intensity of the low-confidence regions.
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @seealso \code{\link{fit_background_model}}, \code{\link{test_binding}}, \code{\link{get_binding_positions}}, \code{get_binding_positions_by_threshold}
#' @importFrom rlang %@%
#' @export
plot_binding_positions <- function(df, bgcolor="lightgrey", ylabels=NA, lowconf=0.2) {
    deps <- purrr::map_lgl(purrr::set_names(c("IRanges", "S4Vectors")), purrr::partial(requireNamespace, ...=, quietly=TRUE))
    if (!all(deps)) {
        pkgs <- paste(names(deps)[!deps], collapse=", ")
        plr <- ifelse(length(pkgs) > 1, "s", "")
        rlang::abort(sprintf("%s package%s not found, please install the %s package%s to enable this functionality", pkgs, plr, pkgs, plr))
    }

    dfname <- as.character(rlang::enexpr(df))
    ref <- df %@% "ref"
    if (is.null(ref))
        rlang::abort(sprintf("'%s' does not have a 'ref' attribute.", dfname))

    nrep <- dplyr::distinct(df, exp, rep) %>%
        dplyr::count(exp, name="nrep")
    if (!rlang::has_name(df, "binding_class"))
        df$binding_class <- 2L
    df <- dplyr::group_by(df, exp, gene) %>%
        dplyr::summarize(bound=list(IRanges::IRanges(start=start, end=end, binding_class=binding_class))) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(ref, by='gene') %>%
        purrr::pmap_dfr(function(exp, gene, bound, cds_length, ...) {
            bclass <- mcols(bound)$binding_class
            weight <- dplyr::if_else(bclass == 2L, 1, dplyr::if_else(bclass == 1L, lowconf, 0))
            cov <- IRanges::coverage(bound, width=cds_length, weight=weight)
            tibble::tibble(exp=exp, gene=gene, xmid=0.5 * (S4Vectors::start(cov) + S4Vectors::end(cov)), width=S4Vectors::width(cov), bound=S4Vectors::runValue(cov), length=cds_length)
        }) %>%
        dplyr::filter(bound > .Machine$double.eps) %>%
        dplyr::inner_join(nrep) %>%
        dplyr::mutate(bound=bound / nrep, gene=factor(gene, levels=unique(gene[order(length, decreasing=TRUE)]), ordered=TRUE))
    p <- ggplot2::ggplot(df, ggplot2::aes(y=gene, fill=exp)) +
        ggplot2::geom_tile(ggplot2::aes(x=0.5 * length + 0.5, width=length), fill=bgcolor) +
        ggplot2::geom_tile(ggplot2::aes(x=xmid, width=width, height=1, alpha=bound)) +
        ggplot2::scale_x_continuous(expand=ggplot2::expand_scale(), name="position / codons") +
        ggplot2::scale_alpha_identity() +
        ggplot2::guides(fill=FALSE) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank())
    if (isFALSE(ylabels) || is.na(ylabels) && max(dplyr::count(df, exp)$n) > 10)
        p <- p + ggplot2::theme(axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank())
    if (length(unique(df$exp)) > 1)
        p <- p + ggplot2::facet_wrap(~exp, scales='free_y')
    p
}
