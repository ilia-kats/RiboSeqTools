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
#' @param exclude Genes to exclude.
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
#'          \eqn{\frac{\sum_{i=1}^n \textrm{sample1}_i}{n}}{sum(sample1)/lengh}.  Note that the actual column
#'          name is the value of \code{sample1}.}
#'      \item{sample2_avg_read_density}{Average read density in sample 2, calculated as
#'          \eqn{\frac{\sum_{i=1}^n \textrm{sample2}_i}{n}}{sum(sample2)/lengh}.  Note that the actual column
#'          name is the value of \code{sample2}.}
#'      \item{rank}{Ranking of the gene within the experiment and replicate. Genes are ranked by \code{lo_CI}
#'          in descending order.}
#'}
#' @seealso \link{defaults}
#' @export
binding_scores <- function(data, sample1, sample2, bin, window_size, skip_5prime=0, skip_3prime=0, exclude=c(), conf.level=0.95, bpparam=BiocParallel::bpparam()) {
    check_serp_class(data)
    stopifnot(!is_normalized(data))

    bin <- get_default_param(data, bin)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)
    window_size <- get_default_param(data, window_size)

    fref <- dplyr::filter(get_reference(data), length > (skip_5prime + window_size + skip_3prime))
    lencol <- 'length'
    if (bin == 'byaa') {
       skip_5prime <- skip_5prime %/% 3
       skip_3prime <- skip_3prime %/% 3
       lencol <- 'cds_length'
    }

    scores <- BiocParallel::bplapply(rlang::set_names(fref$gene),
                                     function(gene, ...)binom_ci_profile(data, gene, ...),
                                     sample1=sample1,
                                     sample2=sample2,
                                     bin=bin,
                                     window_size=window_size,
                                     conf.level=conf.level,
                                     BPPARAM=bpparam) %>%
        dplyr::bind_rows(.id='gene') %>%
        {suppressWarnings(dplyr::inner_join(., get_reference(data), by='gene'))} %>%
        dplyr::group_by(exp, rep, gene, !!rlang::sym(lencol)) %>%
        dplyr::group_map(function(.x, .y) {
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
        dplyr::filter(!(gene %in% exclude)) %>%
        dplyr::ungroup()
    avgscores <- dplyr::group_by(scores, exp, gene) %>%
        dplyr::filter_if(is.numeric, dplyr::all_vars(is.finite(.))) %>%
        dplyr::filter(n() > 1) %>%
        dplyr::group_trim() %>%
        dplyr::summarize_if(is.numeric, mean) %>%
        dplyr::mutate(rep='avg')
    dplyr::bind_rows(scores, avgscores) %>%
        dplyr::group_by(exp, rep) %>%
        dplyr::mutate(rank=dplyr::dense_rank(desc(lo_CI))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(gene=as.factor(gene))
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
#' @seealso \code{\link{get_background_model}}, \code{\link{test_binding}}, \code{\link{get_binding_positions}}
#' @export
fit_background_model <- function(data, sample1, sample2, bin, tunnelcoords=18:90, use_quantile=0.5) {
    check_serp_class(data)
    if(is_normalized(data))
        rlang::abort("Need count data")
    bin <- get_default_param(data, bin)
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)

    coords <- tunnelcoords

    if (bin == 'byaa') {
        coords <- coords %/% 3
    }

    ret <- purrr::map(get_data(data), function(exp) {
        purrr::map(exp, function(rep) {
            s1 <- Matrix::rowSums(rep[[sample1]][[bin]][,coords], na.rm=TRUE)
            s2 <- Matrix::rowSums(rep[[sample2]][[bin]][,coords], na.rm=TRUE)
            genes <- intersect(names(s1), names(s2)[s2 >= quantile(s2, use_quantile)])
            opt <- optim(c('s'=2,'m'=0.5), betabinom_ll, gr=betabinom_gradient, x=s1[genes], n=s1[genes] + s2[genes], method='L-BFGS-B', control=list(fnscale=-1), lower=rep(.Machine$double.eps, 2), upper=c(Inf, 1 - .Machine$double.eps))
            ret <- list(samples=list())
            ret$samples[[sample1]] <- s1
            ret$samples[[sample2]] <- s2
            ret$quantile <- use_quantile
            ret$used_gens <- genes
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
            ggplot2::geom_histogram(aes(!!rlang::sym(s1) / (!!rlang::sym(s1) + !!rlang::sym(s2)), stat(density)), binwidth=0.01) +
            geom_line(aes(x, density), data=densities, color='red') +
            facet_grid(rep~exp)
    } else if (type == 'ratio') {
        counts <- mutate(counts, denom=!!rlang::sym(s1)+!!rlang::sym(s2))
        rcounts <- purrr::map_dfr(bgmodel$model, function(exp) {
            purrr::map_dfr(exp, function(rep) {
                rb <- rbetabinom(length(rep$samples[[s1]]),(rep$samples[[s1]] + rep$samples[[s2]]), m=rep$fit$par['m'], s=rep$fit$par['s'])
                tibble::enframe(rb, name='gene', value=s1) %>%
                    mutate(gene=names(rep$samples[[s1]]), !!s2 := rep$samples[[s2]], denom=rep$samples[[s1]]+rep$samples[[s2]], type='random')
            }, .id='rep')
        }, .id='exp')
        bind_rows(counts, rcounts) %>%
            ggplot(aes(!!rlang::sym(s2), !!rlang::sym(s1)/denom, shape=type, color=type)) +
            geom_point(alpha=0.5) +
            scale_x_log10() +
            scale_shape_manual(values=c(observed=1, random=20)) +
            scale_color_manual(values=c(observed='black', random='red')) +
            facet_grid(rep~exp)
    }
}

windowed_readcounts <- function(mat, windowlen) {
    t(apply(mat, 1, function(x) {
        len <- sum(!is.na(x))
        total <- rep(NA_integer_, length(x))
        selectstart <- floor(0.5 * windowlen + 1)
        selectstop <- windowlen - selectstart + 1
        total[1:len] <- as.integer(round(convolve(x[1:len], rep(1, windowlen), type='open')[selectstart:(len+selectstop)]))
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
#' @seealso \code{\link{get_binding_pvalues}}, \code{\link{fit_background_model}}, \code{\link{get_binding_positions}}
#' @export
#' @importFrom rmutil pbetabinom dbetabinom
test_binding <- function(data, window_size, bpparam=BiocParallel::bpparam()) {
    check_serp_class(data)
    window_size <- get_default_param(data, window_size)

    bgmodel <- get_background_model(data)
    if (is.null(bgmodel))
        rlang::abort("No background model present. Run fit_background_model first.")

    ref <- get_reference(data)

    rawdata <- get_data(data)

    winsize <- window_size
    maxtun <- max(bgmodel$tunnelcoords)
    if (bgmodel$bin == 'byaa') {
        winsize <- winsize %/% 3
        maxtun <- maxtun %/% 3
        ref$length <- ref$length %/% 3
    }

    bgdata <- bgmodel$model
    if (!all(names(rawdata) == names(bgmodel$data))) {
        rlang::warn("Sample mismatch between data and model. Restricting p-value calculation to intersection of samples.")
        smpls <- intersect(names(data), names(bgmodel$data))
        rawdata <- rawdata[smpls]
        bgdata <- bgdata[smpls]
    }

    pvals <- BiocParallel::bpmapply(function(dexp, mexp) {
        purrr::map2_dfr(dexp, mexp, function(drep, mrep) {
            s1 <- windowed_readcounts(drep[[bgmodel$sample1]][[bgmodel$bin]], winsize)[, -(1:maxtun)]
            s2 <- windowed_readcounts(drep[[bgmodel$sample2]][[bgmodel$bin]], winsize)[, -(1:maxtun)]

            genes <- intersect(rownames(s1), rownames(s2))
            fref <- dplyr::filter(ref, gene %in% genes)

            purrr::map2_dfr(rlang::set_names(as.character(fref$gene)), fref$length, function(g, l) {
                l <- l - maxtun
                if (l > 0) {
                    pvals <- 1 - pbetabinom(s1[g,1:l], s1[g,1:l] + s2[g,1:l], m=mrep$fit$par['m'], s=mrep$fit$par['s']) + dbetabinom(s1[g,1:l], s1[g,1:l] + s2[g,1:l], m=mrep$fit$par['m'], s=mrep$fit$par['s'])
                    tibble(pos=1:l + maxtun, pval=pvals)
                } else {
                    tibble()
                }
            }, .id='gene')
        }, .id='rep')
    }, rawdata, bgdata, SIMPLIFY=FALSE, BPPARAM=bpparam) %>%
        dplyr::bind_rows(.id='exp') %>%
        group_by(exp, rep) %>%
        mutate(p.adj=p.adjust(pval, method='BY')) %>%
        ungroup()
    set_binding_pvalues(data, pvals)
}

#' Extract statistically significant binding regions
#'
#' Calculates continuous binding regions at given FDR.
#'
#' @param data A \code{serp_data} object. \code{\link{test_binding}} must have been run on the data.
#' @param fdr False discovery rate.
#' @return  A \link[tibble]{tibble} with the following columns: \describe{
#'      \item{exp}{Experiment.}
#'      \item{rep}{Replicate.}
#'      \item{gene}{Gene.}
#'      \item{start}{Start of continuous binding region. If the binning mode for \code{\link{fit_background_model}}
#'           was \code{byaa}, this will be in codons, otherwise in nucleotides.}
#'      \item{end}{End of continuous binding region.}
#'      \item{width}{Width of continuous binding region.}
#'      \item{sample1}{Sum of counts of sample1 within the binding region. Note that the actual column name
#'              is the value of \code{sample1} given to \code{\link{fit_background_model}}.}
#'      \item{sample2}{Sum of counts of sample2 within the binding region. Note that the actual column name
#'              is the value of \code{sample1} given to \code{\link{fit_background_model}}.}
#'}
#' @seealso \code{\link{fit_background_model}}, \code{\link{test_binding}}
#' @export
get_binding_positions <- function(data, fdr=0.01) {
    check_serp_class(data)
    bgmodel <- get_background_model(data)
    pvals <- get_binding_pvalues(data)
    if (is.null(pvals))
        rlang::abort("No p-values present. Run test_binding first.")
    dplyr::filter(pvals, p.adj < fdr) %>%
        dplyr::group_by(exp, rep, gene) %>%
        dplyr::group_map(function(.x, .y) {
            if (nrow(.x) > 1) {
                diffs <- rle(diff(.x$pos))
                starts <- cumsum(c(1,diffs$lengths[1:(length(diffs$lengths) - 1)]))
                ends <- cumsum(diffs$lengths)
                pidx <- which(diffs$values == 1)

                start <- .x$pos[starts[pidx]]
                end <- .x$pos[ends[pidx]]
            } else {
                start <- end <- .x$pos
            }

            s1 <- as.integer(get_data(data)[[.y$exp]][[.y$rep]][[bgmodel$sample1]][[bgmodel$bin]][.y$gene,])
            s2 <- as.integer(get_data(data)[[.y$exp]][[.y$rep]][[bgmodel$sample2]][[bgmodel$bin]][.y$gene,])

            ret <- tibble::tibble(start=start, end=end, width=end - start + 1) %>%
                mutate(!!bgmodel$sample1 := purrr::map2_int(start, end, function(s,e)sum(s1[s:e])), !!bgmodel$sample2 := purrr::map2_int(start, end, function(s,e)sum(s2[s:e])))
        }) %>%
        dplyr::ungroup()
}
