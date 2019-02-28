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
#' @param wndow_size Neighborhood size for the confidence interval calculation in nucleotides. If missing, the default
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
            gene <- .y$gene
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
        dplyr::ungroup() %>%
        dplyr::mutate(gene=as.factor(gene))
    avgscores <- dplyr::group_by(scores, exp, gene) %>%
        dplyr::filter_if(is.numeric, dplyr::all_vars(is.finite(.))) %>%
        dplyr::filter(n() > 1) %>%
        dplyr::summarize_if(is.numeric, mean) %>%
        dplyr::mutate(rep='avg')
    dplyr::bind_rows(scores, avgscores) %>%
        dplyr::group_by(exp, rep) %>%
        dplyr::mutate(rank=dplyr::dense_rank(desc(lo_CI))) %>%
        dplyr::ungroup()
}
