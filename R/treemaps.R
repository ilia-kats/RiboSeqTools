#' Treemap plot of read counts per gene with an additional grouping variable
#'
#' This function tries really hard to make sure that the ordering of class groups in the plot
#' corresponds to the ordering of the factor levels. The treemap is plotted using \code{\link[treemap]{treemap}}.
#'
#' @param data A \code{serp_data} object.
#' @param exp Experiment name.
#' @param rep Replicate name. If missing or NULL, read counts will be averaged over all replicates.
#' @param sample Sample type.
#' @param geneclass Data frame with columns \code{gene} and \code{class}. \code{class} must be an
#'      ordered factor.
#' @param title Plot title.
#' @param palette RColorBrewer palette to color-code the class variable
#' @param exclude Character vector of genes to exclude from the plot
#' @return An invisible list from \code{\link[treemap]{treemap}}
#' @seealso \code{\link[treemap]{treemap}}
#' @importFrom magrittr %$%
#' @export
plot_treemap <- function(data, exp, rep, sample, geneclass, title='', palette="Set2", exclude=c()) {
    check_serp_class(data)

    if (ncol(geneclass) == 2 && !('class' %in% colnames(geneclass)) && 'gene' %in% colnames(geneclass))
        geneclass <- dplyr::rename(geneclass, class=-dplyr::matches('gene'))
    classes <- levels(geneclass$class)
    if (!('unknown' %in% classes))
        classes <- c(classes, 'unknown')

    reads_per_gene_sample <- if (missing(rep) || is.null(rep)) {
        if (!is_normalized(data)) {
            rlang::warn("Unnormalized data given, normalizing")
            data <- normalize(data)
        }
        allcounts <- purrr::map(get_data(data)[[exp]], function(rep)Matrix::rowSums(rep[[sample]][[1]], na.rm=TRUE))
        genes <- purrr::reduce(allcounts, function(x, y)intersect(names(x), names(y)))
        purrr::reduce(allcounts, function(x, y)x[genes] + y[genes]) / length(allcounts)
    } else {
        Matrix::rowSums(get_data(data)[[exp]][[rep]][[sample]][[1]], na.rm=TRUE)
    }
    reads_per_gene_sample <- tibble::enframe(reads_per_gene_sample, name='gene', value='read_sum') %>%
        dplyr::left_join(geneclass, by='gene') %>%
        tidyr::replace_na(list(class='unknown')) %>%
        dplyr::mutate(class=factor(class, classes, ordered=TRUE)) %>%
        dplyr::filter(!(gene %in% union(exclude, excluded(data))))

    maxreadsums <- dplyr::group_by(reads_per_gene_sample, class) %>%
        dplyr::summarize(s=sum(1/(read_sum + 1))) %$%
        max(s)
    fraction_reads <- dplyr::group_by(reads_per_gene_sample, class) %>%
        dplyr::summarize(n=sum(read_sum)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(perc = n / sum(reads_per_gene_sample$read_sum) * 100) %>%
        tibble::column_to_rownames('class')
    reads_per_gene_sample <- dplyr::group_by(reads_per_gene_sample, class) %>%
        dplyr::mutate(sort=(as.integer(class) * maxreadsums + 1) / n() + 1/(read_sum + 1))
    levels(reads_per_gene_sample$class) <- sprintf("%s\n%.1f%%", levels(reads_per_gene_sample$class), fraction_reads[levels(reads_per_gene_sample$class), ]$perc)

    grid::grid.newpage()
    vp <- grid::viewport(gp=grid::gpar(lineheight=0.8))
    treemap::treemap(reads_per_gene_sample,
            index = c("class", "gene"),
            sortID="sort",
            vSize = "read_sum",
            type = "index",
            palette = palette,
            border.col=c("black","white"),
            border.lwds=c(1,0.5),
            #fontsize.labels=c(24,18),
            fontcolor.labels=c("black", "white"),
            fontface.labels=c(2,1),
            lowerbound.cex.labels=0.2,
            bg.labels=c("transparent"),
            align.labels=list(
                c("right", "bottom"),
                c("center", "center")
            ),
            title = title,
            aspRatio = 10/10,
            vp=vp
    )
}
