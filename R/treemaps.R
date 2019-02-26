#' @importFrom magrittr %$%
#' @export
plot_treemap <- function(data, exp, rep, sample, title, geneclass, palette="Set2", exclude=c()) {
    check_serp_class(data)

    if (ncol(geneclass) == 2 && !('class' %in% colnames(geneclass)) && 'gene' %in% colnames(geneclass))
        geneclass <- dplyr::rename(geneclass, class=dplyr::matches('gene'))
    classes <- c(levels(geneclass$class), 'unknown')

    reads_per_gene_sample <- tibble::enframe(rowSums(get_data(data)[[exp]][[rep]][[sample]][[1]], na.rm=TRUE), name='gene', value='read_sum') %>%
        dplyr::left_join(geneclass) %>%
        tidyr::replace_na(list(class='unknown')) %>%
        dplyr::mutate(class=factor(class, classes, ordered=TRUE)) %>%
        dplyr::filter(!(gene %in% exclude))

    maxreadsums <- dplyr::group_by(reads_per_gene_sample, class) %>%
        dplyr::summarize(s=sum(1/(read_sum + 1))) %$%
        max(s)
    fraction_reads <- dplyr::group_by(reads_per_gene_sample, class) %>%
        dplyr::summarize(n=sum(read_sum)) %>%
        dplyr::ungroup %>%
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
