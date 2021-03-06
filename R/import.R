import_csv <- function(path) {
    m = as.matrix(read.table(path , fill=TRUE, header=TRUE, sep=",", row.names=1))
    acols <- ncol(m) %% 3
    if (acols > 0) {
        acols <- ncol(m):(ncol(m) - acols + 1)
        if (all(is.na(m[,acols])))
        m <- m[,-acols]
    }
    m[is.na(m)] <- 0
    Matrix::Matrix(m, sparse=TRUE)
}

test_hdf5 <- function(path) {
    h5 <- file(path, "rb", raw=TRUE)
    magicnr <- readBin(h5, "raw", n=8)
    close(h5)
    magicnr[1] == 0x89 &&
    magicnr[2] == 0x48 &&
    magicnr[3] == 0x44 &&
    magicnr[4] == 0x46 &&
    magicnr[5] == 0x0d &&
    magicnr[6] == 0x0a &&
    magicnr[7] == 0x1a &&
    magicnr[8] == 0x0a
}

require_h5 <- function() {
    if (!requireNamespace("rhdf5", quietly=TRUE))
        rlang::abort("rhdf5 package not found, please install the rhdf5 package to read HDF5 files")
}

get_hdf5_genes <- function(f) {
    require_h5()
    rhdf5::h5ls(f, recursive=FALSE, all=FALSE, datasetinfo=FALSE, order="H5_ITER_NATIVE") %>%
        dplyr::filter(otype == "H5I_DATASET") %>%
        dplyr::pull(name)
}

make_ref_from_hdf5 <- function(paths) {
    require_h5()
    purrr::map_dfr(paths, function(path) {
        f <- rhdf5::H5Fopen(path)
        genes <- get_hdf5_genes(f)
        genes <- purrr::map_dfr(rlang::set_names(genes), function(g) {
            attrs <- rhdf5::h5readAttributes(f, g)
            if ("gene" %in% names(attrs)) {
                ga <- attrs[["gene"]]
                if (ga != g) {
                    attrs$gene_alt <- ga
                    attrs$gene <- NULL
                }
            }
            tibble::as_tibble(attrs) %>%
                dplyr::rename(length=cds_length)
        }, .id='gene')
        rhdf5::H5Fclose(f)
        genes
    }) %>%
        dplyr::distinct()
}

import_hdf5 <- function(path, ref) {
    require_h5()
    dsets <- rhdf5::h5dump(path, recursive=FALSE, load=TRUE, order="H5_ITER_NATIVE")
    lengths <- sapply(dsets, function(x)ifelse(is.matrix(x), nrow(x), 0))
    idx <- which(lengths > 0)
    dsets <- dsets[idx]
    lengths < lengths[idx]
    m <- do.call(rbind, dsets)
    i <- rep(1:length(idx), times=lengths)

    lengths <- ref$length[match(names(dsets), ref$gene)]
    Matrix::sparseMatrix(i=i, j=m[,1], x=m[,2], dims=c(length(dsets), max(lengths)), dimnames=list(names(dsets), NULL))
}

load_experiment <- function(..., .ref, .bin=c('bynuc', 'byaa'), .exclude=NULL) {
    bin <- match.arg(.bin, several.ok=TRUE)
    ret <- mapply(function(...) {
        paths <- rlang::list2(...)
        ret <- sapply(paths, function(path) {
            ret <- list()

            if (test_hdf5(path))
                m <- import_hdf5(path, .ref)
            else
                m <- import_csv(path)

            libsize <- sum(m, na.rm=TRUE)
            if ('bynuc' %in% .bin)
                ret$bynuc <- m
            if ('byaa' %in% .bin) {
                mod <- ncol(m) %% 3
                if (mod) {
                    rlang::warn(sprintf("Length of at least one CDS in %s is not divisible by 3", path))
                    m <- m[,1:(ncol(m) - mod)]
                }
                ret$byaa <- m[,seq(1, ncol(m), by=3), drop=FALSE] + m[,seq(2, ncol(m), by=3), drop=FALSE] + m[,seq(3, ncol(m), by=3), drop=FALSE]
            }
            if (!is.null(.exclude) && length(.exclude) > 0)
                ret <- lapply(ret, function(x) x[!(rownames(x) %in% .exclude),])
            ret
        }, simplify=FALSE)
        ret
    },  ..., SIMPLIFY=FALSE)
    names(ret) <- 1:unique(purrr::map_int(list(...), length))
    ret
}

#' Import ribosome profiling data
#'
#' Import read count tables. Currently, CSV and HDF5 files are accepted.
#'
#' CSV files are expected to have one row per ORF with the first column containing the ORF names.
#' Other columns represent positions from the 5' end of of the ORF in nucleotides and must contain
#' integer-valued read counts. A header must be present.
#'
#' HDF5 files are assumed to contain one data set per gene at the top level. Each data set must be
#' a two-column matrix, the first column containing the position from the 5' end of the ORF in nucleotides
#' and the second column containing the associated integer-valued read counts. Data set names are assumed
#' to be gene names.
#'
#' @param ... Name-value pairs of lists. The name of each argument will be the name of an experiment.
#'      The name of each element will be the sample type (e.g. TT for total translatome), the value
#'      of each element must be a character vector of file paths, where each file is a read count
#'      table of one replicate experiment. Replicate order must match between sample types.
#' @param ref Reference data frame containing at least the following columns:
#'      \describe{
#'          \item{gene}{Gene/ORF name. Must match the names given in the read count tables.}
#'          \item{length}{ORF length in nucleotides.}
#'      }
#'      If all input files are HDF5 files, this argument can be missing, in which case a refence is
#'      created from the union of all input files. For this to work, each HDF5 data set must have
#'      an attribute \code{cds_length}. Other HDF5 attributes will be included as additional columns. If
#'      HDF5 datasets have an attribute \code{gene}, the corresponding column will be named \code{gene_alt}
#'      to avoid conflicts with the \code{gene} column created from dataset names.
#' @param normalize Normalize the read counts to library size? Output will then be in RPM.
#' @param bin Bin the data. \code{bynuc}: No binning (i.e. counts per nucleotide). \code{byaa}: Bin by residue.
#' @param exclude Genes to exclude in all future analyses. This genes will also be excluded from total read count
#'      calculation. Note that the raw count tables will not be modified. Named list with names corresponding to
#'      experiments. If a character vector of gene names is given, these genes will be excluded from all
#'      experiments.
#' @param defaults Default parameters of the data set.
#' @return An object of class \code{serp_data}
#' @seealso \code{\link{serp_data_accessors}}, \code{\link{defaults}}
#' @examples \dontrun{
#'      data <- load_serp(DnaK=list(ip=c('data/dnak1_ip.csv', 'data/dnak2_ip.csv'),
#'                                  tt=c('data/dnak1_tt.csv', 'data/dnak2_tt.csv')),
#'                        TF=list(ip=c('data/tf1_ip.csv', 'data/tf2_ip.csv'),
#'                                tt=c('data/tf1_tt.csv', 'data/tf2_tt.csv')),
#'                        ref=reference_df,
#'                        bin='byaa')
#'      }
#' @export
load_serp <- function(..., ref, normalize=FALSE, bin=c('bynuc', 'byaa'), exclude=list(), defaults=list()) {
    experiments <- rlang::list2(...)
    if (is.null(names(experiments)) || sum(nchar(names(experiments)) == 0) > 0)
        rlang::abort("all experiments must be named")

    if (missing(ref) && all(sapply(unlist(experiments), test_hdf5))) {
        ref <- make_ref_from_hdf5(unlist(experiments))
    } else if (!('gene' %in% colnames(ref)) || !('length' %in% colnames(ref))) {
        rlang::abort("reference must have columns 'gene'  and 'length'")
    }

    bin <- match.arg(bin)
    if (!is.list(exclude) && is.character(exclude))
        exclude <- purrr::map(experiments, function(...)exclude)
    data <- sapply(experiments, purrr::lift_dl(load_experiment), .bin=bin, .ref=ref, simplify=FALSE)
    what <- ifelse('bynuc' %in% bin, 'bynuc', 'byaa')

    ref$cds_length <- ref$length %/% 3
    ret <- list(ref=ref, data=data, normalized=FALSE)
    ret <- structure(ret, class="serp_data")
    ret <- set_excluded(ret, exclude)
    ret <- set_total_counts(ret, calc_total_counts(ret, what))
    if (normalize)
        ret <- normalize(ret)
    defaults <- purrr::list_modify(.defaults, !!!defaults)

    if (defaults$bin == 'byaa' && !('byaa' %in% bin))
        defaults$bin <- 'bynuc'
    ret$defaults <- defaults

    set_downsampled(ret, FALSE)
}
