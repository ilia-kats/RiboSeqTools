load_experiment <- function(..., bin=c('bynuc', 'byaa'), exclude=NULL) {
    bin <- match.arg(bin, several.ok=TRUE)
    ret <- mapply(function(...) {
        paths <- rlang::list2(...)
        ret <- sapply(paths, function(path) {
            ret <- list()
            m = as.matrix(read.table(file.path(datadir, path) , fill=TRUE, header=TRUE, sep=",", row.names=1))
            acols <- ncol(m) %% 3
            if (acols > 0) {
                acols <- ncol(m):(ncol(m) - acols + 1)
                if (all(is.na(m[,acols])))
                m <- m[,-acols]
            }
            m[is.na(m)] <- 0
            libsize <- sum(m, na.rm=TRUE)
            if ('bynuc' %in% bin)
                ret$bynuc <- m
            if ('byaa' %in% bin) {
                ret$byaa <- m[,seq(1, ncol(m), by=3)] + m[,seq(2, ncol(m), by=3)] + m[,seq(3, ncol(m), by=3)]
            }
            if (!is.null(exclude) && length(exclude) > 0) {
                ret <- lapply(ret, function(x)Matrix::Matrix(x[!(rownames(x) %in% exclude),], sparse=TRUE))
            }
            else {
                ret <- lapply(ret, function(x)Matrix::Matrix(x, sparse=TRUE))
            }
            ret
        }, simplify=FALSE)
        ret
    },  ..., SIMPLIFY=FALSE)
    names(ret) <- 1:length(path_ip)
    ret
}

#' @export
load_serp <- function(..., ref, normalize=FALSE, bin=c('bynuc', 'byaa'), exclude=NULL, defaults=list()) {
    experiments <- rlang::list2(...)
    data <- sapply(experiments, load_experiment, bin=bin, exclude=exclude, simplify=FALSE)
    what <- ifelse('bynuc' %in% bin, 'bynuc', 'byaa')
    total <- sapply(data, function(exp) {
        sapply(exp, function(rep) {
            sapply(rep, function(sample) {
                sum(sample[[what]])
            }, simplify=FALSE)
        }, simplify=FALSE)
    }, simplify=FALSE)

    ref$cds_length <- ref$length %/% 3
    ret <- list(ref=ref, data=data, total=total)
    ret <- structure(ret, class="serp_data")
    if (normalize)
        ret <- normalize(ret)

    if (is.null(defaults$bin)) {
        if ('byaa' %in% bin)
            autodefaults$bin <- 'byaa'
        else
            defaults$bin <- 'bynuc'
    }
    ret$defaults <- defaults

    ret
}
