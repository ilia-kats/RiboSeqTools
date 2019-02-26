normalize.serp_data <- function(data, exclude=c()) {
    mapply(function(exp, texp) {
        mapply(function(rep, trep) {
            mapply(function(sample, tsample) {
                sapply(sample, function(bin) {
                    toinclude <- rownames(bin)
                    toinclude <- toinclude[!(toinclude %in% exclude)]
                    bin[toinclude,] / tsample * 1e6
                })
            }, rep, trep, SIMPLIFY=FALSE)
        }, exp, texp, SIMPLIFY=FALSE)
    }, get_data(data), get_total(data), SIMPLIFY=FALSE)
}

#' @export
normalize <- function(data, exclude=c()) {
    UseMethod("normalize")
}

#' @export
get_elbow_threshold <- function(xvals, yvals) {
    if (missing(yvals) && is.matrix(xvals) && ncol(xvals) == 2) {
        yvals <- xvals[,2]
        xvals <- xvals[,1]
    } else if (missing(yvals)) {
        stop("y coordinates required")
    }
    v <- c(diff(range(xvals)), -diff(range(yvals)))
    w <- mapply(function(x, y, xmin, ymax)c(x - xmin, y - ymax), xvals, yvals, min(xvals), max(yvals))
    x <- apply(w, 2, function(w)sqrt(sum((w - as.vector(v %*% w) * v / sum(v^2))^2)))
    w[1, which.max(x)]
}

#' @export
get_data <- function(data) {
    check_serp_class(data)
    data$data
}

#' @export
get_reference <- function(data) {
    check_serp_class(data)
    data$ref
}

#' @export
get_total_counts <- function(data) {
    check_serp_class(data)
    data$total
}

#' @export
get_defaults <- function(data) {
    check_serp_class(data)
    data$defaults
}

#' @export
set_defaults <- function(data, defaults) {
    check_serp_class(data)
    #TODO validate_defaults(defaults)
    data$defaults <- defaults
    data
}
