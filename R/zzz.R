`%:%` <- function(d, a) {
    ifelse(is.nan(d), a, d)
}

get_default_param <- function(serp_data, param, error=TRUE) {
    pname <- as.character(rlang::ensym(param))
    if (missing(param)) {
        param <- serp_data$defaults[[pname]]
    }
    if (is.null(param) && error)
        stop(sprintf("invalid %s argument", pname))
    param
}

check_serp_class <- function(data) {
    stopifnot(inherits(data, 'serp_data'))
}
