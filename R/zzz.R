`%:%` <- function(d, a) {
    ifelse(is.nan(d), a, d)
}
check_serp_class <- function(data) {
    stopifnot(inherits(data, 'serp_data'))
}
