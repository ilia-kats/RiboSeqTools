`%:%` <- function(d, a) {
    ifelse(is.nan(d), a, d)
}
check_serp_class <- function(data) {
    stopifnot(inherits(data, 'serp_data'))
}

print_list_name <- function(leader, newchar, name, elem, depth=0) {
    if (depth > 0)
        cat(sprintf("%s  %s\n", leader, name))
    if (inherits(elem, 'list')) {
        nms <- names(elem)
        for (i in 1:length(elem)) {
            print_list_name(paste0(leader, newchar), newchar, nms[i], elem[[i]], depth + 1)
        }
    }
}
