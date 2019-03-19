`%:%` <- function(d, a) {
    ifelse(is.nan(d), a, d)
}

check_serp_class <- function(data) {
    stopifnot(inherits(data, 'serp_data'))
}

check_serp_features_class <- function(data) {
    stopifnot(inherits(data, 'serp_features'))
}

check_background_fit_class <- function(data) {
    stopifnot(inherits(data, 'serp_background_fit'))
}

print_list_name <- function(leader, newchar, name, elem, depth=0) {
    if (depth > 0) {
        if (nchar(leader) > 1)
            cat(sprintf("%s  %s\n", leader, name))
        else
            cat(name, sep='\n')
        leader <- paste0(leader, newchar)
    }
    if (inherits(elem, 'list')) {
        nms <- names(elem)
        for (i in 1:length(elem)) {
            print_list_name(leader, newchar, nms[i], elem[[i]], depth + 1)
        }
    }
}
