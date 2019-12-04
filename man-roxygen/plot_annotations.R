#' @param higlightregion Regions in the plot to be highlighted. Can be either a list containing two-element
#'      numeric vectors or a two-column matrix. Rectangles spanning the entire height of the plot will be
#'      drawn between the given x positions.
#' @param highlightargs Named list of additional parameters to pass to \code{\link[ggplot2]{annotate}}.
#'      Each list element will be recycled to match the number of annotations.
