#' @return  A \link[tibble]{tibble} with the following columns: \describe{
#'      \item{exp}{Experiment.}
#'      \item{rep}{Replicate.}
#'      \item{gene}{Gene.}
#'      \item{start}{Start of continuous binding region. If the binning mode for \code{\link{fit_background_model}}
#'           was \code{byaa}, this will be in codons, otherwise in nucleotides.}
#'      \item{end}{End of continuous binding region.}
#'      \item{width}{Width of continuous binding region.}
#'      \item{sample1}{Sum of counts of sample1 within the binding region. Note that the actual column name
#'              is the value of \code{sample1}<%= print_var_if_exists(sample_colname_suffix) %>.}
#'      \item{sample2}{Sum of counts of sample2 within the binding region. Note that the actual column name
#'              is the value of \code{sample1}<%= print_var_if_exists(sample_colname_suffix) %>.}
#'      <%= print_var_if_exists(additional_columns) %>
#' }
