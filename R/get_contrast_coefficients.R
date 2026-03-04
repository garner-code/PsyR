#' Get the contrast coefficient vectors from a list of contrast tables output by
#' emmeans `contrast` function
#'
#' @param contrast_tables
#'
#' @returns a list of the same size as contrast tables. Each element of the list
#' is a c x k matrix of contrast coefficient vectors
#' @export
#'
#' @examples
#' data(spacing)
#' # fit model
#' mod <- afex::aov_ez("subj","yield", spacing, within = "spacing", between = "group")
#' sum_emm_btwn <- emmeans::emmeans(mod, "group") # get marginal means
#' con <- list( # define some contrasts
#'             "12vs34" = c(0.5, 0.5, -0.5, -0.5),
#'              "1vs2" = c(1, -1, 0, 0),
#'              "3vs4" = c(0, 0, 1, -1)
#'              )
#' btwn_con <- emmeans::contrast(sum_emm_btwn, con) # get contrast table
#' get_contrast_coefficients(btwn_con)
get_contrast_coefficients <- function(contrast_tables){

  cs <- lapply(contrast_tables, stats::coef)
  lapply(cs, function(x) x[, grep("^c\\.[0-9]+$", colnames(x))])
}
