#' Check a list of contrast tables for unwanted contrasts
#'
#' This function will check the vectors of contrast co-efficients from a list of
#' emmeans objects that have been generated using the emmeans 'contrast'
#' function. It returns a list that contains an index of which contrast
#' coefficient vectors sum to > 0. Each element of the list corresponds to the
#' equivalent element of the list of contrast tables
#'
#' @param contrast_tables - a list of contrast tables created using the emmeans
#' 'contrast' function
#'
#' @returns junk - a list of indexes of contrast coefficient vectors that sum to > 0
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
#' junk = junk_check(list(btwn_con))
junk_check <- function(contrast_tables){
  cs = lapply(contrast_tables, stats::coef)
  cs = lapply(cs, function(x) apply(x[,grep("c.*", names(x))], 2, sum))
  junk = lapply(cs, function(x) which(x > 0))
  junk
}
