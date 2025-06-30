#' Calculate CIs for a contrast table (emmeans object) given a critical constant
#'
#' @param contrast_table - a table of contrasts, produced using the 'contrast'
#' function from emmeans
#' @param cc - a critical constrant, to be used to calculate CIs.
#'
#' @returns contrast_table - a summary table of the contrasts with the CI info added
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
#' cc <- cc_ind_t(v_e=12, alpha=0.05) # for example
#' get_cis(btwn_con, cc)
get_cis <- function(contrast_table, cc){
  contrast_table = summary(contrast_table)
  se_cont = contrast_table[,"SE"]
  cis = unlist(lapply(se_cont, compute_contrast_ci, cc))
  contrast_table$cc = cc
  contrast_table$lower = contrast_table$estimate - cis
  contrast_table$upper = contrast_table$estimate + cis
  contrast_table
}
