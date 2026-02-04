#' Update attributes of contrast table with Psy CI info
#'
#' This function adds to the attributes of the contrast table, to let you know
#' the method that has been applied by PsyR to compute confidence intervals.
#' The added attributes tell you the key parameters that were used. Specifically,
#' the method applied, the family of contrasts (if relevant), factors that were
#' interpreted as between or within factors (where relevant), degrees of freedom
#' (between (b), within (w), and b x w, where relevant), the alpha value applied,
#' and the critical constant value.
#'
#' @param contrast_table a table of contrast values of class emmGrid
#' @param method  a single string - method for confidence interval computation -
#' "ind", "bf", or "ph" - see "?psyci" for further details
#' @param family a single string - family of contrasts - "b", "w" or "bw" -
#' see "?psyci" for further details
#' @param between_factors a list of between subject factor names, default = NA
#' @param within_factors a list of within subject factor names, default = NA
#' @param v_e single value. df error
#' @param v_w single value. df within. default = NA
#' @param v_b single value. df between. default = NA
#' @param alpha single value. applied alpha rate. default = .05
#'
#' @returns an emmeans contrast table with attributes updated
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
#' contrast_table <- summary(btwn_con)
#' v_e = contrast_table$df[1]
#' v_b = compute_df(mod, fctrs=list("group"))
#' critical_constant = cc_ph_b(v_b=v_b, v_e=v_e, alpha=0.05)
#' update_attributes(contrast_table, method="ph", family="b",
#'   between_factors=list("group"), v_e=v_e,
#'   v_b=v_b, alpha=.05)
update_attributes <- function(contrast_table, method, family = NA,
                              between_factors = NA, within_factors = NA,
                              v_e, v_w = NA, v_b = NA,
                              alpha = 0.05){

  btwn_msg = FALSE
  wthn_msg = FALSE
  # set full family names for messages
  if (family == "b"){
    family_full = "between subject contrasts"
    btwn_msg = TRUE
  } else if (family == "w"){
    family_full = "within subject contrasts"
    wthn_msg = TRUE
  } else if (family == "bw"){
    family_full = "between x within subject contrasts"
    btwn_msg = TRUE
    wthn_msg = TRUE
  } else {
    family_full = "unknown family" #prob not necessary as earlier messages should catch this
  }

  # set full names for method for messages
  if (method == "ind"){
    method = "independent"
  } else if (method == "bf"){
    method = "Bonferroni"
  } else if (method == "ph"){
    if (family == "b"){
      method = "Scheffe"
    } else if (family == "w"){
      method = "post-hoc within"
    } else if (family == "bw"){
      method = "post-hoc between x within (Roy's GCR)"
    } else
    method = "post-hoc"
  } else {
    method = "unknown method" #prob not necessary as earlier messages should catch this
  }

  # apply message update that will be applied to all contrast tables
  attr(contrast_table, "mesg") <- c(attr(contrast_table, "mesg"),
                                    paste("PsyR CI method:", method, "has been applied"),
                                    paste("Family-wise correction assumes current contrasts are:", family_full),
                                    paste("PsyR used an alpha rate of:", alpha),
                                    paste("PsyR used df error of:", v_e)
                                    )
  # now apply specific attributes depending on currently used family
  if (btwn_msg){
    attr(contrast_table, "mesg") <- c(attr(contrast_table, "mesg"),
                                      paste("PsyR assumed between subject factor(s) are:", between_factors),
                                      paste("PsyR used df between of:", v_b)
    )
  }

  if (wthn_msg){
    attr(contrast_table, "mesg") <- c(attr(contrast_table, "mesg"),
                                      paste("PsyR assumed within subject factor(s) are:", within_factors),
                                      paste("PsyR used df within of:", v_w)
    )
  }
  contrast_table
}
