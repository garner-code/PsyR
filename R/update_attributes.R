#' Update attributes of contrast table with Psy CI info
#'
#' This function adds to the attributes of the contrast table, to let you know
#' the methods that have been applied by PsyR to compute confidence intervals
#' and, when requested by `psyci()`, procedure-specific adjusted p-values.
#' The added attributes tell you the key parameters that were used. Specifically,
#' the method applied, the family of contrasts (if relevant), factors that were
#' interpreted as between or within factors (where relevant), degrees of freedom
#' (between and within post-hoc dimensions, where relevant), the alpha value
#' applied, and the critical constant value. The displayed `F_value` and
#' `p_value` always use an ordinary F distribution with one numerator degree of
#' freedom; individual and Bonferroni procedures do not use between or within
#' post-hoc dimensions. Roy GCR and SMR tables also receive a
#' `psyr_probability_parameters` attribute using the parameter names from
#' Psy's probability calculator.
#'
#' @param contrast_table a table of contrast values of class emmGrid
#' @param method  a single string - method for confidence interval computation -
#' "ind", "bf", or "ph" - see "?psyci" for further details
#' @param family a single string - family of contrasts - "b", "w" or "bw" -
#' see "?psyci" for further details
#' @param between_factors a list of between subject factor names, default = NA
#' @param within_factors a list of within subject factor names, default = NA
#' @param v_e single value. df error
#' @param v_w Single value giving the within dimension used by a post-hoc
#'   procedure. The default is `NA`.
#' @param v_b Single value giving the between dimension used by a post-hoc
#'   procedure. The default is `NA`.
#' @param alpha single value. applied alpha rate. default = .05
#' @param smr_params  a list of parameters used for the Studentized Maximum Root method,
#' if this methods was used. Must have named elements p and q. default = NULL
#' @param p_value_method Optional method code used to compute
#'   `p_adj`. When supplied, messages distinguish the ordinary
#'   `F_value` and `p_value` from the procedure-specific adjusted probability
#'   and confidence interval. Default is `NULL` for direct calls that have not
#'   computed PsyR p-values.
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
                              alpha = 0.05,
                              smr_params = NULL,
                              p_value_method = NULL){

  method_code <- method
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

  # set full names for methods used in messages and attributes
  method = .psyci_method_label(method, family)
  p_adjustment = if (is.null(p_value_method)) {
    NULL
  } else {
    .psyci_adjustment_label(p_value_method, family)
  }

  p_value_message = character()
  if (!is.null(p_adjustment)) {
    ordinary_message <- paste0(
      "PsyR reports `F_value` as the ordinary F(1, df2) contrast statistic."
    )
    adjusted_message <- if (identical(p_adjustment, "none")) {
      paste0(
        "PsyR p-value adjustment: none. `p_adj` is unadjusted, ",
        "and the confidence interval is individual."
      )
    } else {
      paste0(
        "PsyR simultaneous-inference procedure: ", p_adjustment,
        ". `p_adj` and the confidence interval use this procedure."
      )
    }
    p_value_message <- c(ordinary_message, adjusted_message)
  }

  probability_parameters <- .psyci_probability_parameters(
    contrast_table = contrast_table,
    method = method_code,
    family = family,
    v_b = v_b,
    v_w = v_w,
    v_e = v_e,
    smr_params = smr_params
  )
  probability_parameter_message <- character()
  if (!is.null(probability_parameters)) {
    attr(
      contrast_table,
      "psyr_probability_parameters"
    ) <- probability_parameters
    probability_parameter_message <- .psyci_probability_parameter_message(
      method = method_code,
      family = family,
      parameters = probability_parameters
    )
  }

  # apply message update that will be applied to all contrast tables
  attr(contrast_table, "mesg") <- c(attr(contrast_table, "mesg"),
                                    paste("PsyR CI method:", method, "has been applied"),
                                    p_value_message,
                                    probability_parameter_message,
                                    paste("Family-wise correction assumes current contrasts are:", family_full),
                                    paste("PsyR used an alpha rate of:", alpha)
                                    )
  # now apply specific attributes depending on currently used family
  if (btwn_msg){
    attr(contrast_table, "mesg") <- c(
      attr(contrast_table, "mesg"),
      paste(
        "PsyR assumed between subject factor(s) are:",
        paste(between_factors, collapse=", ")
      )
    )
    if (method_code == "ph") {
      attr(contrast_table, "mesg") <- c(
        attr(contrast_table, "mesg"),
        paste("PsyR used post-hoc between dimension of:", v_b)
      )
    } else if (method_code == "smr") {
      attr(contrast_table, "mesg") <- c(
        attr(contrast_table, "mesg"),
        paste("PsyR used SMR parameter p of:", smr_params$p),
        paste("PsyR used SMR parameter q of:", smr_params$q)
      )
    }
  }

  if (wthn_msg){
    attr(contrast_table, "mesg") <- c(
      attr(contrast_table, "mesg"),
      paste(
        "PsyR assumed within subject factor(s) are:",
        paste(within_factors, collapse=", ")
      )
    )
    if (method_code == "ph") {
      attr(contrast_table, "mesg") <- c(
        attr(contrast_table, "mesg"),
        paste("PsyR used post-hoc within dimension of:", v_w)
      )
    }
  }

  attr(contrast_table, "mesg") <- c(attr(contrast_table, "mesg"),
                                    paste("PsyR used df error of:", v_e)
  )
  contrast_table
}

.psyci_probability_parameters <- function(
    contrast_table,
    method,
    family,
    v_b,
    v_w,
    v_e,
    smr_params
) {
  if (identical(method, "ph") && identical(family, "bw")) {
    return(list(
      s = min(v_b, v_w),
      m = (abs(v_b - v_w) - 1) / 2,
      n = (v_e - v_w - 1) / 2
    ))
  }
  if (identical(method, "smr")) {
    k <- if ("F_value" %in% names(contrast_table)) {
      sqrt(unname(contrast_table[["F_value"]]))
    } else {
      numeric()
    }
    return(list(
      k = k,
      p = smr_params$p,
      q = smr_params$q,
      df = v_e
    ))
  }
  NULL
}

.psyci_probability_parameter_message <- function(method, family, parameters) {
  if (identical(method, "ph") && identical(family, "bw")) {
    return(paste0(
      "Psy probability-calculator parameters for Roy GCR: s = ",
      parameters$s, ", m = ", parameters$m, ", n = ", parameters$n,
      ". They are stored in `psyr_probability_parameters`; its row-aligned ",
      "`theta` attribute entry contains the values for Psy's GCR field."
    ))
  }
  if (identical(method, "smr")) {
    return(paste0(
      "Psy probability-calculator parameters for SMR: row-specific ",
      "k = sqrt(F_value), p = ", parameters$p, ", q = ", parameters$q,
      ", df = ", parameters$df,
      ". They are stored in `psyr_probability_parameters`."
    ))
  }
  character()
}

.psyci_adjustment_label <- function(method, family) {
  if (method == "ind") {
    "none"
  } else if (method == "bf") {
    "Bonferroni"
  } else if (method == "ph") {
    if (identical(family, "b")) {
      "Scheffe"
    } else if (identical(family, "w")) {
      "Hotelling T-squared"
    } else if (identical(family, "bw")) {
      "Roy GCR"
    } else {
      "post-hoc"
    }
  } else if (method == "smr") {
    "SMR"
  } else {
    "unknown"
  }
}

.psyci_method_label <- function(method, family) {
  if (method == "ind") {
    "independent"
  } else if (method == "bf") {
    "Bonferroni"
  } else if (method == "ph") {
    if (identical(family, "b")) {
      "Scheffe"
    } else if (identical(family, "w")) {
      "post-hoc within (Hotelling's T-squared)"
    } else if (identical(family, "bw")) {
      "post-hoc between x within (Roy's GCR)"
    } else {
      "post-hoc"
    }
  } else if (method == "smr") {
    "Studentized Maximum Root (SMR)"
  } else {
    "unknown method"
  }
}
