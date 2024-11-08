#' Calculate confidence intervals for contrasts of interest from ANOVA models
#'
#' This function computes confidence intervals for contrasts testing effects
#' from ANOVA models. The method used to compute the confidence intervals
#' is selected by the user (see below). The methods are described in Bird
#' (2002) (see below) https://doi.org/10.1177/0013164402062002001
#'
#' # METHODS
#' ind - individual t ci's - A set of individual 100(1 - α)% confidence
#' intervals on planned contrasts controls the noncoverage error rate at
#' α for each contrast in the set, thereby producing a familywise error
#' rate which is greater than α (see Bird, 2002, bottom of page 204).
#'
#' bf - bonferonni t ci's - only valid if contrasts are defined independently
#' of the data. Controls the familywise non-coverage error rate at alpha, given
#' distributional assumptions are met. (see Bird, 2002, pg. 205)
#'
#' ph - post-hoc ci's - valid for post-hoc contrasts, will control the
#' familywise non-coverage error rate at alpha, given distributional assumptions
#' are met (see Bird, 2002, pg. 217, eqs. 7-9).
#'
#' @param model an anova model made using the afex package (class afex_aov)
#' @param contrast_table a table of contrasts generated using emmeans contrast()
#' @param method a single string - method for confidence interval computation -
#' "ind", "bf", or "ph" - see further details below
#' @param family a single string - family of contrasts - "b", "w" or "bw" - see
#  further details below, default = NA
#' @param between_factors - a list of between subject factor names, default = NA
#' @param within_factors - a list of within subject factor names, default = NA
#' @param alpha - desired type 1 error rate
#'
#' @return an emmeans contrast table with confidence intervals added
#' @export
#'
#' @examples
#' x <- 10
psyci <- function(model, contrast_table, method, family = NA,
                  between_factors = NA, within_factors = NA, alpha = 0.05){

  if(!inherits(model, "afex_aov")){
    stop("Error: model needs to be of class afex_aov")
  }

  if (inherits(contrast_table, "emmGrid") |
        inherits(contrast_table, "summary_emm")){
  } else {
    stop("Error: contrast table needs to be of class emmGrid or summary_emm")
  }
  if (!method %in% c("ind", "bf", "ph")){
    stop("Error: method should be ind, bf, or ph")
  }
  if (method %in% "ph" & is.na(family)){
    stop("Error: define family if using ph method (b, w, or bw)")
  }
  if (method != "ph" & !anyNA(family)){
    warning("Warning: family is defined but not with ph method.
            Ignoring family.")
  }
  if (!anyNA(family) & !family %in% c("b", "w", "bw")){
    stop("Error: family should be b, w, or bw")
  }

  if (method %in% "w" | method %in% "bw"){
    if (!inherits(within_factors, "list")){
      stop("Error: within method required but no list of within factors
           supplied")
    }
  }

  if (method %in% "b" | method %in% "bw"){
    if (!inherits(between_factors, "list")){
      stop("Error: within method required but no list of between factors
           supplied")
    }
  }

  contrast_table = summary(contrast_table)
  v_e = contrast_table$df[1]

  # if required, get v_b and v_w
  if (method %in% "ph"){
    if (family %in% "b" | family %in% "bw"){

      v_b = compute_df(model, fctrs = between_factors)
    }
    if (family %in% "w" | family %in% "bw"){

      v_w = compute_df(model, fctrs = within_factors)
    }
  }

  if (method %in% "ind"){

    critical_constant = cc_ind_t(v_e=v_e)
  } else if (method %in% "bf") {

    nk = nrow(contrast_table)
    critical_constant = cc_bonf_t(v_e=v_e, n_k=nk, alpha=alpha)
  } else if (method %in% "ph"){

    if (family %in% "b"){

      v_b = # this needs to be defined, ideally within the code
      critical_constant = cc_ph_b(v_b=v_b, v_e=v_e, alpha=alpha)

    } else if (family %in% "w") {

      v_w = # this needs to be defined, ideally within the code
      critical_constant = cc_ph_w(v_w=v_w, v_e=v_e, alpha=alpha)

    } else if (family %in% "bw") {

      v_b = v_b
      v_w = v_w
      critical_constant = cc_ph_bw(v_w, v_b, v_e, alpha)
    }
  }

  se_cont = contrast_table[,"SE"]
  cis = unlist(lapply(se_cont, compute_contrast_ci, critical_constant))
  contrast_table$lower = contrast_table$estimate - cis
  contrast_table$upper = contrast_table$estimate + cis
  contrast_table
}
