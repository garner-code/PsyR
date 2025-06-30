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
#' @param contrast_tables a LIST of contrast tables generated using emmeans contrast()
#' @param method a single string - method for confidence interval computation -
#' "ind", "bf", or "ph" - see further details below
#' @param between_factors - a list of between subject factor names, default = NA
#' @param within_factors - a list of within subject factor names, default = NA
#' @param alpha - desired type 1 error rate
#'
#' @return an emmeans contrast table with confidence intervals added
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
#' # now add CI's using post-hoc method
#' psyci(model = mod, contrast_tables = list(btwn_con), method = "ph",
#'         family = "b", between_factors = list("group"))
psyci <- function(model, contrast_tables, method,
                  between_factors = NA, within_factors = NA, alpha = 0.05){

  if(!inherits(model, "afex_aov")){
    stop("Error: model needs to be of class afex_aov")
  }
  if (!method %in% c("ind", "bf", "ph")){
    stop("Error: method should be ind, bf, or ph")
  }

  test_emmGrids <- sum(sapply(contrast_tables, inherits, "emmGrid"))
  test_summary_emm <- sum(sapply(contrast_tables, inherits, "summary_emm"))
  if (test_emmGrids == length(contrast_tables) |
      test_summary_emm == length(contrast_tables)){
  } else {
    stop("Error: contrast table needs to be of class emmGrid or summary_emm")
  }

  contrast_info <- do.call(rbind, contrast_tables)
  all_contrast_tables <- summary(contrast_info) # create one big contrast table that one can index easily
  contrast_tables <- contrast_tables # keeping list
  # get the error df
  v_e = unique(all_contrast_tables$df)
  if (length(v_e) > 1){
    stop("Error: more than one contrast analysis passed in")
  }

  # here I am assuming I get a list, called 'family_list' for now
  # which designates which table from the list of contrast_tables contains
  # which type of contrasts


  # if required, get v_b and v_w
  if (method %in% "ph"){
    if (any(family_list == "b") | any(family_list == "bw")){

      v_b = compute_df(model, fctrs = between_factors)
    }
    if (any(family_list == "w")  | any(family_list == "bw")){

      v_w = compute_df(model, fctrs = within_factors)
    }
  }

  # now get the appropriate critical constant, given requested
  # method
  # first, make a list of critical constants, that is as long
  # as the number of contrast tables that have been passed in as
  # contrast tables
  critical_constant = list()
  length(critical_constant) = length(family_list)
  names(critical_constant) = family_list

  if (method %in% "ind"){

    critical_constant[1:length(critical_constant)] = rep(cc_ind_t(v_e=v_e, alpha=alpha),
                                                         length(critical_constant))

  } else if (method %in% "bf") {

    cs <- lapply(contrast_tables, stats::coef)
    nk = sum(unlist(lapply(cs, function(x) ncol(x[,grep("c.*", names(x))]))))
    critical_constant[1:length(critical_constant)] = rep(cc_bonf_t(v_e=v_e, n_k=nk, alpha=alpha),
                                                         length(critical_constant))

  } else if (method %in% "ph"){

    # here we make a list of critical constants, to be applied to construct confidence
    # intervals for each family of interest


    if (any(family_list == "b")){

      cc_b = cc_ph_b(v_b=v_b, v_e=v_e, alpha=alpha)
      critical_constant[["b"]] = cc_b

    }

    if (any(family_list == "w")) {

      cc_w = cc_ph_w(v_w=v_w, v_e=v_e, alpha=alpha)
      critical_constant[["w"]] = cc_w

    }

    if (any(family_list == "bw")) {

      cc_bw = cc_ph_bw(v_w=v_w, v_b=v_b, v_e=v_e, alpha=alpha)
      critical_constant[["bw"]] = cc_bw
    }
  }

  # now get SE from each contrast table in contrast_tables, compute CIs using appropriate
  # cc, and save results to a concatenated summary table
  contrasts_w_cis <- mapply(get_cis, contrast_tables, critical_constant)


  # # now I want to add information about what happened as an attribute to the
  # # contrast table, before spitting it back out to the user
  # # so if the method is ind or bf, then all other factors should be set to NA
  # if (method %in% c("ind", "bf")){
  #
  #   contrast_table <- update_attributes(contrast_table,
  #                                       method=method,
  #                                       v_e=v_e,
  #                                       alpha=alpha,
  #                                       critical_constant=round(critical_constant,3))
  # } else if (method %in% "ph") {
  #
  #   if (family %in% "b"){
  #
  #     contrast_table <- update_attributes(contrast_table,
  #                                         method=method,
  #                                         family=family,
  #                                         between_factors = between_factors,
  #                                         v_e=v_e,
  #                                         v_b=v_b,
  #                                         alpha=alpha,
  #                                         critical_constant=round(critical_constant,3))
  #
  #   } else if (family %in% "w") {
  #
  #     contrast_table <- update_attributes(contrast_table,
  #                                         method=method,
  #                                         family=family,
  #                                         within_factors = within_factors,
  #                                         v_e=v_e,
  #                                         v_w=v_w,
  #                                         alpha=alpha,
  #                                         critical_constant=round(critical_constant,3))
  #
  #
  #   } else if (family %in% "bw") {
  #
  #     contrast_table <- update_attributes(contrast_table,
  #                                         method=method,
  #                                         family=family,
  #                                         between_factors = between_factors,
  #                                         within_factors = within_factors,
  #                                         v_e=v_e,
  #                                         v_w=v_w,
  #                                         v_b=v_b,
  #                                         alpha=alpha,
  #                                         critical_constant=round(critical_constant,3))
  #
  #   }
  # }

  contrasts_w_cis
}
