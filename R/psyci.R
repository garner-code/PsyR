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
#' Note that this list should contain only orthogonal contrast tables, or non-orthogonal.
#' Do not mix or match orthogonal or non-orthogonal. If there is only one table of contrasts,
#' it is not necesssary to enter it as a list.
#' @param method a single string - method for confidence interval computation -
#' "ind", "bf", or "ph" - see further details below
#' @param family_list - a LIST the same length of contrast tables that specifies
#' if the associated contrast table is of the family "b", "w", or "bw"
#' @param between_factors - a list of between subject factor names, default = NA
#' @param within_factors - a list of within subject factor names, default = NA
#' @param alpha - desired type 1 error rate. Either a single value (e.g. alpha = .05)
#' that will be applied to each family, or a list of alphas of length
#' contrast_tables. Assumes that the alphas are provided in the order as should be
#' applied to contrast_tables. default = .05.
#' @param independent - a boolean value, default = TRUE. If TRUE, correction will
#' be applied to each family independently. If FALSE, all contrasts will be summed
#' together and treated as one family for the purposes of correction.
#' This is only relevant if method = "bf" or method = "ph", and will be ignored if method = "ind".
#' @param nu1 - if independent is set to FALSE, and method = "ph", then supply nu1 (df1) that
#' is used to adjust the critical constant to control for type 1 error (see eqs 7, 8, 9 of
#' Bird (2002) https://doi.org/10.1177/0013164402062002001). default = NA
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
#'         family_list = list("b"), between_factors = list("group"))
psyci <- function(model, contrast_tables, method, family_list,
                  between_factors = NA, within_factors = NA, alpha = 0.05,
                  independent = TRUE, nu1 = NA){

  if(!inherits(model, "afex_aov")){
    stop("Error: model needs to be of class afex_aov")
  }
  if (!method %in% c("ind", "bf", "ph")){
    stop("Error: method should be ind, bf, or ph")
  }

  if (!is.list(family_list)){
    stop("Error: family_list should be a list of family codes") # may update this later
  }

  # first check if a single contrast table has been entered, and if so, convert
  # to be a list
  types <- c("emmGrid", "summary_emm")
  if (sum(sapply(types, inherits, x=contrast_tables))){
    # if this is passed, then someone has entered a single object rather than a
    # list
    contrast_tables <- list(contrast_tables) # convert to list
  }
  test_emmGrids <- sum(sapply(contrast_tables, inherits, "emmGrid"))
  test_summary_emm <- sum(sapply(contrast_tables, inherits, "summary_emm"))
  if (test_emmGrids == length(contrast_tables) |
      test_summary_emm == length(contrast_tables)){
  } else {
    stop("Error: contrast table needs to be of class emmGrid or summary_emm")
  }

  if (!independent){
    if (method %in% "ph"){
      if (is.na(nu1)) {
        stop("Error: post-hoc method with non-orthogonal families requested but nu1 remains undefined. Please supply via `nu1=` argument")
      }
    }
  }

  contrast_info <- do.call(rbind, contrast_tables)
  all_contrast_tables <- summary(contrast_info) # create one big contrast table that one can index easily
  contrast_tables <- contrast_tables # keeping list
  # get the error df
  v_e = unique(all_contrast_tables$df)

  if (length(v_e) > 1){
    stop("Error: more than one contrast analysis passed in")
  }

  # if required, get v_b and v_w
  v_b = 1 # setting these as 1, in case, if not they will be updated in the if statement below
  v_w = 1

  if (method %in% "ph"){
    if (any(family_list == "b") | any(family_list == "bw")){

      if (independent){
        v_b = compute_df(model, fctrs = between_factors)
      } else {
        # here I need to get the factors, and the correct families
        v_b = nu1
      }
    }
    if (any(family_list == "w")  | any(family_list == "bw")){

      if (independent){
        v_w = compute_df(model, fctrs = within_factors)
      } else {
        v_w = nu1
      }
    }
  }

  # now get the appropriate critical constant, given requested
  # method
  # first, make a list for critical constants, that is as long
  # as the number of contrast tables that have been passed in as
  # contrast tables
  nfamilies = length(contrast_tables)
  if (length(family_list) == 1 & nfamilies > 1){
    family_list = rep(family_list, times=nfamilies)
  } else if (length(family_list) == nfamilies) {
    family_list = family_list
  } else {
    stop("Error: length of family_list is incompatible with length of contrast_tables")
  }

  # now make sure provided alphas are compatible with the number of families
  if (length(alpha) == 1){
    alphas = rep(list(alpha), length(contrast_tables))
  } else {
    alphas = alpha
  }
  if (length(alphas) != length(contrast_tables)){
    stop(sprintf("Error! Provided alpha is of length > 1 and does not match
                   length of contrast tables"))
  }
  names(alphas) = family_list

  if (method %in% "ind"){

    critical_constant = lapply(alphas, cc_ind_t, v_e=v_e)

  } else if (method %in% "bf") {

    # first, get the nk per family
    cs = get_contrast_coefficients(contrast_tables)
    nk = unlist(lapply(cs, function(x) ncol(x[,grep("c.*", names(x))])))
   if (!independent){ # if not independent families
     nk = rep(sum(nk), length(nk))
   }

   critical_constant = mapply(cc_bonf_t, n_k=nk, alpha=alphas, MoreArgs=list(v_e=v_e), SIMPLIFY=FALSE)

  } else if (method %in% "ph"){

    # here we make a list of critical constants, to be applied to construct confidence
    # intervals for each family of interest
    critical_constant = vector("list", length(contrast_tables)) # need to make a list for the critical constants
    # that we will collect, and we name them after the families
    names(critical_constant) = family_list

    if (any(family_list == "b")){

      cc_bs = lapply(alphas[names(alphas) %in% "b"], cc_ph_b, v_b=v_b, v_e=v_e)
      critical_constant[names(critical_constant) %in% "b"] = cc_bs
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
  contrasts_w_cis <- mapply(get_cis, contrast_tables, critical_constant, SIMPLIFY=FALSE)

  # now I want to add information about what happened as an attribute to the
  # list of contrast tables
  n_c_tables = length(contrast_tables)
  methods = as.list(rep(method, n_c_tables))
  families = family_list
  btwn_fctrs = wthn_fctrs = v_ws = v_bs = as.list(rep(NA, n_c_tables))
  names(btwn_fctrs) = names(wthn_fctrs) = names(v_ws) = names(v_bs) = families
  btwn_fctrs_idx = names(btwn_fctrs) %in% c("b", "bw")
  btwn_fctrs[btwn_fctrs_idx] = lapply(btwn_fctrs[btwn_fctrs_idx], function(x){
    unlist(between_factors)
  })
  wthn_fctrs_idx = names(wthn_fctrs) %in% c("w", "bw")
   wthn_fctrs[wthn_fctrs_idx] = lapply(wthn_fctrs[wthn_fctrs_idx], function(x){
     unlist(within_factors)
   })

  v_bs[names(v_bs) %in% c("b", "bw")] = v_b
  v_ws[names(v_ws) %in% c("w", "bw")] = v_w
  v_es = as.list(rep(v_e, n_c_tables))
  alphas = alphas

  contrasts_w_cis <- mapply(update_attributes, contrasts_w_cis, method = methods,
                            family=families, between_factors=btwn_fctrs,
                            within_factors=wthn_fctrs, v_b=v_bs, v_w=v_ws,
                            v_e=v_es, alpha=alphas, SIMPLIFY = FALSE)
  names(contrasts_w_cis) = families

  return(contrasts_w_cis)
}
