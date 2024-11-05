psyci <- function(contrast_table, method, family = NA, alpha = 0.05){

  if (!class(contrast_table)[[1]] == "emmGrid" |
        !class(contrast_table)[[1]] == "summary_emm"){
    stop("Error: contrast table needs to be of class emmGrid or summary_emm")
  }
  if (!method %in% c("ind", "bf", "ph")){
    stop("Error: method should be ind, bf, or ph")
  }
  if (method == "ph" & is.na(family)){
    stop("Error: define family if using ph method (b, w, or bw)")
  }
  if (method != "ph" & !is.na(family)){
    warning("Warning: family is defined but not with ph method.
            Ignoring family.")
  }
  if (!is.na(family) & !family %in% c("b", "w", "bw")){
    stop("Error: family should be b, w, or bw")
  }

  contrast_table = summary(contrast_table)
  v_e = contrast_table$df[1]
  if (method %in% "ind"){

    critical_constant = cc_individual(v_e=v_e)
  } else if (method %in% "bf") {

    nk = nrow(contrast_table)
    crtical_constant = cc_bonf_t(v_e=v_e, n_k=nk, alpha=alpha)
  } else if (method %in% "ph"){

    if (family %in% "b"){

      v_b = # this needs to be defined, ideally within the code
      critical_constant = cc_ph_b(v_b=v_b, v_e=v_e, alpha=alpha)

    } else if (family %in% "w") {

      v_e = # this needs to be defined, ideally within the code
      critical_constant = cc_ph_w(v_w=v_w, v_e=v_e, alpha=alpha)

    } else if (family %in% "bw") {

      v_b = #
      v_e = #
      critical_constant = cc_ph_bw(v_w, v_b, v_e, alpha)
    }
  }

  se_cont = contrast_table[,"SE"]
  cis = unlist(lapply(se_cont, contrast_ci, critical_constant))
  contrast_table$lower = contrast_table$estimate - cis
  contrast_table$upper = contrast_table$estimate + cis
  contrast_table
}
