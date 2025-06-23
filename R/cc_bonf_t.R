#' Compute critical constant using bonferroni t procedure
#'
#' This provides a critical constant for use in the computation of
#' simultaneous confidence intervals over raw effect sizes
#' that uses the Bonferroni-t approach for FWER. This approach provides
#' more conservative confidence intervals relative to other approaches
#' when planned contrasts k <= v_b, or the between degrees of freedom.
#' If k > v_b, then this approach should be preferred if it produces smaller CIs
#' relative to competitors.
#'
#' Simultaneous confidence intervals controlling for FWER
#' as per Bird (2002) p. 205. https://doi.org/10.1177/0013164402062002001
#'
#' @param v_e a single numeric value that is the df residual
#' @param n_k a single numeric value, total number of contrasts
#' for that family (between, within, between x within are each a family)
#' @param alpha a single numeric value - the acceptable type 1 error rate
#' for each family of contrasts. Default = 0.05.
#'
#' @return a single value that is the critical constant for that family,
#' using the Bonferroni-t method
#' @export
#'
#' @examples
#' N = 16 # total sample size
#' J = 4 # number of groups
#' v_e = N - J
#' n_k = 3 # number of contrasts in that family
#' alpha = 0.05
#' cc_bonf_t(v_e, n_k, alpha)
cc_bonf_t <- function(v_e, n_k, alpha){

  stats::qt(alpha/(2*n_k), v_e, lower.tail = FALSE)
}
