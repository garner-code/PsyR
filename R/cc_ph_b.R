#' Critical constant for post-hoc simultaneous between CIs
#'
#' Computes the critical constant for use in calculating post-hoc simultaneous
#' confidence intervals for between groups contrast effects. For further
#' details see equation 7 (p. 217) from Bird (2002),
#' https://doi.org/10.1177/0013164402062002001.
#'
#' @param v_b - a single numeric value - df between
#' @param v_e - a single numeric value - df residual
#' @param alpha - a single numeric value - acceptable type 1 error rate for
#' the family of contrasts
#'
#' @return - a single numeric value that is the critical constant for use in the
#' construction of post-hoc simultaneous confidence intervals for between
#' subject effects
#' @export
#'
#' @examples
#' N = 16 # total number of subjects
#' J = 4 # number of groups
#' v_b = J-1
#' v_e = N - J
#' alpha = 0.05
#' cc_between(v_b, v_e, alpha)
cc_between <- function(v_b, v_e, alpha = 0.05){

  crit_F <- stat::qf(alpha, v_b, v_e, lower.tail = FALSE)
  sqrt(v_b * crit_F)
}
