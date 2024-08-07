#' Critical constant for simultaneous post-hoc intervals for within-subject
#' contrasts
#'
#' Produce the critical constant for use in calculating simultaneous post-hoc
#' within-subject contrast effects. Note that this method controls the FWER at
#' alpha.
#'
#' For further details, see eq 8 from Bird (2002)
#' DOI: https://doi.org/10.1177/0013164402062002001
#'
#' @param v_w single numeric value - within subjects df: e.g. for a
#' single within subjects factor with 3 levels, v_w = 3-1, for a 2 x 2
#' within subjects design, with 4 levels for the first factor, and 2 levels for
#' the second factor, v_e = (4-1)(2-1) = 6
#' @param v_e single numeric value - residual degrees of freedom - N - J,
#' where N is the total sample size, and J = the number of between subjects
#' groups
#' @param alpha single numeric value - acceptable type 1 error rate at the
#' family level (aka all within-subjects contrasts)
#'
#' @return single numeric value - critical constant to be used to construct
#' simultaneous, within subject, post-hoc confidence intervals for a given
#' contrast
#' @export
#'
#' @examples
#' N = 16
#' J = 4
#' n_within = 3 # e.g. 1 factor with 3 levels
#' v_w = n_within - 1
#' v_e = N - J
#' alpha = 0.05
#' cc_ph_w(v_w, v_e, alpha)
cc_ph_w <- function(v_w, v_e, alpha){

  crit_F <- stats::qf(alpha, v_w, v_e-v_w+1, lower.tail = FALSE)
  sqrt((v_w*v_e)/(v_e-v_w+1) * crit_F)
}
