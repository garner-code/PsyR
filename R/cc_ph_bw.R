#' Critical constant for post-hoc simultaneous between x within CIs
#'
#' Computes the critical constant for use in calculating post-hoc simultaneous
#' confidence intervals for between x within contrast effects. For further
#' details see equation 9 (p. 217) from Bird (2002),
#' https://doi.org/10.1177/0013164402062002001.
#'
#' @param v_w a single numeric value - df within
#' @param v_b a single numeric value - df between
#' @param v_e a single numeric value - df error
#' @param alpha desired type 1 error rate e.g. p=.05
#'
#' @return a single numeric value that is the critical constant for use in the
#' construction of post-hoc simultaneous confidence intervals for between x
#' within subject effects
#' @export
#'
#' @examples
#' N = 16 # total number of subjects
#' W = 3 # number of levels within
#' J = 4 # number of groups
#' v_w = W - 1
#' v_b = J - 1
#' v_e = N - J
#' alpha = 0.05
#' cc_ph_bw(v_w, v_b, v_e, alpha)
cc_ph_bw <- function(v_w, v_b, v_e, alpha){

  s = min(v_b, v_w)
  m = (abs(v_b - v_w)-1)/2
  n = (v_e - v_w - 1)/2

  theta = gcr_crit(alpha, s, m, n)
  sqrt((v_e*theta)/(1-theta))
}
