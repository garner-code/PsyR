#' Probability calulcator for Roy's (1953) Greatest Characteristic Root
#' distribution
#'
#' Compute the probability of observing root of a covariance
#' matrix of size theta by chance, given s, m, n.
#'
#' See Chiani et al (2016) for in relation to GCR:
#' http://dx.doi.org/10.1016/j.jmva.2015.10.007
#'
#' and See Bird (2002) for definitions of s, m, & n
#' https://doi.org/10.1177/0013164402062002001
#'
#' @param v_w a single numeric value - df within
#' @param v_b a single numeric value - df between
#' @param v_e a single numeric value - df error
#' @param theta Roy's GCR statistic, for which the probability of is desired
#'
#' @return a single numeric value that is p(theta)
#' @export
#'
#' @examples
#' # a 2 x 2 factorial design where factor A is between and B is within, w
#' # 100 participants
#' v_w <- 1
#' v_b <- 1
#' v_e <- 96
#' theta <- 0.4
#' pgcr(theta, v_w, v_b, v_e)
pgcr <- function(theta, v_w, v_b, v_e){

  s = min(v_b, v_w)
  m = (abs(v_b - v_w)-1)/2
  n = (v_e - v_w - 1)/2

  1-gcr_cdf(s, m, n, theta)
}
