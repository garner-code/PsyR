#' Find the critical value theta (GCR), given alpha, s, m, & n
#'
#' This function takes a given alpha rate (e.g. p=.05), and finds the
#' corresponding critical value of theta - i.e. the greatest characteristic root
#' that would be expected p=0.05 of the time, if the null hypothesis is true,
#' and given the parameter values of s,m, and n (see below).
#' The greatest characteristic root (GCR) is the largest linear combination of
#' effects from a covariance matrix, and can be used as the critical value
#' used to compute confidence intervals for between x within ANOVA designs.
#'
#' For further info see:
#' Bird (2002) https://doi.org/10.1177/0013164402062002001
#' Bird & Hadzi-Pavlovic (1983) https://doi.org/10.1037/0033-2909.93.1.167
#' Pillai (1965) https://www.jstor.org/stable/2333693
#'
#' @param alpha alpha/type 1 error rate - e.g. .05
#' @param s min(v_b, v_w) where v_b = df between, and v_w = df within (Bird,
#' 2002, p. 217)
#' @param m (|v_b - v_w|-1)/2
#' @param n (v_e - v_w - 1)/2 where v_e = N-J, where N is the total sample
#' size and J is the number of groups
#'
#' @return a single value, giving the critical theta/GCR expected under the
#' null hypothesis (given the covariance matrix)
#' @export
#'
#' @examples
#' gcr_crit(.05, 1, -0.5, 96)
#'
gcr_crit <- function(alpha, s, m, n){

  search_gcr <- function(s, m, n, alpha, x) 1-gcr_cdf(s, m, n, x)-alpha
  stats::uniroot(search_gcr, interval = 0:1, s=s, m=m, n=n, alpha=alpha)$root
}
