#' Convert the critical constant for the SMR test into a critical constant for the
#' computation of confidence intervals. This takes the square root of the SMR parameter
#' output by smr_crit, as per Bird & Hadzi-Pavlovic (2005), page 360.
#'
#' @param alpha Numeric. The significance level (e.g., 0.05 for a 95%
#' confidence level).
#' @param p Integer. The matrix dimension (number of variables).
#' @param q Integer. The degrees of freedom for the Wishart distribution.
#' @param n Integer. The degrees of freedom for the test statistic.
#' @param n_sim Integer. The number of simulations to run for estimating
#' the moments. Default is 100000.
#' @param seed Integer or NULL. The random seed for reproducibility.
#' Default is NULL.
#'
#' @return Numeric. The critical value for the confidence interval computation.
#' @export
#'
#' @examples
#' critical_constant <- cc_smr(2, 2, 96, .05)
#' print(critical_constant)
cc_smr <- function(p, q, v_e, alpha, n_sim = 100000, seed = NULL){

  cc <- smr_crit(p=p, q=q, n=v_e, alpha=alpha,
                 n_sim = n_sim, seed = seed)
  cc <- sqrt(cc)
  cc
}
