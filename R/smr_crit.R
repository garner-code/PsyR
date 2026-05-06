#' Compute the critical value for the Studentized Maxmimum Root (SMR)
#' distribution
#'
#' The `smr_crit` function calculates the critical value for a
#' statistical test using the moments of the largest roots of Wishart
#' matrices and an approximation to the F-distribution.
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
#' @return Numeric. The critical value for the statistical test.
#'
#' @details
#' The function first simulates the largest roots of Wishart matrices
#' using the `wishartlr_sample` function to estimate the first three
#' moments (mean, variance, and skewness). These moments are then used
#' to approximate the parameters of an F-distribution, which is used to
#' compute the critical value for the given significance level
#' (`alpha`).
#'
#' @examples
#' # Compute critical value for a test
#' crit_value <- smr_crit(alpha = 0.05, p = 5, q = 10, n = 20, n_sim = 100000, seed = 123)
#' print(crit_value)
#'
#' @seealso wishartlr_sample
#' @export
smr_crit <- function(alpha, p, q, n, n_sim = 100000, seed = NULL) {

  calc_f_params <- function(mu1, mu2, mu3, tau) {

    # Generate the parameters for the F-distribution based on the moments
    # mu1: first moment (mean)
    # mu2: second moment (variance)
    # mu3: third moment (skewness)
    # tau: degrees of freedom

    k <- (mu2 + mu1^2) * (tau - 2) / ((tau - 4) * mu1^2)
    u <- ((mu3 + 3 * mu2 * mu1 + mu1^3)
          * (tau - 2)^2 / ((mu1^3 * (tau - 4) * (tau - 6))))

    f_pi <- 2 * (k + 3 * u - 4 * k^2) / (k + u - 2 * k^2)
    f_gamma <- 2 * (f_pi - 2) / ((f_pi - 2) * (k - 1) - 2 * k)
    f_omega <- mu1 * tau * (f_pi - 2) / ((tau - 2) * f_pi)

    return(list(f_pi = f_pi, f_gamma = f_gamma, f_omega = f_omega))
  }

  moments <- wishartlr_sample(p = p, q = q, n_sim = n_sim, seed = seed)
  f_params <- calc_f_params(moments$mu1, moments$mu2, moments$mu3, n)
  return(qf(1 - alpha, df1 = f_params$f_gamma, df2 = f_params$f_pi)
         * f_params$f_omega)
}