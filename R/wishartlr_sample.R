#' Simulate largest roots of Wishart matrices
#'
#' The `wishartlr_sample` function simulates the largest roots of
#' Wishart matrices given the matrix dimensions and number of
#' simulations and calculates the first three moments of the resulting
#' distribution.
#'
#' @param p Integer. The first matrix dimension.
#' @param q Integer. The second matrix dimension.
#' @param n_sim Integer. The number of simulations to run. Default is
#' 10000.
#' @param seed Integer or NULL. The random seed for reproducibility.
#' Default is NULL.
#'
#' @return A list containing the following elements:
#' \item{mu1}{The mean of the largest roots.}
#' \item{mu2}{The second moment (variance) of the largest roots. }
#' \item{mu3}{The third moment (skewness) of the largest roots.}
#' \item{largest_roots}{A numeric vector of the largest roots simulated.
#' Only returned if `return_samples` is TRUE.}
#'
#' @details The function generates Wishart matrices with an identity
#' scale matrix and computes the largest eigenvalue (root) of each
#' matrix. The Wishart distribution is commonly used in multivariate
#' statistics, particularly in covariance matrix estimation.
#'
#' @examples
#' # Simulate largest roots of Wishart matrices
#' wishart_moments <- wishartlr_sample(p = 5, q = 10, n_sim = 1000, seed =
#' 123)
#'
#' @export
wishartlr_sample <- function(p, q, n_sim = 100000, return_samples = FALSE,
                             seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  largest_roots <- numeric(n_sim)

  for (i in 1:n_sim) {
    # Generate Wishart matrix with identity scale matrix
    # W ~ W_p(q, I)
    X <- matrix(rnorm(p * q), nrow = q, ncol = p)
    W <- crossprod(X)  # W = X'X ~ Wishart(q, I)

    # Extract eigenvalues
    evals <- eigen(W, symmetric = TRUE, only.values = TRUE)$values

    # Store largest root
    largest_roots[i] <- max(evals)
  }
  mu1 <- mean(largest_roots)
  mu2 <- mean((largest_roots - mu1)^2)
  mu3 <- mean((largest_roots - mu1)^3)
  if (return_samples) {
    result <- list(mu1 = mu1,
                   mu2 = mu2,
                   mu3 = mu3,
                   largest_roots = largest_roots)
  } else {
    result <- list(mu1 = mu1,
                   mu2 = mu2,
                   mu3 = mu3)
  }
  return(result)
}
