#' Compute the standard error of a contrast
#'
#' This computes the estimated standard error of
#' a contrast, given the mean square error (MSE),
#' the contrast vector, and the relevant n.
#'
#' See the definition of sigma from equation 3 of
#' Bird (2002), p. 203
#'
#' @param MSE mean square error: a single numeric value
#' @param contrast_vector a vector defining the contrast
#' @param n n observations: a single numeric value
#'
#' @return SE: a single numeric value
#' @export
#'
#' @examples
#' x <- 2.2
#' y <- c(1, -1)
#' z <- 50
#' se_cont(x, y, z)
se_cont <- function(MSE, contrast_vector, n){

  sqrt(MSE * sum(contrast_vector^2/n))
}
