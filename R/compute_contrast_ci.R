#' Compute the range over which the average effect would likely vary, to
#' be added/subtracted from the main effect to get the confidence intervals
#' for the contrast
#'
#' @param cc The critical constant, given the effect (Bird (2002), eq 3, term CC)
#' @param se_contrast The standard error of the contrast (Bird (2002), eq 3, sigma)
#'
#' @return the range of confidence, to be added to/subtracted from the effect
#' @export
#'
#' @examples
#' x <- 2.5
#' y <- 1.2
#' compute_contrast_ci(x, y)
#'
#' #' see Bird (2002) Confidence intervals for effect sizes in analysis of variance.
#' Educational and Psychological Measurement, 62(2).
#' https://doi.org/10.1177/00131644020620020
compute_contrast_ci <- function(cc, se_contrast){

  cc * se_contrast
}
