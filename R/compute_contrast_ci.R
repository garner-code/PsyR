compute_contrast_ci <- function(cc, se_contrast){
  # compute the confidence interval using the
  # critical constant and the standard error of the
  # contrast
  # --
  # Kwargs
  # cc [int] - appropriate critical constant
  # se_contrast

  cc * se_contrast
}
