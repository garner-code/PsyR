#' compute critical constant for individual CI
#'
#' Compute the critical constant for an individual planned
#' between subjects contrasts. This can be applied to either single
#' factor or between subjects designs, with v_e degrees of freedom
#' for error.
#'
#' This will control the noncoverage
#' error rate at the defined alpha, for each contrast in the set.
#'
#' WARNING: if you are using this method for more than one contrast
#' then your family-wise error rate will be greater than alpha
#'
#' See Bird (2002) https://doi.org/10.1177/0013164402062002001 p. 205.
#'
#' @param v_e single value
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
cc_individual <- function(v_e, alpha = 0.05){

  qt(alpha/2, v_e, lower.tail = FALSE)
}
