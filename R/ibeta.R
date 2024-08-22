#' Incomplete beta function
#'
#' Compute the incomplete beta function
#' Required to compute the greatest characteristic root
#' that is expected under the null hypothesis.
#'
#' See Chiani et al (2016) for in relation to GCR:
#' http://dx.doi.org/10.1016/j.jmva.2015.10.007
#'
#' And Boik & Robinson-Cox (1998) https://www.jstatsoft.org/article/view/v003i01
#' for more info on the incomplete beta function.
#'
#' @param x a single numeric value that is a given theta between 0 and 1
#' @param a a single numeric value, the alpha parameter from the beta distribution
#' @param b a single numeric value, the beta parameter as above
#'
#' @return a single numeric value that is p(theta) given the incomplete beta
#' @export
#'
#' @examples
#' x <- 0.4 # a theta to be evaluated
#' a <- 5 # alpha parameter (i.e. N hits)
#' b <- 10 # beta parameter (i.e. N not hits)
ibeta <- function(x,a,b){
  pbeta(x,a,b)*beta(a,b)
}
