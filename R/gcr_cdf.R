#' Cumulative density function for Roy's (1953) Greatest Characteristic Root
#' distribution
#'
#' Compute the cumulative density function for the largest root of a covariance
#' matrix, given s, m, n, and theta.
#'
#' See Chiani et al (2016) for in relation to GCR:
#' http://dx.doi.org/10.1016/j.jmva.2015.10.007
#'
#' and See Bird (2002) for definitions of s, m, & n
#' https://doi.org/10.1177/0013164402062002001
#'
#' @param s min(v_b, v_w) where v_b = df between, and v_w = df within (Bird,
#' 2002, p. 217)
#' @param m (|v_b - v_w|-1)/2
#' @param n (v_e - v_w - 1)/2 where v_e = N-J, where N is the total sample
#' size and J is the number of groups
#' @param theta the critical GCR value, for which the probability of is desired
#'
#' @return a single numeric value that is p(theta)
#' @export
#'
#' @examples
#' # a 2 x 2 factorial design where factor A is between and B is within, w
#' # 100 participants
#' s <- 1
#' m <- -0.5
#' n <- 96
#' theta <- 0.4
#' gcr_cdf(s,m,n,theta)
#'
gcr_cdf <- function(s, m, n, theta){
  if(s%%2 == 0){
    A <- matrix(0, nrow = s, ncol = s)
  } else {
    A <- matrix(0, nrow = s+1, ncol = s+1)
  }

  b <- c()
  for (i in 1:s){
    b[i] <- (ibeta(theta, m+i, n+1)^2)/2
    if((s-1)>=i) {
      for (j in i:(s-1)){
        b[j+1] <- (((m+j)/(m+j+n+1))*b[j])-(ibeta(theta, 2*m+i+j, 2*n+2)/(m+j+n+1))
        A[i, j+1] <- ibeta(theta, m+i, n+1)*ibeta(theta, m + j + 1, n + 1) - 2*b[j+1]
      }
    }
  }

  if (s%%2 != 0 ){
    for (i in 1:s){
      A[i, s+1] <- ibeta(theta, m+i, n+1)
    }

  }

  A = A - t(A)
  c_fun <- function(s, m, n) {
    gamma_product <- c()
    for (i in 1:s) {
      gamma_product[i] <- gamma((i+2*m+2*n+s+2)/2)/(gamma(i/2)*gamma((i+2*m+1)/2)*gamma((i+2*n+1)/2))
    }
    gamma_product <- prod(gamma_product)
    gamma_product*(pi^(s/2))
  }

  c_fun(s, m, n)*sqrt(det(A))
}
