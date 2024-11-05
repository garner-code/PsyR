#' get number of levels in a factor
#'
#' from a long form data frame, get the number of levels
#' in the factor specified by factor
#' @param data a longform dataframe. Predictor variables should be factors
#' @param fctr a single string which matches the name of the factor in the
#' dataframe
#'
#' @return a single numeric value which is the number of levels in factor fctr
#' @export
#'
#' @examples
#' data = data.frame(A = c("X", "X", "X", "Y", "Y", "Y"))
#' data$A <- as.factor(data$A)
#' fctr = "A"
#' get_levels(data, fctr)
#'
get_levels <- function(data, fctr){
  if (!is.factor(data[,fctr][[1]])){
    stop(sprintf("Error: '%s' is not of type factor", fctr))
  }

  length(levels(data[,fctr][[1]]))
}
