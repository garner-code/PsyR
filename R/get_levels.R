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
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.", call. = FALSE)
    }
    if (!is.character(fctr) || length(fctr) != 1L || is.na(fctr) ||
        !nzchar(fctr)) {
        stop("`fctr` must be one non-empty column name.", call. = FALSE)
    }
    if (!fctr %in% names(data)) {
        stop(sprintf("Column '%s' was not found in `data`.", fctr), call. = FALSE)
    }

    factor_column <- data[[fctr]]
    if (!is.factor(factor_column)) {
        stop(sprintf("Column '%s' is not a factor.", fctr), call. = FALSE)
    }

    nlevels(factor_column)
}
