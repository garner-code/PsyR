#' Compute degrees of freedom for the within or between family of contrasts
#'
#' Given a list of within or between factors, and an object of the afex_aov
#' model class, this function will compute the degrees of freedom for that
#' family For further details see p. 216-217 from Bird (2002),
#' https://doi.org/10.1177/0013164402062002001.
#'
#' @param model an anova model object of the afex_aov class
#' @param fctrs a list of characters or strings - list of relevant factor names
#'
#' @return a single value which is the degrees of freedom for that family
#' @export
#'
#' @examples
#' data(spacing)
#' mod <- aov_ez("subj", "yield", spacing,
#'                within = "spacing", between = "group")
#' compute_df(mod, list("group")) # get the degrees of freedom for the
#' between family
#' compute_df(mod, list("spacing")) # get df for the within family
#'
#' or aritrary example to show you can add more than one factor to the list:
#' mod <- aov_ez(some data and inputs)
#' fctrs <- list("groupA", "groupB", "groupC")
#' compute_df(mod, fctrs)
compute_df <- function(model, fctrs){

  if(!class(model) == "afex_aov"){
    stop("Error: model needs to be of class afex_aov")
  }
  if(!is.list(fctrs)){
    stop("Error: factor names should be provided as a list")
  }

  fctrs = unlist(fctrs)
  data = model$data$long
  N = unlist(lapply(fctrs, get_levels, data = data))
  df = N-1
  prod(df)
}
