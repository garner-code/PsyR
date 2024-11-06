#' Compute degrees of freedom for the within or between family of contrasts
#'
#' Given a list of within or between factors, and an object of the afex_aov
#' model class, this function will compute the degrees of freedom for that
#' family For further details see p. 216-217 from Bird (2002),
#' https://doi.org/10.1177/0013164402062002001.
#'
#' @param model an anova model object of the afex_aov class
#' @param data instead of the anova model you can use a long form data frame
#' @param fctrs a list of characters or strings - list of relevant factor names
#'
#' @return a single value which is the degrees of freedom for that family
#' @export
#'
#' @examples
#' data("spacing")
#' spacing$group <- as.factor(spacing$group)
#' spacing$spacing <- as.factor(spacing$spacing)
#'  # get the degrees of freedom for the between family
#'  data("spacing")
#'  spacing$group <- as.factor(spacing$group)
#'  spacing$spacing <- as.factor(spacing$spacing)
#'  compute_df(model=NA, data=spacing, fctrs=list("group"))
compute_df <- function(model, data = NA, fctrs){

  if (!is.na(model)){
    if(!inherits(model, "afex_aov")){
      stop("Error: model needs to be of class afex_aov")
    }
  }
  if (is.na(model)){
    if(!inherits(data, "data.frame")){
      stop("Error: data needs to be of class data.frame")
    }
  }
  if(!is.list(fctrs)){
    stop("Error: factor names should be provided as a list")
  }

  fctrs = unlist(fctrs)
  if (!is.na(model)) {
      data = model$data$long
  } else if (is.na(model)){
      data = data
  }
  N = unlist(lapply(fctrs, get_levels, data = data))
  df = N-1
  prod(df)
}
