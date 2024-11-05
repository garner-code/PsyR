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
