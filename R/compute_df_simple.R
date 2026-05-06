compute_df_simple <- function(model, data = NA, simp_fx){

  if (!anyNA(model)){
    if(!inherits(model, "afex_aov")){
      stop("Error: model needs to be of class afex_aov")
    }
  }
  if (anyNA(model)){
    if(!inherits(data, "data.frame")){
      stop("Error: data needs to be of class data.frame")
    }
  }
  if(!is.list(fctrs)){
    stop("Error: factor names should be provided as a list")
  }

  fctrs = unlist(fctrs)
  if (!anyNA(model)) {
    data = model$data$long
  } else if (anyNA(model)){
    data = data
  }

  length(attributes(data[,fctr])$levels)

}
