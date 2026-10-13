# Deprecated internal compatibility helper. `psyci()` now identifies and
# removes non-contrast rows automatically.
junk_check <- function(contrast_tables){
  .Deprecated(
    new = "psyci",
    package = "PsyR",
    msg = paste(
      "`junk_check()` is deprecated; `psyci()` now identifies and removes",
      "non-contrast rows automatically."
    )
  )

  coefficients <- get_contrast_coefficients(contrast_tables)
  lapply(coefficients, function(values) {
    which(abs(colSums(values)) > sqrt(.Machine$double.eps))
  })
}
