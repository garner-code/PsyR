.validate_alpha_scalar <- function(alpha, argument = "alpha") {
    valid <- is.numeric(alpha) && length(alpha) == 1L &&
        !is.na(alpha) && is.finite(alpha) && alpha > 0 && alpha < 1
    if (!valid) {
        stop(
            sprintf(
                "`%s` must be one finite number strictly between 0 and 1.",
                argument
            ),
            call. = FALSE
        )
    }
    invisible(alpha)
}

.validate_alpha <- function(alpha) {
    if (is.list(alpha)) {
        if (length(alpha) == 0L) {
            stop("`alpha` must not be empty.", call. = FALSE)
        }
        invisible(lapply(alpha, .validate_alpha_scalar))
    } else {
        if (!is.numeric(alpha) || length(alpha) == 0L) {
            stop(
                paste(
                    "`alpha` must be a finite number, or a vector or list of",
                    "finite numbers, strictly between 0 and 1."
                ),
                call. = FALSE
            )
        }
        invisible(lapply(as.list(alpha), .validate_alpha_scalar))
    }
}
