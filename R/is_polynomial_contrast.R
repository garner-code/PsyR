#' Identify a standard polynomial contrast
#'
#' Matches a single equally spaced polynomial trend over at least three levels,
#' allowing any nonzero scaling or sign reversal. Two-level comparisons and
#' all-zero vectors return `FALSE`. Custom contrasts with equivalent coefficients
#' also match; unequal level spacing is not inferred from numeric level labels.
#'
#' @param coefficients A non-empty finite numeric coefficient vector, in level order.
#' @param tol Positive finite relative comparison tolerance.
#' @return A single logical value.
#' @export
#' @examples
#' is_polynomial_contrast(c(-1, 0, 1))
#' is_polynomial_contrast(c(1, -1))
is_polynomial_contrast <- function(coefficients, tol = sqrt(.Machine$double.eps)) {
    if (!is.numeric(coefficients) || !is.null(dim(coefficients)) ||
        !length(coefficients) || any(!is.finite(coefficients))) {
        stop("`coefficients` must be a non-empty finite numeric vector.", call. = FALSE)
    }
    if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
        stop("`tol` must be a positive finite number.", call. = FALSE)
    }
    .polynomial_degree(coefficients, tol) > 0L
}

.polynomial_degree <- function(coefficients, tol) {
    magnitude <- max(abs(coefficients))
    if (length(coefficients) < 3L || magnitude == 0) {
        return(0L)
    }
    coefficients <- coefficients / magnitude
    basis <- stats::contr.poly(length(coefficients))
    basis <- sweep(basis, 2L, apply(abs(basis), 2L, max), "/")
    matches <- apply(basis, 2L, function(column) {
        min(max(abs(coefficients - column)), max(abs(coefficients + column))) <= tol
    })
    if (any(matches)) which(matches)[[1L]] else 0L
}

.polynomial_contrast_rows <- function(extraction, tol = sqrt(.Machine$double.eps)) {
    grid <- extraction$grid
    vapply(seq_len(ncol(extraction$coefficients)), function(index) {
        coefficients <- extraction$coefficients[, index]
        threshold <- tol * max(abs(coefficients))
        any(vapply(names(grid), function(factor_name) {
            values <- as.character(grid[[factor_name]])
            level_order <- if (is.factor(grid[[factor_name]])) {
                levels(grid[[factor_name]])[levels(grid[[factor_name]]) %in% values]
            } else {
                extraction$level_order[[factor_name]]
            }
            if (length(level_order) < 3L) return(FALSE)
            others <- setdiff(names(grid), factor_name)
            groups <- if (length(others)) {
                split(seq_along(values), interaction(grid[others], drop = TRUE))
            } else {
                list(seq_along(values))
            }
            degrees <- vapply(groups, function(rows) {
                if (all(abs(coefficients[rows]) <= threshold)) return(NA_integer_)
                positions <- match(level_order, values[rows])
                if (anyNA(positions) || anyDuplicated(values[rows])) return(0L)
                .polynomial_degree(coefficients[rows][positions], tol)
            }, integer(1))
            degrees <- degrees[!is.na(degrees)]
            length(degrees) > 0L && all(degrees > 0L) && length(unique(degrees)) == 1L
        }, logical(1)))
    }, logical(1))
}
