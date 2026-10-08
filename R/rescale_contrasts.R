#' Rescale contrast coefficients using Psy conventions
#'
#' Rescales a named list of contrast coefficient vectors before it is supplied
#' to [emmeans::contrast()]. The available modes implement the two coefficient
#' transformations from the original Psy program: mean-difference contrasts
#' and interaction contrasts.
#'
#' @param contrasts A non-empty ordinary named list of equal-length numeric
#'   coefficient vectors. The coefficient order must match the order of the
#'   estimated marginal means to which the contrasts will be applied.
#' @param mode Scaling mode. `"mean_difference"` makes the positive
#'   coefficients sum to one. `"interaction"` scales coefficients for the
#'   order given by `interaction_order`.
#' @param interaction_order A whole number from 0 to 9. Order 0 is
#'   mean-difference scaling, order 1 is appropriate for a two-factor product
#'   interaction, and order 2 is appropriate for a three-factor product
#'   interaction. This must be 0 unless `mode = "interaction"`.
#' @param zero_sum How a non-zero-sum vector is handled in a rescaling mode.
#'   `"error"` stops before rescaling. `"warn"` reproduces Psy's permissive
#'   behaviour by warning and applying the legacy arithmetic when the positive
#'   coefficient sum is defined.
#' @param tol A finite relative tolerance greater than or equal to zero and less
#'   than one, used to determine whether a coefficient vector sums to zero.
#'
#' @details
#' For a coefficient vector `x`, let `P` be the sum of its positive
#' coefficients. Mean-difference mode returns `x / P`, while interaction mode
#' returns `2^interaction_order * x / P`. At least one positive coefficient is
#' therefore required, and all-zero or all-nonpositive vectors are rejected.
#' No coefficients are rounded.
#'
#' The interaction order describes each supplied coefficient vector, not the
#' table that will eventually be created by `emmeans`. Two separately supplied
#' nonzero, zero-sum order-0 contrasts form a correctly scaled order-1
#' interaction when crossed using `interaction = list(...)`. More generally,
#' crossing nonzero, zero-sum component orders `b` and `w` produces order
#' `b + w + 1`. An already formed two-factor product vector must instead be
#' supplied using interaction order 1.
#'
#' Zero-sum balance is assessed after dividing a vector by its largest absolute
#' coefficient. If `z = x / max(abs(x))`, its relative imbalance is
#' `abs(sum(z)) / sum(abs(z))`, and it is treated as zero-sum when that value is
#' no greater than `tol`.
#'
#' A single mode and interaction order apply to every vector in `contrasts`.
#' Lists containing contrasts of different orders should be split, rescaled in
#' separate calls, and recombined. Selectors, averaging vectors, and contrasts
#' that should retain their supplied scale should bypass this function. That
#' bypass is the R-workflow equivalent of Psy's original No Rescaling choice;
#' it does not require a function call.
#'
#' @return A named list with the same structure, names, and ordering as
#'   `contrasts`, containing the rescaled coefficient vectors.
#' @export
#'
#' @examples
#' raw <- list(
#'     "Groups 1 and 2 versus 3 and 4" = c(1, 1, -1, -1)
#' )
#'
#' rescale_contrasts(raw, mode = "mean_difference")
#'
#' product <- list(
#'     "Two-factor interaction" = c(1, -1, -1, 1)
#' )
#' rescale_contrasts(
#'     product,
#'     mode = "interaction",
#'     interaction_order = 1
#' )
rescale_contrasts <- function(
    contrasts,
    mode = c("mean_difference", "interaction"),
    interaction_order = 0L,
    zero_sum = c("error", "warn"),
    tol = sqrt(.Machine$double.eps)
) {
    mode <- match.arg(mode)
    zero_sum <- match.arg(zero_sum)

    .validate_rescale_tolerance(tol)
    .validate_rescale_interaction_order(interaction_order)
    .validate_rescale_contrast_list(contrasts)

    if (mode != "interaction" && interaction_order != 0) {
        stop(
            "`interaction_order` must be 0 unless `mode = \"interaction\"`.",
            call. = FALSE
        )
    }

    has_positive_coefficient <- vapply(
        contrasts,
        function(x) any(x > 0),
        logical(1)
    )
    if (!all(has_positive_coefficient)) {
        invalid <- names(contrasts)[!has_positive_coefficient]
        stop(
            sprintf(
                paste0(
                    "Cannot rescale contrast coefficient vectors without a ",
                    "positive coefficient: %s."
                ),
                .format_rescale_names(invalid)
            ),
            call. = FALSE
        )
    }

    relative_imbalance <- vapply(
        contrasts,
        .relative_contrast_imbalance,
        numeric(1)
    )
    non_zero_sum <- relative_imbalance > tol
    condition_message <- NULL

    if (any(non_zero_sum)) {
        invalid <- names(contrasts)[non_zero_sum]
        condition_message <- sprintf(
            paste0(
                "Contrast coefficient vectors must sum to zero within `tol`: ",
                "%s."
            ),
            .format_rescale_names(invalid)
        )

        if (zero_sum == "error") {
            stop(condition_message, call. = FALSE)
        }
    }

    target_positive_sum <- if (mode == "mean_difference") {
        1
    } else {
        2^interaction_order
    }

    rescaled <- lapply(
        contrasts,
        .rescale_contrast_vector,
        target_positive_sum = target_positive_sum
    )
    finite_output <- vapply(rescaled, function(x) all(is.finite(x)), logical(1))
    if (!all(finite_output)) {
        invalid <- names(contrasts)[!finite_output]
        stop(
            sprintf(
                paste0(
                    "Rescaling would produce non-finite coefficients for: ",
                    "%s. The coefficient magnitudes may be too disparate."
                ),
                .format_rescale_names(invalid)
            ),
            call. = FALSE
        )
    }

    if (!is.null(condition_message)) {
        warning(condition_message, call. = FALSE)
    }

    rescaled
}

.validate_rescale_tolerance <- function(tol) {
    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol < 0 || tol >= 1) {
        stop(
            "`tol` must be one finite number greater than or equal to 0 and less than 1.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

.validate_rescale_interaction_order <- function(interaction_order) {
    if (!is.numeric(interaction_order) || length(interaction_order) != 1L ||
        is.na(interaction_order) || !is.finite(interaction_order) ||
        interaction_order != floor(interaction_order) ||
        interaction_order < 0 || interaction_order > 9) {
        stop(
            "`interaction_order` must be one whole number from 0 to 9.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

.validate_rescale_contrast_list <- function(contrasts) {
    if (!is.vector(contrasts, mode = "list") ||
        length(contrasts) == 0L) {
        stop(
            "`contrasts` must be a non-empty ordinary named list.",
            call. = FALSE
        )
    }

    contrast_names <- names(contrasts)
    if (is.null(contrast_names) || anyNA(contrast_names) ||
        any(!nzchar(trimws(contrast_names))) ||
        anyDuplicated(contrast_names)) {
        stop(
            "`contrasts` must have non-empty, unique names.",
            call. = FALSE
        )
    }

    valid_vectors <- vapply(
        contrasts,
        function(x) is.numeric(x) && is.null(dim(x)) && length(x) > 0L,
        logical(1)
    )
    if (!all(valid_vectors)) {
        invalid <- contrast_names[!valid_vectors]
        stop(
            sprintf(
                "Each contrast must be a non-empty numeric vector: %s.",
                .format_rescale_names(invalid)
            ),
            call. = FALSE
        )
    }

    finite_vectors <- vapply(
        contrasts,
        function(x) all(is.finite(x)),
        logical(1)
    )
    if (!all(finite_vectors)) {
        invalid <- contrast_names[!finite_vectors]
        stop(
            sprintf(
                "Contrast coefficients must all be finite: %s.",
                .format_rescale_names(invalid)
            ),
            call. = FALSE
        )
    }

    vector_lengths <- lengths(contrasts)
    if (length(unique(vector_lengths)) != 1L) {
        stop(
            "All contrast coefficient vectors must have the same length.",
            call. = FALSE
        )
    }

    all_zero <- vapply(
        contrasts,
        function(x) all(x == 0),
        logical(1)
    )
    if (any(all_zero)) {
        invalid <- contrast_names[all_zero]
        stop(
            sprintf(
                "Contrast coefficient vectors cannot be all zero: %s.",
                .format_rescale_names(invalid)
            ),
            call. = FALSE
        )
    }

    invisible(NULL)
}

.relative_contrast_imbalance <- function(coefficients) {
    normalized <- coefficients / max(abs(coefficients))
    abs(sum(normalized)) / sum(abs(normalized))
}

.rescale_contrast_vector <- function(coefficients, target_positive_sum) {
    normalized <- coefficients / max(abs(coefficients))
    positive_sum <- sum(normalized[normalized > 0])
    normalized * target_positive_sum / positive_sum
}

.format_rescale_names <- function(names) {
    paste(sprintf("`%s`", names), collapse = ", ")
}
