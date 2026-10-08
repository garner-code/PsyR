#' Compute a PsyR-compatible p-value for an ANOVA contrast
#'
#' Computes two-sided p-values that use the same reference distributions as
#' PsyR's confidence-interval procedures. The returned value is the smallest
#' family-wise alpha at which the corresponding PsyR confidence interval would
#' reject the supplied null value. It is therefore a contrast-specific p-value,
#' not an omnibus test of the complete contrast family.
#'
#' @param estimate A numeric vector of contrast estimates.
#' @param se_contrast A positive numeric vector containing the standard errors
#'   of the contrast estimates.
#' @param method A single method code: `"ind"` for an individual t procedure,
#'   `"bf"` for a Bonferroni-t procedure, `"ph"` for a post-hoc procedure, or
#'   `"smr"` for a Studentized Maximum Root procedure.
#' @param v_e A single positive number giving the error degrees of freedom.
#' @param family For `method = "ph"`, a single family code: `"b"` for a
#'   between-subjects family, `"w"` for a within-subjects family, or `"bw"`
#'   for a between-by-within family.
#' @param n_k For `method = "bf"`, the positive integer number of contrasts in
#'   the family.
#' @param v_b For post-hoc between-subjects calculations, the positive integer
#'   between-subjects degrees of freedom.
#' @param v_w For post-hoc within-subjects calculations, the positive integer
#'   within-subjects degrees of freedom.
#' @param smr_params For `method = "smr"`, a named list containing positive
#'   integers `p` and `q`. It may also contain `n_sim`, the number of Wishart
#'   simulations (an integer of at least 3; default 100000), and `seed`, an
#'   optional non-negative integer random seed.
#' @param null A numeric null value. A scalar is recycled across contrasts.
#'   Default is zero.
#'
#' @return A numeric vector of p-values with the same length as the recycled
#'   `estimate`, `se_contrast`, and `null` inputs. Names are retained from the
#'   first named, non-scalar input. Missing values in any of those inputs
#'   produce a missing p-value.
#'
#' @details
#' `method = "ph"` selects the reference distribution from `family`. Between-
#' subjects contrasts use Scheffe's F procedure, within-subjects contrasts use
#' the corresponding Hotelling procedure, and between-by-within contrasts use
#' Roy's Greatest Characteristic Root distribution through `pgcr()`.
#'
#' The SMR calculation uses the same simulated moment-based scaled-F
#' approximation as `smr_crit()`. Set `smr_params$seed` when results must be
#' reproducible across separate calls.
#'
#' @export
#'
#' @examples
#' compute_contrast_p_value(
#'     estimate = 2.5,
#'     se_contrast = 1.2,
#'     method = "ind",
#'     v_e = 24
#' )
#'
#' compute_contrast_p_value(
#'     estimate = c(2.5, -1.1),
#'     se_contrast = c(1.2, 0.8),
#'     method = "ph",
#'     family = "b",
#'     v_e = 24,
#'     v_b = 3
#' )
compute_contrast_p_value <- function(
    estimate,
    se_contrast,
    method,
    v_e,
    family = NULL,
    n_k = NULL,
    v_b = NULL,
    v_w = NULL,
    smr_params = NULL,
    null = 0
) {
    method <- .match_contrast_p_choice(
        method,
        choices = c("ind", "bf", "ph", "smr"),
        argument = "method"
    )
    .validate_contrast_p_scalar(v_e, "v_e", lower = 0)

    inputs <- .recycle_contrast_p_inputs(estimate, se_contrast, null)
    estimate <- inputs$estimate
    se_contrast <- inputs$se_contrast
    null <- inputs$null

    missing_input <- is.na(estimate) | is.na(se_contrast) | is.na(null)
    observed <- !missing_input
    if (any(!is.na(estimate) & !is.finite(estimate))) {
        stop("`estimate` must contain only finite values or `NA`.", call. = FALSE)
    }
    if (any(!is.na(null) & !is.finite(null))) {
        stop("`null` must contain only finite values or `NA`.", call. = FALSE)
    }
    available_se <- !is.na(se_contrast)
    if (any(!is.finite(se_contrast[available_se])) ||
        any(se_contrast[available_se] <= 0)) {
        stop(
            "`se_contrast` must contain only positive finite values or `NA`.",
            call. = FALSE
        )
    }

    statistic <- rep(NA_real_, length(estimate))
    statistic[observed] <- (
        (estimate[observed] - null[observed]) / se_contrast[observed]
    )^2

    if (method == "ind") {
        p_value <- stats::pf(statistic, 1, v_e, lower.tail = FALSE)
    } else if (method == "bf") {
        .validate_contrast_p_scalar(
            n_k,
            "n_k",
            lower = 0,
            integer = TRUE
        )
        p_value <- n_k * stats::pf(statistic, 1, v_e, lower.tail = FALSE)
    } else if (method == "ph") {
        family <- .match_contrast_p_choice(
            family,
            choices = c("b", "w", "bw"),
            argument = "family"
        )
        if (family == "b") {
            .validate_contrast_p_scalar(
                v_b,
                "v_b",
                lower = 0,
                integer = TRUE
            )
            p_value <- stats::pf(
                statistic / v_b,
                df1 = v_b,
                df2 = v_e,
                lower.tail = FALSE
            )
        } else if (family == "w") {
            .validate_contrast_p_scalar(
                v_w,
                "v_w",
                lower = 0,
                integer = TRUE
            )
            denominator_df <- v_e - v_w + 1
            if (denominator_df <= 0) {
                stop(
                    "`v_e - v_w + 1` must be positive.",
                    call. = FALSE
                )
            }
            f_statistic <- statistic * denominator_df / (v_w * v_e)
            p_value <- stats::pf(
                f_statistic,
                df1 = v_w,
                df2 = denominator_df,
                lower.tail = FALSE
            )
        } else {
            .validate_contrast_p_scalar(
                v_b,
                "v_b",
                lower = 0,
                integer = TRUE
            )
            .validate_contrast_p_scalar(
                v_w,
                "v_w",
                lower = 0,
                integer = TRUE
            )
            if (v_e - v_w + 1 <= 0) {
                stop(
                    "`v_e - v_w + 1` must be positive.",
                    call. = FALSE
                )
            }
            theta <- statistic / (v_e + statistic)
            theta[is.infinite(statistic)] <- 1
            p_value <- vapply(theta, function(value) {
                if (is.na(value)) {
                    return(NA_real_)
                }
                pgcr(value, v_w = v_w, v_b = v_b, v_e = v_e)
            }, numeric(1))
            if (any(observed & !is.finite(p_value))) {
                stop(
                    paste(
                        "Roy's GCR calculation produced a non-finite",
                        "probability for the supplied parameters."
                    ),
                    call. = FALSE
                )
            }
        }
    } else {
        smr_params <- .validate_smr_p_value_params(smr_params, v_e)
        distribution <- .fit_smr_reference_distribution(
            p = smr_params$p,
            q = smr_params$q,
            n = v_e,
            n_sim = smr_params$n_sim,
            seed = smr_params$seed
        )
        p_value <- .smr_reference_p_value(statistic, distribution)
    }

    p_value <- pmin(1, pmax(0, p_value))
    names(p_value) <- inputs$output_names
    p_value
}

.match_contrast_p_choice <- function(value, choices, argument) {
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
        stop(
            sprintf("`%s` must be one of %s.", argument, .format_p_choices(choices)),
            call. = FALSE
        )
    }
    if (!value %in% choices) {
        stop(
            sprintf("`%s` must be one of %s.", argument, .format_p_choices(choices)),
            call. = FALSE
        )
    }
    value
}

.format_p_choices <- function(choices) {
    paste(sprintf("'%s'", choices), collapse = ", ")
}

.validate_contrast_p_scalar <- function(
    value,
    argument,
    lower = -Inf,
    integer = FALSE
) {
    valid <- is.numeric(value) && length(value) == 1L &&
        !is.na(value) && is.finite(value) && value > lower
    if (valid && integer) {
        valid <- value == floor(value)
    }
    if (!valid) {
        requirement <- if (integer) {
            sprintf(
                "one integer greater than or equal to %s",
                format(floor(lower) + 1, scientific = FALSE)
            )
        } else if (lower >= 0) {
            "one positive finite number"
        } else {
            "one finite number"
        }
        stop(
            sprintf("`%s` must be %s.", argument, requirement),
            call. = FALSE
        )
    }
    invisible(value)
}

.recycle_contrast_p_inputs <- function(estimate, se_contrast, null) {
    values <- list(
        estimate = estimate,
        se_contrast = se_contrast,
        null = null
    )
    numeric_input <- vapply(values, is.numeric, logical(1))
    non_empty <- vapply(values, length, integer(1)) > 0L
    if (!all(numeric_input & non_empty)) {
        stop(
            "`estimate`, `se_contrast`, and `null` must be non-empty numeric vectors.",
            call. = FALSE
        )
    }

    lengths <- vapply(values, length, integer(1))
    output_length <- max(lengths)
    if (any(!lengths %in% c(1L, output_length))) {
        stop(
            paste(
                "`estimate`, `se_contrast`, and `null` must have a common",
                "length or length one."
            ),
            call. = FALSE
        )
    }

    output_names <- NULL
    for (value in values) {
        if (length(value) == output_length && !is.null(names(value))) {
            output_names <- names(value)
            break
        }
    }
    values <- lapply(values, rep_len, length.out = output_length)
    values$output_names <- output_names
    values
}

.validate_smr_p_value_params <- function(smr_params, v_e) {
    if (!is.list(smr_params)) {
        stop(
            "`smr_params` must be a named list containing `p` and `q`.",
            call. = FALSE
        )
    }
    if (is.null(names(smr_params)) ||
        !all(c("p", "q") %in% names(smr_params))) {
        stop(
            "`smr_params` must be a named list containing `p` and `q`.",
            call. = FALSE
        )
    }
    .validate_contrast_p_scalar(
        smr_params$p,
        "smr_params$p",
        lower = 0,
        integer = TRUE
    )
    .validate_contrast_p_scalar(
        smr_params$q,
        "smr_params$q",
        lower = 0,
        integer = TRUE
    )
    if (v_e <= 6) {
        stop("`v_e` must be greater than 6 for the SMR method.", call. = FALSE)
    }

    if (is.null(smr_params$n_sim)) {
        smr_params$n_sim <- 100000
    }
    .validate_contrast_p_scalar(
        smr_params$n_sim,
        "smr_params$n_sim",
        lower = 2,
        integer = TRUE
    )

    if (is.null(smr_params$seed)) {
        smr_params$seed <- NULL
    } else {
        .validate_contrast_p_scalar(
            smr_params$seed,
            "smr_params$seed",
            lower = -1,
            integer = TRUE
        )
        if (smr_params$seed > .Machine$integer.max) {
            stop(
                "`smr_params$seed` must not exceed the maximum integer value.",
                call. = FALSE
            )
        }
    }

    smr_params
}
