#' Algina-Keselman confidence intervals for standardized effects
#'
#' Computes the approximate confidence interval proposed by Algina and
#' Keselman (2003) for a standardized contrast of correlated means. The method
#' first inverts the noncentral t distribution for the observed t statistic,
#' then multiplies the noncentrality-parameter limits by the contrast standard
#' error divided by the sample standard deviation used to standardize the
#' effect.
#'
#' @param estimate A numeric vector of estimated contrasts in
#'   dependent-variable units.
#' @param se_contrast A numeric vector of positive contrast standard errors.
#' @param df A numeric vector of positive error degrees of freedom.
#' @param scale A numeric vector of positive sample standard deviations used to
#'   standardize the contrasts. For the pairwise method in Equation 9 of Algina
#'   and Keselman (2003), this is the square root of the mean of the two sample
#'   variances involved in the contrast.
#' @param alpha A single type-I error rate. Default is `0.05`.
#' @param adjust A single positive integer giving the number of comparisons in
#'   a Bonferroni family. Use `1` for an individual interval. Default is `1`.
#'
#' @return A data frame containing the standardized effect estimate and
#'   standard error, the observed t statistic, the estimated noncentrality
#'   parameter, its lower and upper confidence limits, and the corresponding
#'   standardized-effect confidence limits. A missing value in any row of the
#'   four numeric inputs produces missing calculated values for that row.
#'
#' @details
#' The confidence limits are generally asymmetric. `adjust` implements the
#' Bonferroni construction used in the paper by using tail probability
#' `alpha / (2 * adjust)`. This method was developed for contrasts of
#' correlated (repeated-measures) means; it should not be used for independent
#' groups.
#'
#' Algina, J., & Keselman, H. J. (2003). Approximate confidence intervals for
#' effect sizes. *Educational and Psychological Measurement, 63*(4), 537-553.
#' \doi{10.1177/0013164403256358}
#'
#' @export
#'
#' @examples
#' algina_keselman_ci(
#'     estimate = 32.500 - 22.467,
#'     se_contrast = sqrt((58.121 + 49.016 - 2 * 35.517) / 30),
#'     df = 29,
#'     scale = sqrt(mean(c(58.121, 49.016))),
#'     adjust = 3
#' )
algina_keselman_ci <- function(
    estimate,
    se_contrast,
    df,
    scale,
    alpha = 0.05,
    adjust = 1
) {
    .validate_ak_probability(alpha, "alpha")
    .validate_ak_adjust(adjust)

    inputs <- .recycle_ak_inputs(estimate, se_contrast, df, scale)
    estimate <- inputs$estimate
    se_contrast <- inputs$se_contrast
    df <- inputs$df
    scale <- inputs$scale

    .validate_ak_vector(estimate, "estimate", positive = FALSE)
    .validate_ak_vector(se_contrast, "se_contrast", positive = TRUE)
    .validate_ak_vector(df, "df", positive = TRUE)
    .validate_ak_vector(scale, "scale", positive = TRUE)

    complete <- stats::complete.cases(estimate, se_contrast, df, scale)
    statistic <- estimate / se_contrast
    tail_probability <- alpha / (2 * adjust)
    ncp_lower <- rep(NA_real_, length(statistic))
    ncp_upper <- rep(NA_real_, length(statistic))
    ncp_lower[complete] <- mapply(
        .invert_noncentral_t,
        statistic[complete],
        df[complete],
        MoreArgs = list(probability = 1 - tail_probability)
    )
    ncp_upper[complete] <- mapply(
        .invert_noncentral_t,
        statistic[complete],
        df[complete],
        MoreArgs = list(probability = tail_probability)
    )
    multiplier <- se_contrast / scale

    result <- data.frame(
        effect_size = estimate / scale,
        standard_error = multiplier,
        statistic = statistic,
        df = df,
        ncp = statistic,
        ncp_lower = unname(ncp_lower),
        ncp_upper = unname(ncp_upper),
        lower = unname(ncp_lower) * multiplier,
        upper = unname(ncp_upper) * multiplier,
        scale = scale,
        row.names = NULL
    )
    if (!is.null(inputs$output_names)) {
        rownames(result) <- inputs$output_names
    }
    result
}

.invert_noncentral_t <- function(statistic, df, probability) {
    objective <- function(ncp) {
        suppressWarnings(
            stats::pt(statistic, df = df, ncp = ncp) - probability
        )
    }

    lower <- -1
    upper <- 1
    lower_value <- objective(lower)
    upper_value <- objective(upper)
    limit <- 2^20

    while (lower_value < 0 && abs(lower) < limit) {
        lower <- lower * 2
        lower_value <- objective(lower)
    }
    while (upper_value > 0 && upper < limit) {
        upper <- upper * 2
        upper_value <- objective(upper)
    }
    if (!is.finite(lower_value) || !is.finite(upper_value) ||
        lower_value < 0 || upper_value > 0) {
        stop(
            "Could not bracket a noncentral-t confidence limit for the supplied values.",
            call. = FALSE
        )
    }

    stats::uniroot(
        objective,
        interval = c(lower, upper),
        f.lower = lower_value,
        f.upper = upper_value,
        tol = .Machine$double.eps^0.75
    )$root
}

.recycle_ak_inputs <- function(estimate, se_contrast, df, scale) {
    values <- list(
        estimate = estimate,
        se_contrast = se_contrast,
        df = df,
        scale = scale
    )
    if (!all(vapply(values, is.numeric, logical(1))) ||
        any(lengths(values) == 0L)) {
        stop(
            "`estimate`, `se_contrast`, `df`, and `scale` must be non-empty numeric vectors.",
            call. = FALSE
        )
    }
    output_length <- max(lengths(values))
    if (any(!lengths(values) %in% c(1L, output_length))) {
        stop(
            paste(
                "`estimate`, `se_contrast`, `df`, and `scale` must have a",
                "common length or length one."
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

.validate_ak_vector <- function(value, argument, positive) {
    valid <- is.na(value) | is.finite(value)
    if (positive) {
        valid <- valid & (is.na(value) | value > 0)
    }
    if (!all(valid)) {
        qualifier <- if (positive) {
            "positive finite values or `NA`"
        } else {
            "finite values or `NA`"
        }
        stop(
            sprintf("`%s` must contain only %s.", argument, qualifier),
            call. = FALSE
        )
    }
    invisible(value)
}

.validate_ak_probability <- function(value, argument) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value <= 0 || value >= 1) {
        stop(
            sprintf("`%s` must be one finite number strictly between 0 and 1.", argument),
            call. = FALSE
        )
    }
    invisible(value)
}

.validate_ak_adjust <- function(adjust) {
    if (!is.numeric(adjust) || length(adjust) != 1L || is.na(adjust) ||
        !is.finite(adjust) || adjust < 1 || adjust != floor(adjust)) {
        stop("`adjust` must be one positive whole number.", call. = FALSE)
    }
    invisible(adjust)
}
