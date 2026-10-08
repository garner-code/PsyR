test_that("the Psy rescaling modes reproduce the legacy arithmetic", {
    raw <- list(
        "Balanced four-level contrast" = c(-1, -1, 1, 1)
    )

    mean_scaled <- rescale_contrasts(raw, mode = "mean_difference")
    expect_equal(
        mean_scaled[[1L]],
        c(-0.5, -0.5, 0.5, 0.5)
    )

    order_zero <- rescale_contrasts(
        raw,
        mode = "interaction",
        interaction_order = 0
    )
    expect_equal(order_zero, mean_scaled)

    order_one <- rescale_contrasts(
        raw,
        mode = "interaction",
        interaction_order = 1
    )
    expect_equal(order_one[[1L]], raw[[1L]])

    order_two <- rescale_contrasts(
        raw,
        mode = "interaction",
        interaction_order = 2
    )
    expect_equal(order_two[[1L]], c(-2, -2, 2, 2))
})

test_that("rescaling has the intended positive and negative coefficient mass", {
    raw <- c(-2, -1, 1, 2)

    for (interaction_order in 0:2) {
        scaled <- rescale_contrasts(
            list("Interaction" = raw),
            mode = "interaction",
            interaction_order = interaction_order
        )[[1L]]
        target <- 2^interaction_order

        expect_equal(sum(scaled[scaled > 0]), target)
        expect_equal(sum(scaled[scaled < 0]), -target)
    }
})

test_that("rescaling is invariant to positive input scale and preserves sign", {
    raw <- c(-1, -1, 2)
    scaled <- rescale_contrasts(list("Contrast" = raw))[[1L]]
    scaled_larger <- rescale_contrasts(list("Contrast" = 10 * raw))[[1L]]
    reversed <- rescale_contrasts(list("Contrast" = -raw))[[1L]]

    expect_equal(scaled, c(-0.5, -0.5, 1))
    expect_equal(scaled_larger, scaled)
    expect_equal(reversed, -scaled)
})

test_that("list and coefficient names and ordering are preserved", {
    raw <- list(
        "Second" = c(left = 1, middle = -1, right = 0),
        "First" = c(left = 0, middle = 1, right = -1)
    )

    scaled <- rescale_contrasts(raw)

    expect_identical(names(scaled), names(raw))
    expect_identical(names(scaled[[1L]]), names(raw[[1L]]))
    expect_identical(names(scaled[[2L]]), names(raw[[2L]]))
})

test_that("separately scaled components form the correctly scaled product", {
    between <- rescale_contrasts(
        list("Between main effect" = c(1, -1)),
        mode = "mean_difference"
    )[[1L]]
    within <- rescale_contrasts(
        list("Within main effect" = c(1, -1)),
        mode = "mean_difference"
    )[[1L]]

    product <- as.vector(outer(between, within))
    expect_equal(product, c(1, -1, -1, 1))
    expect_equal(sum(product[product > 0]), 2)

    direct_product <- rescale_contrasts(
        list("Between by within" = product),
        mode = "interaction",
        interaction_order = 1
    )[[1L]]
    expect_equal(direct_product, product)

    incorrectly_mean_scaled <- rescale_contrasts(
        list("Between by within" = product),
        mode = "mean_difference"
    )[[1L]]
    expect_equal(incorrectly_mean_scaled, product / 2)

    between_order_one <- rescale_contrasts(
        list("Between product" = c(1, -1)),
        mode = "interaction",
        interaction_order = 1
    )[[1L]]
    order_two_product <- as.vector(outer(between_order_one, within))
    expect_equal(sum(order_two_product[order_two_product > 0]), 4)
})

test_that("selectors are rejected and should bypass the function", {
    selector <- list("At the second level" = c(0, 1, 0))

    expect_error(
        rescale_contrasts(selector, mode = "mean_difference"),
        "At the second level"
    )
})

test_that("non-zero-sum handling is explicit", {
    invalid <- list("Unbalanced" = c(2, -1))

    expect_error(
        rescale_contrasts(invalid),
        "Unbalanced"
    )

    legacy_result <- NULL
    expect_warning(
        legacy_result <- rescale_contrasts(invalid, zero_sum = "warn"),
        "Unbalanced"
    )
    expect_equal(legacy_result[[1L]], c(1, -0.5))

    nearly_balanced <- list("Rounding residual" = c(1, -1 + 1e-12))
    expect_silent(
        rescale_contrasts(nearly_balanced, tol = 1e-10)
    )
})

test_that("invalid contrast inputs are rejected before rescaling", {
    expect_error(
        rescale_contrasts(c(1, -1)),
        "ordinary named list"
    )
    expect_error(
        rescale_contrasts(list()),
        "ordinary named list"
    )
    expect_error(
        rescale_contrasts(data.frame(A = c(1, -1), B = c(-1, 1))),
        "ordinary named list"
    )
    expect_error(
        rescale_contrasts(as.pairlist(list(A = c(1, -1)))),
        "ordinary named list"
    )
    expect_error(
        rescale_contrasts(
            structure(
                list(c(1, -1)),
                dim = 1L,
                dimnames = list("A")
            )
        ),
        "ordinary named list"
    )
    expect_error(
        rescale_contrasts(list(c(1, -1))),
        "non-empty, unique names"
    )
    expect_error(
        rescale_contrasts(
            stats::setNames(list(c(1, -1)), "")
        ),
        "non-empty, unique names"
    )
    expect_error(
        rescale_contrasts(list(A = c(1, -1), A = c(1, -1))),
        "non-empty, unique names"
    )
    expect_error(
        rescale_contrasts(list(A = character())),
        "non-empty numeric vector"
    )
    expect_error(
        rescale_contrasts(list(A = matrix(c(1, -1), nrow = 1))),
        "non-empty numeric vector"
    )
    expect_error(
        rescale_contrasts(list(A = c(1, -1), B = c(1, 0, -1))),
        "same length"
    )
    expect_error(
        rescale_contrasts(list(A = c(1, NA_real_))),
        "must all be finite"
    )
    expect_error(
        rescale_contrasts(list(A = c(1, Inf))),
        "must all be finite"
    )
    expect_error(
        rescale_contrasts(list(A = c(0, 0))),
        "cannot be all zero"
    )
    expect_error(
        rescale_contrasts(
            list(A = c(-1, -2)),
            mode = "interaction",
            zero_sum = "warn"
        ),
        "without a positive coefficient"
    )
    expect_error(
        rescale_contrasts(
            list(A = c(-1e308, 5e-324)),
            zero_sum = "error"
        ),
        "must sum to zero"
    )
    expect_error(
        rescale_contrasts(
            list(A = c(-1e308, 5e-324)),
            zero_sum = "warn"
        ),
        "non-finite coefficients"
    )
})

test_that("interaction order and tolerance are validated", {
    raw <- list(A = c(1, -1))

    expect_error(
        rescale_contrasts(raw, mode = "no_rescaling"),
        "one of"
    )
    expect_error(
        rescale_contrasts(raw, zero_sum = "ignore"),
        "one of"
    )

    order_nine <- rescale_contrasts(
        raw,
        mode = "interaction",
        interaction_order = 9
    )[[1L]]
    expect_equal(sum(order_nine[order_nine > 0]), 512)

    for (invalid_order in list(-1, 10, 0.5, c(0, 1), NA_real_, Inf)) {
        expect_error(
            rescale_contrasts(
                raw,
                mode = "interaction",
                interaction_order = invalid_order
            ),
            "whole number from 0 to 9"
        )
    }

    expect_error(
        rescale_contrasts(
            raw,
            mode = "mean_difference",
            interaction_order = 1
        ),
        "must be 0"
    )
    for (invalid_tol in list(-1, 1, c(0, 1), NA_real_, Inf, "small")) {
        expect_error(
            rescale_contrasts(raw, tol = invalid_tol),
            "greater than or equal to 0"
        )
    }
})

test_that("zero-sum diagnostics aggregate names and honor the tolerance", {
    invalid <- list(
        "First invalid" = c(2, -1),
        "Second invalid" = c(3, -1)
    )
    warned_result <- NULL
    expect_warning(
        warned_result <- rescale_contrasts(invalid, zero_sum = "warn"),
        "First invalid.*Second invalid"
    )
    expect_length(warned_result, 2L)

    boundary_vector <- c(1, -0.8)
    boundary <- abs(sum(boundary_vector)) / sum(abs(boundary_vector))
    expect_silent(
        rescale_contrasts(
            list("At boundary" = boundary_vector),
            tol = boundary
        )
    )
    expect_error(
        rescale_contrasts(
            list("Below boundary" = boundary_vector),
            tol = boundary / 2
        ),
        "Below boundary"
    )
})

test_that("rescaled lists work directly in emmeans interactions", {
    data <- expand.grid(
        between = factor(c("B1", "B2")),
        within = factor(c("W1", "W2")),
        replicate = 1:6
    )
    data$response <- with(
        data,
        10 +
            1 * (between == "B2") +
            2 * (within == "W2") +
            3 * (between == "B2" & within == "W2") +
            c(-2, -1, 0, 1, 2, 0)[replicate]
    )

    model <- stats::lm(response ~ between * within, data = data)
    between_means <- suppressMessages(
        emmeans::emmeans(model, "between")
    )
    cell_means <- emmeans::emmeans(model, c("between", "within"))

    raw_between <- list("Between" = c(2, -2))
    raw_within <- list("Within" = c(3, -3))
    scaled_between <- rescale_contrasts(raw_between)
    scaled_within <- rescale_contrasts(raw_within)

    raw_between_table <- emmeans::contrast(
        between_means,
        method = raw_between,
        adjust = "none"
    )
    scaled_between_table <- emmeans::contrast(
        between_means,
        method = scaled_between,
        adjust = "none"
    )
    raw_between_summary <- as.data.frame(summary(raw_between_table))
    scaled_between_summary <- as.data.frame(summary(scaled_between_table))

    expect_equal(
        raw_between_summary$estimate,
        2 * scaled_between_summary$estimate
    )
    expect_equal(raw_between_summary$SE, 2 * scaled_between_summary$SE)
    expect_equal(
        raw_between_summary$t.ratio,
        scaled_between_summary$t.ratio
    )
    expect_equal(
        raw_between_summary$p.value,
        scaled_between_summary$p.value
    )

    raw_table <- emmeans::contrast(
        cell_means,
        interaction = list(raw_between, raw_within),
        adjust = "none"
    )
    scaled_table <- emmeans::contrast(
        cell_means,
        interaction = list(scaled_between, scaled_within),
        adjust = "none"
    )

    raw_summary <- as.data.frame(summary(raw_table))
    scaled_summary <- as.data.frame(summary(scaled_table))

    expect_equal(raw_summary$estimate, 6 * scaled_summary$estimate)
    expect_equal(raw_summary$SE, 6 * scaled_summary$SE)
    expect_equal(raw_summary$t.ratio, scaled_summary$t.ratio)
    expect_equal(raw_summary$p.value, scaled_summary$p.value)

    raw_interval <- as.data.frame(stats::confint(raw_table, adjust = "none"))
    scaled_interval <- as.data.frame(
        stats::confint(scaled_table, adjust = "none")
    )
    expect_equal(raw_interval$lower.CL, 6 * scaled_interval$lower.CL)
    expect_equal(raw_interval$upper.CL, 6 * scaled_interval$upper.CL)

    coefficients <- stats::coef(scaled_table)
    coefficient_column <- grep("^c\\.", names(coefficients), value = TRUE)
    expect_length(coefficient_column, 1L)
    expect_equal(
        coefficients[[coefficient_column]],
        c(1, -1, -1, 1)
    )
})
