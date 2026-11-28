make_polynomial_fixture <- function() {
    old <- afex::afex_options("emmeans_model")
    on.exit(afex::afex_options(emmeans_model = old))
    fixture <- make_mixed_fixture()
    fixture$emmeans <- fixture$emmeans$within
    fixture$tables <- list(emmeans::contrast(fixture$emmeans, list(
        trend = c(1, 0, -1), ordinary = c(0, 1, -1)
    ), adjust = "none"))
    fixture
}

test_that("standard polynomial detection is invariant to sign and scaling", {
    for (n in 3:8) {
        basis <- stats::contr.poly(n)
        for (degree in seq_len(n - 1L)) {
            for (scale in c(-1e100, -2, 1e-100, 3)) {
                expect_true(is_polynomial_contrast(scale * basis[, degree]))
            }
        }
    }
    expect_false(is_polynomial_contrast(c(1, -1)))
    expect_false(is_polynomial_contrast(c(1, 0, 0, -1)))
    expect_false(is_polynomial_contrast(c(0, 0, 0)))
    expect_false(is_polynomial_contrast(1))
    expect_true(is_polynomial_contrast(c(-1, 1e-10, 1)))
    expect_false(is_polynomial_contrast(c(-1, 1e-10, 1), tol = 1e-12))
    for (invalid in list(NULL, numeric(), "a", TRUE, c(NA, 1), c(Inf, 1),
                         matrix(1:3), 1i)) {
        expect_error(is_polynomial_contrast(invalid), "coefficients")
    }
    for (tol in list(NA, Inf, 0, -1, "a", numeric(), c(1, 2))) {
        expect_error(is_polynomial_contrast(c(-1, 0, 1), tol), "tol")
    }
})

test_that("factor profiles detect simple effects and polynomial interactions", {
    grid <- expand.grid(A = factor(1:3), B = factor(1:3))
    linear <- c(-1, 0, 1)
    quadratic <- c(1, -2, 1)
    ordinary <- c(1, -1, 0)
    coefficients <- cbind(
        main = rep(linear, 3),
        simple = c(linear, rep(0, 6)),
        product = as.vector(outer(linear, ordinary)),
        two_trends = as.vector(outer(quadratic, linear)),
        ordinary = as.vector(outer(ordinary, ordinary)),
        mixed_degrees = c(linear, quadratic, ordinary),
        zero = rep(0, 9)
    )
    extraction <- list(grid = grid, coefficients = coefficients)
    expected <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
    expect_identical(.polynomial_contrast_rows(extraction), expected)
    order <- c(9, 2, 4, 1, 8, 6, 3, 7, 5)
    extraction$grid <- grid[order, ]
    extraction$coefficients <- coefficients[order, ]
    expect_identical(.polynomial_contrast_rows(extraction), expected)
})

test_that("get_cis skips polynomial widths and rejects summaries", {
    fixture <- make_polynomial_fixture()
    table <- fixture$tables[[1L]]
    expected <- as.data.frame(summary(table))
    local_mocked_bindings(compute_contrast_ci = function(cc, se_contrast) {
        expect_length(se_contrast, 1L)
        cc * se_contrast
    })
    observed <- get_cis(table, 2)
    expect_identical(observed$is_polynomial, c(TRUE, FALSE))
    expect_equal(observed$estimate, expected$estimate)
    expect_true(is.na(observed$lower[1]))
    expect_true(is.na(observed$upper[1]))
    expect_equal(observed$lower[2], expected$estimate[2] - 2 * expected$SE[2])
    expect_error(get_cis(summary(table), 2), "original emmGrid")
    table <- update(table, infer = c(TRUE, TRUE))
    expect_false(any(c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL") %in%
        names(get_cis(table, 2))))
})

test_that("psyci retains polynomial tests and correction counts in mixed tables", {
    fixture <- make_polynomial_fixture()
    observed <- psyci(fixture$model, fixture$tables, method = "bf", legacy_summary = TRUE)
    raw <- observed[[1L]]
    scaled <- attr(observed, "scaled_tables")[[1L]]
    expect_identical(raw$is_polynomial, c(TRUE, FALSE))
    expect_equal(raw$p_value_adjusted, pmin(1, 2 * raw$p_value))
    expect_equal(raw$cc, rep(stats::qt(1 - .05 / 4, raw$df_error[1]), 2))
    expect_true(all(is.finite(raw$estimate)))
    expect_true(all(is.finite(raw$F_value)))
    expect_equal(sum(attr(observed, "legacy_anova_table")$row_type == "contrast"), 2L)
    for (table in list(raw, scaled)) {
        expect_true(is.na(table$lower[1]))
        expect_true(is.na(table$upper[1]))
        expect_true(is.finite(table$lower[2]))
        expect_true(is.finite(table$upper[2]))
        expect_length(grep("suppressed confidence", attr(table, "mesg")), 1L)
        output <- capture.output(print(table))
        expect_false(any(grepl("\\bNA\\b", output)))
    }
})

test_that("all-polynomial tables never call interval calculators", {
    fixture <- make_polynomial_fixture()
    table <- update(emmeans::contrast(fixture$emmeans, "poly"), infer = c(TRUE, TRUE))
    local_mocked_bindings(
        compute_contrast_ci = function(...) stop("width must not be calculated"),
        algina_keselman_ci = function(...) stop("inversion must not be calculated")
    )
    direct <- get_cis(table, 2)
    expect_true(all(direct$is_polynomial))
    expect_true(all(is.na(direct$lower)))
    for (method in c("bird", "algina_keselman")) {
        result <- psyci(fixture$model, table, method = "ind", standardized_ci = method)
        expect_true(all(result[[1L]]$is_polynomial))
        expect_true(all(is.na(result[[1L]]$lower)))
        scaled <- attr(result, "scaled_tables")[[1L]]
        expect_true(all(is.na(scaled$upper)))
        if (method == "algina_keselman") {
            expect_true(is.finite(scaled$estimate[1]))
            expect_true(is.na(scaled$estimate[2]))
        }
    }
})

test_that("Algina-Keselman skips polynomial inversion and retains ordinary intervals", {
    fixture <- make_polynomial_fixture()
    original <- algina_keselman_ci
    local_mocked_bindings(algina_keselman_ci = function(estimate, ...) {
        expect_length(estimate, 1L)
        original(estimate, ...)
    })
    result <- psyci(fixture$model, fixture$tables, method = "bf",
                    standardized_ci = "algina_keselman")
    raw <- result[[1L]]
    scaled <- attr(result, "scaled_tables")[[1L]]
    expected <- original(raw$estimate[2], raw$SE[2], raw$df_error[2],
                         scaled$scale[2], adjust = 2)
    expect_equal(scaled$lower[2], expected$lower)
    expect_equal(scaled$upper[2], expected$upper)
    expect_equal(scaled$estimate[1], raw$estimate[1] / scaled$scale[1])
    expect_true(is.na(scaled$ncp_lower[1]))
    expect_true(is.na(scaled$lower[1]))
})

test_that("emmeans interactions and conditional trends retain row alignment", {
    old <- afex::afex_options("emmeans_model")
    on.exit(afex::afex_options(emmeans_model = old))
    fixture <- make_mixed_fixture()
    interaction <- emmeans::contrast(fixture$emmeans$interaction,
        interaction = list(list(pair = c(1, -1, 0, 0)),
                           list(trend = c(-1, 0, 1), ordinary = c(1, -1, 0))))
    result <- psyci(fixture$model, interaction, method = "ind")[[1L]]
    expect_identical(result$is_polynomial, c(TRUE, FALSE))
    expect_identical(is.na(result$lower), c(TRUE, FALSE))
    conditional <- emmeans::contrast(
        emmeans::emmeans(fixture$model, ~ spacing | group), "poly")
    result <- psyci(fixture$model, conditional, method = "bf")[[1L]]
    expect_equal(nrow(result), 8L)
    expect_true(all(result$is_polynomial))
    expect_true(all(is.na(result$upper)))
})
