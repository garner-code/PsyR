test_that("Algina-Keselman intervals reproduce the paper's worked example", {
    estimate <- 32.500 - 22.467
    covariance <- matrix(
        c(
            58.121, 35.517, 36.690,
            35.517, 49.016, 40.384,
            36.690, 40.384, 49.430
        ),
        nrow = 3,
        byrow = TRUE
    )
    se <- sqrt(c(1, -1, 0) %*% covariance %*% c(1, -1, 0) / 30)
    paper_scale <- sqrt(mean(diag(covariance)))

    observed <- algina_keselman_ci(
        estimate = estimate,
        se_contrast = se,
        df = 29,
        scale = paper_scale,
        alpha = 0.05,
        adjust = 3
    )

    expect_equal(observed$effect_size, 1.39, tolerance = 0.005)
    expect_equal(observed$lower, 0.82, tolerance = 0.005)
    expect_equal(observed$upper, 1.95, tolerance = 0.005)
})

test_that("noncentral-t limits have the requested tail probabilities", {
    observed <- algina_keselman_ci(
        estimate = c(first = -2.5, second = 4),
        se_contrast = c(1.2, 0.8),
        df = 18,
        scale = c(2, 3),
        alpha = 0.04,
        adjust = 2
    )
    tail_probability <- 0.04 / 4

    expect_equal(
        stats::pt(observed$statistic, observed$df, observed$ncp_lower),
        rep(1 - tail_probability, 2),
        tolerance = 1e-8
    )
    expect_equal(
        stats::pt(observed$statistic, observed$df, observed$ncp_upper),
        rep(tail_probability, 2),
        tolerance = 1e-8
    )
    expect_identical(rownames(observed), c("first", "second"))
    expect_true(all(observed$lower < observed$effect_size))
    expect_true(all(observed$upper > observed$effect_size))
})

test_that("Algina-Keselman inputs are validated", {
    expect_error(
        algina_keselman_ci(1, 0, 10, 2),
        "`se_contrast`"
    )
    expect_error(
        algina_keselman_ci(1, 1, -1, 2),
        "`df`"
    )
    expect_error(
        algina_keselman_ci(1, 1, 10, 0),
        "`scale`"
    )
    expect_error(
        algina_keselman_ci(1, 1, 10, 2, alpha = 1),
        "`alpha`"
    )
    expect_error(
        algina_keselman_ci(1, 1, 10, 2, adjust = 1.5),
        "`adjust`"
    )
    expect_error(
        algina_keselman_ci(1:2, 1:3, 10, 2),
        "common length"
    )
})

test_that("missing inputs produce missing interval rows", {
    observed <- algina_keselman_ci(
        estimate = c(1, NA_real_),
        se_contrast = 1,
        df = 10,
        scale = 2
    )

    expect_true(all(is.finite(unlist(observed[1, ], use.names = FALSE))))
    expect_true(all(is.na(observed[2, c(
        "effect_size", "statistic", "ncp", "ncp_lower", "ncp_upper",
        "lower", "upper"
    )])))
    expect_equal(observed$standard_error[2], 0.5)
})

test_that("psyci adds Algina-Keselman scaling without changing raw results", {
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    fixture <- make_mixed_fixture()
    within_table <- fixture$tables[2]

    existing <- psyci(
        fixture$model,
        within_table,
        method = "bf"
    )
    explicit_existing <- psyci(
        fixture$model,
        within_table,
        method = "bf",
        standardized_ci = "bird"
    )
    observed <- psyci(
        fixture$model,
        within_table,
        method = "bf",
        standardized_ci = "algina_keselman"
    )

    expect_identical(explicit_existing, existing)
    expect_identical(observed[[1L]], existing[[1L]])
    expect_identical(
        attr(observed, "standardized_ci_method"),
        "algina_keselman"
    )

    scaled <- attr(observed, "scaled_tables")[[1L]]
    raw <- observed[[1L]]
    expect_equal(scaled$estimate, raw$estimate / scaled$scale)
    expect_equal(scaled$SE, raw$SE / scaled$scale)
    expect_true(all(scaled$lower <= scaled$estimate))
    expect_true(all(scaled$upper >= scaled$estimate))
    expect_false(isTRUE(all.equal(
        scaled$upper - scaled$estimate,
        scaled$estimate - scaled$lower
    )))
    expect_equal(
        stats::pt(raw$t.ratio, raw$df, scaled$ncp_lower),
        rep(1 - 0.05 / (2 * nrow(raw)), nrow(raw)),
        tolerance = 1e-8
    )
})

test_that("psyci enforces the published scope of Algina-Keselman intervals", {
    between <- make_between_fixture()
    expect_error(
        psyci(
            between$model,
            between$tables,
            method = "ind",
            standardized_ci = "algina_keselman"
        ),
        "within-subject family"
    )
    expect_error(
        psyci(
            between$model,
            between$tables,
            method = "ph",
            standardized_ci = "algina_keselman"
        ),
        'method = "ind" or method = "bf"',
        fixed = TRUE
    )

    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    mixed <- make_mixed_fixture()
    nonpairwise <- emmeans::contrast(
        mixed$emmeans$within,
        list("Quadratic" = c(1, -2, 1)),
        adjust = "none"
    )
    expect_error(
        psyci(
            mixed$model,
            nonpairwise,
            method = "ind",
            standardized_ci = "algina_keselman"
        ),
        "exactly two"
    )
})
