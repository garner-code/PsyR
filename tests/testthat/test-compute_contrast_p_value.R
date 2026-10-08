test_that("individual p-values use the two-sided t distribution", {
    estimate <- c(lower = -2.5, null = 0.25, upper = 4)
    se_contrast <- c(0.5, 1, 2)
    null <- 0.25
    v_e <- 17

    expected <- 2 * stats::pt(
        abs((estimate - null) / se_contrast),
        df = v_e,
        lower.tail = FALSE
    )

    observed <- compute_contrast_p_value(
        estimate = estimate,
        se_contrast = se_contrast,
        method = "ind",
        v_e = v_e,
        null = null
    )

    expect_equal(observed, expected)
    expect_identical(names(observed), names(estimate))
    expect_equal(
        compute_contrast_p_value(estimate, se_contrast, "ind", v_e),
        compute_contrast_p_value(-estimate, se_contrast, "ind", v_e)
    )
})

test_that("Bonferroni p-values use the requested family size", {
    estimate <- c(-3, 0, 1.5, 8)
    se_contrast <- c(1, 2, 0.5, 1)
    v_e <- 24
    n_k <- 4
    t_statistic <- abs(estimate / se_contrast)
    expected <- pmin(
        1,
        n_k * 2 * stats::pt(t_statistic, df = v_e, lower.tail = FALSE)
    )

    observed <- compute_contrast_p_value(
        estimate = estimate,
        se_contrast = se_contrast,
        method = "bf",
        v_e = v_e,
        n_k = n_k
    )

    expect_equal(observed, expected)
    expect_equal(
        compute_contrast_p_value(estimate, se_contrast, "bf", v_e, n_k = 1),
        compute_contrast_p_value(estimate, se_contrast, "ind", v_e)
    )
})

test_that("post-hoc between and within p-values invert their F statistics", {
    estimate <- c(-2.5, 0.75, 4)
    se_contrast <- c(1, 0.5, 2)
    v_e <- 20
    v_b <- 3
    v_w <- 2
    statistic_squared <- ((estimate - 0.25) / se_contrast)^2

    expected_between <- stats::pf(
        statistic_squared / v_b,
        df1 = v_b,
        df2 = v_e,
        lower.tail = FALSE
    )
    expected_within <- stats::pf(
        statistic_squared * (v_e - v_w + 1) / (v_w * v_e),
        df1 = v_w,
        df2 = v_e - v_w + 1,
        lower.tail = FALSE
    )

    expect_equal(
        compute_contrast_p_value(
            estimate,
            se_contrast,
            "ph",
            v_e,
            family = "b",
            v_b = v_b,
            null = 0.25
        ),
        expected_between
    )
    expect_equal(
        compute_contrast_p_value(
            estimate,
            se_contrast,
            "ph",
            v_e,
            family = "w",
            v_w = v_w,
            null = 0.25
        ),
        expected_within
    )
})

test_that("post-hoc mixed p-values use the greatest-root upper tail", {
    estimate <- c(-2.5, 0.75, 4)
    se_contrast <- c(1, 0.5, 2)
    v_e <- 20
    v_b <- 3
    v_w <- 2
    statistic_squared <- (estimate / se_contrast)^2
    theta <- statistic_squared / (v_e + statistic_squared)
    expected <- vapply(
        theta,
        pgcr,
        numeric(1),
        v_w = v_w,
        v_b = v_b,
        v_e = v_e
    )

    observed <- compute_contrast_p_value(
        estimate,
        se_contrast,
        "ph",
        v_e,
        family = "bw",
        v_b = v_b,
        v_w = v_w
    )

    expect_equal(observed, expected)

    one_df_mixed <- compute_contrast_p_value(
        estimate,
        se_contrast,
        "ph",
        v_e,
        family = "bw",
        v_b = 1,
        v_w = 1
    )
    individual <- compute_contrast_p_value(
        estimate,
        se_contrast,
        "ind",
        v_e
    )
    expect_equal(one_df_mixed, individual, tolerance = 1e-12)
})

test_that("p-values reproduce the existing confidence-interval boundaries", {
    alpha <- 0.05
    se_contrast <- 1.7
    v_e <- 20

    individual_boundary <- cc_ind_t(v_e, alpha) * se_contrast
    expect_equal(
        compute_contrast_p_value(
            c(individual_boundary, -individual_boundary),
            se_contrast,
            "ind",
            v_e
        ),
        rep(alpha, 2),
        tolerance = 1e-12
    )

    bonferroni_boundary <- cc_bonf_t(v_e, n_k = 4, alpha) * se_contrast
    expect_equal(
        compute_contrast_p_value(
            c(bonferroni_boundary, -bonferroni_boundary),
            se_contrast,
            "bf",
            v_e,
            n_k = 4
        ),
        rep(alpha, 2),
        tolerance = 1e-12
    )

    between_boundary <- cc_ph_b(v_b = 3, v_e, alpha) * se_contrast
    expect_equal(
        compute_contrast_p_value(
            c(between_boundary, -between_boundary),
            se_contrast,
            "ph",
            v_e,
            family = "b",
            v_b = 3
        ),
        rep(alpha, 2),
        tolerance = 1e-12
    )

    within_boundary <- cc_ph_w(v_w = 2, v_e, alpha) * se_contrast
    expect_equal(
        compute_contrast_p_value(
            c(within_boundary, -within_boundary),
            se_contrast,
            "ph",
            v_e,
            family = "w",
            v_w = 2
        ),
        rep(alpha, 2),
        tolerance = 1e-12
    )

    mixed_boundary <- cc_ph_bw(v_w = 2, v_b = 5, v_e = 5, alpha) *
        se_contrast
    mixed_p_value <- compute_contrast_p_value(
        c(mixed_boundary, -mixed_boundary),
        se_contrast,
        "ph",
        v_e = 5,
        family = "bw",
        v_b = 5,
        v_w = 2
    )
    expect_lt(max(abs(mixed_p_value - alpha)), 5e-5)
})

test_that("core numeric inputs are vectorized with scalar recycling", {
    estimate <- c(first = 1, second = -2, third = 4)
    observed <- compute_contrast_p_value(
        estimate,
        se_contrast = 2,
        method = "ind",
        v_e = 30,
        null = c(0, -1, 1)
    )
    expected <- 2 * stats::pt(
        abs((estimate - c(0, -1, 1)) / 2),
        df = 30,
        lower.tail = FALSE
    )

    expect_type(observed, "double")
    expect_length(observed, length(estimate))
    expect_identical(names(observed), names(estimate))
    expect_equal(observed, expected)

    expect_length(
        compute_contrast_p_value(
            estimate = 2,
            se_contrast = c(1, 2, 4),
            method = "ind",
            v_e = 30
        ),
        3
    )
    expect_error(
        compute_contrast_p_value(
            estimate = 1:2,
            se_contrast = 1:3,
            method = "ind",
            v_e = 30
        )
    )
})

test_that("missing core inputs propagate row by row", {
    observed <- compute_contrast_p_value(
        estimate = c(NA_real_, 1, 1, 2),
        se_contrast = c(1, NA_real_, 1, 1),
        method = "ind",
        v_e = 12,
        null = c(0, 0, NA_real_, 0)
    )

    expect_true(all(is.na(observed[1:3])))
    expect_false(is.na(observed[4]))
})

test_that("returned p-values are bounded despite numerical round-off", {
    observed <- compute_contrast_p_value(
        estimate = c(0, 1e8),
        se_contrast = 1,
        method = "ph",
        v_e = 20,
        family = "bw",
        v_b = 3,
        v_w = 2
    )

    expect_equal(observed[1], 1)
    expect_true(all(observed >= 0 & observed <= 1))

    overflowed_statistic <- compute_contrast_p_value(
        estimate = .Machine$double.xmax,
        se_contrast = .Machine$double.xmin,
        method = "ph",
        v_e = 20,
        family = "bw",
        v_b = 3,
        v_w = 2
    )
    expect_true(is.finite(overflowed_statistic))
    expect_true(overflowed_statistic >= 0 && overflowed_statistic <= 1)

    expect_error(
        compute_contrast_p_value(
            estimate = 1,
            se_contrast = 1,
            method = "ph",
            v_e = 500,
            family = "bw",
            v_b = 2,
            v_w = 2
        ),
        "non-finite probability"
    )
})

test_that("core arguments are validated", {
    base_call <- function(...) {
        compute_contrast_p_value(
            estimate = 1,
            se_contrast = 1,
            method = "ind",
            v_e = 12,
            ...
        )
    }

    expect_error(compute_contrast_p_value("one", 1, "ind", 12))
    expect_error(compute_contrast_p_value(Inf, 1, "ind", 12))
    expect_error(compute_contrast_p_value(1, "one", "ind", 12))
    expect_error(compute_contrast_p_value(1, 0, "ind", 12))
    expect_error(compute_contrast_p_value(1, -1, "ind", 12))
    expect_error(compute_contrast_p_value(1, Inf, "ind", 12))
    expect_error(base_call(null = Inf))
    expect_error(compute_contrast_p_value(1, 1, "unknown", 12))
    expect_error(compute_contrast_p_value(1, 1, "i", 12))
    expect_error(compute_contrast_p_value(1, 1, "ind", 0))
    expect_error(compute_contrast_p_value(1, 1, "ind", Inf))
    expect_error(compute_contrast_p_value(1, 1, "ind", c(10, 12)))
})

test_that("method-specific arguments are validated", {
    expect_error(compute_contrast_p_value(1, 1, "bf", 12))
    expect_error(compute_contrast_p_value(1, 1, "bf", 12, n_k = 0))
    expect_error(compute_contrast_p_value(1, 1, "bf", 12, n_k = 1.5))
    expect_error(
        compute_contrast_p_value(1, 1, "bf", 12, n_k = 100000000.5)
    )
    expect_error(compute_contrast_p_value(1, 1, "bf", 12, n_k = c(2, 3)))

    expect_error(compute_contrast_p_value(1, 1, "ph", 12))
    expect_error(
        compute_contrast_p_value(1, 1, "ph", 12, family = "unknown")
    )
    expect_error(
        compute_contrast_p_value(1, 1, "ph", 12, family = "bet", v_b = 2)
    )
    expect_error(
        compute_contrast_p_value(1, 1, "ph", 12, family = "b")
    )
    expect_error(
        compute_contrast_p_value(1, 1, "ph", 12, family = "b", v_b = 0)
    )
    expect_error(
        compute_contrast_p_value(1, 1, "ph", 12, family = "w")
    )
    expect_error(
        compute_contrast_p_value(1, 1, "ph", 12, family = "w", v_w = 0)
    )
    expect_error(
        compute_contrast_p_value(1, 1, "ph", 2, family = "w", v_w = 3)
    )
    expect_error(
        compute_contrast_p_value(
            1,
            1,
            "ph",
            12,
            family = "bw",
            v_b = 2
        )
    )
    expect_error(
        compute_contrast_p_value(
            1,
            1,
            "ph",
            12,
            family = "bw",
            v_w = 2
        )
    )
})

test_that("SMR parameters are validated before simulation", {
    expect_error(compute_contrast_p_value(1, 1, "smr", 20))
    expect_error(
        compute_contrast_p_value(
            1,
            1,
            "smr",
            20,
            smr_params = list(p = 2)
        )
    )
    expect_error(
        compute_contrast_p_value(
            1,
            1,
            "smr",
            20,
            smr_params = list(p = 0, q = 2, n_sim = 100, seed = 1)
        )
    )
    expect_error(
        compute_contrast_p_value(
            1,
            1,
            "smr",
            20,
            smr_params = list(p = 2, q = 2.5, n_sim = 100, seed = 1)
        )
    )
    expect_error(
        compute_contrast_p_value(
            1,
            1,
            "smr",
            20,
            smr_params = list(p = 2, q = 2, n_sim = 0, seed = 1)
        )
    )
    expect_error(
        compute_contrast_p_value(
            1,
            1,
            "smr",
            20,
            smr_params = list(p = 2, q = 2, n_sim = 100, seed = c(1, 2))
        )
    )
})

test_that("seeded SMR p-values match the existing critical-value boundary", {
    critical_constant <- cc_smr(
        p = 2,
        q = 2,
        v_e = 20,
        alpha = 0.05,
        n_sim = 5000,
        seed = 123
    )
    se_contrast <- 1.5
    smr_params <- list(p = 2, q = 2, n_sim = 5000, seed = 123)

    observed <- compute_contrast_p_value(
        estimate = c(
            critical_constant * se_contrast,
            -critical_constant * se_contrast,
            2.5 * se_contrast
        ),
        se_contrast = se_contrast,
        method = "smr",
        v_e = 20,
        smr_params = smr_params
    )

    expect_equal(observed[1:2], rep(0.05, 2), tolerance = 1e-8)
    expect_equal(observed[3], 0.181347444212378, tolerance = 1e-6)

    distribution <- .fit_smr_reference_distribution(
        p = 2,
        q = 2,
        n = 20,
        n_sim = 5000,
        seed = 123
    )
    expect_equal(
        sqrt(.smr_reference_quantile(0.05, distribution)),
        critical_constant,
        tolerance = 1e-12
    )
})

test_that("seeded SMR p-values preserve the caller's RNG state", {
    set.seed(765)
    seed_before <- .Random.seed

    compute_contrast_p_value(
        estimate = 2,
        se_contrast = 1,
        method = "smr",
        v_e = 20,
        smr_params = list(p = 2, q = 2, n_sim = 100, seed = 123)
    )

    expect_identical(.Random.seed, seed_before)
})
