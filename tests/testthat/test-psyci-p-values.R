test_that("psyci adds individual and Bonferroni p-values", {
    fixture <- make_between_fixture()
    before <- as.data.frame(summary(fixture$tables[[1L]]))
    v_e <- unique(before$df)

    individual <- psyci(
        model = fixture$model,
        contrast_tables = fixture$tables,
        method = "ind",
        family_list = fixture$families,
        between_factors = fixture$between_factors
    )[[1L]]
    expected_individual <- compute_contrast_p_value(
        estimate = before$estimate,
        se_contrast = before$SE,
        method = "ind",
        v_e = v_e
    )

    expect_equal(individual$p.value, before$p.value)
    expect_identical(individual$contrast, before$contrast)
    expect_equal(individual$psyr_p_value, expected_individual)
    expect_true(all(individual$psyr_p_value >= 0 &
        individual$psyr_p_value <= 1))
    expect_true(any(abs(individual$psyr_p_value - individual$p.value) > 1e-8))
    expect_identical(
        attr(individual, "psyr_p_value_method"),
        "independent"
    )
    expect_true(any(grepl(
        "PsyR p-value method: independent",
        attr(individual, "mesg"),
        fixed = TRUE
    )))

    bonferroni <- psyci(
        model = fixture$model,
        contrast_tables = fixture$tables,
        method = "bf",
        family_list = fixture$families,
        between_factors = fixture$between_factors
    )[[1L]]
    expected_bonferroni <- compute_contrast_p_value(
        estimate = before$estimate,
        se_contrast = before$SE,
        method = "bf",
        v_e = v_e,
        n_k = nrow(before)
    )

    expect_equal(bonferroni$p.value, before$p.value)
    expect_equal(bonferroni$psyr_p_value, expected_bonferroni)
    expect_identical(
        attr(bonferroni, "psyr_p_value_method"),
        "Bonferroni"
    )
    expect_true(any(grepl(
        "PsyR p-value method: Bonferroni",
        attr(bonferroni, "mesg"),
        fixed = TRUE
    )))
})

test_that("psyci rejects summary tables without coefficient provenance", {
    fixture <- make_between_fixture()
    contrast_table <- summary(fixture$tables[[1L]])
    estimable_df <- contrast_table$df[[1L]]
    contrast_table$estimate[[2L]] <- NA_real_
    contrast_table$SE[[2L]] <- NA_real_
    contrast_table$df[[2L]] <- NA_real_

    expect_error(
        psyci(
            model = fixture$model,
            contrast_tables = list(contrast_table),
            method = "ind",
            family_list = fixture$families,
            between_factors = fixture$between_factors,
            family_check = "none"
        ),
        regexp = "original emmGrid|genuine contrasts|orthogonality",
        ignore.case = TRUE
    )
})

test_that("psyci uses the pooled Bonferroni family size", {
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    fixture <- make_mixed_fixture()
    before <- lapply(
        fixture$tables,
        function(table) as.data.frame(summary(table))
    )
    pooled_n_k <- sum(vapply(before, nrow, integer(1)))
    v_e <- unique(unlist(lapply(before, function(table) table$df)))

    observed <- psyci(
        model = fixture$model,
        contrast_tables = fixture$tables,
        method = "bf",
        family_list = fixture$families,
        between_factors = fixture$between_factors,
        within_factors = fixture$within_factors,
        independent = FALSE
    )

    for (index in seq_along(observed)) {
        expected <- compute_contrast_p_value(
            estimate = before[[index]]$estimate,
            se_contrast = before[[index]]$SE,
            method = "bf",
            v_e = v_e,
            n_k = pooled_n_k
        )
        expect_equal(observed[[index]]$p.value, before[[index]]$p.value)
        expect_equal(observed[[index]]$psyr_p_value, expected)
    }
})

test_that("psyci maps post-hoc p-values to each contrast family", {
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    fixture <- make_mixed_fixture()
    before <- lapply(
        fixture$tables,
        function(table) as.data.frame(summary(table))
    )
    v_e <- unique(unlist(lapply(before, function(table) table$df)))
    v_b <- compute_df(fixture$model, fctrs = fixture$between_factors)
    v_w <- compute_df(fixture$model, fctrs = fixture$within_factors)

    observed <- psyci(
        model = fixture$model,
        contrast_tables = fixture$tables,
        method = "ph",
        family_list = fixture$families,
        between_factors = fixture$between_factors,
        within_factors = fixture$within_factors
    )
    expected_labels <- c(
        "Scheffe",
        "post-hoc within",
        "post-hoc between x within (Roy's GCR)"
    )

    for (index in seq_along(observed)) {
        family <- fixture$families[[index]]
        expected <- compute_contrast_p_value(
            estimate = before[[index]]$estimate,
            se_contrast = before[[index]]$SE,
            method = "ph",
            v_e = v_e,
            family = family,
            v_b = if (family %in% c("b", "bw")) v_b else NULL,
            v_w = if (family %in% c("w", "bw")) v_w else NULL
        )

        expect_equal(observed[[index]]$p.value, before[[index]]$p.value)
        expect_equal(observed[[index]]$psyr_p_value, expected)
        expect_identical(
            attr(observed[[index]], "psyr_p_value_method"),
            expected_labels[[index]]
        )
        expect_true(any(grepl(
            paste0("PsyR p-value method: ", expected_labels[[index]]),
            attr(observed[[index]], "mesg"),
            fixed = TRUE
        )))
    }

    expect_identical(observed[[1L]]$contrast, before[[1L]]$contrast)
    expect_identical(observed[[2L]]$contrast, before[[2L]]$contrast)
    expect_identical(
        observed[[3L]][c("group_custom", "spacing_custom")],
        before[[3L]][c("group_custom", "spacing_custom")]
    )
})

test_that("psyci reuses one seeded SMR distribution for CIs and p-values", {
    fixture <- make_between_fixture()
    before <- as.data.frame(summary(fixture$tables[[1L]]))
    v_e <- unique(before$df)
    alpha <- 0.05
    smr_params <- list(p = 2, q = 2, n_sim = 1000, seed = 123)
    expected_distribution <- .fit_smr_reference_distribution(
        p = smr_params$p,
        q = smr_params$q,
        n = v_e,
        n_sim = smr_params$n_sim,
        seed = smr_params$seed
    )
    expected_cc <- sqrt(
        .smr_reference_quantile(alpha, expected_distribution)
    )
    expected_p_value <- .smr_reference_p_value(
        (before$estimate / before$SE)^2,
        expected_distribution
    )

    original_sampler <- wishartlr_sample
    sampler_calls <- 0L
    testthat::local_mocked_bindings(
        wishartlr_sample = function(...) {
            sampler_calls <<- sampler_calls + 1L
            original_sampler(...)
        },
        .package = "PsyR"
    )
    set.seed(765)
    seed_before <- .Random.seed

    observed <- psyci(
        model = fixture$model,
        contrast_tables = list(
            fixture$tables[[1L]],
            fixture$tables[[1L]]
        ),
        method = "smr",
        family_list = list("b", "b"),
        between_factors = fixture$between_factors,
        alpha = alpha,
        smr_params = smr_params
    )

    expect_identical(sampler_calls, 1L)
    expect_identical(.Random.seed, seed_before)
    for (table in observed) {
        expect_equal(table$p.value, before$p.value)
        expect_equal(table$cc, rep(expected_cc, nrow(before)))
        expect_equal(table$psyr_p_value, expected_p_value)
        expect_identical(
            attr(table, "psyr_p_value_method"),
            "Studentized Maximum Root (SMR)"
        )
        expect_true(any(grepl(
            "PsyR p-value method: Studentized Maximum Root (SMR)",
            attr(table, "mesg"),
            fixed = TRUE
        )))
    }

    observed_with_top_level_seed <- psyci(
        model = fixture$model,
        contrast_tables = fixture$tables,
        method = "smr",
        family_list = fixture$families,
        between_factors = fixture$between_factors,
        alpha = alpha,
        smr_params = smr_params[c("p", "q", "n_sim")],
        seed = smr_params$seed
    )[[1L]]

    expect_identical(sampler_calls, 2L)
    expect_identical(.Random.seed, seed_before)
    expect_equal(observed_with_top_level_seed$cc, rep(expected_cc, nrow(before)))
    expect_equal(observed_with_top_level_seed$psyr_p_value, expected_p_value)
})
