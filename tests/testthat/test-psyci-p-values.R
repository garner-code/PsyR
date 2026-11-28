test_that("psyci separates ordinary and adjusted p-values", {
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

    expect_false(any(c("t.ratio", "p.value", "df") %in% names(individual)))
    expect_equal(individual$df1, rep(1L, nrow(before)))
    expect_equal(individual$df2, before$df)
    expect_equal(individual$F_value, (before$estimate / before$SE)^2)
    expect_identical(individual$contrast, before$contrast)
    expect_equal(individual$p_value, expected_individual)
    expect_equal(individual$p_adj, expected_individual)
    expect_identical(individual$method, rep("none", nrow(before)))
    expect_true(all(individual$p_adj >= 0 &
        individual$p_adj <= 1))
    expect_true(any(grepl(
        "PsyR p-value adjustment: none",
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

    expect_equal(bonferroni$F_value, individual$F_value)
    expect_equal(bonferroni$p_value, individual$p_value)
    expect_equal(bonferroni$p_adj, expected_bonferroni)
    expect_equal(
        bonferroni$p_adj,
        pmin(1, nrow(before) * bonferroni$p_value)
    )
    comparable <- stats::complete.cases(
        bonferroni$p_adj,
        bonferroni$lower,
        bonferroni$upper
    )
    expect_identical(
        bonferroni$p_adj[comparable] <= 0.05,
        bonferroni$lower[comparable] > 0 | bonferroni$upper[comparable] < 0
    )
    expect_identical(
        bonferroni$method,
        rep("Bonferroni", nrow(before))
    )
    expect_false(any(grepl(
        "PsyR simultaneous-inference procedure: Bonferroni",
        attr(bonferroni, "mesg"),
        fixed = TRUE
    )))
    expect_false(any(grepl(
        "dimension of: NA",
        attr(bonferroni, "mesg"),
        fixed = TRUE
    )))
    expect_false(any(grepl(
        "df between of:",
        attr(bonferroni, "mesg"),
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "ordinary F(1, df2) contrast statistic",
        attr(bonferroni, "mesg"),
        fixed = TRUE
    )))
    expect_false(any(grepl(
        paste("PsyR used df error of:", v_e),
        attr(bonferroni, "mesg"),
        fixed = TRUE
    )))
    expect_false("psyr_p_value" %in% names(bonferroni))
    expect_false(any(grepl(
        "adjusted `F_value`|Bonferroni-adjusted F",
        attr(bonferroni, "mesg"),
        ignore.case = TRUE
    )))
})

test_that("Bonferroni probabilities handle one contrast and the unit cap", {
    statistic <- c(0, 10)
    raw <- .compute_unadjusted_contrast_p_value(statistic, df_error = 20)
    one_contrast <- .compute_contrast_p_value_from_f(
        statistic = statistic,
        method = "bf",
        v_e = 20,
        n_k = 1
    )
    four_contrasts <- .compute_contrast_p_value_from_f(
        statistic = statistic,
        method = "bf",
        v_e = 20,
        n_k = 4
    )

    expect_equal(one_contrast, raw)
    expect_equal(four_contrasts, pmin(1, 4 * raw))
    expect_identical(four_contrasts[[1L]], 1)
    expect_lt(four_contrasts[[2L]], 1)
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
        regexp = "original emmGrid|interpretable contrasts|orthogonality",
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
        expect_equal(
            observed[[index]]$F_value,
            (before[[index]]$estimate / before[[index]]$SE)^2
        )
        expect_equal(
            observed[[index]]$p_value,
            stats::pf(
                observed[[index]]$F_value,
                df1 = 1,
                df2 = observed[[index]]$df2,
                lower.tail = FALSE
            )
        )
        expect_equal(observed[[index]]$p_adj, expected)
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
    expected_labels <- c("Scheffe", "Hotelling T-squared", "Roy GCR")

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

        expect_equal(
            observed[[index]]$F_value,
            (before[[index]]$estimate / before[[index]]$SE)^2
        )
        expect_equal(
            observed[[index]]$p_value,
            stats::pf(
                observed[[index]]$F_value,
                df1 = 1,
                df2 = observed[[index]]$df2,
                lower.tail = FALSE
            )
        )
        expect_equal(observed[[index]]$p_adj, expected)
        if (family == "bw") {
            expect_false("theta" %in% names(observed[[index]]))
            expect_equal(
                attr(observed[[index]], "psyr_probability_parameters")$theta,
                observed[[index]]$F_value /
                    (observed[[index]]$df2 + observed[[index]]$F_value)
            )
        } else {
            expect_false("theta" %in% names(observed[[index]]))
        }
        expect_identical(
            observed[[index]]$method,
            rep(expected_labels[[index]], nrow(observed[[index]]))
        )
        expect_false(any(grepl(
            paste0(
                "PsyR simultaneous-inference procedure: ",
                expected_labels[[index]]
            ),
            attr(observed[[index]], "mesg"),
            fixed = TRUE
        )))
    }

    expect_null(attr(
        observed[[1L]],
        "psyr_probability_parameters",
        exact = TRUE
    ))
    expect_null(attr(
        observed[[2L]],
        "psyr_probability_parameters",
        exact = TRUE
    ))
    expect_equal(
        attr(observed[[3L]], "psyr_probability_parameters"),
        list(
            s = min(v_b, v_w),
            m = (abs(v_b - v_w) - 1) / 2,
            n = (v_e - v_w - 1) / 2,
            theta = observed[[3L]]$F_value /
                (observed[[3L]]$df2 + observed[[3L]]$F_value)
        )
    )
    expect_true(any(grepl(
        paste0(
            "Psy probability-calculator parameters for Roy GCR: s = ",
            min(v_b, v_w)
        ),
        attr(observed[[3L]], "mesg"),
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "`theta` attribute entry contains the values for Psy's GCR field",
        attr(observed[[3L]], "mesg"),
        fixed = TRUE
    )))
    expect_equal(
        attr(attr(observed, "scaled_tables")[[3L]], "psyr_probability_parameters")$theta,
        attr(observed[[3L]], "psyr_probability_parameters")$theta
    )

    between_messages <- attr(observed[[1L]], "mesg")
    within_messages <- attr(observed[[2L]], "mesg")
    interaction_messages <- attr(observed[[3L]], "mesg")
    expect_true(any(grepl(
        paste("post-hoc between dimension of:", v_b),
        between_messages,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        paste("post-hoc within dimension of:", v_w),
        within_messages,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        paste("post-hoc between dimension of:", v_b),
        interaction_messages,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        paste("post-hoc within dimension of:", v_w),
        interaction_messages,
        fixed = TRUE
    )))

    expect_identical(observed[[1L]]$contrast, before[[1L]]$contrast)
    expect_identical(observed[[2L]]$contrast, before[[2L]]$contrast)
    expect_identical(
        as.data.frame(observed[[3L]][c("group_custom", "spacing_custom")]),
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
        expect_equal(table$F_value, (before$estimate / before$SE)^2)
        expect_equal(sqrt(table$F_crit), rep(expected_cc, nrow(before)))
        expect_equal(
            table$p_value,
            stats::pf(table$F_value, 1, table$df2, lower.tail = FALSE)
        )
        expect_equal(table$p_adj, expected_p_value)
        expect_identical(table$method, rep("SMR", nrow(table)))
        probability_parameters <- attr(
            table,
            "psyr_probability_parameters"
        )
        expect_named(probability_parameters, c("k", "p", "q", "df"))
        expect_equal(probability_parameters$k, sqrt(table$F_value))
        expect_identical(probability_parameters$p, smr_params$p)
        expect_identical(probability_parameters$q, smr_params$q)
        expect_equal(probability_parameters$df, v_e)
        expect_false(any(grepl(
            "PsyR simultaneous-inference procedure: SMR",
            attr(table, "mesg"),
            fixed = TRUE
        )))
        expect_true(any(grepl(
            "row-specific k = sqrt(F_value)",
            attr(table, "mesg"),
            fixed = TRUE
        )))
    }
    scaled_tables <- attr(observed, "scaled_tables")
    output <- paste(capture.output(print(observed)), collapse = "\n")
    expect_match(output, "parameters for SMR:", fixed = TRUE)
    expect_false(grepl("ordinary F(1, df2)", output, fixed = TRUE))
    for (index in seq_along(observed)) {
        expect_equal(
            attr(
                scaled_tables[[index]],
                "psyr_probability_parameters"
            ),
            attr(observed[[index]], "psyr_probability_parameters")
        )
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
    expect_equal(sqrt(observed_with_top_level_seed$F_crit), rep(expected_cc, nrow(before)))
    expect_equal(
        observed_with_top_level_seed$p_adj,
        expected_p_value
    )
})

test_that("ordinary p-values reproduce the Psy probability calculator", {
    observed <- .compute_unadjusted_contrast_p_value(
        F_value = 7.910922,
        df_error = 87
    )

    expect_equal(round(observed, 9), 0.006072727)
})

test_that("Roy GCR theta reproduces the Psy calculator input", {
    observed <- .compute_gcr_theta(F_value = 7.910922, df_error = 87)

    expect_equal(observed, 7.910922 / (87 + 7.910922))
    expect_equal(round(observed, 6), 0.083351)
    expect_identical(.compute_gcr_theta(Inf, 87), 1)
})

test_that("printed tables show only adjusted probabilities", {
    fixture <- make_between_fixture()
    table <- psyci(fixture$model, fixture$tables, method = "bf")[[1L]]
    output <- paste(capture.output(print(table, verbose = TRUE)), collapse = "\n")

    expect_false(grepl("p_value", output, fixed = TRUE))
    expect_match(output, "p_adj", fixed = TRUE)
    expect_match(output, "method", fixed = TRUE)
    expect_false(grepl("psyr_p_value", output, fixed = TRUE))
    expect_match(
        output,
        "ordinary F(1, df2) contrast statistic",
        fixed = TRUE
    )
})
