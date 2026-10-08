make_tutorial_family_fixture <- function() {
    data("depression", package = "PsyR", envir = environment())
    depression$Group <- factor(
        depression$Group,
        levels = c("Ctrl", "Tmt1", "Tmt2")
    )
    depression$Time <- factor(
        depression$Time,
        levels = c("Pre", "Post", "FU")
    )
    old_emmeans_model <- afex::afex_options("emmeans_model")
    afex::afex_options(emmeans_model = "multivariate")
    model <- afex::aov_ez(
        "Subject", "Happiness", depression,
        within = "Time", between = "Group"
    )
    emmeans <- emmeans::emmeans(model, c("Group", "Time"))
    between <- list(
        Ts_v_Ctrl = c(-1, 0.5, 0.5),
        T1_v_T2 = c(0, 1, -1)
    )
    between_with_selectors <- c(
        between,
        list(`T1*` = c(0, 1, 0), `T2*` = c(0, 0, 1), `Ctrl*` = c(1, 0, 0))
    )
    within <- list(
        Post_vs_Pre = c(-1, 1, 0),
        FU_vs_Post = c(0, -1, 1)
    )
    within_with_selectors <- c(
        within,
        list(`Post*` = c(0, 1, 0), `FU*` = c(0, 0, 1))
    )
    tables <- list(
        emmeans::contrast(
            emmeans,
            interaction = list(between_with_selectors, within_with_selectors)
        ),
        emmeans::contrast(emmeans::emmeans(model, "Group"), between),
        emmeans::contrast(emmeans::emmeans(model, "Time"), within)
    )
    list(
        model = model,
        tables = tables,
        old_emmeans_model = old_emmeans_model
    )
}

test_that("family_group pools genuine contrasts after selector rows are dropped", {
    fixture <- make_tutorial_family_fixture()
    on.exit(
        afex::afex_options(emmeans_model = fixture$old_emmeans_model),
        add = TRUE
    )

    observed <- expect_no_warning(psyci(
        model = fixture$model,
        contrast_tables = fixture$tables,
        method = "bf",
        alpha = 0.15,
        family_group = rep("all_factorial", 3)
    ))

    expect_equal(
        unname(vapply(observed, nrow, integer(1))),
        c(14L, 2L, 2L)
    )
    expected_cc <- sqrt(stats::qf(1 - 0.15 / 18, 1, 12))
    expect_equal(unique(unlist(lapply(observed, `[[`, "cc"))), expected_cc)
    expect_equal(expected_cc, 3.15268131217, tolerance = 1e-10)

    dropped <- attr(observed[[1L]], "psyr_dropped_noncontrasts")
    expect_equal(nrow(dropped), 6L)
    expect_setequal(
        dropped$label,
        c(
            "T1* x Post*", "T2* x Post*", "Ctrl* x Post*",
            "T1* x FU*", "T2* x FU*", "Ctrl* x FU*"
        )
    )
    expect_identical(
        attr(observed[[1L]], "psyr_family_group"),
        "all_factorial"
    )
    expect_identical(attr(observed[[1L]], "psyr_family_size"), 18L)
    expect_true(any(grepl("dropped 6 non-contrast", attr(observed[[1L]], "mesg"))))
    expect_true(any(grepl(
        "family_group 'all_factorial' contains 18",
        attr(observed[[1L]], "mesg"),
        fixed = TRUE
    )))
})

test_that("the grouped tutorial intervals reproduce the published values", {
    fixture <- make_tutorial_family_fixture()
    on.exit(
        afex::afex_options(emmeans_model = fixture$old_emmeans_model),
        add = TRUE
    )
    observed <- psyci(
        fixture$model,
        fixture$tables,
        method = "bf",
        alpha = 0.15,
        family_group = rep("all", 3)
    )

    between <- observed[[2L]]
    within <- observed[[3L]]
    interaction <- observed[[1L]]
    expect_equal(round(between$lower, 3), c(-0.462, -5.997))
    expect_equal(round(between$upper, 3), c(6.462, 1.997))
    expect_equal(round(within$lower, 3), c(3.674, -5.021))
    expect_equal(round(within$upper, 3), c(8.326, -0.979))
    expect_equal(round(interaction$lower[1:5], 3),
        c(1.065, -1.698, 5.971, 1.971, -2.029))
    expect_equal(round(interaction$upper[1:5], 3),
        c(10.935, 9.698, 14.029, 10.029, 6.029))
})

test_that("separate correction warns for cross-table nonorthogonality", {
    fixture <- make_tutorial_family_fixture()
    on.exit(
        afex::afex_options(emmeans_model = fixture$old_emmeans_model),
        add = TRUE
    )
    observed <- NULL
    expect_warning(
        observed <- psyci(
            fixture$model,
            fixture$tables,
            method = "bf",
            alpha = 0.15
        ),
        regexp = "Non-orthogonality.*separately corrected",
        ignore.case = TRUE
    )
    expect_equal(
        unique(observed[[1L]]$cc),
        cc_bonf_t(v_e = 12, n_k = 14, alpha = 0.15)
    )
    expect_equal(
        unique(observed[[2L]]$cc),
        cc_bonf_t(v_e = 12, n_k = 2, alpha = 0.15)
    )
    expect_true(length(
        attr(observed[[1L]], "psyr_orthogonality")$not_orthogonal_to
    ) > 0L)
})

test_that("family_group validates method, length, and alpha agreement", {
    fixture <- make_tutorial_family_fixture()
    on.exit(
        afex::afex_options(emmeans_model = fixture$old_emmeans_model),
        add = TRUE
    )
    expect_error(
        psyci(fixture$model, fixture$tables, "ind", family_group = rep("a", 3)),
        regexp = "only.*method.*bf",
        ignore.case = TRUE
    )
    expect_error(
        psyci(fixture$model, fixture$tables, "bf", family_group = "a"),
        regexp = "one identifier per contrast table",
        fixed = TRUE
    )
    expect_error(
        psyci(
            fixture$model, fixture$tables, "bf",
            alpha = list(0.15, 0.05, 0.15),
            family_group = rep("a", 3)
        ),
        regexp = "same alpha",
        fixed = TRUE
    )
})

test_that("different family_group identifiers receive different pooled sizes", {
    fixture <- make_tutorial_family_fixture()
    on.exit(
        afex::afex_options(emmeans_model = fixture$old_emmeans_model),
        add = TRUE
    )
    observed <- suppressWarnings(psyci(
        fixture$model,
        fixture$tables,
        method = "bf",
        alpha = 0.15,
        family_group = c("primary", "primary", "secondary")
    ))
    expect_equal(
        unique(observed[[1L]]$cc),
        cc_bonf_t(v_e = 12, n_k = 16, alpha = 0.15)
    )
    expect_equal(
        unique(observed[[2L]]$cc),
        cc_bonf_t(v_e = 12, n_k = 16, alpha = 0.15)
    )
    expect_equal(
        unique(observed[[3L]]$cc),
        cc_bonf_t(v_e = 12, n_k = 2, alpha = 0.15)
    )
})

test_that("orthogonal tables do not trigger the separate-family warning", {
    fixture <- make_between_fixture()
    first <- emmeans::contrast(
        fixture$emmeans,
        list("Alert - Mild" = c(1, -1, 0, 0))
    )
    second <- emmeans::contrast(
        fixture$emmeans,
        list("Moderate - Extreme" = c(0, 0, 1, -1))
    )
    observed <- expect_no_warning(psyci(
        fixture$model,
        list(first, second),
        method = "bf"
    ))
    expect_identical(
        attr(observed[[1L]], "psyr_orthogonality")$orthogonal_to,
        "contrast_tables[[2]]"
    )
})

test_that("a table containing only selectors is rejected", {
    fixture <- make_between_fixture()
    selector <- emmeans::contrast(
        fixture$emmeans,
        list(Alert = c(1, 0, 0, 0))
    )
    expect_error(
        psyci(fixture$model, selector, "bf"),
        regexp = "no genuine contrasts",
        fixed = TRUE
    )
})
