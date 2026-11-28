make_emmeans_basis_fixture <- function() {
    data("spacing", package = "PsyR", envir = environment())
    spacing$group <- factor(spacing$group)
    model <- afex::aov_ez(
        "subj", "yield", spacing, within = "spacing", between = "group",
        include_aov = TRUE
    )
    tables <- lapply(c("univariate", "multivariate"), function(basis) {
        cell_means <- emmeans::emmeans(
            model, c("group", "spacing"), model = basis
        )
        list(
            b = emmeans::contrast(
                emmeans::emmeans(model, "group", model = basis),
                list("Group 1 - Group 2" = c(1, -1, 0, 0))
            ),
            w = emmeans::contrast(
                emmeans::emmeans(model, "spacing", model = basis),
                list("Twenty - Forty" = c(1, -1, 0))
            ),
            bw = emmeans::contrast(
                cell_means,
                interaction = list(
                    list("Group 1 - Group 2" = c(1, -1, 0, 0)),
                    list("Twenty - Forty" = c(1, -1, 0))
                )
            )
        )
    })
    list(model = model, univariate = tables[[1]], multivariate = tables[[2]])
}

test_that("actual table bases determine validity independently of options", {
    fixture <- make_emmeans_basis_fixture()
    old <- afex::afex_options("emmeans_model")
    on.exit(afex::afex_options(emmeans_model = old))
    afex::afex_options(emmeans_model = "univariate")
    checked <- withVisible(check_emmeans_model(
        fixture$model, fixture$multivariate
    ))
    expect_false(checked$visible)
    expect_identical(checked$value$family, c("b", "w", "bw"))
    expect_identical(checked$value$detected_basis, rep("multivariate", 3))
    expect_identical(checked$value$multivariate_required, c(FALSE, TRUE, TRUE))
    expect_true(all(checked$value$valid))
    expect_identical(checked$value$table, paste0("contrast_tables[['", c("b", "w", "bw"), "']]"))
    for (family in c("w", "bw")) {
        expected <- psyci(fixture$model, fixture$multivariate[[family]], method = "ind")
        afex::afex_options(emmeans_model = "multivariate")
        actual <- psyci(fixture$model, fixture$multivariate[[family]], method = "ind")
        expect_equal(actual, expected, tolerance = 1e-12)
        expect_error(
            check_emmeans_model(fixture$model, fixture$univariate[[family]]),
            "contrast_tables\\[\\[1\\]\\].*basis: univariate"
        )
        afex::afex_options(emmeans_model = "univariate")
    }
    expect_identical(afex::afex_options("emmeans_model"), "univariate")
})

test_that("psyci enforces the basis requirement in every family check mode", {
    fixture <- make_emmeans_basis_fixture()
    for (mode in c("error", "warn", "none")) {
        for (family in c("w", "bw")) {
            expect_error(
                psyci(fixture$model, fixture$univariate[[family]],
                    method = "ind", family_check = mode),
                "requires multivariate"
            )
        }
    }
    for (mode in c("warn", "none")) {
        expect_error(
            suppressWarnings(psyci(
                fixture$model, fixture$univariate$w, method = "ind",
                family_list = "b", family_check = mode
            )),
            "requires multivariate"
        )
    }
    expect_error(
        check_emmeans_model(fixture$model, fixture$univariate$w, "b"),
        "requires multivariate"
    )
    expect_error(
        check_emmeans_model(fixture$model, fixture$univariate$b, "w"),
        "requires multivariate"
    )
    expect_silent(check_contrast_families(fixture$model, fixture$univariate$w))
})

test_that("unknown required bases and mixed-basis lists identify affected tables", {
    fixture <- make_emmeans_basis_fixture()
    unknown <- fixture$multivariate$w
    unknown@model.info$call <- NULL
    expect_error(
        check_emmeans_model(fixture$model, list(mystery = unknown), "w"),
        "mystery.*basis: unknown"
    )
    expect_error(
        psyci(fixture$model, list(mystery = unknown), method = "ind", family_list = "w",
            family_check = "none"),
        "mystery.*basis: unknown"
    )
    expect_error(
        check_emmeans_model(fixture$model, list(
            good = fixture$multivariate$w, bad = fixture$univariate$bw
        )),
        "bad.*basis: univariate"
    )
    checked <- check_emmeans_model(fixture$model, list(
        between = fixture$univariate$b, within = fixture$multivariate$w
    ))
    expect_true(all(checked$valid))
    expect_identical(checked$detected_basis, c("univariate", "multivariate"))
})

test_that("failed checks preserve inputs and global options", {
    fixture <- make_emmeans_basis_fixture()
    before_model <- serialize(fixture$model, NULL)
    before_tables <- unserialize(serialize(fixture$univariate, NULL))
    before_options <- afex::afex_options()
    expect_error(check_emmeans_model(fixture$model, fixture$univariate), "requires multivariate")
    expect_error(psyci(fixture$model, fixture$univariate$w, method = "ind"), "requires multivariate")
    expect_identical(serialize(fixture$model, NULL), before_model)
    # R may compile the retained df functions during evaluation; compare their
    # behavior and contents rather than the bytecode representation.
    expect_equal(fixture$univariate, before_tables)
    expect_identical(afex::afex_options(), before_options)
})

test_that("pure within designs are checked and between designs remain accepted", {
    old <- afex::afex_options("emmeans_model")
    on.exit(afex::afex_options(emmeans_model = old))
    data("spacing", package = "PsyR", envir = environment())
    model <- afex::aov_ez("subj", "yield", spacing, within = "spacing")
    table <- emmeans::contrast(
        emmeans::emmeans(model, "spacing", model = "multivariate"),
        list("Twenty - Forty" = c(1, -1, 0))
    )
    expect_true(check_emmeans_model(model, table)$valid)
    between <- make_between_fixture()
    checked <- check_emmeans_model(between$model, between$tables)
    expect_false(checked$multivariate_required)
    expect_true(checked$valid)
})

test_that("checker validates original objects and family arguments", {
    fixture <- make_emmeans_basis_fixture()
    for (bad in list(NULL, list(), 1, summary(fixture$multivariate$w))) {
        expect_error(check_emmeans_model(fixture$model, bad), "original emmGrid")
    }
    expect_error(check_emmeans_model(NULL, fixture$multivariate$w), "afex_aov")
    expect_error(check_emmeans_model(fixture$model, fixture$multivariate$w, "bad"), "family code")
    expect_error(check_emmeans_model(fixture$model, fixture$multivariate, c("w", "bw")), "length")
    checked <- check_emmeans_model(fixture$model, rep(list(fixture$multivariate$w), 2), "w")
    expect_true(all(checked$valid))
})
