make_legacy_anxiety_fixture <- function() {
    data("social_anxiety", package = "PsyR", envir = environment())
    social_anxiety$Session <- factor(
        social_anxiety$Session,
        levels = c("Pre", "Post", "FU1", "FU2")
    )
    social_anxiety$Group <- factor(
        social_anxiety$Group,
        levels = c(1, 2, 3),
        labels = c("New", "Standard", "MinContact")
    )
    model <- afex::aov_ez(
        id = "Subject",
        dv = "Score",
        data = social_anxiety,
        between = "Group",
        within = "Session"
    )
    between <- rescale_contrasts(list(
        "Ts - C" = c(1, 1, -2),
        "NT - ST" = c(1, -1, 0)
    ))
    within <- rescale_contrasts(list(
        "Pre - Rest" = c(3, -1, -1, -1),
        "Post - FUs" = c(0, 2, -1, -1),
        "FU1 - FU2" = c(0, 0, 1, -1)
    ))
    list(
        model = model,
        tables = list(
            emmeans::contrast(emmeans::emmeans(model, "Session"), within),
            emmeans::contrast(emmeans::emmeans(model, "Group"), between),
            emmeans::contrast(
                emmeans::emmeans(model, c("Group", "Session")),
                interaction = list(between, within)
            )
        )
    )
}

test_that("psyci exposes the new inferential schema", {
    fixture <- make_between_fixture()
    result <- psyci(fixture$model, fixture$tables, method = "ind")
    table <- result[[1L]]

    expect_named(
        table,
        c(
            "contrast", "is_polynomial", "estimate", "SE", "df_num", "df_error", "F_value",
            "p_value", "p_value_adjusted", "p_adjustment", "cc", "lower",
            "upper"
        )
    )
    expect_null(attr(result, "legacy_anova_table", exact = TRUE))
    expect_equal(table$df_num, rep(1L, nrow(table)))
    expect_equal(table$F_value, (table$estimate / table$SE)^2)
    expect_equal(
        table$p_value,
        stats::pf(table$F_value, 1, table$df_error, lower.tail = FALSE)
    )
    expect_equal(table$df_error, rep(40, nrow(table)))
})

test_that("legacy_anova is validated", {
    fixture <- make_between_fixture()
    for (invalid in list(NA, 1, "yes", logical(), c(TRUE, FALSE))) {
        expect_error(
            psyci(
                fixture$model,
                fixture$tables,
                method = "ind",
                legacy_anova = invalid
            ),
            "one non-missing logical"
        )
    }
})

test_that("legacy table is attached and printed only on request", {
    fixture <- make_mixed_fixture()
    without <- psyci(fixture$model, fixture$tables, method = "ind")
    with <- psyci(
        fixture$model,
        fixture$tables,
        method = "ind",
        legacy_anova = TRUE
    )
    legacy <- attr(with, "legacy_anova_table", exact = TRUE)

    expect_null(attr(without, "legacy_anova_table", exact = TRUE))
    expect_s3_class(legacy, "psyr_legacy_anova")
    expect_named(
        legacy,
        c("section", "label", "source", "row_type", "SS", "df", "MS", "F_value")
    )
    output_without <- paste(capture.output(print(without)), collapse = "\n")
    output_with <- paste(capture.output(print(with)), collapse = "\n")
    expect_equal(
        lengths(regmatches(
            output_without,
            gregexpr("Analysis of Variance Summary Table", output_without, fixed = TRUE)
        )),
        0L
    )
    expect_equal(
        lengths(regmatches(
            output_with,
            gregexpr("Analysis of Variance Summary Table", output_with, fixed = TRUE)
        )),
        1L
    )
    expect_match(output_with, "Raw contrast effects", fixed = TRUE)
    expect_lt(
        regexpr("Analysis of Variance Summary Table", output_with, fixed = TRUE),
        regexpr("Raw contrast effects", output_with, fixed = TRUE)
    )
})

test_that("legacy table obeys Psy SS, MS, df, and error invariants", {
    fixture <- make_mixed_fixture()
    result <- psyci(
        fixture$model,
        fixture$tables,
        method = "ind",
        legacy_anova = TRUE
    )
    legacy <- attr(result, "legacy_anova_table")
    contrasts <- legacy$row_type == "contrast"
    errors <- legacy$row_type == "error"

    expect_true(all(legacy$df[contrasts] == 1L))
    expect_true(all(legacy$df[errors] == 12L))
    expect_equal(legacy$MS[contrasts], legacy$SS[contrasts])
    expect_true(all(is.na(legacy$F_value[errors])))

    between_error <- legacy$MS[legacy$section == "Between" & errors]
    between_rows <- legacy$section == "Between" & contrasts
    expect_equal(
        legacy$F_value[between_rows],
        legacy$SS[between_rows] / between_error
    )
    for (within_id in c("W1", "W2")) {
        group_rows <- legacy$section == "Within" &
            grepl(paste0(within_id, "$"), legacy$source) & contrasts
        last_group_row <- max(which(group_rows))
        error_ms <- legacy$MS[[last_group_row + 1L]]
        expect_equal(
            legacy$F_value[group_rows],
            legacy$SS[group_rows] / error_ms
        )
    }

    standard_f <- unlist(lapply(result, `[[`, "F_value"), use.names = FALSE)
    legacy_f <- legacy$F_value[contrasts]
    # The standard list is table ordered (between, within, interaction), while
    # the legacy display is between first and then nested by within direction.
    expect_setequal(round(standard_f, 12), round(legacy_f, 12))
})

test_that("legacy output does not synthesize absent interaction rows", {
    fixture <- make_mixed_fixture()
    result <- psyci(
        fixture$model,
        fixture$tables[1:2],
        method = "ind",
        legacy_anova = TRUE
    )
    legacy <- attr(result, "legacy_anova_table")
    expect_false(any(grepl("^B[0-9]+W[0-9]+$", legacy$source)))
    expect_equal(sum(legacy$row_type == "contrast"), 4L)
})

test_that("Bonferroni p-values are two-tailed probabilities from F", {
    fixture <- make_between_fixture()
    result <- psyci(fixture$model, fixture$tables, method = "bf")[[1L]]
    individual <- psyci(fixture$model, fixture$tables, method = "ind")[[1L]]
    family_size <- nrow(result)
    t_value <- result$estimate / result$SE
    expected_t <- pmin(
        1,
        family_size * 2 * stats::pt(
            -abs(t_value),
            df = result$df_error
        )
    )
    expected_f <- pmin(
        1,
        family_size * stats::pf(
            result$F_value,
            df1 = 1,
            df2 = result$df_error,
            lower.tail = FALSE
        )
    )
    expect_equal(result$p_value_adjusted, expected_t)
    expect_equal(result$p_value_adjusted, expected_f)
    expect_equal(
        result$p_value,
        stats::pf(
            result$F_value,
            df1 = 1,
            df2 = result$df_error,
            lower.tail = FALSE
        )
    )
    expect_equal(result$F_value, individual$F_value)
    expect_equal(result$p_value, individual$p_value)
    expect_true(any(
        result$p_value_adjusted != individual$p_value_adjusted
    ))
})

test_that("legacy anxiety statistics reproduce the original Psy table", {
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    afex::afex_options(emmeans_model = "multivariate")
    fixture <- make_legacy_anxiety_fixture()
    result <- psyci(
        fixture$model,
        fixture$tables,
        method = "ph",
        legacy_anova = TRUE
    )
    legacy <- attr(result, "legacy_anova_table")
    contrast_rows <- legacy$row_type == "contrast"
    error_rows <- legacy$row_type == "error"

    expect_equal(
        legacy$SS[contrast_rows],
        c(
            931.612, 1228.538,
            1106.156, 205.967, 148.513,
            1294.252, 470.712, 198.025,
            176.022, 186.336, 114.075
        ),
        tolerance = 1e-3
    )
    expect_equal(
        round(legacy$F_value[contrast_rows], 3),
        c(
            5.192, 6.847,
            42.486, 7.911, 5.704,
            122.790, 44.658, 18.787,
            30.778, 32.581, 19.946
        )
    )
    expect_equal(
        legacy$SS[error_rows],
        c(15610.075, 2265.114, 917.011, 497.567),
        tolerance = 1e-3
    )
    expect_true(all(legacy$df[error_rows] == 87L))
})
