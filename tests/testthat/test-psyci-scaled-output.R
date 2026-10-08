expected_psy_sample_sd <- function(model) {
    data <- model$data$long
    response <- attr(model, "dv")
    between <- names(attr(model, "between"))
    within <- names(attr(model, "within"))
    between_cell <- if (length(between) == 0L) {
        factor(rep("1", nrow(data)))
    } else {
        interaction(data[between], drop = TRUE, lex.order = TRUE)
    }
    within_cell <- if (length(within) == 0L) {
        factor(rep("1", nrow(data)))
    } else {
        interaction(data[within], drop = TRUE, lex.order = TRUE)
    }

    variances <- vapply(split(seq_len(nrow(data)), within_cell), function(rows) {
        groups <- droplevels(between_cell[rows])
        residuals <- data[[response]][rows] -
            ave(data[[response]][rows], groups, FUN = mean)
        sum(residuals^2) / (length(rows) - nlevels(groups))
    }, numeric(1))
    sqrt(mean(variances))
}

test_that("psyci scales between-subject results using Psy's pooled sample SD", {
    fixture <- make_between_fixture()
    observed <- psyci(
        fixture$model,
        fixture$tables,
        method = "ind"
    )
    scale_info <- attr(observed, "scale_info")
    scaled <- attr(observed, "scaled_tables")
    expected_divisor <- expected_psy_sample_sd(fixture$model)

    expect_s3_class(observed, "psyr_ci")
    expect_true(isTRUE(scale_info$available))
    expect_equal(scale_info$divisor, expected_divisor, tolerance = 1e-12)
    expect_named(scaled, names(observed))

    raw <- observed[[1L]]
    standardized <- scaled[[1L]]
    for (column in c("estimate", "SE", "lower", "upper")) {
        expect_equal(
            standardized[[column]],
            raw[[column]] / expected_divisor,
            tolerance = 1e-12
        )
    }
    for (column in c("df", "t.ratio", "p.value", "cc", "psyr_p_value")) {
        expect_equal(standardized[[column]], raw[[column]])
    }
    expect_equal(attr(standardized, "psyr_scale_divisor"), expected_divisor)
    expect_identical(attr(standardized, "psyr_scale_units"), "sample SD")
})

test_that("psyci pools SDs over repeated-measure cells for mixed designs", {
    old_emmeans_model <- afex::afex_options("emmeans_model")
    on.exit(
        afex::afex_options(emmeans_model = old_emmeans_model),
        add = TRUE
    )
    fixture <- make_mixed_fixture()
    observed <- psyci(
        fixture$model,
        fixture$tables,
        method = "ind"
    )
    scale_info <- attr(observed, "scale_info")
    expected_divisor <- expected_psy_sample_sd(fixture$model)

    expect_equal(length(scale_info$cell_sd), 3L)
    expect_equal(scale_info$divisor, expected_divisor, tolerance = 1e-12)
    expect_equal(
        scale_info$divisor,
        sqrt(mean(scale_info$cell_sd^2)),
        tolerance = 1e-12
    )
    expect_length(attr(observed, "scaled_tables"), 3L)
})

test_that("printing a psyci result shows raw and scaled tables", {
    fixture <- make_between_fixture()
    observed <- psyci(fixture$model, fixture$tables, method = "ind")
    output <- paste(capture.output(print(observed)), collapse = "\n")

    expect_match(
        output,
        "Raw contrast effects and confidence intervals",
        fixed = TRUE
    )
    expect_match(
        output,
        "Scaled contrast effects and confidence intervals",
        fixed = TRUE
    )
    expect_match(output, "Scaling divisor (pooled sample SD)", fixed = TRUE)
    invisible(capture.output(returned <- print(observed)))
    expect_identical(returned, observed)
})

test_that("zero sample variance leaves raw results available", {
    fixture <- make_between_fixture()
    zero_variance_model <- fixture$model
    response <- attr(zero_variance_model, "dv")
    zero_variance_model$data$long[[response]] <- 1

    observed <- psyci(
        zero_variance_model,
        fixture$tables,
        method = "ind"
    )
    scale_info <- attr(observed, "scale_info")
    output <- paste(capture.output(print(observed)), collapse = "\n")

    expect_false(scale_info$available)
    expect_null(attr(observed, "scaled_tables"))
    expect_match(output, "Scaled contrast effects", fixed = TRUE)
    expect_match(output, "Unavailable:", fixed = TRUE)
    expect_s3_class(observed[[1L]], "summary_emm")
})
