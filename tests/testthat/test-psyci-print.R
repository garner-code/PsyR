test_that("printing combines families without changing stored results", {
    fixture <- make_mixed_fixture()
    result <- psyci(fixture$model, fixture$tables, "ind")
    before <- serialize(result, NULL)
    old_options <- options(width = 300)
    on.exit(options(old_options), add = TRUE)
    output <- capture.output(returned <- withVisible(print(result)))
    text <- paste(output, collapse = "\n")

    expect_identical(serialize(result, NULL), before)
    expect_identical(returned$value, result)
    expect_false(returned$visible)
    expect_equal(sum(grepl("family.*contrast", output)), 2L)
    expect_equal(sum(grepl("^Table [123] ", output)), 6L)
    expect_false(any(grepl("^\\$", output)))
    expect_match(text, "df1 df2", fixed = TRUE)
    expect_match(text, "F_crit", fixed = TRUE)
    expect_match(text, "p_adj", fixed = TRUE)
    expect_match(text, "method", fixed = TRUE)
    expect_false(grepl("p_value|p\\.value|df_num|df_error|is[._]polynomial|\\bCC\\b|\\bcc\\b", text))
    labels <- c(b = "between", w = "within", bw = "mixed")
    for (family in names(labels)) {
        expect_equal(sum(grepl(paste0("^ *[0-9]+ +", labels[[family]], " +"), output)),
                     2L * nrow(result[[family]]))
        expect_match(text, paste0("(", labels[[family]], ")"), fixed = TRUE)
    }
})

test_that("post-hoc mixed labels print together without altering stored tables", {
    fixture <- make_mixed_fixture()
    result <- psyci(fixture$model, fixture$tables, "ph")
    old_options <- options(width = 500)
    on.exit(options(old_options), add = TRUE)
    old_names <- c("df_num", "df_error", "p_value_adjusted", "p_adjustment")
    new_names <- c("df1", "df2", "p_adj", "method")
    for (tables in list(result, attr(result, "scaled_tables"))) {
        for (table in tables) {
            expect_true(all(new_names %in% names(table)))
            expect_false(any(old_names %in% names(table)))
            expect_true(all(c("is_polynomial", "p_value") %in% names(table)))
        }
        table <- tables[["bw"]]
        before <- serialize(table, NULL)
        components <- attr(table, "pri.vars", exact = TRUE)
        expected <- paste0("(", table[[components[1]]], ") \u00d7 (",
                           table[[components[2]]], ")")
        display <- PsyR:::.psyci_display_table(table)
        expect_identical(names(display)[1], "contrast")
        expect_identical(display$contrast, expected)
        expect_false(any(components %in% names(display)))
        output <- paste(capture.output(returned <- withVisible(print(table))),
                        collapse = "\n")
        # R escapes the multiplication sign when the session locale is ASCII.
        output <- gsub("<U+00D7>", "\u00d7", output, fixed = TRUE)
        for (label in expected) expect_match(output, label, fixed = TRUE)
        combined <- paste(capture.output(PsyR:::.print_psyci_tables(tables)),
                          collapse = "\n")
        combined <- gsub("<U+00D7>", "\u00d7", combined, fixed = TRUE)
        for (label in expected) expect_match(combined, label, fixed = TRUE)
        expect_identical(serialize(table, NULL), before)
        expect_identical(returned$value, table)
        expect_false(returned$visible)
    }
})

test_that("label combination preserves grouping and falls back without metadata", {
    table <- structure(data.frame(
        first_custom = c("A - B", "A - C"), second_custom = c("X - Y", "X - Z"),
        group = c("one", "two"), estimate = c(1, 2),
        is_polynomial = FALSE, is.polynomial = FALSE, p_value = 0.1
    ), class = c("psyr_contrast_table", "data.frame"),
    pri.vars = c("first_custom", "second_custom"), by.vars = "group")
    display <- PsyR:::.psyci_display_table(table)
    expect_named(display, c("contrast", "group", "estimate"))
    expect_identical(display$group, table$group)
    table$contrast <- c("custom one", "custom two")
    display <- PsyR:::.psyci_display_table(table)
    expect_identical(display$contrast, table$contrast)
    expect_true(all(c("first_custom", "second_custom") %in% names(display)))
    table$contrast <- NULL
    for (metadata in list(NULL, c("first_custom", "missing"))) {
        attr(table, "pri.vars") <- metadata
        display <- PsyR:::.psyci_display_table(table)
        expect_false("contrast" %in% names(display))
        expect_true(all(c("first_custom", "second_custom", "group") %in% names(display)))
        expect_false(any(c("is_polynomial", "is.polynomial", "p_value") %in% names(display)))
    }
})

test_that("printing formats values and fills absent grouping columns", {
    first <- structure(data.frame(
        contrast = "first", group = "A", estimate = 1.234567,
        SE = NA_real_, df1 = 1L, df2 = 20,
        F_value = Inf, F_crit = 4.123456, p_value = 0.876543,
        p_adj = 0.049999, method = "none", sig = "*",
        lower = NA_real_, upper = NaN
    ), class = c("psyr_contrast_table", "data.frame"))
    second <- first
    second$contrast <- "second"
    second$group <- NULL
    second$sig <- NA_character_
    old_options <- options(width = 300)
    on.exit(options(old_options), add = TRUE)
    output <- capture.output(PsyR:::.print_psyci_tables(list(b = first, w = second)))
    text <- paste(output, collapse = "\n")

    expect_length(output, 3L)
    expect_match(output[2], "1.235", fixed = TRUE)
    expect_match(output[2], "4.123", fixed = TRUE)
    expect_match(output[2], "0.050", fixed = TRUE)
    expect_match(output[2], "1 +20 +Inf")
    expect_match(output[2], "A +1.235")
    expect_match(output[3], "second +1.235")
    expect_false(grepl("NA|NaN|0.876543|1.234567", text))
    expect_match(output[2], "*", fixed = TRUE)
    expect_false(grepl("*", output[3], fixed = TRUE))
})

test_that("significance uses full precision and a strict alpha boundary", {
    fixture <- make_between_fixture()
    baseline <- psyci(fixture$model, fixture$tables, "ind")[[1L]]
    boundary <- baseline$p_adj[4L]
    for (alpha in c(boundary - 1e-8, boundary, boundary + 1e-8)) {
        result <- psyci(fixture$model, fixture$tables, "ind", alpha = alpha)
        expect_identical(result[[1L]]$sig[4L], if (alpha > boundary) "*" else "")
        expect_identical(attr(result, "scaled_tables")[[1L]]$sig, result[[1L]]$sig)
        expect_equal(result[[1L]]$F_crit,
                     rep(stats::qf(1 - alpha, 1, baseline$df2[1]), nrow(baseline)))
        expect_false(any(c("cc", "CC") %in% names(result[[1L]])))
    }
})

test_that("legacy inference follows source rows with duplicate labels", {
    fixture <- make_mixed_fixture()
    # Repeat the between table under a different alpha: labels cannot identify it.
    tables <- c(fixture$tables, fixture$tables[1L])
    alphas <- list(0.05, 0.1, 0.15, 0.2)
    result <- psyci(fixture$model, tables, "ind", alpha = alphas, legacy_summary = TRUE)
    legacy <- attr(result, "legacy_anova_table")
    inference <- attr(legacy, "psyr_inference")
    between <- which(legacy$section == "Between" & legacy$row_type == "contrast")
    expect_equal(inference$F_crit[between], c(result[[1]]$F_crit, result[[4]]$F_crit))
    expect_equal(inference$p_adj[between], c(result[[1]]$p_adj,
                                             result[[4]]$p_adj))
    # Within rows interleave main effects and interactions by direction.
    within <- which(legacy$section == "Within" & legacy$row_type == "contrast")
    expect_equal(inference$F_crit[within], c(
        result[[2]]$F_crit[1], result[[3]]$F_crit[1:2],
        result[[2]]$F_crit[2], result[[3]]$F_crit[3:4]
    ))
    expect_true(all(is.na(inference$F_crit[legacy$row_type == "error"])))
    expect_true(all(is.na(inference$p_adj[legacy$row_type == "error"])))
    for (index in seq_along(result)) {
        expect_identical(result[[index]]$sig,
                         ifelse(result[[index]]$p_adj < alphas[[index]], "*", ""))
    }
    output <- paste(capture.output(print(legacy)), collapse = "\n")
    expect_match(output, "F_crit", fixed = TRUE)
    expect_match(output, "p_adj", fixed = TRUE)
})

test_that("legacy ANOVA sections share column widths and preserve results", {
    table <- structure(data.frame(
        section = c("Between", "Between", "Within", "Within"),
        label = c("Short", "", "Much_longer_interaction_label", ""),
        source = c("B1", "Error", "B12W34", "Error"),
        SS = c(123456.789, 100, 2, 10), df = c(1, 120, 1, 20),
        MS = c(123456.789, 0.833, 2345678.901, 0.5),
        F_value = c(Inf, NA, 12.345, NA)
    ), class = c("psyr_legacy_anova", "data.frame"))
    old_options <- options(width = 300)
    on.exit(options(old_options), add = TRUE)

    for (with_inference in c(FALSE, TRUE)) {
        if (with_inference) {
            attr(table, "psyr_inference") <- data.frame(
                F_crit = c(4.123, NA, 1234.567, NA),
                p_adj = c(0.125, NA, 0.001, NA)
            )
        }
        before <- serialize(table, NULL)
        output <- capture.output(returned <- withVisible(print(table)))
        headers <- grep("Label.*Source.*SS.*df.*MS.*F", output)
        expect_length(headers, 2L)
        expect_identical(output[headers[1]], output[headers[2]])
        expect_identical(output[c(3, 8)], c("Between", "Within"))
        expect_identical(serialize(table, NULL), before)
        expect_identical(returned$value, table)
        expect_false(returned$visible)
        expect_false(any(grepl("NA|NaN", output)))

        columns <- c("Label", "Source", "SS", "df", "MS", "F")
        values <- list(
            c("Short", "B1", "123456.789", "1", "123456.789", "Infinite"),
            c("", "Error", "100.000", "120", "0.833", ""),
            c("Much_longer_interaction_label", "B12W34", "2.000", "1",
              "2345678.901", "12.345"),
            c("", "Error", "10.000", "20", "0.500", "")
        )
        if (with_inference) {
            columns <- c(columns, "F_crit", "p_adj")
            inference_values <- list(c("4.123", "0.125"), c("", ""),
                                     c("1234.567", "0.001"), c("", ""))
            values <- Map(c, values, inference_values)
        }
        # Each cell ends at its header's right edge, including blank error cells.
        header_tokens <- gregexpr("\\S+", output[headers[1]])[[1]]
        ends <- header_tokens + attr(header_tokens, "match.length") - 1L
        expect_length(ends, length(columns))
        starts <- c(1L, head(ends, -1L) + 1L)
        lines <- output[c(headers[1] + 1:2, headers[2] + 1:2)]
        for (row in seq_along(lines)) {
            cells <- trimws(substring(lines[row], starts, ends))
            expect_identical(cells, values[[row]])
        }
    }
})

test_that("legacy ANOVA printing handles single and empty sections", {
    old_options <- options(width = 300)
    on.exit(options(old_options), add = TRUE)
    for (section in c("Between", "Within")) {
        table <- structure(data.frame(
            section = section, label = "Only", source = "B1", SS = 2,
            df = 1, MS = 2, F_value = 4
        ), class = c("psyr_legacy_anova", "data.frame"))
        output <- capture.output(print(table))
        expected <- data.frame(Label = "Only", Source = "B1", SS = "2.000",
                               df = "1", MS = "2.000", F = "4.000")
        expect_identical(output, c(
            "Analysis of Variance Summary Table", "", section,
            capture.output(print(expected, row.names = FALSE, right = TRUE,
                                 quote = FALSE))
        ))
        empty <- table[FALSE, , drop = FALSE]
        output <- capture.output(returned <- withVisible(print(empty)))
        expect_identical(output, c("Analysis of Variance Summary Table",
                                  "No supplied contrasts were available."))
        expect_identical(returned$value, empty)
        expect_false(returned$visible)
    }
})

test_that("Bonferroni notes show k and interpretable terminology", {
    fixture <- make_between_fixture()
    result <- psyci(fixture$model, fixture$tables, "bf")
    notes <- paste(attr(result[[1]], "mesg"), collapse = "\n")
    expect_match(notes, "k=4 interpretable contrast(s)", fixed = TRUE)
    expect_false(grepl("genuine", notes, fixed = TRUE))
    expect_equal(result[[1]]$F_crit,
                 rep(stats::qt(1 - .05 / 8, result[[1]]$df2[1])^2, 4))
})

test_that("non-estimable contrasts have unavailable significance", {
    fixture <- make_between_fixture()
    table <- fixture$tables[[1]]
    table@nbasis <- matrix(table@linfct[1, ], ncol = 1)
    result <- psyci(fixture$model, table, "ind")
    unavailable <- is.na(result[[1]]$p_adj)
    expect_true(any(unavailable))
    expect_true(all(is.na(result[[1]]$sig[unavailable])))
    expect_identical(attr(result, "scaled_tables")[[1]]$sig, result[[1]]$sig)
    output <- paste(capture.output(print(result)), collapse = "\n")
    expect_false(grepl("\\bNA\\b|NaN|nonEst", output))
})
