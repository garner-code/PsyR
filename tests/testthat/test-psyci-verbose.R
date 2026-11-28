test_that("verbosity controls retained notes without changing calculations", {
    fixture <- make_mixed_fixture()
    quiet <- psyci(fixture$model, fixture$tables, "ph")
    detailed <- psyci(fixture$model, fixture$tables, "ph", verbose = TRUE)
    quiet_text <- paste(capture.output(print(quiet)), collapse = "\n")
    detailed_text <- paste(capture.output(print(detailed)), collapse = "\n")
    expect_false(attr(quiet, "verbose"))
    expect_identical(capture.output(print(quiet, verbose = TRUE)),
                     capture.output(print(detailed)))
    expect_identical(capture.output(print(detailed, verbose = FALSE)),
                     capture.output(print(quiet)))
    for (text in c(quiet_text, detailed_text)) {
        expect_match(text, "PsyR used an alpha rate of:", fixed = TRUE)
        expect_match(text, "parameters for Roy GCR:", fixed = TRUE)
        expect_false(grepl("PsyR CI method:|PsyR simultaneous-inference procedure:|PsyR used df error of:", text))
    }
    for (note in c("ordinary F(1, df2)", "Family-wise correction assumes",
                   "PsyR assumed", "non-contrast rows", "Orthogonality",
                   "post-hoc between dimension", "Scaling divisor",
                   "PsyR scaled estimate")) {
        expect_false(grepl(note, quiet_text, fixed = TRUE), info = note)
        expect_true(grepl(note, detailed_text, fixed = TRUE), info = note)
    }
    for (index in seq_along(quiet)) {
        expect_equal(as.data.frame(quiet[[index]]), as.data.frame(detailed[[index]]),
                     ignore_attr = TRUE)
        expect_identical(attr(quiet[[index]], "mesg"), attr(detailed[[index]], "mesg"))
        expect_identical(capture.output(print(quiet[[index]], verbose = TRUE)),
                         capture.output(print(detailed[[index]])))
    }
})

test_that("Bonferroni family size stays visible and verbose is validated", {
    fixture <- make_between_fixture()
    result <- psyci(fixture$model, fixture$tables, "bf")
    text <- paste(capture.output(print(result)), collapse = "\n")
    expect_match(text, "containing k=", fixed = TRUE)
    expect_false(grepl("ordinary F(1, df2)", text, fixed = TRUE))
    for (invalid in list(NA, NULL, c(TRUE, FALSE), 1, "TRUE")) {
        expect_error(psyci(fixture$model, fixture$tables, "bf", verbose = invalid),
                     "`verbose` must be one non-missing logical value.", fixed = TRUE)
        expect_error(print(result, verbose = invalid), "`verbose`", fixed = TRUE)
    }
})

test_that("message filtering keeps only required parameters by default", {
    messages <- c("Results have been averaged over: group",
                  "PsyR used an alpha rate of: 0.05",
                  "Psy probability-calculator parameters for SMR: k, p, q, df",
                  "PsyR used SMR parameter p of: 2",
                  "PsyR used SMR parameter q of: 3",
                  "PsyR used Algina-Keselman noncentral-t confidence limits")
    expect_identical(PsyR:::.psyci_visible_messages(messages, FALSE), messages[2:5])
    expect_identical(PsyR:::.psyci_visible_messages(messages, TRUE), messages)
})
