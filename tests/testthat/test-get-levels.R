test_that("get_levels supports data frames and tibble-backed package data", {
    ordinary <- data.frame(group = factor(c("a", "b", "a")))
    expect_identical(get_levels(ordinary, "group"), 2L)

    data("spacing", package = "PsyR", envir = environment())
    spacing$group <- factor(spacing$group)
    expect_identical(get_levels(spacing, "group"), 4L)
    expect_identical(
        compute_df(model = NA, data = spacing, fctrs = list("group")),
        3
    )
})

test_that("get_levels validates its inputs", {
    data <- data.frame(group = c("a", "b"))
    expect_error(get_levels(data, "missing"), "was not found")
    expect_error(get_levels(data, "group"), "is not a factor")
    expect_error(get_levels(data, character()), "one non-empty column name")
})
