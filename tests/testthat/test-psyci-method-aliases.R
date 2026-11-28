test_that("full method names produce identical results to short aliases", {
    fixture <- make_between_fixture()
    aliases <- c(individual = "ind", bonferroni = "bf", "post-hoc" = "ph")

    for (full_name in names(aliases)) {
        full_result <- psyci(
            model = fixture$model,
            contrast_tables = fixture$tables,
            method = full_name
        )
        short_result <- psyci(
            model = fixture$model,
            contrast_tables = fixture$tables,
            method = unname(aliases[[full_name]])
        )
        expect_identical(full_result, short_result, info = full_name)
    }
})

test_that("invalid methods report full names and aliases", {
    fixture <- make_between_fixture()
    expect_error(
        psyci(fixture$model, fixture$tables, method = "unknown"),
        '"individual" \\("ind"\\), "bonferroni" \\("bf"\\), "post-hoc" \\("ph"\\), or "smr"'
    )
})
