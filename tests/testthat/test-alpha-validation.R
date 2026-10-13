test_that("critical-constant helpers reject invalid alpha values", {
    invalid <- list(0, 1, -0.1, 1.1, NA_real_, Inf, c(0.05, 0.1), "0.05")
    helpers <- list(
        function(alpha) cc_ind_t(12, alpha),
        function(alpha) cc_bonf_t(12, 3, alpha),
        function(alpha) cc_ph_b(3, 12, alpha),
        function(alpha) cc_ph_w(2, 12, alpha)
    )

    for (helper in helpers) {
        for (alpha in invalid) {
            expect_error(helper(alpha), "strictly between 0 and 1", fixed = TRUE)
        }
    }
})

test_that("psyci validates scalar and family-specific alpha values early", {
    fixture <- make_between_fixture()

    expect_error(
        psyci(fixture$model, fixture$tables, method = "ind", alpha = 0),
        "strictly between 0 and 1",
        fixed = TRUE
    )
    expect_error(
        psyci(fixture$model, fixture$tables, method = "ind", alpha = list(NA_real_)),
        "strictly between 0 and 1",
        fixed = TRUE
    )
})
