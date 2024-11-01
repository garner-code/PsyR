test_that("ibeta works", {
  expect_equal(ibeta(0.4, 5, 10), stats::pbeta(0.4,5,10)*beta(5,10))
})
