test_that("ibeta works", {
  expect_equal(ibeta(0.4, 5, 10), stats::pbeta(x,a,b)*beta(a,b))
})
