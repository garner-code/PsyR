test_that("gcr_crit output matches Psy output", {
  expect_equal(round(gcr_crit(.05, 2, 1, 1),2), 0.93)
})
