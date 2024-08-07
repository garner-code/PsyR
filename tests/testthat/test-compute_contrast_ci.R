test_that("compute_contrast_ci multiplies two values", {
  expect_equal(compute_contrast_ci(2.2,3.5), 2.2 * 3.5)
})
