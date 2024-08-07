test_that("cc_bonf_t does as it should", {
  expect_equal(round(cc_bonf_t(12,3,0.05),3),
               round((8.839 - 2.958)/2.116, 3))
})
