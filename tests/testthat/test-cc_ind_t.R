test_that("t_indiv provides the correct answer, according
          to Psy output", {
  expect_equal(round(cc_ind_t(12.0, 0.05),3),
               round((7.568 - 2.958)/2.116,3))
})
