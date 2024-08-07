test_that("t_indiv provides the correct answer, according
          to Psy output", {
  expect_equal(round(t_indiv(12.0, 0.05),3),
               round((7.568 - 2.958)/2.116,3))
})
