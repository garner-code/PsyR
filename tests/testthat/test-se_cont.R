test_that("se_cont works", {
  expect_equal(se_cont(2.2, c(1,-1), 50),
               sqrt(2.2 * sum(c(1,-1)^2/50)))
})
