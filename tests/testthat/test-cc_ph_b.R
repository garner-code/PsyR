test_that("cc_ph_b replicates Psy output", {
  expect_equal(round(cc_ph_b(3.0,12.0,0.05),3),
               round((9.805 - 2.958)/2.116 ,3))
})
