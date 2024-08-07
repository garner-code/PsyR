test_that("cc_ph_w matches Psy outputs", {
  expect_equal(round(cc_ph_w(2.0,12.0,0.05),3),
               round((0.703--1.750)/0.832,3))
})
