test_that("coef.rhabit", {
  data(tracks)
  fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
  ref <-  5.417143
  expect_equal(unname(coef(fitted_langevin)), 5.417143)
})
