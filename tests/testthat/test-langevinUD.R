test_that("fit_langevin_ud", {
  data("tracks")
  expect_error(fit_langevin_ud(x~grad_c1+grad_c2, data = tracks), 
               "the LHS of the formula should have the form")
  expect_error(fit_langevin_ud(x~grad_c3, data = tracks), 
               "the LHS of the formula should have the")
  tracks_no_tims <- tracks %>% dplyr::select(-t)
  expect_error(fit_langevin_ud(cbind(x,y) ~ grad_c1 + grad_c2, data = tracks_no_tims),
               "No t column for specifying times in the dataset!")
  expect_error(fit_langevin_ud(cbind(x,y) ~ grad_c3, data = tracks),
               "Missing variable")
  expect_error(fit_langevin_ud(cbind(x,y) ~ grad_c3, data = tracks),
               "Missing variable")
  
})
