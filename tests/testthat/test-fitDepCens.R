



test_that('Check fit dependent censoring', {


  ## Check the fit of dependent censoring model

  #Check summary table
  expect_message(summary(fit))

  # check the scale and dependency parameter are positive

  expect_true(all(parE[(l-1):l]>0))
  expect_equal(fit$copula, "Frank")
  expect_equal(fit$censoringDistribution, "Weibull")
  expect_equal(length(parE), 5)
  expect_true(inherits(fit, "depFit"))
})
