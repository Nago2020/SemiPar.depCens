


test_that('Check fit independent censoring', {


  ## Check the fit of dependent censoring model

  # Check summary table
  expect_message(summary(fitI))

  # check the scale and dependency parameter are positive

  expect_true(inherits(fitI, "indepFit"))
  expect_true(parI[m]>0)
  expect_equal(fitI$censoringDistribution, "Weibull")
  expect_equal(length(fitI$parameterEstimates), 4)
})
