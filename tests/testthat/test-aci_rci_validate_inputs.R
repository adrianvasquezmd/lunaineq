test_that("aci_ineqeco structural checks work", {
  df <- generate_aci_rci_mock_data(50)
  
  expect_error(
    aci_ineqeco(
      data = NULL,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000
    ),
    "`data` must be a data frame."
  )
})

test_that("aci_ineqeco ensures rate scaling factor is provided for rates", {
  df <- generate_aci_rci_mock_data(50)
  
  expect_error(
    aci_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE
      # missing rate_scaling_factor
    ),
    "`rate_scaling_factor` is required"
  )
})

test_that("rci_ineqeco bounded corrections checks", {
  df <- generate_aci_rci_mock_data(50)
  
  expect_error(
    rci_ineqeco(
      data = df,
      health_indicator_type = "proportion",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      higher_ineq_is_favorable = FALSE,
      bounded_corrections = "invalid_correction"
    ),
    "should be one of"
  )
})

test_that("rci_ineqeco handles missing stratifier variables with a warning", {
  df <- generate_aci_rci_mock_data(50)
  df$equity_stratifier[1] <- NA
  
  expect_warning({
    res <- rci_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000,
      verbose = TRUE
    )
  })
})
