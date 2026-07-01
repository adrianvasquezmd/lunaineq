test_that("sii_ineqeco structural checks work", {
  df <- generate_sii_rii_mock_data(50)
  
  expect_error(
    sii_ineqeco(
      data = NULL,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000
    ),
    "`data` must be a data.frame or tibble."
  )
})

test_that("sii_ineqeco unsupported arguments are rejected", {
  df <- generate_sii_rii_mock_data(50)
  
  expect_error(
    sii_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000,
      health_se_var = "se_column" # Not supported for SII/RII
    ),
    "not part of this version"
  )
  
  expect_error(
    sii_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000,
      made_up_argument = 123
    ),
    "Unsupported argument"
  )
})

test_that("sii_ineqeco handles missing stratifier variables with warning", {
  df <- generate_sii_rii_mock_data(50, "nas_in_stratifier")
  
  expect_error({
    res <- sii_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000
    )
  }, "must be numeric, finite, and non-missing")
})

test_that("rii_ineqeco ensures rate scaling factor is provided for rates", {
  df <- generate_sii_rii_mock_data(50)
  
  expect_error(
    rii_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE
      # missing rate_scaling_factor
    ),
    "`rate_scaling_factor` must be a single numeric value > 0"
  )
})
