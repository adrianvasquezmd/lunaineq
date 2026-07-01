test_that("rii_ineqeco runs successfully with normal data", {
  df <- generate_sii_rii_mock_data(50, "normal")
  
  res <- rii_ineqeco(
    data = df,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    ci_method = "wald",
    vcov_type = "HC1"
  )
  
  expect_s3_class(res, "rii_ineqeco_result")
  expect_true(is.numeric(res$results_overall$Estimate))
  expect_true(is.numeric(res$results_overall$Standard_Error))
})

test_that("rii_ineqeco behaves securely under perfect separation or flat relationships", {
  df_sep <- generate_sii_rii_mock_data(50, "perfect_separation")
  df_flat <- generate_sii_rii_mock_data(50, "flat_relationship")
  
  # Flat relationship (RII should be very close to 1)
  res_flat <- rii_ineqeco(
    data = df_flat,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    ci_method = "wald"
  )
  
  expect_s3_class(res_flat, "rii_ineqeco_result")
  
  res_sep <- rii_ineqeco(
    data = df_sep,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    ci_method = "wald"
  )
  
  expect_s3_class(res_sep, "rii_ineqeco_result")
})

test_that("rii_ineqeco handles missing outcome values correctly", {
  df <- generate_sii_rii_mock_data(50, "nas_in_outcome")
  
  expect_error({
    res <- rii_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000,
      ci_method = "wald"
    )
  }, "must be numeric, finite, and non-missing")
})
