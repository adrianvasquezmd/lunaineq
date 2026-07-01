test_that("sii_ineqeco runs successfully with normal data", {
  df <- generate_sii_rii_mock_data(50, "normal")
  
  res <- sii_ineqeco(
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
  
  expect_s3_class(res, "sii_ineqeco_result")
  expect_true(is.numeric(res$results_overall$Estimate))
  expect_true(is.numeric(res$results_overall$Standard_Error))
})

test_that("sii_ineqeco handles model selection gracefully", {
  df <- generate_sii_rii_mock_data(50, "normal")
  
  # Forcing multiple models to be evaluated
  res <- sii_ineqeco(
    data = df,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    models = c("poisson", "negative_binomial"),
    ci_method = "wald"
  )
  
  expect_s3_class(res, "sii_ineqeco_result")
  expect_true(nrow(res$results_models) > 0)
  expect_true(!is.na(res$results_models$AIC[1]))
})

test_that("sii_ineqeco behaves securely under perfect separation or flat relationships", {
  df_sep <- generate_sii_rii_mock_data(50, "perfect_separation")
  df_flat <- generate_sii_rii_mock_data(50, "flat_relationship")
  
  # Flat relationship (slope should be very close to 0)
  res_flat <- sii_ineqeco(
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
  
  expect_s3_class(res_flat, "sii_ineqeco_result")
  
  # Perfect separation could throw warnings or just fit. We don't force expect_warning.
  res_sep <- sii_ineqeco(
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
  
  expect_s3_class(res_sep, "sii_ineqeco_result")
})

test_that("sii_ineqeco robust variance fallback when HC1 fails", {
  df <- generate_sii_rii_mock_data(10, "extreme_overdispersion") # Low N to force unstable covariance
  
  res <- sii_ineqeco(
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
  
  expect_s3_class(res, "sii_ineqeco_result")
})
