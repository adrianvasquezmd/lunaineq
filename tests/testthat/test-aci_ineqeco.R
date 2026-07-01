test_that("aci_ineqeco runs successfully with normal data", {
  df <- generate_aci_rci_mock_data(50, "normal")
  
  res <- aci_ineqeco(
    data = df,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    ci_method = "analytic"
  )
  
  expect_s3_class(res, "aci_ineqeco_result")
  expect_true(is.numeric(res$results_overall$Estimate))
  expect_true(is.numeric(res$results_overall$Standard_Error))
  # ACI returns an estimate > 0 or < 0 but not NA
  expect_true(!is.na(res$results_overall$Estimate[1]))
})

test_that("aci_ineqeco computes bounded corrections when outcome is bounded", {
  df <- generate_aci_rci_mock_data(50, "normal")
  
  res <- aci_ineqeco(
    data = df,
    health_indicator_type = "proportion", # Proportions are bounded [0, 1]
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    ci_method = "analytic"
  )
  
  expect_s3_class(res, "aci_ineqeco_result")
})

test_that("aci_ineqeco behaves securely under perfect equality", {
  df_eq <- generate_aci_rci_mock_data(50, "perfect_equality")
  
  res_eq <- aci_ineqeco(
    data = df_eq,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    ci_method = "analytic"
  )
  
  expect_s3_class(res_eq, "aci_ineqeco_result")
  
  # ACI should be very close to 0 if everyone has the same rate
  base_aci <- res_eq$results_overall$Estimate[1]
  expect_true(abs(base_aci) < 0.01)
})

test_that("aci_ineqeco handles missing outcome values with a warning", {
  df <- generate_aci_rci_mock_data(50, "nas_in_outcome")
  
  expect_warning({
    res <- aci_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      equity_stratifier_var = "equity_stratifier",
      analysis_unit_var = "id",
      higher_ineq_is_favorable = FALSE,
      rate_scaling_factor = 1000,
      ci_method = "analytic",
      verbose = TRUE
    )
  })
})
