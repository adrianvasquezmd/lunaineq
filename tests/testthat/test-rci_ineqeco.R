test_that("rci_ineqeco runs successfully with normal data", {
  df <- generate_aci_rci_mock_data(50, "normal")
  
  res <- rci_ineqeco(
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
  
  expect_s3_class(res, "rci_ineqeco_result")
  expect_true(is.numeric(res$results_overall$Estimate))
  expect_true(is.numeric(res$results_overall$Standard_Error))
})

test_that("rci_ineqeco behaves securely under perfect equality", {
  df_eq <- generate_aci_rci_mock_data(50, "perfect_equality")
  
  res_eq <- rci_ineqeco(
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
  
  expect_s3_class(res_eq, "rci_ineqeco_result")
  
  # RCI should be very close to 0 if everyone has the same rate
  base_rci <- res_eq$results_overall$Estimate[1]
  expect_true(abs(base_rci) < 0.01)
})

test_that("rci_ineqeco calculates bounds for bounded outcome types", {
  df <- generate_aci_rci_mock_data(50, "normal")
  
  res <- rci_ineqeco(
    data = df,
    health_indicator_type = "proportion", # Proportions are bounded [0, 1]
    health_numerator_var = "cases",
    health_denominator_var = "population",
    equity_stratifier_var = "equity_stratifier",
    analysis_unit_var = "id",
    higher_ineq_is_favorable = FALSE,
    bounded_corrections = "wagstaff",
    ci_method = "analytic"
  )
  
  expect_s3_class(res, "rci_ineqeco_result")
  expect_true(any(grepl("Wagstaff", res$results_overall$Metric)))
})

test_that("rci_ineqeco handles perfect inequality gracefully", {
  df <- generate_aci_rci_mock_data(50, "perfect_inequality")
  
  res <- rci_ineqeco(
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
  
  expect_s3_class(res, "rci_ineqeco_result")
  expect_true(is.numeric(res$results_overall$Estimate))
})
