test_that("rg_ineqeco runs successfully with normal data", {
  df <- generate_mock_ecological_data(50, "normal")
  
  res <- rg_ineqeco(
    data = df,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "delta"
  )
  
  expect_s3_class(res, "rg_ineqeco")
  expect_true(is.numeric(res$results_overall$Estimate))
  expect_true(res$results_overall$Estimate >= 0)
})

test_that("rg_ineqeco handles missing values correctly (nas_in_indicator)", {
  df <- generate_mock_ecological_data(50, "nas_in_indicator")
  
  expect_warning({
    res <- rg_ineqeco(
      data = df,
      health_indicator_type = "rate",
      health_numerator_var = "cases",
      health_denominator_var = "population",
      analysis_unit_var = "id",
      equity_stratifier_var = "equity_stratifier",
      rate_scaling_factor = 1000,
      grouping_approach = "territorial",
      territorial_method = "quantile",
      n_groups = 5,
      ci_method = "delta"
    )
  }, "ecological unit\\(s\\) were removed")
  
  expect_s3_class(res, "rg_ineqeco")
})

test_that("rg_ineqeco handles zero cases appropriately without breaking", {
  # For relative gaps, a zero in the advantaged group denominator (reference) 
  # can cause Inf or NaN. We test the behavior on zero cases.
  df <- generate_mock_ecological_data(50, "zero_cases")
  
  res <- rg_ineqeco(
    data = df,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "delta"
  )
  
  expect_s3_class(res, "rg_ineqeco")
  # Depending on mathematical handling, we just want to ensure it doesn't crash.
  # the estimate could be Inf, NaN, or numeric depending on which group got zeroes.
})

test_that("rg_ineqeco bootstrap estimation works with edge data", {
  df <- generate_mock_ecological_data(50, "normal")
  
  res <- rg_ineqeco(
    data = df,
    health_indicator_type = "rate",
    health_numerator_var = "cases",
    health_denominator_var = "population",
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "bootstrap",
    n_boot = 100 # Low boot for fast test execution
  )
  
  expect_s3_class(res, "rg_ineqeco")
  expect_true(is.numeric(res$results_overall$Standard_Error))
})
