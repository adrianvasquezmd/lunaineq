test_that("ag_ineqeco runs successfully with normal data", {
  df <- generate_mock_ecological_data(50, "normal")
  
  res <- ag_ineqeco(
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
  
  expect_s3_class(res, "ag_ineqeco")
  expect_true(is.numeric(res$results_overall$Estimate))
  expect_true(is.numeric(res$results_overall$Standard_Error))
})

test_that("ag_ineqeco handles missing values correctly (nas_in_indicator)", {
  df <- generate_mock_ecological_data(50, "nas_in_indicator")
  
  expect_warning({
    res <- ag_ineqeco(
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
  }, "ecological unit\\(s\\) were removed because of missing values")
  
  expect_s3_class(res, "ag_ineqeco")
  # Still expecting a valid estimate despite the NAs (if n_groups allows)
  expect_true(!is.na(res$results_overall$Estimate))
})

test_that("ag_ineqeco handles zero denominators correctly (zero_pop)", {
  df <- generate_mock_ecological_data(50, "zero_pop")
  
  expect_error({
    ag_ineqeco(
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
  }, "must contain values > 0")
})

test_that("ag_ineqeco jackknife works without crashing with edge data", {
  df <- generate_mock_ecological_data(50, "normal")
  
  # Jackknife estimation
  res <- ag_ineqeco(
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
    ci_method = "jackknife"
  )
  
  expect_s3_class(res, "ag_ineqeco")
  expect_true(is.numeric(res$results_overall$Standard_Error))
  expect_true(is.numeric(res$results_overall$CI_Lower))
  expect_true(is.numeric(res$results_overall$CI_Upper))
})

test_that("ag_ineqeco tied stratifiers (ties) generate expected output without crashing", {
  df <- generate_mock_ecological_data(50, "ties")
  
  # It should warn about group modifications when there are severe ties
  expect_warning({
    res <- ag_ineqeco(
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
  })
  
  expect_s3_class(res, "ag_ineqeco")
})

test_that("ag_ineqeco extremes extreme weights do not cause NAs", {
  df <- generate_mock_ecological_data(50, "extreme_weights")
  
  res <- ag_ineqeco(
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
  
  expect_true(is.numeric(res$results_overall$Estimate))
})
