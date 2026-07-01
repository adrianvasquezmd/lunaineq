test_that(".ag_rg_validate_inputs structural checks work", {
  df <- generate_mock_ecological_data(50)
  
  # Requires data frame
  expect_error(.ag_rg_validate_inputs(
    data = NULL,
    health_indicator_type = "rate",
    health_indicator_var = "indicator",
    population_weights_var = "population",
    health_numerator_var = NULL,
    health_denominator_var = NULL,
    health_se_var = NULL,
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "delta",
    n_boot = 2000,
    conf_level = 0.95,
    verbose = FALSE
  ), "`data` must be supplied")
  
  # Requires at least 2 rows
  expect_error(.ag_rg_validate_inputs(
    data = df[1, , drop=FALSE],
    health_indicator_type = "rate",
    health_indicator_var = "indicator",
    population_weights_var = "population",
    health_numerator_var = NULL,
    health_denominator_var = NULL,
    health_se_var = NULL,
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "delta",
    n_boot = 2000,
    conf_level = 0.95,
    verbose = FALSE
  ), "at least two ecological units")
})

test_that(".ag_rg_validate_inputs missing and mismatched columns throw errors", {
  df <- generate_mock_ecological_data(50)
  
  # Mismatched numerator/denominator
  expect_error(.ag_rg_validate_inputs(
    data = df,
    health_indicator_type = "rate",
    health_indicator_var = NULL,
    population_weights_var = NULL,
    health_numerator_var = "cases",
    health_denominator_var = NULL,
    health_se_var = NULL,
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "delta",
    n_boot = 2000,
    conf_level = 0.95,
    verbose = FALSE
  ), "must be supplied together")
})

test_that(".ag_rg_validate_inputs valid scenarios return standardized df", {
  df <- generate_mock_ecological_data(50)
  
  # counts_no_se
  res <- .ag_rg_validate_inputs(
    data = df,
    health_indicator_type = "rate",
    health_indicator_var = NULL,
    population_weights_var = NULL,
    health_numerator_var = "cases",
    health_denominator_var = "population",
    health_se_var = NULL,
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "delta",
    n_boot = 2000,
    conf_level = 0.95,
    verbose = FALSE
  )
  
  expect_equal(res$scenario, "counts_no_se")
  expect_true(is.data.frame(res$df))
  expect_true("numerator" %in% names(res$df))
})

test_that(".ag_rg_validate_inputs handles edge cases properly", {
  df_neg_denom <- generate_mock_ecological_data(50)
  df_neg_denom$population[1] <- -100
  
  expect_error(.ag_rg_validate_inputs(
    data = df_neg_denom,
    health_indicator_type = "rate",
    health_indicator_var = NULL,
    population_weights_var = NULL,
    health_numerator_var = "cases",
    health_denominator_var = "population",
    health_se_var = NULL,
    analysis_unit_var = "id",
    equity_stratifier_var = "equity_stratifier",
    higher_ineq_is_favorable = FALSE,
    rate_scaling_factor = 1000,
    grouping_approach = "territorial",
    territorial_method = "quantile",
    n_groups = 5,
    ci_method = "delta",
    n_boot = 2000,
    conf_level = 0.95,
    verbose = FALSE
  ), "must contain values > 0")
})
