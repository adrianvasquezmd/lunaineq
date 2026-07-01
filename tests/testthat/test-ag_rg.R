test_that("AG and RG execute correctly under stress scenarios", {
  
  # Basic data with no missing values
  set.seed(123)
  df <- data.frame(
    id = 1:50,
    health = runif(50, 0, 1),
    pop = runif(50, 1000, 5000),
    equity = rnorm(50, 50000, 10000)
  )
  
  # Absolute Gap (AG) with different grouping methods
  expect_no_error(ag_ineqeco(df, "rate", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", analysis_unit_var = "id", rate_scaling_factor = 1000, grouping_approach = "territorial", territorial_method = "quantile", n_groups = 5, ci_method = "delta"))
  expect_no_error(ag_ineqeco(df, "rate", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", rate_scaling_factor = 1000, grouping_approach = "territorial", territorial_method = "fisher", n_groups = 5, ci_method = "delta"))
  expect_no_error(ag_ineqeco(df, "rate", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", rate_scaling_factor = 1000, grouping_approach = "territorial", territorial_method = "kmeans", n_groups = 5, ci_method = "delta"))
  
  # Relative Gap (RG) with bootstrap
  expect_no_error(rg_ineqeco(df, "rate", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", rate_scaling_factor = 1000, grouping_approach = "territorial", territorial_method = "quantile", n_groups = 5, ci_method = "bootstrap", n_boot = 100))
  
  # Missing data handling (should silently remove NAs if complete.cases is used inside)
  df_na <- df
  df_na$health[c(1, 10)] <- NA
  expect_no_error(ag_ineqeco(df_na, "rate", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", rate_scaling_factor = 1000, grouping_approach = "territorial", territorial_method = "quantile", n_groups = 5, ci_method = "delta"))
})
