test_that("ACI and RCI execute correctly under stress scenarios", {
  
  set.seed(456)
  df <- data.frame(
    id = 1:50,
    health = runif(50, 0.1, 0.9),
    pop = runif(50, 1000, 5000),
    equity = rnorm(50, 50000, 10000)
  )
  
  # ACI with analytic
  expect_no_error(aci_ineqeco(df, "proportion", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", analysis_unit_var = "id", ci_method = "analytic"))
  
  # RCI with bootstrap
  expect_no_error(rci_ineqeco(df, "proportion", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", ci_method = "bootstrap", n_boot = 100))
  
  # RCI with jackknife
  expect_no_error(rci_ineqeco(df, "proportion", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", ci_method = "jackknife"))
  
  # Test with bounded corrections auto
  expect_no_error(rci_ineqeco(df, "proportion", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", ci_method = "analytic", bounded_corrections = "auto"))
})
