test_that("SII and RII execute correctly under stress scenarios", {
  
  set.seed(789)
  df <- data.frame(
    id = 1:50,
    health = runif(50, 0.1, 0.9),
    pop = runif(50, 1000, 5000),
    equity = rnorm(50, 50000, 10000),
    num = rpois(50, lambda=50),
    denom = rpois(50, lambda=1000)
  )
  
  # SII with wald
  expect_no_error(sii_ineqeco(df, "proportion", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", analysis_unit_var = "id", higher_ineq_is_favorable = TRUE, ci_method = "wald"))
  
  # RII with wald and counts
  expect_no_error(rii_ineqeco(df, "rate", health_numerator_var = "num", health_denominator_var = "denom", equity_stratifier_var = "equity", analysis_unit_var = "id", higher_ineq_is_favorable = TRUE, ci_method = "wald", rate_scaling_factor = 1000))
  
  # SII with bounded social position
  expect_no_error(sii_ineqeco(df, "proportion", health_indicator_var = "health", population_weights_var = "pop", equity_stratifier_var = "equity", analysis_unit_var = "id", higher_ineq_is_favorable = TRUE, ci_method = "wald", social_position_method = "bounded"))
})
