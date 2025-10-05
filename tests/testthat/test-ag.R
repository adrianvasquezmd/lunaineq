test_that("ag() returns correct class and structure", {
  df <- tibble::tibble(
    health_numerator = c(10, 20, 30, 40),
    health_denominator = c(100, 100, 100, 100),
    population_weights = c(1000, 1000, 1000, 1000),
    equity_stratifier = c(1, 2, 3, 4)
  )

  result <- ag(
    health_numerator = df$health_numerator,
    health_denominator = df$health_denominator,
    equity_stratifier = df$equity_stratifier,
    population_weights = df$population_weights,
    health_indicator_type = "rate"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("inequality_metric", "valor", "ic_inf", "ic_sup") %in% names(result)))
})
