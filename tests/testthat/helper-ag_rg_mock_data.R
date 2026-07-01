# Helper functions for generating extreme ecological mock data

#' Generate mock ecological data for testing
#'
#' @param n Number of territorial units
#' @param scenario "normal", "zero_pop", "zero_cases", "nas_in_stratifier", "nas_in_indicator", "outliers", "ties", "negative_indicator"
#' @return A data.frame
generate_mock_ecological_data <- function(n = 50, scenario = "normal") {
  set.seed(42) # Fixed seed for reproducibility in tests
  
  df <- data.frame(
    id = 1:n,
    population = round(runif(n, 1000, 100000)),
    equity_stratifier = rnorm(n, 50, 15),
    cases = NA_integer_,
    indicator = NA_real_,
    se = NA_real_
  )
  
  # Base cases and indicator
  df$cases <- round(df$population * runif(n, 0.001, 0.05))
  df$indicator <- df$cases / df$population
  df$se <- sqrt((df$indicator * (1 - df$indicator)) / df$population)
  
  if (scenario == "zero_pop") {
    # Introduct zero population in some regions
    idx <- sample(1:n, min(3, n))
    df$population[idx] <- 0
    df$cases[idx] <- 0
    df$indicator[idx] <- 0
    df$se[idx] <- 0
  } else if (scenario == "zero_cases") {
    # Introduce zero cases (indicator = 0)
    idx <- sample(1:n, min(5, n))
    df$cases[idx] <- 0
    df$indicator[idx] <- 0
    df$se[idx] <- 0
  } else if (scenario == "nas_in_stratifier") {
    idx <- sample(1:n, min(5, n))
    df$equity_stratifier[idx] <- NA
  } else if (scenario == "nas_in_indicator") {
    idx <- sample(1:n, min(5, n))
    df$indicator[idx] <- NA
    df$cases[idx] <- NA
  } else if (scenario == "outliers") {
    # Extreme outlier in stratifier
    df$equity_stratifier[1] <- 999999
    # Extreme outlier in indicator
    df$indicator[2] <- 0.99
    df$cases[2] <- round(df$population[2] * 0.99)
  } else if (scenario == "ties") {
    # Many units have exactly the same equity stratifier
    df$equity_stratifier <- rep(c(10, 20, 30, 40, 50), length.out = n)
  } else if (scenario == "negative_indicator") {
    df$indicator[1:3] <- -0.05
  } else if (scenario == "extreme_weights") {
    df$population[1] <- 10000000 # 10 Million
    df$population[2:n] <- 1      # 1 person
    df$cases[1] <- 500000
    df$cases[2:n] <- 0
    df$indicator <- df$cases / df$population
  }
  
  return(df)
}
