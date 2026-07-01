# Helper functions for generating extreme regression mock data

#' Generate mock ecological data for testing SII/RII regressions
#'
#' @param n Number of territorial units
#' @param scenario "normal", "perfect_separation", "flat_relationship", "extreme_overdispersion", "nas_in_stratifier", "zero_variance_stratifier"
#' @return A data.frame
generate_sii_rii_mock_data <- function(n = 50, scenario = "normal") {
  set.seed(123)
  
  df <- data.frame(
    id = 1:n,
    population = round(runif(n, 5000, 150000)),
    equity_stratifier = rnorm(n, 50, 15),
    cases = NA_integer_,
    indicator = NA_real_
  )
  
  # Base cases
  df$cases <- round(df$population * runif(n, 0.01, 0.10))
  df$indicator <- df$cases / df$population
  
  if (scenario == "perfect_separation") {
    # All cases in the top half of the stratifier, 0 cases in bottom half
    median_strat <- median(df$equity_stratifier)
    df$cases[df$equity_stratifier < median_strat] <- 0
    df$indicator <- df$cases / df$population
  } else if (scenario == "flat_relationship") {
    # Cases exactly proportional to population
    df$cases <- round(df$population * 0.05)
    df$indicator <- df$cases / df$population
  } else if (scenario == "extreme_overdispersion") {
    # Tiny population, huge cases (impossible proportions but useful for testing rate model robustness)
    df$population <- 1:n
    df$cases <- round(runif(n, 10000, 50000))
    df$indicator <- df$cases / df$population
  } else if (scenario == "zero_variance_stratifier") {
    df$equity_stratifier <- 100
  } else if (scenario == "nas_in_stratifier") {
    df$equity_stratifier[sample(1:n, 5)] <- NA
  } else if (scenario == "nas_in_outcome") {
    df$cases[sample(1:n, 5)] <- NA
    df$indicator <- df$cases / df$population
  }
  
  return(df)
}
