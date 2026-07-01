# Helper functions for generating extreme ACI/RCI mock data

#' Generate mock ecological data for testing ACI/RCI computations
#'
#' @param n Number of territorial units
#' @param scenario "normal", "perfect_equality", "perfect_inequality", "nas_in_outcome"
#' @return A data.frame
generate_aci_rci_mock_data <- function(n = 50, scenario = "normal") {
  set.seed(123)
  
  df <- data.frame(
    id = 1:n,
    population = round(runif(n, 5000, 150000)),
    equity_stratifier = rnorm(n, 50, 15),
    cases = NA_integer_,
    indicator = NA_real_
  )
  
  if (scenario == "normal") {
    df$cases <- round(df$population * runif(n, 0.01, 0.10))
    df$indicator <- df$cases / df$population
  } else if (scenario == "perfect_equality") {
    # Indicator is completely equal across all strata
    df$indicator <- 0.05
    df$cases <- round(df$population * df$indicator)
  } else if (scenario == "perfect_inequality") {
    # All cases concentrated in one single stratifier value (the most advantaged or disadvantaged)
    df$cases <- 0
    top_index <- which.max(df$equity_stratifier)
    df$cases[top_index] <- df$population[top_index] # 100% prevalence in the richest
    df$indicator <- df$cases / df$population
  } else if (scenario == "nas_in_outcome") {
    df$cases <- round(df$population * runif(n, 0.01, 0.10))
    df$indicator <- df$cases / df$population
    df$cases[sample(1:n, 5)] <- NA
    df$indicator[sample(1:n, 5)] <- NA
  }
  
  return(df)
}
