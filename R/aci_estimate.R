#' Core Analytical Engine for the Absolute Concentration Index (ACI)
#'
#' @description
#' `aci_estimate` performs the principal statistical computation of the Absolute
#' Concentration Index (ACI) using pre-processed fractional ranks (ridits).
#'
#' The function computes the empirical covariance-based Relative Concentration Index (RCI)
#' first, scaling it by the global weighted mean to yield the ACI:
#' \eqn{ACI = \mu * RCI}.
#'
#' Variance is estimated analytically using a grouped approximation of the Kakwani
#' covariance method, or non-parametrically through delegated resampling engines
#' (`.aci_rci_resample`).
#'
#' @importFrom stats qnorm var
#' @noRd

.aci_estimate <- function(prepared, config, seed = NULL, tolerance = 1e-10) {
  warnings_out <- character()
  diagnostics <- list()

  add_warning <- function(x) {
    warnings_out <<- c(warnings_out, x)
    warning(x, call. = FALSE)
    invisible(NULL)
  }

  stop_estimate <- function(x) {
    stop(x, call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # 1. Base Estimation
  # -------------------------------------------------------------------------
  data_units <- prepared$data_units
  y <- data_units$.health
  w <- data_units$.population_weight
  f <- data_units$.population_share
  rank <- data_units$.fractional_rank
  q <- data_units$.cum_health
  q_lag <- data_units$.cum_health_lag

  total_weight <- sum(w)
  weighted_health_total <- sum(w * y)
  mean_health <- weighted_health_total / total_weight

  rci <- (2 / mean_health) * sum(f * y * rank) - 1
  aci <- mean_health * rci

  base_metrics <- data.frame(
    Metric = "ACI",
    Metric_Label = "Absolute Concentration Index",
    Estimate = aci,
    stringsAsFactors = FALSE
  )

  estimate_diag <- list(
    n_units = nrow(data_units),
    total_population_weight = total_weight,
    weighted_health_total = weighted_health_total,
    mean_health = mean_health,
    rci_covariance = rci,
    aci = aci,
    tolerance = tolerance
  )

  estimate_out <- list(
    aci = aci,
    rci = rci,
    mean_health = mean_health,
    total_population_weight = total_weight,
    weighted_health_total = weighted_health_total,
    base_metrics = base_metrics,
    warnings = unique(warnings_out),
    diagnostics = estimate_diag
  )

  # -------------------------------------------------------------------------
  # 2. Corrections (None for ACI)
  # -------------------------------------------------------------------------
  corrections_out <- list(
    metrics = base_metrics,
    correction_factors = c(ACI = mean_health),
    correction_metadata = data.frame(
      Metric = "ACI",
      Correction_Method = "absolute_generalized",
      Transformation_Factor = mean_health,
      Formula = "mu * C",
      Requires_Bounded_Outcome = FALSE,
      Computed = TRUE,
      Failure_Reason = NA_character_,
      stringsAsFactors = FALSE
    ),
    warnings = character(),
    diagnostics = list(bounded_outcome = FALSE)
  )

  # -------------------------------------------------------------------------
  # 3. Confidence Intervals
  # -------------------------------------------------------------------------
  ci_method <- config$ci_method
  conf_level <- config$conf_level
  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha / 2)

  make_results <- function(standard_errors, ci_lower, ci_upper, ci_method_label, se_method_label, interval_type, n_replicates = NA_integer_) {
    out <- base_metrics
    out$Standard_Error <- as.numeric(standard_errors[out$Metric])
    out$CI_Lower <- as.numeric(ci_lower[out$Metric])
    out$CI_Upper <- as.numeric(ci_upper[out$Metric])
    out$CI_Method <- ci_method_label
    out$SE_Method <- se_method_label
    out$Interval_Type <- interval_type
    out$Confidence_Level <- conf_level
    out$N_Replicates <- n_replicates
    out
  }

  ci_diag <- list(
    ci_method = ci_method,
    z_value = z_value,
    health_se_supplied = !is.null(config$health_se_var)
  )

  if (ci_method == "analytic") {
    n_groups <- nrow(data_units)
    a_i <- (y / mean_health) * (2 * rank - 1 - rci) + 2 - q_lag - q
    fa2 <- f * (a_i^2)

    has_health_se <- !is.null(config$health_se_var)
    if (!has_health_se) {
      raw_var_rci <- (sum(fa2) - (1 + rci)^2) / n_groups
      se_method <- "KWvD_grouped_no_external_se"
    } else {
      health_se <- data_units$.health_se
      first_component <- (sum(fa2) - (1 + rci)^2) / total_weight
      second_component_numerator <- sum(f * (health_se^2) * ((2 * rank - 1 - rci)^2))
      second_component <- second_component_numerator / (total_weight * mean_health^2)
      raw_var_rci <- first_component + second_component
      se_method <- "KWvD_grouped_with_ecological_health_se"
    }

    if (raw_var_rci < 0) {
      var_rci <- 0
    } else {
      var_rci <- raw_var_rci
    }

    se_rci <- sqrt(var_rci)
    se_aci <- mean_health * se_rci

    standard_errors <- c(ACI = se_aci)
    ci_lower <- c(ACI = aci - z_value * se_aci)
    ci_upper <- c(ACI = aci + z_value * se_aci)

    results_overall <- make_results(standard_errors, ci_lower, ci_upper, "analytic", se_method, "Wald_normal")
    results_replicates <- NULL

    ci_diag$se_method <- se_method
    ci_diag$interval_type <- "Wald_normal"
    ci_diag$n_groups <- n_groups
    ci_diag$raw_var_rci <- raw_var_rci
    ci_diag$se_rci <- se_rci

  } else {
    resample_out <- .aci_rci_resample(
      prepared = prepared,
      config = config,
      method = ci_method,
      n_boot = config$n_boot,
      seed = seed,
      verbose = isTRUE(config$verbose)
    )
    
    if (length(resample_out$warnings) > 0L) {
      warnings_out <<- c(warnings_out, resample_out$warnings)
    }
    
    results_replicates <- resample_out$replicates
    results_replicates <- results_replicates[, c("Replicate_ID", "Replicate_Type", "ACI"), drop = FALSE]
    
    valid_vals <- results_replicates$ACI
    valid_vals <- valid_vals[is.finite(valid_vals)]
    n_valid <- length(valid_vals)

    if (ci_method == "bootstrap") {
      se_aci <- stats::sd(valid_vals)
      q <- stats::quantile(valid_vals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE, type = 7)
      lower <- q[[1L]]
      upper <- q[[2L]]
      se_method <- "bootstrap_ecological_units"
      interval_type <- "bootstrap_percentile"
    } else {
      theta_bar <- mean(valid_vals)
      se_aci <- sqrt((n_valid - 1) / n_valid * sum((valid_vals - theta_bar)^2))
      lower <- aci - z_value * se_aci
      upper <- aci + z_value * se_aci
      se_method <- "jackknife_ecological_units"
      interval_type <- "jackknife_Wald_normal"
    }

    standard_errors <- c(ACI = se_aci)
    ci_lower <- c(ACI = lower)
    ci_upper <- c(ACI = upper)

    results_overall <- make_results(standard_errors, ci_lower, ci_upper, ci_method, se_method, interval_type, n_valid)

    ci_diag$se_method <- se_method
    ci_diag$interval_type <- interval_type
    ci_diag$n_replicates_valid <- n_valid
    ci_diag$resample_diagnostics <- resample_out$diagnostics
  }

  ci_out <- list(
    results_overall = results_overall,
    results_replicates = results_replicates,
    warnings = unique(warnings_out),
    diagnostics = ci_diag
  )

  list(
    estimate = estimate_out,
    corrections = corrections_out,
    ci = ci_out
  )
}
