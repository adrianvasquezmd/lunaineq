#' Core Analytical Engine for the Relative Concentration Index (RCI)
#'
#' @description
#' `rci_estimate` performs the principal statistical computation of the Relative
#' Concentration Index (RCI) using pre-processed fractional ranks (ridits).
#'
#' The function computes the empirical covariance-based Relative Concentration Index (RCI).
#' When specified and valid for bounded outcome spaces, it systematically applies Wagstaff
#' normalization and Erreygers corrections to account for structural metric bounds.
#'
#' Variance is estimated analytically using a grouped approximation of the Kakwani
#' covariance method, or non-parametrically through delegated resampling engines
#' (`.aci_rci_resample`).
#'
#' @importFrom stats qnorm var
#' @noRd

.rci_estimate <- function(prepared, config, seed = NULL, tolerance = 1e-10) {
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

  base_metrics <- data.frame(
    Metric = "RCI",
    Metric_Label = "Relative Concentration Index",
    Estimate = rci,
    stringsAsFactors = FALSE
  )

  estimate_diag <- list(
    n_units = nrow(data_units),
    total_population_weight = total_weight,
    weighted_health_total = weighted_health_total,
    mean_health = mean_health,
    rci_covariance = rci,
    tolerance = tolerance
  )

  estimate_out <- list(
    rci = rci,
    mean_health = mean_health,
    total_population_weight = total_weight,
    weighted_health_total = weighted_health_total,
    base_metrics = base_metrics,
    warnings = unique(warnings_out),
    diagnostics = estimate_diag
  )

  # -------------------------------------------------------------------------
  # 2. Corrections
  # -------------------------------------------------------------------------
  metrics <- base_metrics
  correction_factors <- c(RCI = 1)
  correction_metadata <- data.frame(
    Metric = "RCI",
    Correction_Method = "none",
    Transformation_Factor = 1,
    Formula = "C",
    Requires_Bounded_Outcome = FALSE,
    Computed = TRUE,
    Failure_Reason = NA_character_,
    stringsAsFactors = FALSE
  )

  bounded_outcome <- isTRUE(config$bounded_outcome)
  lower_bound <- config$health_lower_bound
  upper_bound <- config$health_upper_bound
  correction_methods <- config$correction_methods
  if (is.null(correction_methods)) correction_methods <- character()
  correction_methods <- unique(tolower(correction_methods))

  computed_methods <- character()

  if (bounded_outcome && length(correction_methods) > 0L) {
    if (is.null(lower_bound) || is.null(upper_bound)) {
      add_warning("Bounded correction methods were requested, but lower and/or upper bounds are missing.")
    } else {
      bounds_range <- as.numeric(upper_bound) - as.numeric(lower_bound)

      # Wagstaff
      if ("wagstaff" %in% correction_methods) {
        wagstaff_denominator <- (upper_bound - mean_health) * (mean_health - lower_bound)
        if (!is.finite(wagstaff_denominator) || abs(wagstaff_denominator) <= tolerance) {
          correction_metadata <- rbind(correction_metadata, data.frame(
            Metric = "Wagstaff_Corrected_RCI", Correction_Method = "wagstaff", Transformation_Factor = NA_real_,
            Formula = "k_W * C", Requires_Bounded_Outcome = TRUE, Computed = FALSE,
            Failure_Reason = "Mean at or too close to bounds.", stringsAsFactors = FALSE
          ))
        } else {
          k_wagstaff <- mean_health * bounds_range / wagstaff_denominator
          wagstaff_estimate <- k_wagstaff * rci
          metrics <- rbind(metrics, data.frame(Metric = "Wagstaff_Corrected_RCI", Metric_Label = "Wagstaff-normalized Relative Concentration Index", Estimate = wagstaff_estimate, stringsAsFactors = FALSE))
          correction_factors <- c(correction_factors, Wagstaff_Corrected_RCI = k_wagstaff)
          correction_metadata <- rbind(correction_metadata, data.frame(
            Metric = "Wagstaff_Corrected_RCI", Correction_Method = "wagstaff", Transformation_Factor = k_wagstaff,
            Formula = "k_W * C", Requires_Bounded_Outcome = TRUE, Computed = TRUE, Failure_Reason = NA_character_, stringsAsFactors = FALSE
          ))
          computed_methods <- c(computed_methods, "wagstaff")
        }
      }

      # Erreygers
      if ("erreygers" %in% correction_methods) {
        if (!is.finite(bounds_range) || bounds_range <= tolerance) {
          correction_metadata <- rbind(correction_metadata, data.frame(
            Metric = "Erreygers_Corrected_RCI", Correction_Method = "erreygers", Transformation_Factor = NA_real_,
            Formula = "k_E * C", Requires_Bounded_Outcome = TRUE, Computed = FALSE,
            Failure_Reason = "Invalid bounds.", stringsAsFactors = FALSE
          ))
        } else {
          k_erreygers <- 4 * mean_health / bounds_range
          erreygers_estimate <- k_erreygers * rci
          metrics <- rbind(metrics, data.frame(Metric = "Erreygers_Corrected_RCI", Metric_Label = "Erreygers-corrected Relative Concentration Index", Estimate = erreygers_estimate, stringsAsFactors = FALSE))
          correction_factors <- c(correction_factors, Erreygers_Corrected_RCI = k_erreygers)
          correction_metadata <- rbind(correction_metadata, data.frame(
            Metric = "Erreygers_Corrected_RCI", Correction_Method = "erreygers", Transformation_Factor = k_erreygers,
            Formula = "k_E * C", Requires_Bounded_Outcome = TRUE, Computed = TRUE, Failure_Reason = NA_character_, stringsAsFactors = FALSE
          ))
          computed_methods <- c(computed_methods, "erreygers")
        }
      }
    }
  }

  corrections_out <- list(
    metrics = metrics,
    correction_factors = correction_factors,
    correction_metadata = correction_metadata,
    warnings = unique(warnings_out),
    diagnostics = list(bounded_outcome = bounded_outcome, correction_methods_computed = computed_methods)
  )

  # -------------------------------------------------------------------------
  # 3. Confidence Intervals
  # -------------------------------------------------------------------------
  ci_method <- config$ci_method
  conf_level <- config$conf_level
  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha / 2)

  make_results <- function(standard_errors, ci_lower, ci_upper, ci_method_label, se_method_label, interval_type, n_replicates = NA_integer_) {
    out <- metrics
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

  ci_diag <- list(ci_method = ci_method, z_value = z_value, health_se_supplied = !is.null(config$health_se_var))
  metric_names <- metrics$Metric
  point_estimates <- stats::setNames(metrics$Estimate, metric_names)

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

    standard_errors <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    ci_lower <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    ci_upper <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)

    for (metric in metric_names) {
      factor_value <- as.numeric(correction_factors[[metric]])
      if (!is.finite(factor_value)) next
      standard_errors[[metric]] <- abs(factor_value) * se_rci
      ci_lower[[metric]] <- point_estimates[[metric]] - z_value * standard_errors[[metric]]
      ci_upper[[metric]] <- point_estimates[[metric]] + z_value * standard_errors[[metric]]
    }

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
    keep_cols <- intersect(c("Replicate_ID", "Replicate_Type", metric_names), names(results_replicates))
    results_replicates <- results_replicates[, keep_cols, drop = FALSE]

    standard_errors <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    ci_lower <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)
    ci_upper <- stats::setNames(rep(NA_real_, length(metric_names)), metric_names)

    valid_count <- resample_out$diagnostics$n_replicates_valid

    if (ci_method == "bootstrap") {
      se_method <- "bootstrap_ecological_units"
      interval_type <- "bootstrap_percentile"
      for (metric in metric_names) {
        if (!metric %in% names(results_replicates)) next
        vals <- results_replicates[[metric]]
        vals <- vals[is.finite(vals)]
        if (length(vals) < 2L) next
        standard_errors[[metric]] <- stats::sd(vals)
        q <- stats::quantile(vals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE, type = 7)
        ci_lower[[metric]] <- q[[1L]]
        ci_upper[[metric]] <- q[[2L]]
      }
    } else {
      se_method <- "jackknife_ecological_units"
      interval_type <- "jackknife_Wald_normal"
      for (metric in metric_names) {
        if (!metric %in% names(results_replicates)) next
        vals <- results_replicates[[metric]]
        vals <- vals[is.finite(vals)]
        n_jack <- length(vals)
        if (n_jack < 2L) next
        theta_bar <- mean(vals)
        standard_errors[[metric]] <- sqrt((n_jack - 1) / n_jack * sum((vals - theta_bar)^2))
        ci_lower[[metric]] <- point_estimates[[metric]] - z_value * standard_errors[[metric]]
        ci_upper[[metric]] <- point_estimates[[metric]] + z_value * standard_errors[[metric]]
      }
    }

    results_overall <- make_results(standard_errors, ci_lower, ci_upper, ci_method, se_method, interval_type, valid_count)
    ci_diag$se_method <- se_method
    ci_diag$interval_type <- interval_type
    ci_diag$n_replicates_valid <- valid_count
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
