#' Structural Compilation Engine for ACI/RCI Results
#'
#' @description
#' `aci_rci_format_outputs` rigorously structures the empirical metrics, variance arrays,
#' topological datasets, and geometric curves derived during the ACI/RCI algorithmic lifecycle.
#'
#' It maps internal operational metadata onto public-facing methodological typologies,
#' consolidates analytical flags, and returns a unified composite object formatted
#' specifically for high-level ecological summaries. This function applies structural
#' coercion without mutating foundational point estimates or their probabilistic bounds.
#'
#'
#' Numeric outputs are intentionally not rounded. Users who want to display more
#' digits in the R console can use `options(digits = 17)`.
#'
#' @param prepared A list returned by `.aci_rci_prepare_data()`.
#' @param estimate A list returned by `.aci_rci_estimate_index()`.
#' @param corrections A list returned by `.aci_rci_apply_corrections()`.
#' @param ci A list returned by `.aci_rci_compute_ci()`.
#' @param smooth_curve A list returned by `.aci_rci_smooth_curve()`.
#' @param config A validated configuration list returned by
#'   `.aci_rci_validate_inputs()`.
#' @param call The matched call from the public wrapper, usually
#'   `match.call()`.
#' @param result_type Character. Either `"rci"` or `"aci"`. Determines the
#'   primary result class and primary metric metadata.
#' @param include_replicates Logical. Whether to include bootstrap/jackknife
#'   replicate estimates in the returned object when available.
#'
#' @return A list with class `"rci_ineqeco_result"` or `"aci_ineqeco_result"`
#'   plus shared class `"aci_rci_ineqeco_result"`. The list contains:
#' \describe{
#'   \item{summary_stats}{Scalar summary statistics.}
#'   \item{results_overall}{Metric estimates, standard errors, and confidence
#'     intervals.}
#'   \item{results_units}{Prepared ecological-unit analytic data.}
#'   \item{results_curve}{Empirical and optional smoothed concentration-curve
#'     data.}
#'   \item{results_replicates}{Bootstrap/jackknife replicate estimates, or
#'     NULL.}
#'   \item{methodology}{Methodological metadata.}
#'   \item{diagnostics}{Nested diagnostic information.}
#'   \item{warnings}{Unique warnings generated across the family workflow.}
#'   \item{call}{Original matched call.}
#' }
#'
#' @keywords internal
#' @noRd
#'
.aci_rci_format_outputs <- function(
  prepared,
  estimate,
  corrections,
  ci,
  smooth_curve,
  config,
  call = NULL,
  result_type = c("rci", "aci"),
  include_replicates = TRUE
) {
  warnings_out <- character()
  diagnostics <- list()

  add_warning <- function(message) {
    if (isTRUE(config$verbose)) warning(message, call. = FALSE, immediate. = TRUE)
    warnings_out <<- c(warnings_out, message)
    invisible(NULL)
  }

  add_message <- function(message) {
    if (isTRUE(config$verbose)) message(message)
    invisible(NULL)
  }

  stop_format <- function(x) {
    stop(x, call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Basic checks
  # -------------------------------------------------------------------------

  result_type <- match.arg(result_type)

  if (!is.list(prepared)) {
    stop_format("`prepared` must be the list returned by `.aci_rci_prepare_data()`.")
  }

  if (!is.list(estimate)) {
    stop_format("`estimate` must be the list returned by `.aci_rci_estimate_index()`.")
  }

  if (!is.list(corrections)) {
    stop_format("`corrections` must be the list returned by `.aci_rci_apply_corrections()`.")
  }

  if (!is.list(ci)) {
    stop_format("`ci` must be the list returned by `.aci_rci_compute_ci()`.")
  }

  if (!is.list(smooth_curve)) {
    stop_format("`smooth_curve` must be the list returned by `.aci_rci_smooth_curve()`.")
  }

  if (!is.list(config)) {
    stop_format("`config` must be the list returned by `.aci_rci_validate_inputs()`.")
  }

  if (!is.logical(include_replicates) || length(include_replicates) != 1L || is.na(include_replicates)) {
    stop_format("`include_replicates` must be a single TRUE/FALSE value.")
  }

  if (is.null(prepared$data_units) || !is.data.frame(prepared$data_units)) {
    stop_format("`prepared$data_units` must be a data frame.")
  }

  if (is.null(smooth_curve$curve_data) || !is.data.frame(smooth_curve$curve_data)) {
    stop_format("`smooth_curve$curve_data` must be a data frame.")
  }

  if (is.null(ci$results_overall) || !is.data.frame(ci$results_overall)) {
    stop_format("`ci$results_overall` must be a data frame.")
  }

  # -------------------------------------------------------------------------
  # Results tables
  # -------------------------------------------------------------------------

  results_overall <- ci$results_overall

  required_overall_cols <- c(
    "Metric",
    "Metric_Label",
    "Estimate",
    "Standard_Error",
    "CI_Lower",
    "CI_Upper",
    "CI_Method",
    "SE_Method",
    "Interval_Type",
    "Confidence_Level"
  )

  missing_overall_cols <- setdiff(required_overall_cols, names(results_overall))

  if (length(missing_overall_cols) > 0L) {
    stop_format(
      paste0(
        "`ci$results_overall` is missing required column(s): ",
        paste(missing_overall_cols, collapse = ", "),
        "."
      )
    )
  }

  # Add a primary metric flag for easier downstream printing/use.
  primary_metric <- if (result_type == "rci") {
    "RCI"
  } else {
    "ACI"
  }

  results_overall$Primary_Metric <- results_overall$Metric == primary_metric

  if (!any(results_overall$Primary_Metric)) {
    add_warning(
      paste0(
        "The expected primary metric `",
        primary_metric,
        "` was not found in `results_overall`."
      )
    )
  }

  results_units <- prepared$data_units
  results_curve <- smooth_curve$curve_data

  # Keep replicates optional because they can be large.
  results_replicates <- NULL

  if (include_replicates) {
    results_replicates <- ci$results_replicates
  } else if (!is.null(ci$results_replicates)) {
    add_warning(
      "`include_replicates = FALSE`; bootstrap/jackknife replicate estimates were not included in the final object."
    )
  }

  # -------------------------------------------------------------------------
  # Summary statistics
  # -------------------------------------------------------------------------

  mean_health <- if (!is.null(estimate$mean_health)) {
    estimate$mean_health
  } else if (!is.null(prepared$summary$mean_health)) {
    prepared$summary$mean_health
  } else {
    NA_real_
  }

  total_population_weight <- if (!is.null(estimate$total_population_weight)) {
    estimate$total_population_weight
  } else if (!is.null(prepared$summary$total_population_weight)) {
    prepared$summary$total_population_weight
  } else {
    NA_real_
  }

  weighted_health_total <- if (!is.null(estimate$weighted_health_total)) {
    estimate$weighted_health_total
  } else if (!is.null(prepared$summary$weighted_health_total)) {
    prepared$summary$weighted_health_total
  } else {
    NA_real_
  }

  primary_estimate <- if (primary_metric %in% results_overall$Metric) {
    results_overall$Estimate[match(primary_metric, results_overall$Metric)]
  } else {
    NA_real_
  }

  primary_standard_error <- if (primary_metric %in% results_overall$Metric) {
    results_overall$Standard_Error[match(primary_metric, results_overall$Metric)]
  } else {
    NA_real_
  }

  primary_ci_lower <- if (primary_metric %in% results_overall$Metric) {
    results_overall$CI_Lower[match(primary_metric, results_overall$Metric)]
  } else {
    NA_real_
  }

  primary_ci_upper <- if (primary_metric %in% results_overall$Metric) {
    results_overall$CI_Upper[match(primary_metric, results_overall$Metric)]
  } else {
    NA_real_
  }

  summary_stats <- data.frame(
    Result_Type = result_type,
    Primary_Metric = primary_metric,
    Primary_Estimate = primary_estimate,
    Primary_Standard_Error = primary_standard_error,
    Primary_CI_Lower = primary_ci_lower,
    Primary_CI_Upper = primary_ci_upper,
    Confidence_Level = config$conf_level,
    CI_Method = config$ci_method,
    N_Units = nrow(results_units),
    Mean_Health = mean_health,
    Total_Population_Weight = total_population_weight,
    Weighted_Health_Total = weighted_health_total,
    Bounded_Outcome = isTRUE(config$bounded_outcome),
    Health_Lower_Bound = if (!is.null(config$health_lower_bound)) config$health_lower_bound else NA_real_,
    Health_Upper_Bound = if (!is.null(config$health_upper_bound)) config$health_upper_bound else NA_real_,
    Health_SE_Supplied = !is.null(config$health_se_var),
    Rate_Scaling_Factor = if (!is.null(config$rate_scaling_factor)) config$rate_scaling_factor else NA_real_,
    Smoothed_Curves_Added = isTRUE(config$add_smoothed_curves),
    stringsAsFactors = FALSE
  )

  # -------------------------------------------------------------------------
  # Methodology metadata
  # -------------------------------------------------------------------------

  correction_methods <- if (!is.null(config$correction_methods) &&
                            length(config$correction_methods) > 0L) {
    paste(config$correction_methods, collapse = ";")
  } else {
    NA_character_
  }

  smoothing_methods_requested <- if (!is.null(config$smoothing_methods) &&
                                     length(config$smoothing_methods) > 0L) {
    paste(config$smoothing_methods, collapse = ";")
  } else {
    NA_character_
  }

  smoothing_methods_computed <- if (!is.null(smooth_curve$diagnostics$smoothing_methods_computed) &&
                                    length(smooth_curve$diagnostics$smoothing_methods_computed) > 0L) {
    paste(smooth_curve$diagnostics$smoothing_methods_computed, collapse = ";")
  } else {
    NA_character_
  }

  methodology <- data.frame(
    Result_Type = result_type,
    Index_Estimation_Method = "grouped_empirical_concentration_index",
    RCI_Formula = "(2 / mu) * sum(f_i * y_i * R_i) - 1",
    ACI_Formula = "mu * RCI",
    Socioeconomic_Order = "most_disadvantaged_to_most_advantaged",
    Equity_Order_Convention = if (isTRUE(config$higher_ineq_is_favorable)) {
      "higher equity stratifier values represent greater advantage"
    } else {
      "higher equity stratifier values represent greater disadvantage; values were internally reversed for ordering"
    },
    Input_Source = config$input_source,
    Health_Indicator_Type = config$health_indicator_type,
    Rate_Scaling_Factor = if (!is.null(config$rate_scaling_factor)) config$rate_scaling_factor else NA_real_,
    Health_SE_Interpretation = if (!is.null(config$health_se_var)) {
      "standard error of the health indicator for each ecological unit/group"
    } else {
      NA_character_
    },
    Bounded_Outcome = isTRUE(config$bounded_outcome),
    Bounded_Corrections = correction_methods,
    CI_Method = config$ci_method,
    SE_Method = if (!is.null(ci$diagnostics$se_method)) ci$diagnostics$se_method else NA_character_,
    Interval_Type = if (!is.null(ci$diagnostics$interval_type)) ci$diagnostics$interval_type else NA_character_,
    Confidence_Level = config$conf_level,
    Smoothed_Curves_Requested = isTRUE(config$add_smoothed_curves),
    Smoothing_Methods_Requested = smoothing_methods_requested,
    Smoothing_Methods_Computed = smoothing_methods_computed,
    Smoothing_Used_For_Estimation = FALSE,
    Rounding_Applied = FALSE,
    stringsAsFactors = FALSE
  )

  # -------------------------------------------------------------------------
  # Harmonize warnings
  # -------------------------------------------------------------------------

  collect_warnings <- function(x) {
    if (is.list(x) && !is.null(x$warnings)) {
      return(as.character(x$warnings))
    }
    character()
  }

  all_warnings <- unique(c(
    collect_warnings(config),
    collect_warnings(prepared),
    collect_warnings(estimate),
    collect_warnings(corrections),
    collect_warnings(ci),
    collect_warnings(smooth_curve),
    warnings_out
  ))

  if (length(all_warnings) == 0L) {
    all_warnings <- character()
  }

  add_message(
    "Output numeric values are not rounded. Use `options(digits = 17)` in R to display more digits in the console."
  )

  all_warnings <- unique(c(all_warnings, warnings_out))

  # -------------------------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------------------------

  diagnostics$result_type <- result_type
  diagnostics$primary_metric <- primary_metric
  diagnostics$primary_metric_found <- any(results_overall$Primary_Metric)
  diagnostics$n_units <- nrow(results_units)
  diagnostics$n_curve_points <- nrow(results_curve)
  diagnostics$n_overall_metrics <- nrow(results_overall)
  diagnostics$include_replicates <- include_replicates
  diagnostics$n_replicates <- if (!is.null(results_replicates) && is.data.frame(results_replicates)) {
    nrow(results_replicates)
  } else {
    0L
  }
  diagnostics$workflow <- list(
    validation = if (!is.null(config$diagnostics)) config$diagnostics else NULL,
    preparation = if (!is.null(prepared$diagnostics)) prepared$diagnostics else NULL,
    estimation = if (!is.null(estimate$diagnostics)) estimate$diagnostics else NULL,
    corrections = if (!is.null(corrections$diagnostics)) corrections$diagnostics else NULL,
    confidence_intervals = if (!is.null(ci$diagnostics)) ci$diagnostics else NULL,
    smoothing = if (!is.null(smooth_curve$diagnostics)) smooth_curve$diagnostics else NULL
  )

  diagnostics$warnings_count <- length(all_warnings)
  diagnostics$rounding_applied <- FALSE
  diagnostics$smoothing_used_for_estimation <- FALSE

  # -------------------------------------------------------------------------
  # Final object
  # -------------------------------------------------------------------------

  result <- list(
    summary_stats = summary_stats,
    results_overall = results_overall,
    results_units = results_units,
    results_curve = results_curve,
    results_replicates = results_replicates,
    methodology = methodology,
    diagnostics = diagnostics,
    warnings = all_warnings,
    call = call
  )

  class(result) <- c(
    paste0(result_type, "_ineqeco_result"),
    "aci_rci_ineqeco_result",
    "ineqeco_result",
    "list"
  )

  result
}
