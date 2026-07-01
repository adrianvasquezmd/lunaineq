#' Absolute concentration index for ecological data
#'
#' @description
#' Estimates the Absolute Concentration Index (ACI) using grouped ecological
#' data. Each row of `data` is treated as an ecological unit or group, such as
#' a department, region, district, province, country, quintile, or other
#' aggregated unit.
#'
#' The function first estimates the empirical grouped Relative Concentration
#' Index (RCI):
#'
#' \deqn{C = (2 / \mu) \sum_i f_i y_i R_i - 1}
#'
#' and then estimates the ACI as:
#'
#' \deqn{ACI = \mu C}
#'
#' where \eqn{\mu} is the weighted mean of the health indicator.
#'
#' The base RCI is also returned in `results_overall` for transparency.
#'
#' Smoothed concentration curves, if requested, are computed for visualization
#' only and are not used to estimate the RCI, ACI, standard errors, or 
#' confidence intervals.
#'
#' @param data A data frame containing ecological-unit data.
#' @param health_indicator_type Character. One of `"proportion"`,
#'   `"percentage"`, `"rate"`, or `"ratio"`.
#' @param health_indicator_var Character or NULL. Name of the health indicator
#'   variable if the indicator is supplied directly.
#' @param health_numerator_var Character or NULL. Name of the health numerator
#'   variable if the health indicator should be constructed from counts.
#' @param health_denominator_var Character or NULL. Name of the health
#'   denominator variable if the health indicator should be constructed from
#'   counts.
#' @param health_se_var Character or NULL. Name of the standard error variable
#'   for the ecological-unit/group health indicator. When supplied and
#'   `ci_method = "analytic"`, it is used as group-level uncertainty in the
#'   grouped analytic variance approximation.
#' @param rate_scaling_factor Numeric or NULL. Scaling factor applied when
#'   `health_indicator_type` is `"rate"` or `"ratio"` and the indicator is
#'   constructed from numerator/denominator data. For example, use 100000 for
#'   maternal mortality ratio per 100,000 live births. Directly supplied rate
#'   indicators are assumed to already be on the desired scale.
#' @param population_weights_var Character or NULL. Name of the population or
#'   exposure weight variable. Required when `health_indicator_var` is supplied
#'   directly. Optional when numerator/denominator variables are supplied; if
#'   omitted, the denominator is used as the population/exposure weight.
#' @param equity_stratifier_var Character. Name of the socioeconomic/equity
#'   stratifier variable.
#' @param analysis_unit_var Character or NULL. Name of the ecological-unit
#'   identifier variable.
#' @param higher_ineq_is_favorable Logical. If TRUE, higher values of the equity
#'   stratifier represent greater advantage. If FALSE, higher values represent
#'   greater disadvantage and are internally reversed for ordering.
#' @param health_lower_bound Numeric or NULL. Lower bound of the health indicator.
#'   Included for API compatibility with `rci_ineqeco`, but ignored for ACI.
#' @param health_upper_bound Numeric or NULL. Upper bound of the health indicator.
#'   Included for API compatibility with `rci_ineqeco`, but ignored for ACI.
#' @param bounded_corrections Character. Included for API compatibility with
#'   `rci_ineqeco`, but ignored for ACI (always "none").
#' @param ci_method Character. One of `"analytic"`, `"bootstrap"`, or
#'   `"jackknife"`.
#' @param n_boot Integer. Number of bootstrap replicates when
#'   `ci_method = "bootstrap"`.
#' @param conf_level Numeric. Confidence level. Default is 0.95.
#' @param add_smoothed_curves Logical. Whether to add smoothed concentration
#'   curves to `results_curve`.
#' @param smoothing_methods Character vector. Any combination of `"ml"` and
#'   `"loess"`. Used only when `add_smoothed_curves = TRUE`.
#' @param loess_span Numeric. Span parameter used for LOESS smoothing.
#' @param seed Integer or NULL. Optional random seed for bootstrap
#'   reproducibility.
#' @param include_replicates Logical. Whether to include bootstrap/jackknife
#'   replicate estimates in the final object.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return An object of class `"aci_ineqeco_result"` and
#'   `"aci_rci_ineqeco_result"`, containing:
#' \describe{
#'   \item{summary_stats}{Scalar summary statistics and the primary ACI result.}
#'   \item{results_overall}{Metric estimates, standard errors, and confidence
#'     intervals.}
#'   \item{results_units}{Prepared ecological-unit analytic data.}
#'   \item{results_curve}{Empirical and optional smoothed concentration-curve
#'     data.}
#'   \item{results_replicates}{Bootstrap/jackknife replicate estimates, or
#'     NULL.}
#'   \item{methodology}{Methodological metadata.}
#'   \item{diagnostics}{Nested diagnostics from each workflow step.}
#'   \item{warnings}{Unique warnings generated during the workflow.}
#'   \item{call}{Original matched call.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Use the built-in paho_data
#' data(paho_data)
#' df_2020 <- subset(paho_data, year == 2020)
#' 
#' # Calculate Absolute Concentration Index using analytical variance
#' aci <- aci_ineqeco(
#'   data = df_2020,
#'   health_indicator_type = "rate",
#'   health_numerator_var = "maternal_death",
#'   health_denominator_var = "live_births",
#'   rate_scaling_factor = 100000,
#'   analysis_unit_var = "state",
#'   equity_stratifier_var = "ubn",
#'   ci_method = "analytic"
#' )
#' print(aci)
#' }
aci_ineqeco <- function(
  data,
  health_indicator_type = c("proportion", "percentage", "rate", "ratio"),
  health_indicator_var = NULL,
  health_numerator_var = NULL,
  health_denominator_var = NULL,
  health_se_var = NULL,
  rate_scaling_factor = NULL,
  population_weights_var = NULL,
  equity_stratifier_var,
  analysis_unit_var = NULL,
  higher_ineq_is_favorable = TRUE,
  health_lower_bound = NULL,
  health_upper_bound = NULL,
  bounded_corrections = "none",
  ci_method = c("analytic", "bootstrap", "jackknife"),
  n_boot = 1000L,
  conf_level = 0.95,
  add_smoothed_curves = TRUE,
  smoothing_methods = c("ml", "loess"),
  loess_span = 0.75,
  seed = NULL,
  include_replicates = TRUE,
  verbose = FALSE
) {
  call <- match.call()

  config <- .aci_rci_validate_inputs(
    data = data,
    measure = "aci",
    health_indicator_type = health_indicator_type,
    health_indicator_var = health_indicator_var,
    health_numerator_var = health_numerator_var,
    health_denominator_var = health_denominator_var,
    health_se_var = health_se_var,
    rate_scaling_factor = rate_scaling_factor,
    population_weights_var = population_weights_var,
    equity_stratifier_var = equity_stratifier_var,
    analysis_unit_var = analysis_unit_var,
    higher_ineq_is_favorable = higher_ineq_is_favorable,
    health_lower_bound = health_lower_bound,
    health_upper_bound = health_upper_bound,
    bounded_corrections = bounded_corrections,
    ci_method = ci_method,
    n_boot = n_boot,
    conf_level = conf_level,
    add_smoothed_curves = add_smoothed_curves,
    smoothing_methods = smoothing_methods,
    verbose = verbose
  )

  prepared <- .aci_rci_prepare_data(
    data = data,
    config = config
  )

  estimation_results <- .aci_estimate(
    prepared = prepared,
    config = config,
    seed = seed
  )

  smooth_curve <- .aci_rci_smooth_curve(
    prepared = prepared,
    config = config,
    loess_span = loess_span
  )

  .aci_rci_format_outputs(
    prepared = prepared,
    estimate = estimation_results$estimate,
    corrections = estimation_results$corrections,
    ci = estimation_results$ci,
    smooth_curve = smooth_curve,
    config = config,
    call = call,
    result_type = "aci",
    include_replicates = include_replicates
  )
}
