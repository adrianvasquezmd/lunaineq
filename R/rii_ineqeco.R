#' Relative Index of Inequality for Ecological Data
#'
#' @description
#' Estimates the relative index of inequality for ecological data from ecological data using a
#' population-weighted ridit/social-position variable and model-based prediction
#' at the disadvantaged and advantaged social-position extremes.
#'
#' @details
#' The estimand is `RII = Y_advantaged / Y_disadvantaged`. For binomial and quasibinomial models, RII is computed as a ratio of predicted probabilities on the response scale, not as an odds ratio.
#'
#' The function is intentionally a user-facing orchestrator. Shared internal
#' helpers with the `.sii_rii_*` prefix validate inputs, prepare data, fit and
#' select models, resample when requested, and format auditable outputs.
#'
#' @param data A data frame with one row per ecological unit.
#' @param health_indicator_type Indicator type: `"proportion"`, `"percentage"`,
#'   `"rate"`, or `"ratio"`.
#' @param health_indicator_var Optional column name with a precomputed indicator.
#' @param population_weights_var Optional column name with population/exposure
#'   weights when using a precomputed indicator.
#' @param health_numerator_var Optional column name with numerator counts.
#' @param health_denominator_var Optional column name with denominator/exposure
#'   counts.
#' @param equity_stratifier_var Column name of the equity stratifier used to
#'   order ecological units.
#' @param analysis_unit_var Column name uniquely identifying ecological units.
#' @param higher_ineq_is_favorable Logical. `TRUE` if higher values of the equity
#'   stratifier indicate greater advantage; `FALSE` otherwise.
#' @param rate_scaling_factor Scaling factor for rates or ratios from counts.
#' @param models Optional character vector of candidate models. If `NULL`, all
#'   allowed models for the indicator type and input source are considered.
#' @param model_selection_metric Metric used to select the model.
#' @param social_position_method `"classic"` evaluates extremes at 0 and 1;
#'   `"bounded"` evaluates at the observed minimum and maximum social position.
#' @param ci_method Confidence interval method: `"wald"`, `"bootstrap"`, or
#'   `"jackknife"`.
#' @param n_boot Number of bootstrap replicates when `ci_method = "bootstrap"`.
#' @param conf_level Confidence level.
#' @param seed Optional random seed for resampling.
#' @param verbose Logical; passed to internal validation metadata.
#' @param vcov_type Robust covariance type passed to `sandwich::vcovHC()` when
#'   available for Wald intervals.
#' @param ... Reserved. Unsupported arguments are rejected with an actionable
#'   message.
#'
#' @return An auditable list with overall estimates, model-selection results,
#'   unit-level predictions, ridit construction, optional resampling replicates,
#'   methodology metadata, and diagnostics.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Use the built-in paho_data
#' data(paho_data)
#' df_2020 <- subset(paho_data, year == 2020)
#' 
#' # Calculate Relative Index of Inequality
#' rii <- rii_ineqeco(
#'   data = df_2020,
#'   health_indicator_type = "rate",
#'   health_numerator_var = "maternal_death",
#'   health_denominator_var = "live_births",
#'   rate_scaling_factor = 100000,
#'   analysis_unit_var = "state",
#'   equity_stratifier_var = "ubn",
#'   higher_ineq_is_favorable = FALSE,
#'   ci_method = "wald"
#' )
#' print(rii)
#' }
rii_ineqeco <- function(
  data,
  health_indicator_type = c("proportion", "percentage", "rate", "ratio"),
  health_indicator_var = NULL,
  population_weights_var = NULL,
  health_numerator_var = NULL,
  health_denominator_var = NULL,
  equity_stratifier_var,
  analysis_unit_var,
  higher_ineq_is_favorable,
  rate_scaling_factor = NULL,
  models = NULL,
  model_selection_metric = c("MSE", "MAE", "RMSE", "MAPE", "AIC", "BIC"),
  social_position_method = c("classic", "bounded"),
  ci_method = c("wald", "bootstrap", "jackknife"),
  n_boot = 2000,
  conf_level = 0.95,
  seed = NULL,
  verbose = TRUE,
  vcov_type = "HC1",
  ...
) {
  extra_args <- list(...)
  if (length(extra_args) > 0) {
    unsupported <- names(extra_args)
    if (is.null(unsupported)) {
      unsupported <- rep("", length(extra_args))
    }

    if ("health_se_var" %in% unsupported) {
      stop(
        "`health_se_var` is not part of this version. ",
        "External standard errors are not incorporated into SII/RII estimation ",
        "or uncertainty intervals. ",
        "Action: omit `health_se_var` and use `ci_method = 'wald'`, ",
        "`ci_method = 'bootstrap'`, or `ci_method = 'jackknife'`.",
        call. = FALSE
      )
    }

    stop(
      "Unsupported argument(s) supplied to `rii_ineqeco()`: ",
      paste(unsupported, collapse = ", "), ". ",
      "Action: remove these argument(s).",
      call. = FALSE
    )
  }


  validation <- .sii_rii_validate_inputs(
    data = data,
    measure = "rii",
    health_indicator_type = health_indicator_type,
    health_indicator_var = health_indicator_var,
    population_weights_var = population_weights_var,
    health_numerator_var = health_numerator_var,
    health_denominator_var = health_denominator_var,
    equity_stratifier_var = equity_stratifier_var,
    analysis_unit_var = analysis_unit_var,
    higher_ineq_is_favorable = higher_ineq_is_favorable,
    rate_scaling_factor = rate_scaling_factor,
    models = models,
    model_selection_metric = model_selection_metric,
    social_position_method = social_position_method,
    ci_method = ci_method,
    n_boot = n_boot,
    conf_level = conf_level,
    seed = seed,
    verbose = verbose
  )

  prepared <- .sii_rii_prepare_data(
    data = data,
    validation = validation
  )

  fitted <- .sii_rii_fit_models(
    prepared = prepared,
    validation = validation
  )

  selected <- .sii_rii_select_model(
    fitted = fitted,
    validation = validation
  )

  estimate <- .sii_rii_estimate_rii(
    selected = selected,
    validation = validation,
    conf_level = validation$conf_level,
    vcov_type = vcov_type
  )

  resample <- NULL
  if (validation$ci_method %in% c("bootstrap", "jackknife")) {
    resample <- .sii_rii_resample(
      selected = selected,
      validation = validation,
      point_estimate = estimate$results_overall$Estimate[1L],
      ci_method = validation$ci_method,
      n_boot = validation$n_boot,
      conf_level = validation$conf_level,
      seed = validation$seed
    )
  }

  out <- .sii_rii_format_outputs(
    validation = validation,
    prepared = prepared,
    fitted = fitted,
    selected = selected,
    estimate = estimate,
    resample = resample
  )

  out$call <- match.call()
  out
}
