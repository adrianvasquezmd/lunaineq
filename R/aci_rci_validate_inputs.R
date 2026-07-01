#' Input Validation and Standardization Engine for ACI/RCI Functions
#'
#' @description
#' `aci_rci_validate_inputs` enforces strict epidemiological and statistical typing constraints
#' on all user-supplied data prior to concentration index estimation.
#'
#' It performs systematic checks on input dimensions, numeric bounds, logical flags, and missing data.
#' It identifies the source of the health indicator, specifies and restricts applicable bounded-outcome
#' corrections, validates resampling topologies, and compiles a comprehensive execution configuration.
#'
#' @importFrom stats na.omit complete.cases
#' @noRd

#' @param data A data frame containing ecological-unit data.
#' @param measure Character. Either `"rci"` or `"aci"`. Used only for metadata.
#' @param health_indicator_type Character. One of `"proportion"`,
#'   `"percentage"`, `"rate"`, or `"ratio"`.
#' @param health_indicator_var Character or NULL. Name of the health indicator
#'   variable if the indicator is supplied directly.
#' @param health_numerator_var Character or NULL. Name of the health numerator
#'   variable if the health indicator is constructed from counts.
#' @param health_denominator_var Character or NULL. Name of the health
#'   denominator variable if the health indicator is constructed from counts.
#' @param health_se_var Character or NULL. Name of the ecological-unit/group
#'   standard error variable for the health indicator.
#' @param rate_scaling_factor Numeric or NULL. Scaling factor applied when
#'   `health_indicator_type` is `"rate"` or `"ratio"` and the indicator is
#'   constructed from numerator/denominator data. For example, use 100000 for
#'   maternal mortality ratio per 100,000 live births. Directly supplied rate
#'   indicators are assumed to already be on the desired scale.
#' @param population_weights_var Character or NULL. Name of the population
#'   weight variable. Required when `health_indicator_var` is used directly.
#'   Optional when numerator/denominator are supplied; in that case the
#'   denominator is used as the default population weight.
#' @param equity_stratifier_var Character. Name of the socioeconomic/equity
#'   stratifier variable.
#' @param analysis_unit_var Character or NULL. Name of the ecological-unit
#'   identifier variable.
#' @param higher_ineq_is_favorable Logical. If TRUE, higher values of the
#'   equity stratifier represent greater advantage. If FALSE, higher values
#'   represent greater disadvantage.
#' @param health_lower_bound Numeric or NULL. Lower bound for bounded outcomes.
#' @param health_upper_bound Numeric or NULL. Upper bound for bounded outcomes.
#' @param bounded_corrections Character. One of `"auto"`, `"none"`,
#'   `"wagstaff"`, `"erreygers"`, or `"both"`.
#' @param ci_method Character. One of `"analytic"`, `"bootstrap"`, or
#'   `"jackknife"`.
#' @param n_boot Integer. Number of bootstrap replicates when
#'   `ci_method = "bootstrap"`.
#' @param conf_level Numeric. Confidence level between 0 and 1.
#' @param add_smoothed_curves Logical. Whether to add smoothed concentration
#'   curves to the curve output.
#' @param smoothing_methods Character vector. Any of `"ml"` and `"loess"`.
#'   Used only when `add_smoothed_curves = TRUE`.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list with validated and standardized inputs, including
#'   `input_source`, `bounded_outcome`, `correction_methods`, `warnings`,
#'   and `diagnostics`.
#'
#' @keywords internal
#' @noRd
.aci_rci_validate_inputs <- function(
  data,
  measure = c("rci", "aci"),
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
  bounded_corrections = c("auto", "none", "wagstaff", "erreygers", "both"),
  ci_method = c("analytic", "bootstrap", "jackknife"),
  n_boot = 1000L,
  conf_level = 0.95,
  add_smoothed_curves = TRUE,
  smoothing_methods = c("ml", "loess"),
  verbose = FALSE
) {
  warnings_out <- character()
  diagnostics <- list()

  add_warning <- function(x) {
    warnings_out <<- c(warnings_out, x)
    if (isTRUE(verbose)) warning(x, call. = FALSE, immediate. = TRUE)
    invisible(NULL)
  }

  add_message <- function(x) {
    if (isTRUE(verbose)) message(x)
    invisible(NULL)
  }

  stop_input <- function(x) {
    stop(x, call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Basic object validation
  # -------------------------------------------------------------------------

  if (!is.data.frame(data)) {
    stop_input("`data` must be a data frame.")
  }

  if (nrow(data) == 0L) {
    stop_input("`data` has zero rows. At least two ecological units/groups are required.")
  }

  if (nrow(data) < 2L) {
    stop_input("At least two ecological units/groups are required to compute ACI/RCI.")
  }

  measure <- match.arg(measure)
  health_indicator_type <- match.arg(health_indicator_type)
  bounded_corrections <- match.arg(bounded_corrections)
  ci_method <- match.arg(ci_method)

  if (!is.logical(higher_ineq_is_favorable) ||
      length(higher_ineq_is_favorable) != 1L ||
      is.na(higher_ineq_is_favorable)) {
    stop_input("`higher_ineq_is_favorable` must be a single TRUE/FALSE value.")
  }

  if (!is.logical(add_smoothed_curves) ||
      length(add_smoothed_curves) != 1L ||
      is.na(add_smoothed_curves)) {
    stop_input("`add_smoothed_curves` must be a single TRUE/FALSE value.")
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop_input("`verbose` must be a single TRUE/FALSE value.")
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level)) {
    stop_input("`conf_level` must be a single numeric value between 0 and 1.")
  }

  if (conf_level <= 0 || conf_level >= 1) {
    stop_input("`conf_level` must be greater than 0 and less than 1.")
  }

  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot)) {
    stop_input("`n_boot` must be a single numeric value.")
  }

  n_boot <- as.integer(n_boot)

  if (n_boot < 1L) {
    stop_input("`n_boot` must be at least 1.")
  }

  if (ci_method == "bootstrap" && n_boot < 100L) {
    add_warning(
      "`n_boot` is less than 100. Bootstrap confidence intervals may be unstable."
    )
  }

  # -------------------------------------------------------------------------
  # Rate/ratio scaling validation
  # -------------------------------------------------------------------------

  if (health_indicator_type %in% c("proportion", "percentage")) {
    if (!is.null(rate_scaling_factor)) {
      add_warning(
        "`rate_scaling_factor` was supplied but will be ignored for proportion/percentage indicators."
      )
    }
    rate_scaling_factor <- 1
  } else {
    if (is.null(rate_scaling_factor)) {
      stop_input(
        paste0(
          "`rate_scaling_factor` is required when `health_indicator_type = \"",
          health_indicator_type,
          "\"`. Use 1 only if the rate/ratio should remain unscaled."
        )
      )
    }

    if (!is.numeric(rate_scaling_factor) || length(rate_scaling_factor) != 1L ||
        is.na(rate_scaling_factor) || rate_scaling_factor <= 0) {
      stop_input("`rate_scaling_factor` must be a single positive numeric value for rate/ratio indicators.")
    }
  }

  # -------------------------------------------------------------------------
  # Variable name helpers
  # -------------------------------------------------------------------------

  is_null_or_string <- function(x) {
    is.null(x) || (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x))
  }

  require_null_or_string <- function(x, arg_name) {
    if (!is_null_or_string(x)) {
      stop_input(paste0("`", arg_name, "` must be NULL or a single non-empty character string."))
    }
  }

  require_var_exists <- function(var, arg_name, required = TRUE) {
    if (is.null(var)) {
      if (required) {
        stop_input(paste0("`", arg_name, "` is required."))
      }
      return(invisible(FALSE))
    }

    if (!var %in% names(data)) {
      stop_input(
        paste0(
          "Variable `", var, "` supplied in `", arg_name,
          "` was not found in `data`."
        )
      )
    }

    invisible(TRUE)
  }

  require_null_or_string(health_indicator_var, "health_indicator_var")
  require_null_or_string(health_numerator_var, "health_numerator_var")
  require_null_or_string(health_denominator_var, "health_denominator_var")
  require_null_or_string(health_se_var, "health_se_var")
  require_null_or_string(population_weights_var, "population_weights_var")
  require_null_or_string(equity_stratifier_var, "equity_stratifier_var")
  require_null_or_string(analysis_unit_var, "analysis_unit_var")

  require_var_exists(equity_stratifier_var, "equity_stratifier_var", required = TRUE)
  require_var_exists(analysis_unit_var, "analysis_unit_var", required = FALSE)
  require_var_exists(health_indicator_var, "health_indicator_var", required = FALSE)
  require_var_exists(health_numerator_var, "health_numerator_var", required = FALSE)
  require_var_exists(health_denominator_var, "health_denominator_var", required = FALSE)
  require_var_exists(health_se_var, "health_se_var", required = FALSE)
  require_var_exists(population_weights_var, "population_weights_var", required = FALSE)

  # -------------------------------------------------------------------------
  # Determine input source
  # -------------------------------------------------------------------------

  has_indicator <- !is.null(health_indicator_var)
  has_numerator <- !is.null(health_numerator_var)
  has_denominator <- !is.null(health_denominator_var)
  has_counts <- has_numerator && has_denominator

  if (has_numerator && !has_denominator) {
    stop_input(
      "`health_denominator_var` is required when `health_numerator_var` is supplied."
    )
  }

  if (!has_numerator && has_denominator) {
    stop_input(
      "`health_numerator_var` is required when `health_denominator_var` is supplied."
    )
  }

  if (!has_indicator && !has_counts) {
    stop_input(
      paste0(
        "You must supply either `health_indicator_var` with `population_weights_var`, ",
        "or both `health_numerator_var` and `health_denominator_var`."
      )
    )
  }

  if (has_indicator && is.null(population_weights_var) && !has_counts) {
    stop_input(
      "`population_weights_var` is required when `health_indicator_var` is supplied directly."
    )
  }

  if (has_indicator && has_counts) {
    add_warning(
      paste0(
        "Both `health_indicator_var` and numerator/denominator variables were supplied. ",
        "The numerator/denominator source will be used to construct the health indicator; ",
        "`health_indicator_var` will be ignored for estimation."
      )
    )
  }

  input_source <- if (has_counts) {
    "counts"
  } else {
    "indicator_population"
  }

  diagnostics$input_source <- input_source

  if (input_source == "counts" && health_indicator_type %in% c("rate", "ratio")) {
    add_message(
      paste0(
        "For numerator/denominator rate or ratio inputs, the health indicator will be calculated as ",
        "(numerator / denominator) * rate_scaling_factor. Current rate_scaling_factor = ",
        rate_scaling_factor, ". `health_se_var`, if supplied, must be on this same scaled metric."
      )
    )
  }

  if (input_source == "indicator_population" && health_indicator_type %in% c("rate", "ratio")) {
    add_message(
      paste0(
        "For directly supplied rate or ratio indicators, `health_indicator_var` is assumed to already be on the desired scale. ",
        "`rate_scaling_factor` is kept as metadata and is not reapplied to `health_indicator_var`."
      )
    )
  }

  if (input_source == "counts" && is.null(population_weights_var)) {
    add_message(
      paste0(
        "`population_weights_var` was not supplied. ",
        "`health_denominator_var` will be used as the population/exposure weight."
      )
    )
  }

  # -------------------------------------------------------------------------
  # Numeric variable validation
  # -------------------------------------------------------------------------

  numeric_vars <- unique(stats::na.omit(c(
    health_indicator_var,
    health_numerator_var,
    health_denominator_var,
    health_se_var,
    population_weights_var,
    equity_stratifier_var
  )))

  for (var in numeric_vars) {
    if (!is.numeric(data[[var]]) && !is.integer(data[[var]])) {
      stop_input(
        paste0(
          "Variable `", var, "` must be numeric. ",
          "Please convert it before calling the function."
        )
      )
    }
  }

  if (!is.null(analysis_unit_var)) {
    if (any(is.na(data[[analysis_unit_var]]))) {
      add_warning(
        paste0(
          "Missing values were detected in `analysis_unit_var` (`",
          analysis_unit_var, "`). These rows may be removed during data preparation."
        )
      )
    }

    duplicated_units <- sum(duplicated(data[[analysis_unit_var]]), na.rm = TRUE)

    if (duplicated_units > 0L) {
      add_warning(
        paste0(
          "Duplicated ecological-unit/group identifiers were detected in `analysis_unit_var` (`",
          analysis_unit_var, "`). The function will treat each row as one ecological unit/group; ",
          "aggregate the data beforehand if duplicate rows should represent one unit/group."
        )
      )
    }
  }

  # -------------------------------------------------------------------------
  # Value checks for input variables
  # -------------------------------------------------------------------------

  if (input_source == "counts") {
    if (any(data[[health_numerator_var]] < 0, na.rm = TRUE)) {
      stop_input("`health_numerator_var` contains negative values.")
    }

    if (any(data[[health_denominator_var]] <= 0, na.rm = TRUE)) {
      stop_input("`health_denominator_var` must contain only positive values.")
    }

    if (any(data[[health_numerator_var]] > data[[health_denominator_var]], na.rm = TRUE) &&
        health_indicator_type %in% c("proportion", "percentage")) {
      stop_input(
        paste0(
          "For `health_indicator_type = \"", health_indicator_type,
          "\"`, numerator values cannot exceed denominator values."
        )
      )
    }
  }

  if (!is.null(population_weights_var)) {
    if (any(data[[population_weights_var]] <= 0, na.rm = TRUE)) {
      add_warning(
        paste0(
          "`population_weights_var` contains zero or negative values. ",
          "Rows with non-positive weights will be removed during data preparation."
        )
      )
    }
  }

  if (!is.null(health_indicator_var)) {
    if (any(data[[health_indicator_var]] < 0, na.rm = TRUE)) {
      add_warning(
        paste0(
          "`health_indicator_var` contains negative values. ",
          "This may be inappropriate for concentration-index calculations unless the indicator is truly allowed to be negative."
        )
      )
    }
  }

  if (!is.null(health_se_var)) {
    if (any(data[[health_se_var]] < 0, na.rm = TRUE)) {
      stop_input("`health_se_var` contains negative values. Standard errors must be non-negative.")
    }

    add_message(
      paste0(
        "`health_se_var` was supplied and will be interpreted as the standard error ",
        "of the health indicator for each ecological unit/group."
      )
    )
  }

  if (any(is.na(data[[equity_stratifier_var]]))) {
    add_warning(
      paste0(
        "Missing values were detected in `equity_stratifier_var` (`",
        equity_stratifier_var, "`). These rows will be removed during data preparation."
      )
    )
  }

  if (length(unique(stats::na.omit(data[[equity_stratifier_var]]))) < 2L) {
    stop_input(
      "`equity_stratifier_var` must contain at least two distinct non-missing values."
    )
  }

  if (any(duplicated(stats::na.omit(data[[equity_stratifier_var]])))) {
    add_warning(
      paste0(
        "Ties were detected in `equity_stratifier_var`. ",
        "The preparation step should use a stable tie-handling rule so tied units/groups are ordered consistently."
      )
    )
  }

  # -------------------------------------------------------------------------
  # Bounds and bounded outcome detection
  # -------------------------------------------------------------------------

  if (!is.null(health_lower_bound)) {
    if (!is.numeric(health_lower_bound) || length(health_lower_bound) != 1L || is.na(health_lower_bound)) {
      stop_input("`health_lower_bound` must be NULL or a single numeric value.")
    }
  }

  if (!is.null(health_upper_bound)) {
    if (!is.numeric(health_upper_bound) || length(health_upper_bound) != 1L || is.na(health_upper_bound)) {
      stop_input("`health_upper_bound` must be NULL or a single numeric value.")
    }
  }

  if (xor(is.null(health_lower_bound), is.null(health_upper_bound))) {
    stop_input(
      "`health_lower_bound` and `health_upper_bound` must be supplied together."
    )
  }

  if (!is.null(health_lower_bound) && !is.null(health_upper_bound)) {
    if (health_upper_bound <= health_lower_bound) {
      stop_input("`health_upper_bound` must be greater than `health_lower_bound`.")
    }
  }

  bounded_outcome <- health_indicator_type %in% c("proportion", "percentage") ||
    (!is.null(health_lower_bound) && !is.null(health_upper_bound))

  if (health_indicator_type == "proportion") {
    if (is.null(health_lower_bound)) health_lower_bound <- 0
    if (is.null(health_upper_bound)) health_upper_bound <- 1
  }

  if (health_indicator_type == "percentage") {
    if (is.null(health_lower_bound)) health_lower_bound <- 0
    if (is.null(health_upper_bound)) health_upper_bound <- 1

    add_message(
      paste0(
        "`health_indicator_type = \"percentage\"` was supplied. ",
        "Directly supplied percentage indicators and their standard errors, if supplied, will be converted ",
        "from percent units to proportions during data preparation."
      )
    )
  }

  if (bounded_outcome && is.null(health_lower_bound)) {
    stop_input(
      paste0(
        "A bounded outcome was detected, but bounds could not be determined. ",
        "Please supply `health_lower_bound` and `health_upper_bound`."
      )
    )
  }

  if (!bounded_outcome && bounded_corrections != "none") {
    if (bounded_corrections == "auto") {
      bounded_corrections_resolved <- character()
    } else {
      add_warning(
        paste0(
          "Bounded corrections were requested, but the outcome does not appear bounded. ",
          "Corrections will not be computed unless valid bounds are supplied."
        )
      )
      bounded_corrections_resolved <- character()
    }
  } else if (bounded_outcome) {
    bounded_corrections_resolved <- switch(
      bounded_corrections,
      auto = c("wagstaff", "erreygers"),
      none = character(),
      wagstaff = "wagstaff",
      erreygers = "erreygers",
      both = c("wagstaff", "erreygers")
    )
  } else {
    bounded_corrections_resolved <- character()
  }

  if (bounded_outcome && length(bounded_corrections_resolved) > 0L) {
    add_message(
      paste0(
        "Bounded-outcome corrections will be computed using lower bound = ",
        health_lower_bound, " and upper bound = ", health_upper_bound, "."
      )
    )
  }

  # -------------------------------------------------------------------------
  # CI method-specific warnings
  # -------------------------------------------------------------------------

  if (ci_method == "analytic") {
    if (is.null(health_se_var)) {
      add_message(
        paste0(
          "`ci_method = \"analytic\"` without `health_se_var` will use the grouped-data ",
          "Kakwani-Wagstaff-van Doorslaer analytic variance approximation with the ",
          "number of ecological units/groups as the effective sample size."
        )
      )
    } else {
      add_message(
        paste0(
          "`ci_method = \"analytic\"` with `health_se_var` will use the supplied ",
          "ecological-unit/group health standard errors in the analytic variance calculation."
        )
      )
    }
  }

  if (ci_method == "bootstrap") {
    add_message(
      paste0(
        "`ci_method = \"bootstrap\"` will resample ecological units/groups with replacement. ",
        "The supplied `health_se_var`, if any, is not used directly by the non-parametric bootstrap."
      )
    )
  }

  if (ci_method == "jackknife") {
    add_message(
      paste0(
        "`ci_method = \"jackknife\"` will use leave-one-ecological-unit/group-out resampling. ",
        "The supplied `health_se_var`, if any, is not used directly by the jackknife."
      )
    )

    if (nrow(data) < 3L) {
      stop_input("`ci_method = \"jackknife\"` requires at least three ecological units/groups.")
    }
  }

  # -------------------------------------------------------------------------
  # Smoothing validation
  # -------------------------------------------------------------------------

  allowed_smoothing_methods <- c("ml", "loess")

  if (add_smoothed_curves) {
    if (is.null(smoothing_methods) || length(smoothing_methods) == 0L) {
      stop_input(
        "`smoothing_methods` must contain at least one method when `add_smoothed_curves = TRUE`."
      )
    }

    if (!is.character(smoothing_methods)) {
      stop_input("`smoothing_methods` must be a character vector.")
    }

    smoothing_methods <- unique(tolower(smoothing_methods))

    invalid_smoothing <- setdiff(smoothing_methods, allowed_smoothing_methods)

    if (length(invalid_smoothing) > 0L) {
      stop_input(
        paste0(
          "Invalid smoothing method(s): ",
          paste(invalid_smoothing, collapse = ", "),
          ". Allowed methods are: ",
          paste(allowed_smoothing_methods, collapse = ", "),
          "."
        )
      )
    }

    add_message(
      paste0(
        "Smoothed concentration curves will be computed for visualization only. ",
        "They will not be used to estimate ACI/RCI or confidence intervals."
      )
    )
  } else {
    smoothing_methods <- character()

    add_message(
      paste0(
        "`add_smoothed_curves = FALSE`; only the empirical concentration curve will be returned."
      )
    )
  }

  # -------------------------------------------------------------------------
  # Missingness diagnostics
  # -------------------------------------------------------------------------

  vars_needed_for_rows <- unique(stats::na.omit(c(
    if (input_source == "counts") {
      c(health_numerator_var, health_denominator_var)
    } else {
      c(health_indicator_var, population_weights_var)
    },
    health_se_var,
    equity_stratifier_var,
    analysis_unit_var
  )))

  complete_rows <- stats::complete.cases(data[, vars_needed_for_rows, drop = FALSE])
  n_missing_rows <- sum(!complete_rows)

  if (n_missing_rows > 0L) {
    add_warning(
      paste0(
        n_missing_rows,
        " row(s) have missing values in variables required for estimation ",
        "and will be removed during data preparation."
      )
    )
  }

  diagnostics$n_rows_input <- nrow(data)
  diagnostics$n_rows_complete_for_required_vars <- sum(complete_rows)
  diagnostics$n_rows_with_missing_required_vars <- n_missing_rows
  diagnostics$bounded_outcome <- bounded_outcome
  diagnostics$correction_methods <- bounded_corrections_resolved
  diagnostics$ci_method <- ci_method
  diagnostics$add_smoothed_curves <- add_smoothed_curves
  diagnostics$smoothing_methods <- smoothing_methods
  diagnostics$health_se_supplied <- !is.null(health_se_var)
  diagnostics$rate_scaling_factor <- rate_scaling_factor

  if (sum(complete_rows) < 2L) {
    stop_input(
      "After excluding rows with missing required values, fewer than two ecological units/groups remain."
    )
  }

  # -------------------------------------------------------------------------
  # Return standardized validated configuration
  # -------------------------------------------------------------------------

  list(
    measure = measure,
    health_indicator_type = health_indicator_type,
    input_source = input_source,

    health_indicator_var = health_indicator_var,
    health_numerator_var = health_numerator_var,
    health_denominator_var = health_denominator_var,
    health_se_var = health_se_var,
    rate_scaling_factor = rate_scaling_factor,
    population_weights_var = population_weights_var,
    equity_stratifier_var = equity_stratifier_var,
    analysis_unit_var = analysis_unit_var,

    higher_ineq_is_favorable = higher_ineq_is_favorable,

    bounded_outcome = bounded_outcome,
    health_lower_bound = health_lower_bound,
    health_upper_bound = health_upper_bound,
    bounded_corrections = bounded_corrections,
    correction_methods = bounded_corrections_resolved,

    ci_method = ci_method,
    n_boot = n_boot,
    conf_level = conf_level,

    add_smoothed_curves = add_smoothed_curves,
    smoothing_methods = smoothing_methods,

    verbose = verbose,

    warnings = unique(warnings_out),
    diagnostics = diagnostics
  )
}
