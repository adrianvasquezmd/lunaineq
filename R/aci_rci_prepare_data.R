#' Pre-processing Engine for ACI/RCI Ecological Data
#'
#' @description
#' `aci_rci_prepare_data` operates the preliminary analytical transformations required for
#' concentration index estimation. It standardizes matrix inputs, enforces directional
#' ordering based on socioeconomic stratifiers, computes ridits (fractional ranks),
#' and aggregates cumulative population and health shares.
#'
#' This algorithmic step converts primary observational cross-sections into an
#' ordered state, strictly prerequisite for the Kakwani covariance formula.
#'
#' @importFrom stats complete.cases
#' @importFrom utils head tail
#' @noRd

#' @param data A data frame containing ecological-unit data.
#' @param config A validated configuration list returned by
#'   `.aci_rci_validate_inputs()`.
#'
#' @return A list with:
#' \describe{
#'   \item{data_units}{Prepared unit-level analytic data.}
#'   \item{summary}{Core scalar summaries used by downstream functions.}
#'   \item{warnings}{Warnings generated during preparation.}
#'   \item{diagnostics}{Preparation diagnostics.}
#' }
#'
#' @keywords internal
#' @noRd
.aci_rci_prepare_data <- function(
  data,
  config
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

  stop_prepare <- function(x) {
    stop(x, call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Basic checks
  # -------------------------------------------------------------------------

  if (!is.data.frame(data)) {
    stop_prepare("`data` must be a data frame.")
  }

  if (!is.list(config)) {
    stop_prepare("`config` must be the list returned by `.aci_rci_validate_inputs()`.")
  }

  required_config <- c(
    "health_indicator_type",
    "input_source",
    "health_indicator_var",
    "health_numerator_var",
    "health_denominator_var",
    "health_se_var",
    "population_weights_var",
    "equity_stratifier_var",
    "analysis_unit_var",
    "higher_ineq_is_favorable",
    "bounded_outcome",
    "health_lower_bound",
    "health_upper_bound"
  )

  missing_config <- setdiff(required_config, names(config))

  if (length(missing_config) > 0L) {
    stop_prepare(
      paste0(
        "`config` is missing required element(s): ",
        paste(missing_config, collapse = ", "),
        "."
      )
    )
  }

  # -------------------------------------------------------------------------
  # Extract config
  # -------------------------------------------------------------------------

  health_indicator_type <- config$health_indicator_type
  input_source <- config$input_source

  health_indicator_var <- config$health_indicator_var
  health_numerator_var <- config$health_numerator_var
  health_denominator_var <- config$health_denominator_var
  health_se_var <- config$health_se_var
  rate_scaling_factor <- if (!is.null(config$rate_scaling_factor)) config$rate_scaling_factor else 1
  population_weights_var <- config$population_weights_var
  equity_stratifier_var <- config$equity_stratifier_var
  analysis_unit_var <- config$analysis_unit_var

  higher_ineq_is_favorable <- config$higher_ineq_is_favorable

  bounded_outcome <- isTRUE(config$bounded_outcome)
  health_lower_bound <- config$health_lower_bound
  health_upper_bound <- config$health_upper_bound

  # -------------------------------------------------------------------------
  # Build standardized data frame
  # -------------------------------------------------------------------------

  n_input <- nrow(data)

  out <- data.frame(
    .row_id = seq_len(n_input),
    stringsAsFactors = FALSE
  )

  if (!is.null(analysis_unit_var)) {
    out$.unit <- data[[analysis_unit_var]]
  } else {
    out$.unit <- out$.row_id
    add_message(
      "`analysis_unit_var` was not supplied. Row numbers will be used as ecological-unit identifiers."
    )
  }

  out$.equity <- as.numeric(data[[equity_stratifier_var]])

  if (input_source == "counts") {
    numerator <- as.numeric(data[[health_numerator_var]])
    denominator <- as.numeric(data[[health_denominator_var]])

    if (health_indicator_type %in% c("proportion", "percentage")) {
      health <- numerator / denominator
    } else if (health_indicator_type %in% c("rate", "ratio")) {
      health <- (numerator / denominator) * rate_scaling_factor
    } else {
      stop_prepare(
        paste0(
          "Unsupported `health_indicator_type`: ",
          health_indicator_type,
          "."
        )
      )
    }

    population_weight <- if (!is.null(population_weights_var)) {
      as.numeric(data[[population_weights_var]])
    } else {
      add_message(
        paste0(
          "`population_weights_var` was not supplied, so `health_denominator_var` (",
          config$health_denominator_var,
          ") is being used as the population/exposure weight."
        )
      )
      denominator
    }

    out$.health_numerator <- numerator
    out$.health_denominator <- denominator
  } else if (input_source == "indicator_population") {
    health <- as.numeric(data[[health_indicator_var]])
    population_weight <- as.numeric(data[[population_weights_var]])
  } else {
    stop_prepare(
      paste0(
        "Unsupported `input_source`: ",
        input_source,
        ". Expected 'counts' or 'indicator_population'."
      )
    )
  }

  health_se <- if (!is.null(health_se_var)) {
    as.numeric(data[[health_se_var]])
  } else {
    rep(NA_real_, n_input)
  }

  # -------------------------------------------------------------------------
  # Scale health indicator and SE consistently
  # -------------------------------------------------------------------------

  scale_factor <- 1

  if (health_indicator_type == "percentage") {
    # Numerator/denominator percentage inputs are first constructed as
    # numerator / denominator, i.e. already as proportions. Directly supplied
    # percentage indicators are converted from percent units to proportions.
    if (input_source == "indicator_population") {
      health <- health / 100
      scale_factor <- 100
    }

    if (!is.null(health_se_var)) {
      health_se <- health_se / 100
    }

    add_message(
      "Percentage inputs are analyzed internally as proportions. Directly supplied percentages were divided by 100; numerator/denominator percentages were already proportions. `health_se_var`, if supplied, was divided by 100 and must originally be in percentage-point units."
    )
  }

  if (health_indicator_type %in% c("rate", "ratio") && input_source == "counts") {
    add_message(
      paste0(
        "Rate/ratio indicator was constructed as (numerator / denominator) * rate_scaling_factor. ",
        "Current rate_scaling_factor = ", rate_scaling_factor, ". ",
        "`health_se_var`, if supplied, is assumed to already be on this same scaled metric."
      )
    )
  }

  if (health_indicator_type %in% c("rate", "ratio") && input_source == "indicator_population") {
    add_warning(
      "Directly supplied rate/ratio indicator is assumed to already be on the desired scale; `rate_scaling_factor` was not reapplied."
    )
  }

  out$.health <- health
  out$.health_se <- health_se
  out$.population_weight <- population_weight

  # -------------------------------------------------------------------------
  # Row exclusion
  # -------------------------------------------------------------------------

  required_cols <- c(".health", ".population_weight", ".equity")

  if (!is.null(health_se_var)) {
    required_cols <- c(required_cols, ".health_se")
  }

  complete_required <- stats::complete.cases(out[, required_cols, drop = FALSE])
  positive_weight <- !is.na(out$.population_weight) & out$.population_weight > 0
  finite_health <- is.finite(out$.health)
  finite_equity <- is.finite(out$.equity)

  keep <- complete_required & positive_weight & finite_health & finite_equity

  n_removed_missing_or_invalid <- sum(!keep)

  if (n_removed_missing_or_invalid > 0L) {
    add_warning(
      paste0(
        n_removed_missing_or_invalid,
        " row(s) were removed because of missing, non-finite, or invalid values in the health indicator, population weight, equity stratifier, or health SE."
      )
    )
  }

  out <- out[keep, , drop = FALSE]

  if (nrow(out) < 2L) {
    stop_prepare(
      "Fewer than two ecological units remain after data preparation. ACI/RCI cannot be computed."
    )
  }

  # -------------------------------------------------------------------------
  # Validate prepared values
  # -------------------------------------------------------------------------

  if (any(out$.population_weight <= 0, na.rm = TRUE)) {
    stop_prepare("Prepared population weights must be strictly positive.")
  }

  if (any(out$.health < 0, na.rm = TRUE)) {
    add_warning(
      "Negative health indicator values remain after preparation. Concentration-index results may be difficult to interpret."
    )
  }

  if (!is.null(health_se_var)) {
    if (any(out$.health_se < 0, na.rm = TRUE)) {
      stop_prepare("Prepared health standard errors must be non-negative.")
    }

    if (all(out$.health_se == 0, na.rm = TRUE)) {
      add_warning(
        "`health_se_var` was supplied, but all prepared standard errors are zero."
      )
    }
  }

  if (bounded_outcome) {
    if (is.null(health_lower_bound) || is.null(health_upper_bound)) {
      stop_prepare(
        "A bounded outcome was requested, but valid lower and upper bounds were not available in `config`."
      )
    }

    tolerance <- sqrt(.Machine$double.eps)

    below_bound <- out$.health < (health_lower_bound - tolerance)
    above_bound <- out$.health > (health_upper_bound + tolerance)

    if (any(below_bound | above_bound, na.rm = TRUE)) {
      add_warning(
        paste0(
          "Some prepared health indicator values fall outside the declared bounds [",
          health_lower_bound, ", ", health_upper_bound,
          "]. Bounded corrections may be invalid unless the bounds or indicator scale are corrected."
        )
      )
    }
  }

  # -------------------------------------------------------------------------
  # Orient equity order
  # -------------------------------------------------------------------------

  # Internal convention:
  #   lower .equity_order = more disadvantaged
  #   higher .equity_order = more advantaged
  if (higher_ineq_is_favorable) {
    out$.equity_order <- out$.equity
  } else {
    out$.equity_order <- -out$.equity
  }

  # Stable ordering. If ties are present, original row order is preserved.
  order_index <- order(out$.equity_order, out$.row_id, na.last = NA)
  out <- out[order_index, , drop = FALSE]
  rownames(out) <- NULL

  tied_equity <- any(duplicated(out$.equity_order))

  if (tied_equity) {
    add_warning(
      "Ties were detected in the oriented equity stratifier. Tied units were assigned the same population-weighted midpoint fractional rank, harmonized with the original luna implementation."
    )
  }

  # -------------------------------------------------------------------------
  # Compute shares, fractional ranks, and cumulative quantities
  # -------------------------------------------------------------------------

  total_weight <- sum(out$.population_weight)

  if (!is.finite(total_weight) || total_weight <= 0) {
    stop_prepare("Total population/exposure weight must be positive and finite.")
  }

  out$.population_share <- out$.population_weight / total_weight
  out$.cum_population <- cumsum(out$.population_share)
  out$.cum_population_lag <- c(0, utils::head(out$.cum_population, -1L))

  # Fractional rank with tie handling harmonized with the original luna
  # implementation: units with the same oriented equity value receive the same
  # midpoint rank over the combined population weight of the tied block.
  out$.fractional_rank <- NA_real_
  tie_groups <- split(seq_len(nrow(out)), out$.equity_order)

  for (idx in tie_groups) {
    tied_pop_start <- min(out$.cum_population_lag[idx])
    tied_pop_end <- max(out$.cum_population[idx])
    out$.fractional_rank[idx] <- tied_pop_start + 0.5 * (tied_pop_end - tied_pop_start)
  }

  weighted_health_total <- sum(out$.population_weight * out$.health)
  weighted_mean_health <- weighted_health_total / total_weight

  if (!is.finite(weighted_mean_health)) {
    stop_prepare("The weighted mean health indicator is not finite.")
  }

  if (weighted_mean_health <= 0) {
    stop_prepare(
      "The weighted mean health indicator is zero or negative. ACI/RCI cannot be computed because the concentration index divides by the mean."
    )
  }

  out$.weighted_health <- out$.population_weight * out$.health
  out$.health_share <- out$.weighted_health / weighted_health_total
  out$.cum_health <- cumsum(out$.health_share)
  out$.cum_health_lag <- c(0, utils::head(out$.cum_health, -1L))

  # Add endpoints useful for concentration-curve plotting/smoothing.
  curve_data <- data.frame(
    .curve_point = seq_len(nrow(out) + 1L),
    .cum_population = c(0, out$.cum_population),
    .cum_health_empirical = c(0, out$.cum_health),
    stringsAsFactors = FALSE
  )

  # Numerical cleanup only for accumulated floating point noise; no rounding of
  # substantive values is performed.
  if (abs(utils::tail(out$.cum_population, 1L) - 1) > 1e-10) {
    add_warning(
      "Final cumulative population share differs from 1 beyond tolerance. Check population weights."
    )
  }

  if (abs(utils::tail(out$.cum_health, 1L) - 1) > 1e-10) {
    add_warning(
      "Final cumulative health share differs from 1 beyond tolerance. Check health indicator and weights."
    )
  }

  # -------------------------------------------------------------------------
  # Diagnostics and summary
  # -------------------------------------------------------------------------

  diagnostics$n_rows_input <- n_input
  diagnostics$n_rows_prepared <- nrow(out)
  diagnostics$n_rows_removed <- n_removed_missing_or_invalid
  diagnostics$tied_equity_order <- tied_equity
  diagnostics$total_population_weight <- total_weight
  diagnostics$weighted_health_total <- weighted_health_total
  diagnostics$weighted_mean_health <- weighted_mean_health
  diagnostics$health_indicator_internal_scale <- if (health_indicator_type == "percentage") {
    "proportion"
  } else {
    health_indicator_type
  }
  diagnostics$percentage_scale_factor <- scale_factor
  diagnostics$rate_scaling_factor <- rate_scaling_factor
  diagnostics$bounded_outcome <- bounded_outcome
  diagnostics$health_lower_bound <- health_lower_bound
  diagnostics$health_upper_bound <- health_upper_bound
  diagnostics$health_se_supplied <- !is.null(health_se_var)

  summary <- list(
    n_units = nrow(out),
    total_population_weight = total_weight,
    weighted_health_total = weighted_health_total,
    mean_health = weighted_mean_health,
    health_indicator_internal_scale = diagnostics$health_indicator_internal_scale,
    rate_scaling_factor = rate_scaling_factor,
    health_lower_bound = health_lower_bound,
    health_upper_bound = health_upper_bound,
    bounded_outcome = bounded_outcome
  )

  # Keep output as base data.frame and numeric values unrounded.
  list(
    data_units = out,
    curve_data = curve_data,
    summary = summary,
    warnings = unique(warnings_out),
    diagnostics = diagnostics
  )
}
