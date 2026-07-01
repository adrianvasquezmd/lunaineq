#' Algorithmic Stratification and Ridit Construction for SII/RII
#'
#' @description
#' `sii_rii_prepare_data` transforms an unordered ecological cross-section into a
#' rigorously ordered and population-weighted fractional rank space (ridits).
#'
#' The function resolves socioeconomic grouping directionality, collapses redundant
#' strata based on continuous vs discrete structures, and anchors the final coordinate
#' vector. The resultant design matrix provides the mandatory theoretical structure for
#' GLM-based estimation at exact lower and upper distributional bounds.
#'
#' @importFrom stats aggregate
#' @noRd
.sii_rii_prepare_data <- function(data, validation) {

  # ---------------------------------------------------------------------------
  # 1. Validate the validation object contract.
  # ---------------------------------------------------------------------------

  required_validation_fields <- c(
    "measure",
    "health_indicator_type",
    "input_source",
    "variable_names",
    "higher_ineq_is_favorable",
    "rate_scaling_factor",
    "candidate_models",
    "model_selection_metric",
    "social_position_method",
    "ci_method",
    "conf_level",
    "estimand"
  )

  missing_validation_fields <- setdiff(required_validation_fields, names(validation))
  if (length(missing_validation_fields) > 0L) {
    stop(
      "`validation` is missing required fields: ",
      paste(missing_validation_fields, collapse = ", "),
      ". Action: call `.sii_rii_validate_inputs()` before `.sii_rii_prepare_data()`.",
      call. = FALSE
    )
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or tibble.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`data` has zero rows.", call. = FALSE)
  }

  measure <- validation$measure
  health_indicator_type <- validation$health_indicator_type
  input_source <- validation$input_source
  vars <- validation$variable_names
  higher_ineq_is_favorable <- validation$higher_ineq_is_favorable
  rate_scaling_factor <- validation$rate_scaling_factor
  social_position_method <- validation$social_position_method

  if (!(measure %in% c("sii", "rii"))) {
    stop("Internal error: unsupported `measure` in `validation`. Expected `sii` or `rii`.", call. = FALSE)
  }
  if (!(health_indicator_type %in% c("proportion", "percentage", "rate", "ratio"))) {
    stop("Internal error: unsupported `health_indicator_type` in `validation`.", call. = FALSE)
  }
  if (!(input_source %in% c("counts", "indicator_population"))) {
    stop("Internal error: unsupported `input_source` in `validation`.", call. = FALSE)
  }
  if (!is.logical(higher_ineq_is_favorable) || length(higher_ineq_is_favorable) != 1L ||
      is.na(higher_ineq_is_favorable)) {
    stop("Internal error: `higher_ineq_is_favorable` must be a single TRUE/FALSE value.",
         call. = FALSE)
  }
  if (!is.numeric(rate_scaling_factor) || length(rate_scaling_factor) != 1L ||
      is.na(rate_scaling_factor) || !is.finite(rate_scaling_factor) || rate_scaling_factor <= 0) {
    stop("Internal error: `rate_scaling_factor` must be a single numeric value > 0 after validation.",
         call. = FALSE)
  }
  if (!(social_position_method %in% c("classic", "bounded"))) {
    stop("Internal error: unsupported `social_position_method` in `validation`.", call. = FALSE)
  }

  required_vars <- c(vars$analysis_unit_var, vars$equity_stratifier_var)
  if (input_source == "counts") {
    required_vars <- c(required_vars, vars$health_numerator_var, vars$health_denominator_var)
  } else {
    required_vars <- c(required_vars, vars$health_indicator_var, vars$population_weights_var)
  }

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      "The following variables required by the validation object were not found in `data`: ",
      paste(missing_vars, collapse = ", "),
      ". Action: use the same `data` object passed to `.sii_rii_validate_inputs()`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Extract common columns and create orientation of the equity stratifier.
  # ---------------------------------------------------------------------------

  n <- nrow(data)

  unit_id <- data[[vars$analysis_unit_var]]
  equity_stratifier <- as.numeric(data[[vars$equity_stratifier_var]])

  # `equity_order` is deliberately oriented from disadvantaged to advantaged.
  # If higher values of the stratifier are favorable, low values should come
  # first. If higher values are unfavorable, high original values should come
  # first, achieved by multiplying by -1.
  equity_order <- if (isTRUE(higher_ineq_is_favorable)) {
    equity_stratifier
  } else {
    -equity_stratifier
  }

  # ---------------------------------------------------------------------------
  # 3. Build outcome, weights, and count/exposure variables by input source.
  # ---------------------------------------------------------------------------

  numerator <- rep(NA_real_, n)
  denominator <- rep(NA_real_, n)
  population_weight <- rep(NA_real_, n)
  health_indicator_internal <- rep(NA_real_, n)
  health_indicator_output <- rep(NA_real_, n)
  model_response <- rep(NA_real_, n)
  count_failures <- rep(NA_real_, n)
  model_offset <- rep(NA_real_, n)

  if (input_source == "counts") {
    numerator <- as.numeric(data[[vars$health_numerator_var]])
    denominator <- as.numeric(data[[vars$health_denominator_var]])
    population_weight <- denominator

    if (health_indicator_type %in% c("proportion", "percentage")) {
      health_indicator_internal <- numerator / denominator
      health_indicator_output <- health_indicator_internal
      if (health_indicator_type == "percentage") {
        health_indicator_output <- 100 * health_indicator_internal
      }
      model_response <- health_indicator_internal
      count_failures <- denominator - numerator
      model_offset <- rep(NA_real_, n)
    }

    if (health_indicator_type %in% c("rate", "ratio")) {
      health_indicator_internal <- (numerator / denominator) * rate_scaling_factor
      health_indicator_output <- health_indicator_internal
      model_response <- numerator
      count_failures <- rep(NA_real_, n)
      model_offset <- log(denominator / rate_scaling_factor)
    }
  }

  if (input_source == "indicator_population") {
    indicator <- as.numeric(data[[vars$health_indicator_var]])
    population_weight <- as.numeric(data[[vars$population_weights_var]])

    if (health_indicator_type == "proportion") {
      health_indicator_internal <- indicator
      health_indicator_output <- indicator
    }

    if (health_indicator_type == "percentage") {
      health_indicator_internal <- indicator / 100
      health_indicator_output <- indicator
    }

    if (health_indicator_type %in% c("rate", "ratio")) {
      health_indicator_internal <- indicator
      health_indicator_output <- indicator
    }

    model_response <- health_indicator_internal
    count_failures <- rep(NA_real_, n)
    model_offset <- rep(NA_real_, n)
  }

  # ---------------------------------------------------------------------------
  # 4. Defensive checks after standardization.
  # ---------------------------------------------------------------------------

  if (any(!is.finite(equity_order))) {
    stop("Internal preparation error: non-finite `equity_order` values were produced.",
         call. = FALSE)
  }
  if (any(!is.finite(population_weight)) || any(population_weight <= 0)) {
    stop("Internal preparation error: `population_weight` must be finite and > 0.",
         call. = FALSE)
  }
  if (any(!is.finite(health_indicator_internal))) {
    stop("Internal preparation error: non-finite `health_indicator_internal` values were produced.",
         call. = FALSE)
  }
  if (any(!is.finite(health_indicator_output))) {
    stop("Internal preparation error: non-finite `health_indicator_output` values were produced.",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 5. Compute ridit/social position with tied stratifier blocks.
  # ---------------------------------------------------------------------------

  total_weight <- sum(population_weight)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("Internal preparation error: total population/exposure weight must be > 0.",
         call. = FALSE)
  }

  block_df <- data.frame(
    equity_order = equity_order,
    population_weight = population_weight,
    stringsAsFactors = FALSE
  )

  block_summary <- stats::aggregate(
    population_weight ~ equity_order,
    data = block_df,
    FUN = sum
  )

  block_summary <- block_summary[order(block_summary$equity_order), , drop = FALSE]
  block_summary$cum_weight_before <- c(0, head(cumsum(block_summary$population_weight), -1L))
  block_summary$cum_weight_after <- cumsum(block_summary$population_weight)
  block_summary$social_position <- (
    block_summary$cum_weight_before + 0.5 * block_summary$population_weight
  ) / total_weight
  block_summary$population_share <- block_summary$population_weight / total_weight

  social_position <- block_summary$social_position[
    match(equity_order, block_summary$equity_order)
  ]

  population_share <- population_weight / total_weight

  # Unit-level cumulative columns are useful for auditing, but the ridit itself
  # is block-based when ties exist. This avoids arbitrary within-tie ordering.
  sort_key <- order(equity_order, as.character(unit_id), seq_len(n))
  cumulative_weight_before_unit <- rep(NA_real_, n)
  cumulative_weight_after_unit <- rep(NA_real_, n)
  cumulative_weight_before_unit[sort_key] <- c(0, head(cumsum(population_weight[sort_key]), -1L))
  cumulative_weight_after_unit[sort_key] <- cumsum(population_weight[sort_key])

  # ---------------------------------------------------------------------------
  # 6. Add transformed variables used by candidate models and diagnostics.
  # ---------------------------------------------------------------------------

  log_response <- rep(NA_real_, n)

  if (all(health_indicator_output > 0)) {
    log_response <- log(health_indicator_output)
  }

  # ---------------------------------------------------------------------------
  # 7. Assemble standardized prepared data.
  # ---------------------------------------------------------------------------

  prepared_data <- data.frame(
    original_row = seq_len(n),
    unit_id = as.character(unit_id),
    equity_stratifier = equity_stratifier,
    equity_order = equity_order,
    population_weight = population_weight,
    population_share = population_share,
    cumulative_weight_before_unit = cumulative_weight_before_unit,
    cumulative_weight_after_unit = cumulative_weight_after_unit,
    social_position = social_position,
    health_indicator_internal = health_indicator_internal,
    health_indicator_output = health_indicator_output,
    model_response = model_response,
    log_response = log_response,
    numerator = numerator,
    denominator = denominator,
    count_failures = count_failures,
    model_offset = model_offset,
    stringsAsFactors = FALSE
  )

  prepared_data <- prepared_data[order(prepared_data$social_position,
                                       prepared_data$equity_order,
                                       prepared_data$unit_id,
                                       prepared_data$original_row), , drop = FALSE]
  rownames(prepared_data) <- NULL

  # ---------------------------------------------------------------------------
  # 8. Define prediction positions for the SII/RII estimand.
  # ---------------------------------------------------------------------------

  if (social_position_method == "classic") {
    disadvantaged_position <- 0
    advantaged_position <- 1
  } else {
    disadvantaged_position <- min(prepared_data$social_position)
    advantaged_position <- max(prepared_data$social_position)
  }

  if (!is.finite(disadvantaged_position) || !is.finite(advantaged_position) ||
      disadvantaged_position >= advantaged_position) {
    stop(
      "Internal preparation error: invalid SII/RII prediction positions. ",
      "The disadvantaged position must be lower than the advantaged position.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 9. Return prepared data and auditable metadata for downstream modules.
  # ---------------------------------------------------------------------------

  list(
    data = prepared_data,
    ridit_blocks = block_summary,
    metadata = list(
      health_indicator_type = health_indicator_type,
      input_source = input_source,
      measure = measure,
      estimand = validation$estimand,
      estimand_type = if (!is.null(validation$estimand_type)) validation$estimand_type else if (measure == "sii") "absolute_difference" else "relative_ratio",
      extreme_contrast_rule = if (!is.null(validation$extreme_contrast_rule)) validation$extreme_contrast_rule else if (measure == "sii") "f(advantaged_position) - f(disadvantaged_position)" else "f(advantaged_position) / f(disadvantaged_position)",
      social_position_direction = "0_most_disadvantaged_to_1_most_advantaged",
      equity_order_rule = if (isTRUE(higher_ineq_is_favorable)) {
        "equity_order = equity_stratifier; lower values are more disadvantaged"
      } else {
        "equity_order = -equity_stratifier; higher original values are more disadvantaged"
      },
      ridit_formula = "social_position = (cumulative population/exposure before tied block + 0.5 * block population/exposure) / total population/exposure",
      tie_handling = "Units with identical oriented equity stratifier values receive the same block-based social position.",
      population_weight_source = if (input_source == "counts") {
        "health_denominator_var"
      } else {
        "population_weights_var"
      },
      outcome_internal_scale = if (health_indicator_type == "percentage") {
        "proportion_0_to_1_for_modeling"
      } else if (health_indicator_type == "proportion") {
        "proportion_0_to_1"
      } else {
        "final_scaled_indicator"
      },
      outcome_output_scale = if (health_indicator_type == "percentage") {
        "percentage_0_to_100"
      } else if (health_indicator_type == "proportion") {
        "proportion_0_to_1"
      } else {
        "original_scaled_indicator"
      },
      rate_scaling_factor = rate_scaling_factor,
      social_position_method = social_position_method,
      disadvantaged_position = disadvantaged_position,
      advantaged_position = advantaged_position,
      total_population_weight = total_weight,
      n_units = nrow(prepared_data),
      n_ridit_blocks = nrow(block_summary)
    ),
    diagnostics = list(
      social_position_min = min(prepared_data$social_position),
      social_position_max = max(prepared_data$social_position),
      social_position_has_ties = any(duplicated(prepared_data$social_position)),
      outcome_min = min(prepared_data$health_indicator_output),
      outcome_max = max(prepared_data$health_indicator_output),
      outcome_has_zero = any(prepared_data$health_indicator_output == 0),
      log_response_available = all(is.finite(prepared_data$log_response)),
      rii_requires_positive_extreme_predictions = if (measure == "rii") TRUE else NA,
      rii_positivity_note = if (measure == "rii") {
        "RII requires strictly positive predicted values at disadvantaged and advantaged positions; this must be checked after model fitting/estimation."
      } else {
        NA_character_
      }
    )
  )
}
