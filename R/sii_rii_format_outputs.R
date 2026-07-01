#' Format outputs for SII/RII estimation
#'
#' Internal helper for `sii_ineqeco()` and `rii_ineqeco()`.
#' Structural Compilation Engine for SII/RII Results
#'
#' @description
#' `sii_rii_format_outputs` rigorously structures the empirical regression metrics,
#' predictions, and variance arrays derived during the SII/RII algorithmic lifecycle.
#'
#' It compiles validation warnings, selection diagnostics, raw unit predictions,
#' and overarching gap estimates into a strictly formatted list object, applying
#' deterministic coercion without altering derived inferential statistics.
#'
#' @importFrom stats weighted.mean
#' @noRd
.sii_rii_format_outputs <- function(validation,
                                prepared,
                                fitted,
                                selected,
                                estimate = NULL,
                                resample = NULL) {

  # ---------------------------------------------------------------------------
  # 1. Local utilities.
  # ---------------------------------------------------------------------------

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  require_fields <- function(x, fields, object_name, action) {
    if (!is.list(x)) {
      stop(
        "`", object_name, "` must be a list produced by the corresponding SII/RII module. ",
        "Action: ", action,
        call. = FALSE
      )
    }
    missing_fields <- setdiff(fields, names(x))
    if (length(missing_fields) > 0L) {
      stop(
        "Invalid `", object_name, "` object: missing field(s): ",
        paste(missing_fields, collapse = ", "), ". ",
        "Action: ", action,
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  add_missing_columns <- function(x, columns) {
    for (nm in names(columns)) {
      if (!nm %in% names(x)) {
        x[[nm]] <- columns[[nm]]
      }
    }
    x
  }

  safe_unique <- function(x) {
    if (is.null(x)) {
      return(character(0))
    }
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    unique(x)
  }

  as_character_or_na <- function(x) {
    if (is.null(x) || length(x) == 0L) {
      return(NA_character_)
    }
    x <- as.character(x[1L])
    if (length(x) == 0L || is.na(x) || !nzchar(x)) NA_character_ else x
  }

  bind_warnings <- function(...) {
    unique(unlist(lapply(list(...), safe_unique), recursive = FALSE, use.names = FALSE))
  }

  # ---------------------------------------------------------------------------
  # 2. Validate upstream contracts.
  # ---------------------------------------------------------------------------

  require_fields(
    validation,
    fields = c(
      "health_indicator_type",
      "input_source",
      "variable_names",
      "higher_ineq_is_favorable",
      "rate_scaling_factor",
      "allowed_models",
      "candidate_models",
      "model_set_source",
      "model_selection_metric",
      "social_position_method",
      "ci_method",
      "conf_level",
      "estimand",
      "diagnostics"
    ),
    object_name = "validation",
    action = "call `.sii_rii_validate_inputs()` before `.sii_rii_format_outputs()`."
  )

  measure <- validation$measure %||% "sii"
  measure <- as.character(measure[1L])
  if (!measure %in% c("sii", "rii")) {
    stop(
      "Invalid `validation$measure`. Expected 'sii' or 'rii'. ",
      "Action: rerun `.sii_rii_validate_inputs()`.",
      call. = FALSE
    )
  }

  metric_label <- if (measure == "sii") "Slope Index of Inequality" else "Relative Index of Inequality"
  result_class <- if (measure == "sii") "sii_ineqeco_result" else "rii_ineqeco_result"
  estimate_helper_name <- if (measure == "sii") ".sii_rii_estimate_sii()" else ".sii_rii_estimate_rii()"
  resample_helper_name <- ".sii_rii_resample()"
  measure_estimand <- validation$estimand
  if (is.null(measure_estimand)) {
    measure_estimand <- if (measure == "sii") {
      "SII = Y_advantaged - Y_disadvantaged"
    } else {
      "RII = Y_advantaged / Y_disadvantaged"
    }
  }
  require_fields(
    prepared,
    fields = c("data", "ridit_blocks", "metadata", "diagnostics"),
    object_name = "prepared",
    action = "call `.sii_rii_prepare_data()` before `.sii_rii_format_outputs()`."
  )

  require_fields(
    fitted,
    fields = c("models", "model_table", "predictions", "prepared_data", "metadata", "diagnostics"),
    object_name = "fitted",
    action = "call `.sii_rii_fit_models()` before `.sii_rii_format_outputs()`."
  )

  require_fields(
    selected,
    fields = c(
      "selected_model_name",
      "selected_model",
      "selected_predictions",
      "model_selection_table",
      "predictions",
      "prepared_data",
      "fitted_models",
      "metadata",
      "diagnostics"
    ),
    object_name = "selected",
    action = "call `.sii_rii_select_model()` before `.sii_rii_format_outputs()`."
  )

  ci_method <- as.character(validation$ci_method)
  if (length(ci_method) != 1L || !ci_method %in% c("wald", "bootstrap", "jackknife")) {
    stop(
      "Invalid `validation$ci_method`. Expected one of: wald, bootstrap, jackknife. ",
      "Action: rerun `.sii_rii_validate_inputs()`.",
      call. = FALSE
    )
  }

  measure_interval_default <- if (measure == "sii") {
    if (ci_method == "wald") "Model-based Wald interval from explicit SII contrast." else NA_character_
  } else {
    if (ci_method == "wald") "Model-based Wald interval for log-RII, exponentiated to the RII scale." else NA_character_
  }

  if (ci_method == "wald") {
    if (is.null(estimate)) {
      stop(
        "`estimate` is required when `ci_method = 'wald'`. ",
        "Action: call the measure-specific estimate helper before `.sii_rii_format_outputs()`.",
        call. = FALSE
      )
    }
    require_fields(
      estimate,
      fields = c(
        "results_overall",
        "extreme_predictions",
        "selected_model_name",
        "model_selection_table",
        "predictions",
        "prepared_data",
        "metadata",
        "diagnostics"
      ),
      object_name = "estimate",
      action = "call the measure-specific estimate helper before `.sii_rii_format_outputs()`."
    )
  } else {
    if (is.null(resample)) {
      stop(
        "`resample` is required when `ci_method = '", ci_method, "'`. ",
        "Action: call `.sii_rii_resample()` before `.sii_rii_format_outputs()`.",
        call. = FALSE
      )
    }
    require_fields(
      resample,
      fields = c("results_overall", "results_replicates", "metadata", "diagnostics"),
      object_name = "resample",
      action = "call `.sii_rii_resample()` before `.sii_rii_format_outputs()`."
    )
    if (!is.null(estimate)) {
      require_fields(
        estimate,
        fields = c(
          "results_overall",
          "extreme_predictions",
          "selected_model_name",
          "metadata",
          "diagnostics"
        ),
        object_name = "estimate",
        action = "pass a valid measure-specific estimate object or leave `estimate = NULL`."
      )
    }
  }

  d <- prepared$data
  if (!is.data.frame(d) || nrow(d) == 0L) {
    stop(
      "`prepared$data` must be a non-empty data.frame. ",
      "Action: inspect `.sii_rii_prepare_data()` output.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 3. Choose the final interval source.
  # ---------------------------------------------------------------------------

  interval_source <- if (ci_method == "wald") estimate else resample
  results_overall <- interval_source$results_overall

  if (!is.data.frame(results_overall) || nrow(results_overall) != 1L) {
    stop(
      "The interval source must contain a one-row `results_overall` data.frame. ",
      "Action: inspect the measure-specific estimate output or `.sii_rii_resample()` output.",
      call. = FALSE
    )
  }

  results_overall <- add_missing_columns(
    results_overall,
    list(
      Interval_Type = if (ci_method == "wald") "model_based_wald_contrast" else NA_character_,
      Model_Selection_Metric = validation$model_selection_metric,
      Model_Set_Source = validation$model_set_source,
      Candidate_Models = paste(validation$candidate_models, collapse = ", "),
      Estimand = measure_estimand
    )
  )

  # If resampling was used, endpoint predictions usually come from the original
  # selected model estimate. They are not redefined by the percentile/jackknife
  # interval itself, so add them only when a point-estimate object is available.
  if (!is.null(estimate) && is.data.frame(estimate$results_overall) && nrow(estimate$results_overall) == 1L) {
    endpoint_cols <- c(
      "Disadvantaged_Position",
      "Advantaged_Position",
      "Disadvantaged_Predicted_Value",
      "Advantaged_Predicted_Value"
    )
    for (nm in endpoint_cols) {
      if (!nm %in% names(results_overall) && nm %in% names(estimate$results_overall)) {
        results_overall[[nm]] <- estimate$results_overall[[nm]][1L]
      }
    }
  }

  # Stable column order without dropping future columns.
  overall_order <- c(
    "Metric",
    "Estimate",
    "Log_Estimate",
    "Standard_Error",
    "Standard_Error_Scale",
    "CI_Lower",
    "CI_Upper",
    "Confidence_Level",
    "CI_Method",
    "Interval_Type",
    "Interval_Scale",
    "Log_CI_Lower",
    "Log_CI_Upper",
    "Selected_Model",
    "Model_Selection_Metric",
    "Model_Set_Source",
    "Candidate_Models",
    "Health_Indicator_Type",
    "Input_Source",
    "Social_Position_Method",
    "Disadvantaged_Position",
    "Advantaged_Position",
    "Disadvantaged_Predicted_Value",
    "Advantaged_Predicted_Value",
    "Estimand"
  )
  existing_order <- intersect(overall_order, names(results_overall))
  results_overall <- results_overall[, c(existing_order, setdiff(names(results_overall), existing_order)), drop = FALSE]

  # ---------------------------------------------------------------------------
  # 4. Build model and unit-level outputs.
  # ---------------------------------------------------------------------------

  results_models <- selected$model_selection_table
  if (!is.data.frame(results_models)) {
    stop(
      "`selected$model_selection_table` must be a data.frame. ",
      "Action: inspect `.sii_rii_select_model()` output.",
      call. = FALSE
    )
  }
  rownames(results_models) <- NULL

  selected_model <- as.character(selected$selected_model_name)
  selected_pred_col <- paste0("Pred_", selected_model)
  pred_table <- selected$predictions
  if (!is.data.frame(pred_table)) {
    stop(
      "`selected$predictions` must be a data.frame. ",
      "Action: inspect `.sii_rii_select_model()` output.",
      call. = FALSE
    )
  }
  if (!selected_pred_col %in% names(pred_table)) {
    stop(
      "Selected prediction column `", selected_pred_col, "` was not found in `selected$predictions`. ",
      "Action: inspect `.sii_rii_fit_models()` and `.sii_rii_select_model()` outputs.",
      call. = FALSE
    )
  }

  # Merge by row order because all upstream modules preserve the prepared-data
  # ordering. Avoid keyed joins to prevent accidental duplicate-row expansion.
  if (nrow(pred_table) != nrow(d)) {
    stop(
      "`selected$predictions` and `prepared$data` have different numbers of rows. ",
      "Action: inspect upstream modules; row order must be preserved across fitting and selection.",
      call. = FALSE
    )
  }

  results_units <- data.frame(
    Unit_ID = d$unit_id,
    Original_Row = d$original_row,
    Equity_Stratifier = d$equity_stratifier,
    Equity_Order = d$equity_order,
    Population_Weight = d$population_weight,
    Population_Share = d$population_share,
    Social_Position = d$social_position,
    Health_Indicator_Internal = d$health_indicator_internal,
    Health_Indicator_Output = d$health_indicator_output,
    Numerator = d$numerator,
    Denominator = d$denominator,
    Count_Failures = d$count_failures,
    Model_Offset = d$model_offset,
    Selected_Model = selected_model,
    Predicted_Value = pred_table[[selected_pred_col]],
    Residual = pred_table$Observed - pred_table[[selected_pred_col]],
    stringsAsFactors = FALSE
  )

  prediction_cols <- grep("^Pred_", names(pred_table), value = TRUE)
  for (nm in prediction_cols) {
    if (!nm %in% names(results_units)) {
      results_units[[nm]] <- pred_table[[nm]]
    }
  }

  unit_order <- c(
    "Unit_ID",
    "Original_Row",
    "Equity_Stratifier",
    "Equity_Order",
    "Population_Weight",
    "Population_Share",
    "Social_Position",
    "Health_Indicator_Internal",
    "Health_Indicator_Output",
    "Numerator",
    "Denominator",
    "Count_Failures",
    "Model_Offset",
    "Selected_Model",
    "Predicted_Value",
    "Residual"
  )
  results_units <- results_units[, c(unit_order, setdiff(names(results_units), unit_order)), drop = FALSE]
  rownames(results_units) <- NULL

  results_ridit <- prepared$ridit_blocks
  if (!is.data.frame(results_ridit)) {
    results_ridit <- data.frame()
  } else {
    rownames(results_ridit) <- NULL
  }

  results_replicates <- NULL
  if (!is.null(resample) && is.data.frame(resample$results_replicates)) {
    results_replicates <- resample$results_replicates
    rownames(results_replicates) <- NULL
  }

  results_extremes <- NULL
  if (!is.null(estimate) && is.data.frame(estimate$extreme_predictions)) {
    results_extremes <- estimate$extreme_predictions
    rownames(results_extremes) <- NULL
  }

  # ---------------------------------------------------------------------------
  # 5. Summary statistics.
  # ---------------------------------------------------------------------------

  summary_stats <- data.frame(
    Metric = c(
      "n_units",
      "n_ridit_blocks",
      "total_population_weight",
      "social_position_min",
      "social_position_max",
      "health_indicator_output_min",
      "health_indicator_output_max",
      "health_indicator_output_mean_weighted",
      "health_indicator_output_mean_unweighted",
      "n_candidate_models",
      "n_viable_models",
      "selected_model",
      "ci_method"
    ),
    Value = c(
      nrow(d),
      nrow(results_ridit),
      sum(d$population_weight),
      min(d$social_position),
      max(d$social_position),
      min(d$health_indicator_output),
      max(d$health_indicator_output),
      stats::weighted.mean(d$health_indicator_output, w = d$population_weight),
      mean(d$health_indicator_output),
      length(validation$candidate_models),
      sum(results_models$Fitted & results_models$Converged & results_models$Prediction_Valid, na.rm = TRUE),
      selected_model,
      ci_method
    ),
    stringsAsFactors = FALSE
  )

  # Keep `Value` as character if mixed numeric/character. This is intentional:
  # summary_stats is a compact audit summary, while numeric results remain in
  # typed technical tables such as results_overall and results_units.
  summary_stats$Value <- as.character(summary_stats$Value)

  # ---------------------------------------------------------------------------
  # 6. Methodology table.
  # ---------------------------------------------------------------------------

  methodology <- data.frame(
    Item = c(
      "Metric",
      "Estimand",
      "Social_Position_Direction",
      "Ridit_Formula",
      "Tie_Handling",
      "Input_Source",
      "Health_Indicator_Type",
      "Outcome_Output_Scale",
      "Population_Weight_Source",
      "Model_Set_Source",
      "Candidate_Models",
      "Model_Selection_Metric",
      "Model_Selection_Rule",
      "Prediction_Comparison_Rule",
      "AIC_BIC_Rule",
      "Selected_Model",
      "CI_Method",
      "Interval_Rule",
      "External_SE_Rule",
      "Rounding_Rule"
    ),
    Description = c(
      metric_label,
      measure_estimand,
      prepared$metadata$social_position_direction %||% "0_most_disadvantaged_to_1_most_advantaged",
      prepared$metadata$ridit_formula %||% NA_character_,
      prepared$metadata$tie_handling %||% NA_character_,
      validation$input_source,
      validation$health_indicator_type,
      prepared$metadata$outcome_output_scale %||% validation$outcome_scale %||% NA_character_,
      prepared$metadata$population_weight_source %||% NA_character_,
      validation$model_set_source,
      paste(validation$candidate_models, collapse = ", "),
      validation$model_selection_metric,
      selected$metadata$selection_rule %||% NA_character_,
      selected$metadata$predictive_metric_rule %||% NA_character_,
      selected$metadata$aic_bic_rule %||% validation$aic_bic_comparability_note %||% NA_character_,
      selected_model,
      ci_method,
      interval_source$metadata$interval_rule %||% measure_interval_default,
      validation$diagnostics$external_standard_errors %||% "External standard errors are not part of this version.",
      "No internal rounding is applied to technical numeric outputs."
    ),
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # 7. Diagnostics.
  # ---------------------------------------------------------------------------

  diagnostic_warnings <- bind_warnings(
    validation$diagnostics$warnings,
    prepared$diagnostics$warnings,
    fitted$diagnostics$warnings,
    fitted$diagnostics$warnings_by_model,
    selected$diagnostics$warnings,
    if (!is.null(estimate)) estimate$diagnostics$estimate_warnings else NULL,
    if (!is.null(resample)) resample$diagnostics$warnings else NULL
  )

  diagnostics <- list(
    measure = measure,
    warnings = diagnostic_warnings,
    ignored_arguments = validation$diagnostics$ignored_arguments %||% character(0),
    external_standard_errors = validation$diagnostics$external_standard_errors %||% "Not part of this version.",
    n_units = nrow(d),
    n_ridit_blocks = nrow(results_ridit),
    n_candidate_models = length(validation$candidate_models),
    n_fitted_models = fitted$metadata$n_fitted_models %||% NA_integer_,
    n_viable_models = fitted$metadata$n_viable_models %||% NA_integer_,
    n_eligible_models = selected$diagnostics$n_eligible_models %||% NA_integer_,
    selected_model = selected_model,
    selected_model_rank = selected$diagnostics$selected_model_rank %||% NA_integer_,
    ci_method = ci_method,
    interval_source = if (ci_method == "wald") paste0(measure, "_estimate") else "sii_rii_resample",
    extreme_check = {
      default_extreme_check <- if (measure == "sii") {
        "Estimate equals Advantaged_Predicted_Value - Disadvantaged_Predicted_Value"
      } else {
        "Estimate equals Advantaged_Predicted_Value / Disadvantaged_Predicted_Value"
      }
      if (!is.null(estimate)) {
        estimate$diagnostics$sign_check %||% estimate$diagnostics$ratio_check %||% default_extreme_check
      } else {
        default_extreme_check
      }
    },
    no_interpretation = TRUE,
    no_rounding = TRUE
  )

  # ---------------------------------------------------------------------------
  # 8. Metadata bundle.
  # ---------------------------------------------------------------------------

  metadata <- list(
    validation = validation[setdiff(names(validation), "diagnostics")],
    preparation = prepared$metadata,
    fitting = fitted$metadata,
    selection = selected$metadata,
    estimation = if (!is.null(estimate)) estimate$metadata else NULL,
    resampling = if (!is.null(resample)) resample$metadata else NULL,
    formatting = list(
      output_module = "sii_rii_format_outputs",
      output_rule = "Technical outputs are assembled without recalculation, interpretation, or rounding.",
      interval_source = if (ci_method == "wald") "estimate" else "resample",
      final_results_overall_source = if (ci_method == "wald") estimate_helper_name else resample_helper_name
    )
  )

  # ---------------------------------------------------------------------------
  # 9. Final structured output.
  # ---------------------------------------------------------------------------

  out <- list(
    summary_stats = summary_stats,
    results_overall = results_overall,
    results_models = results_models,
    results_units = results_units,
    results_ridit = results_ridit,
    results_extremes = results_extremes,
    results_replicates = results_replicates,
    methodology = methodology,
    metadata = metadata,
    diagnostics = diagnostics
  )

  class(out) <- c(result_class, "sii_rii_ineqeco_result", class(out))
  out
}
