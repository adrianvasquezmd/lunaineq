#' Discriminative Model Selection Array
#'
#' @description
#' `sii_rii_select_model` executes comparative predictive accuracy assessments across all
#' fitted GLMs provided by the estimation engine. 
#'
#' It systematically scores models using absolute and quadratic distance metrics (MAE, MSE,
#' RMSE, MAPE) computed explicitly on the un-transformed health outcome space to ensure
#' metric parity across non-nested link functions. If maximum-likelihood paradigms map consistently,
#' it permits classical information criteria (AIC/BIC). Ultimately, the algorithm isolates the
#' optimal model object for the variance phase.
#'
#' The shared SII/RII workflow uses the convention:
#'   social_position = 0 -> most disadvantaged
#'   social_position = 1 -> most advantaged
#'
#' The final estimand is not computed here. This helper only selects the model.
#' The selected model can then be used to estimate:
#'   SII = Y_advantaged - Y_disadvantaged
#'   RII = Y_advantaged / Y_disadvantaged
#'
#' @importFrom stats AIC BIC
#' @keywords internal
#' @noRd
.sii_rii_select_model <- function(fitted, validation) {

  # ---------------------------------------------------------------------------
  # 1. Validate contracts from upstream modules.
  # ---------------------------------------------------------------------------

  required_fitted_fields <- c(
    "models",
    "model_table",
    "predictions",
    "prepared_data",
    "metadata",
    "diagnostics"
  )
  missing_fitted_fields <- setdiff(required_fitted_fields, names(fitted))
  if (length(missing_fitted_fields) > 0L) {
    stop(
      "Invalid `fitted` object: missing field(s): ",
      paste(missing_fitted_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_fit_models()` before `.sii_rii_select_model()`.",
      call. = FALSE
    )
  }

  required_validation_fields <- c(
    "health_indicator_type",
    "input_source",
    "candidate_models",
    "model_selection_metric",
    "model_set_source",
    "aic_bic_comparability_note",
    "diagnostics"
  )
  missing_validation_fields <- setdiff(required_validation_fields, names(validation))
  if (length(missing_validation_fields) > 0L) {
    stop(
      "Invalid `validation` object: missing field(s): ",
      paste(missing_validation_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_validate_inputs()` before `.sii_rii_select_model()`.",
      call. = FALSE
    )
  }

  model_table <- fitted$model_table
  predictions <- fitted$predictions
  model_objects <- fitted$models
  measure <- if (!is.null(validation$measure)) validation$measure else {
    if (!is.null(fitted$metadata$measure)) fitted$metadata$measure else "sii"
  }
  measure <- match.arg(measure, choices = c("sii", "rii"))

  estimand <- if (!is.null(validation$estimand)) validation$estimand else {
    if (identical(measure, "sii")) {
      "SII = Y_advantaged - Y_disadvantaged"
    } else {
      "RII = Y_advantaged / Y_disadvantaged"
    }
  }

  estimand_type <- if (!is.null(validation$estimand_type)) validation$estimand_type else {
    if (identical(measure, "sii")) "absolute_difference" else "relative_ratio"
  }

  extreme_contrast_rule <- if (!is.null(validation$extreme_contrast_rule)) {
    validation$extreme_contrast_rule
  } else if (identical(measure, "sii")) {
    "f(advantaged_position) - f(disadvantaged_position)"
  } else {
    "f(advantaged_position) / f(disadvantaged_position)"
  }


  if (!is.data.frame(model_table) || nrow(model_table) == 0L) {
    stop(
      "Invalid fitted object: `model_table` must be a non-empty data frame. ",
      "Action: inspect `.sii_rii_fit_models()` output.",
      call. = FALSE
    )
  }

  required_model_columns <- c(
    "Model",
    "Fit_Scale",
    "Evaluation_Scale",
    "Outcome_Scale",
    "Status",
    "Fitted",
    "Converged",
    "Prediction_Valid",
    "Reason_Excluded",
    "Message"
  )
  missing_model_columns <- setdiff(required_model_columns, names(model_table))
  if (length(missing_model_columns) > 0L) {
    stop(
      "Invalid fitted object: `model_table` is missing column(s): ",
      paste(missing_model_columns, collapse = ", "), ". ",
      "Action: rebuild fitted models with `.sii_rii_fit_models()`.",
      call. = FALSE
    )
  }

  required_prediction_columns <- c("Unit_ID", "Social_Position", "Observed", "Weight")
  missing_prediction_columns <- setdiff(required_prediction_columns, names(predictions))
  if (length(missing_prediction_columns) > 0L) {
    stop(
      "Invalid fitted object: `predictions` is missing column(s): ",
      paste(missing_prediction_columns, collapse = ", "), ". ",
      "Action: rebuild fitted models with `.sii_rii_fit_models()`.",
      call. = FALSE
    )
  }

  metric <- toupper(as.character(validation$model_selection_metric))
  if (!metric %in% c("MSE", "MAE", "RMSE", "MAPE", "AIC", "BIC")) {
    stop(
      "Unsupported `model_selection_metric`: ", metric, ". ",
      "Action: use one of MSE, MAE, RMSE, MAPE, AIC, or BIC.",
      call. = FALSE
    )
  }

  candidate_models <- unique(as.character(validation$candidate_models))
  if (length(candidate_models) == 0L || anyNA(candidate_models)) {
    stop(
      "No candidate models were available for selection. ",
      "Action: inspect `models`, `health_indicator_type`, and input source.",
      call. = FALSE
    )
  }

  # Keep the model table aligned with validation order where possible.
  model_table$Model <- as.character(model_table$Model)
  model_table <- model_table[match(candidate_models, model_table$Model), , drop = FALSE]
  if (anyNA(model_table$Model)) {
    stop(
      "The fitted model table does not contain all validation candidate models. ",
      "Action: rerun `.sii_rii_fit_models()` using the validation object used for model selection.",
      call. = FALSE
    )
  }
  rownames(model_table) <- NULL

  # ---------------------------------------------------------------------------
  # 2. Small internal helpers.
  # ---------------------------------------------------------------------------

  is_true <- function(x) {
    isTRUE(x) && length(x) == 1L
  }

  safe_numeric <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x
  }

  weighted_mean_safe <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) {
      return(NA_real_)
    }
    sum(w[ok] * x[ok]) / sum(w[ok])
  }

  weighted_mae <- function(y, yhat, w) {
    weighted_mean_safe(abs(y - yhat), w)
  }

  weighted_mse <- function(y, yhat, w) {
    weighted_mean_safe((y - yhat)^2, w)
  }

  weighted_rmse <- function(y, yhat, w) {
    mse <- weighted_mse(y, yhat, w)
    if (!is.finite(mse)) NA_real_ else sqrt(mse)
  }

  weighted_mape <- function(y, yhat, w) {
    if (any(is.finite(y) & y == 0)) {
      return(NA_real_)
    }
    weighted_mean_safe(abs((y - yhat) / y), w) * 100
  }

  is_viable_model_row <- function(row) {
    isTRUE(row$Fitted) && isTRUE(row$Converged) && isTRUE(row$Prediction_Valid)
  }

  aic_bic_set_is_comparable <- function(type, source, models) {
    if (source == "counts" && type %in% c("rate", "ratio") &&
        all(models %in% c("poisson", "negative_binomial"))) {
      return(TRUE)
    }
    if (source == "counts" && type %in% c("proportion", "percentage") &&
        all(models %in% c("binomial"))) {
      return(TRUE)
    }
    FALSE
  }

  get_information_criterion <- function(model_object, criterion) {
    if (is.null(model_object)) {
      return(NA_real_)
    }
    value <- tryCatch(
      {
        if (criterion == "AIC") {
          stats::AIC(model_object)
        } else {
          stats::BIC(model_object)
        }
      },
      error = function(e) NA_real_,
      warning = function(w) NA_real_
    )
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) NA_real_ else value
  }

  # ---------------------------------------------------------------------------
  # 3. Validate observed outcome and weights used for model-selection metrics.
  # ---------------------------------------------------------------------------

  observed <- safe_numeric(predictions$Observed)
  weights <- safe_numeric(predictions$Weight)

  if (length(observed) != nrow(predictions) || any(!is.finite(observed))) {
    stop(
      "Invalid fitted predictions: observed outcome values must be finite numeric values. ",
      "Action: inspect `.sii_rii_prepare_data()` and `.sii_rii_fit_models()` outputs.",
      call. = FALSE
    )
  }

  if (length(weights) != nrow(predictions) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop(
      "Invalid fitted predictions: model-selection weights must be positive finite numeric values. ",
      "Action: inspect denominators or population weights in the input data.",
      call. = FALSE
    )
  }

  if (metric == "MAPE" && any(observed == 0)) {
    stop(
      "`model_selection_metric = 'MAPE'` is not available because the observed outcome contains zero values. ",
      "Action: use MSE, RMSE, or MAE instead. These metrics remain comparable across candidate models when predictions are evaluated on the final outcome scale.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 4. Compute predictive metrics on the final outcome scale for each model.
  # ---------------------------------------------------------------------------

  metrics <- data.frame(
    Model = model_table$Model,
    MAE = rep(NA_real_, nrow(model_table)),
    MSE = rep(NA_real_, nrow(model_table)),
    RMSE = rep(NA_real_, nrow(model_table)),
    MAPE = rep(NA_real_, nrow(model_table)),
    AIC = rep(NA_real_, nrow(model_table)),
    BIC = rep(NA_real_, nrow(model_table)),
    Metric_Available = rep(FALSE, nrow(model_table)),
    Comparable_For_Selected_Metric = rep(FALSE, nrow(model_table)),
    Selected = rep(FALSE, nrow(model_table)),
    Selection_Value = rep(NA_real_, nrow(model_table)),
    Selection_Rank = rep(NA_integer_, nrow(model_table)),
    Selection_Reason = rep(NA_character_, nrow(model_table)),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(model_table))) {
    model_name <- model_table$Model[i]
    pred_col <- paste0("Pred_", model_name)

    if (!is_viable_model_row(model_table[i, , drop = FALSE])) {
      metrics$Selection_Reason[i] <- model_table$Reason_Excluded[i]
      if (is.na(metrics$Selection_Reason[i]) || !nzchar(metrics$Selection_Reason[i])) {
        metrics$Selection_Reason[i] <- model_table$Message[i]
      }
      next
    }

    if (!pred_col %in% names(predictions)) {
      metrics$Selection_Reason[i] <- paste0(
        "Prediction column `", pred_col, "` was not found. Action: rerun `.sii_rii_fit_models()`."
      )
      next
    }

    yhat <- safe_numeric(predictions[[pred_col]])
    if (length(yhat) != length(observed) || any(!is.finite(yhat))) {
      metrics$Selection_Reason[i] <- "Predictions are missing or non-finite on the final outcome scale."
      next
    }

    metrics$MAE[i] <- weighted_mae(observed, yhat, weights)
    metrics$MSE[i] <- weighted_mse(observed, yhat, weights)
    metrics$RMSE[i] <- weighted_rmse(observed, yhat, weights)
    metrics$MAPE[i] <- weighted_mape(observed, yhat, weights)

    if (!is.null(model_objects[[model_name]])) {
      metrics$AIC[i] <- get_information_criterion(model_objects[[model_name]], "AIC")
      metrics$BIC[i] <- get_information_criterion(model_objects[[model_name]], "BIC")
    }
  }

  # ---------------------------------------------------------------------------
  # 5. Apply comparability rules for the selected metric.
  # ---------------------------------------------------------------------------

  if (metric %in% c("MAE", "MSE", "RMSE", "MAPE")) {
    metrics$Comparable_For_Selected_Metric <- model_table$Fitted &
      model_table$Converged &
      model_table$Prediction_Valid

    metrics$Selection_Value <- metrics[[metric]]
    metrics$Metric_Available <- is.finite(metrics$Selection_Value)

    unavailable <- metrics$Comparable_For_Selected_Metric & !metrics$Metric_Available
    metrics$Selection_Reason[unavailable] <- paste0(
      metric,
      " could not be computed on the final outcome scale. Action: inspect predictions and observed outcome values."
    )
  } else {
    viable_for_aic_bic <- model_table$Model[
      model_table$Fitted & model_table$Converged & model_table$Prediction_Valid
    ]

    if (!aic_bic_set_is_comparable(
      type = validation$health_indicator_type,
      source = validation$input_source,
      models = viable_for_aic_bic
    )) {
      stop(
        "`model_selection_metric = '", metric, "'` cannot be used because the viable fitted models are not likelihood-comparable. ",
        "Action: use MSE, RMSE, MAE, or MAPE; or provide a likelihood-comparable model subset. ",
        "For rates/ratios with counts, use `models = c('poisson', 'negative_binomial')`. ",
        "For proportions/percentages with counts, AIC/BIC is only meaningful for binomial-only selection.",
        call. = FALSE
      )
    }

    metrics$Comparable_For_Selected_Metric <- model_table$Fitted &
      model_table$Converged &
      model_table$Prediction_Valid &
      model_table$Model %in% viable_for_aic_bic

    metrics$Selection_Value <- metrics[[metric]]
    metrics$Metric_Available <- metrics$Comparable_For_Selected_Metric &
      is.finite(metrics$Selection_Value)

    unavailable <- metrics$Comparable_For_Selected_Metric & !metrics$Metric_Available
    metrics$Selection_Reason[unavailable] <- paste0(
      metric,
      " was not available for this fitted model. Action: use MSE/RMSE/MAE/MAPE or remove models without a comparable likelihood."
    )
  }

  eligible <- metrics$Comparable_For_Selected_Metric & metrics$Metric_Available

  if (!any(eligible)) {
    attempted <- paste(
      paste0(metrics$Model, ": ", ifelse(is.na(metrics$Selection_Reason), "not eligible", metrics$Selection_Reason)),
      collapse = " | "
    )
    stop(
      "No fitted model was eligible for selection using `model_selection_metric = '", metric, "'`. ",
      "Action: use `models = NULL` to let the function evaluate all allowed models, choose a simpler compatible model, or change the selection metric to MSE, RMSE, or MAE. ",
      "Details: ", attempted,
      call. = FALSE
    )
  }

  # Lower is better for all supported metrics.
  rank_values <- metrics$Selection_Value
  rank_values[!eligible] <- Inf
  rank_order <- order(rank_values, metrics$Model, na.last = TRUE)
  finite_order <- rank_order[is.finite(rank_values[rank_order])]

  metrics$Selection_Rank[finite_order] <- seq_along(finite_order)
  selected_index <- finite_order[1L]
  selected_model_name <- metrics$Model[selected_index]
  metrics$Selected[selected_index] <- TRUE
  metrics$Selection_Reason[selected_index] <- paste0(
    "Selected because it had the lowest valid ", metric,
    " among eligible candidate models."
  )

  # Add clear reasons for eligible non-selected models.
  eligible_not_selected <- which(eligible & !metrics$Selected)
  if (length(eligible_not_selected) > 0L) {
    metrics$Selection_Reason[eligible_not_selected] <- paste0(
      "Eligible but not selected because another model had a lower ", metric, "."
    )
  }

  not_comparable <- which(!metrics$Comparable_For_Selected_Metric &
                            model_table$Fitted & model_table$Converged &
                            model_table$Prediction_Valid)
  if (length(not_comparable) > 0L) {
    metrics$Selection_Reason[not_comparable] <- paste0(
      "Fitted model was not comparable for the selected metric `", metric, "`."
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Assemble auditable model-selection table.
  # ---------------------------------------------------------------------------

  selection_table <- cbind(
    model_table,
    metrics[, setdiff(names(metrics), "Model"), drop = FALSE]
  )

  # Keep metric columns close to model metadata for readability.
  desired_order <- c(
    "Model",
    "Fit_Scale",
    "Evaluation_Scale",
    "Outcome_Scale",
    "Status",
    "Fitted",
    "Converged",
    "Prediction_Valid",
    "Comparable_For_Selected_Metric",
    "Metric_Available",
    "MAE",
    "MSE",
    "RMSE",
    "MAPE",
    "AIC",
    "BIC",
    "Selection_Value",
    "Selection_Rank",
    "Selected",
    "Reason_Excluded",
    "Selection_Reason",
    "Message"
  )
  selection_table <- selection_table[, c(desired_order, setdiff(names(selection_table), desired_order)), drop = FALSE]
  rownames(selection_table) <- NULL

  selected_model_object <- model_objects[[selected_model_name]]
  if (is.null(selected_model_object)) {
    stop(
      "Internal selection error: selected model object is NULL. ",
      "Action: rerun `.sii_rii_fit_models()` and inspect model fitting diagnostics.",
      call. = FALSE
    )
  }

  selected_prediction_col <- paste0("Pred_", selected_model_name)
  selected_predictions <- predictions[[selected_prediction_col]]

  # ---------------------------------------------------------------------------
  # 7. Return selected model and complete audit trail.
  # ---------------------------------------------------------------------------

  list(
    selected_model_name = selected_model_name,
    selected_model = selected_model_object,
    selected_predictions = selected_predictions,
    model_selection_table = selection_table,
    predictions = predictions,
    prepared_data = fitted$prepared_data,
    fitted_models = fitted$models,
    metadata = list(
      health_indicator_type = validation$health_indicator_type,
      input_source = validation$input_source,
      candidate_models = candidate_models,
      model_set_source = validation$model_set_source,
      model_selection_metric = metric,
      selection_rule = paste0("Lowest valid ", metric, " among eligible candidate models."),
      predictive_metric_rule = "MAE, MSE, RMSE, and MAPE are computed on the final outcome scale and weighted by population/exposure weights.",
      aic_bic_rule = "AIC/BIC are used only for likelihood-comparable fitted model sets; otherwise selection stops with an actionable message.",
      transformed_model_rule = "Models fitted on transformed scales are compared only after back-transformation to the final outcome scale.",
      metric_weighting_rule = "Predictive model-selection metrics are always weighted by population/exposure weights.",
      selected_model = selected_model_name,
      selected_metric_value = metrics$Selection_Value[selected_index],
      measure = measure,
      estimand = estimand,
      estimand_type = estimand_type,
      extreme_contrast_rule = extreme_contrast_rule
    ),
    diagnostics = list(
      n_candidate_models = length(candidate_models),
      n_eligible_models = sum(eligible),
      n_metric_available = sum(metrics$Metric_Available),
      selected_model = selected_model_name,
      selected_model_rank = metrics$Selection_Rank[selected_index],
      measure = measure,
      estimand = estimand,
      excluded_or_ineligible_models = selection_table$Model[!selection_table$Metric_Available | !selection_table$Comparable_For_Selected_Metric],
      metric_scale = if (metric %in% c("MAE", "MSE", "RMSE", "MAPE")) {
        "final_outcome_scale_weighted_by_population_or_exposure"
      } else {
        "likelihood_scale_comparable_models_only"
      },
      warnings = character(0)
    )
  )
}
