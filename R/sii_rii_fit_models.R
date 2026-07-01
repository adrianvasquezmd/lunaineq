#' Generalized Linear Modeling Engine for SII/RII Distributions
#'
#' @description
#' `sii_rii_fit_models` programmatically cycles through candidate regression classes
#' (e.g., standard binomial, quasibinomial, Poisson, quasi-Poisson, negative binomial,
#' log-linear, linear) dependent on the underlying distribution specified by the input
#' configuration matrix.
#'
#' The engine iteratively fits models to the empirical ridit space, extracting convergence
#' diagnostics, response-scale likelihood components, and primary model objects required
#' for downstream `emmeans` contrast processing.
#'
#' @importFrom stats glm binomial quasibinomial poisson quasipoisson lm predict setNames
#' @noRd

#'
#' @keywords internal
.sii_rii_fit_models <- function(prepared, validation) {

  # ---------------------------------------------------------------------------
  # 1. Validate contracts from upstream modules.
  # ---------------------------------------------------------------------------

  required_prepared_fields <- c("data", "metadata", "diagnostics")
  missing_prepared_fields <- setdiff(required_prepared_fields, names(prepared))
  if (length(missing_prepared_fields) > 0L) {
    stop(
      "Invalid `prepared` object: missing field(s): ",
      paste(missing_prepared_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_prepare_data()` before `.sii_rii_fit_models()`.",
      call. = FALSE
    )
  }

  required_validation_fields <- c(
    "health_indicator_type",
    "input_source",
    "candidate_models",
    "model_set_source",
    "model_selection_metric",
    "diagnostics"
  )
  missing_validation_fields <- setdiff(required_validation_fields, names(validation))
  if (length(missing_validation_fields) > 0L) {
    stop(
      "Invalid `validation` object: missing field(s): ",
      paste(missing_validation_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_validate_inputs()` before `.sii_rii_fit_models()`.",
      call. = FALSE
    )
  }

  d <- prepared$data
  if (!is.data.frame(d) || nrow(d) < 2L) {
    stop("Invalid prepared data: at least two ecological units are required.",
         call. = FALSE)
  }

  required_data_columns <- c(
    "unit_id",
    "population_weight",
    "social_position",
    "health_indicator_internal",
    "health_indicator_output",
    "model_response",
    "log_response",
    "numerator",
    "denominator",
    "count_failures",
    "model_offset"
  )
  missing_data_columns <- setdiff(required_data_columns, names(d))
  if (length(missing_data_columns) > 0L) {
    stop(
      "Invalid prepared data: missing standardized column(s): ",
      paste(missing_data_columns, collapse = ", "), ". ",
      "Action: rebuild the prepared data with `.sii_rii_prepare_data()`.",
      call. = FALSE
    )
  }

  health_indicator_type <- validation$health_indicator_type
  input_source <- validation$input_source
  candidate_models <- unique(as.character(validation$candidate_models))

  if (length(candidate_models) == 0L || anyNA(candidate_models)) {
    stop("No candidate models were provided by validation. Action: check `models` and indicator type.",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 2. Small internal helpers.
  # ---------------------------------------------------------------------------

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  empty_model_record <- function(model_name,
                                 fit_scale,
                                 evaluation_scale,
                                 outcome_scale,
                                 status = "not_attempted",
                                 fitted = FALSE,
                                 converged = FALSE,
                                 prediction_valid = FALSE,
                                 reason_excluded = NA_character_,
                                 message = NA_character_) {
    data.frame(
      Model = model_name,
      Fit_Scale = fit_scale,
      Evaluation_Scale = evaluation_scale,
      Outcome_Scale = outcome_scale,
      Status = status,
      Fitted = fitted,
      Converged = converged,
      Prediction_Valid = prediction_valid,
      Reason_Excluded = reason_excluded,
      Message = message,
      stringsAsFactors = FALSE
    )
  }

  model_scales <- function(model_name) {
    output_scale <- prepared$metadata$outcome_output_scale %||% "final_outcome_scale"

    if (model_name %in% c("binomial", "quasibinomial")) {
      return(list(
        fit_scale = "binomial_logit_link",
        evaluation_scale = "final_outcome_scale_after_inverse_logit",
        outcome_scale = output_scale
      ))
    }

    if (model_name %in% c("poisson", "quasipoisson", "negative_binomial")) {
      return(list(
        fit_scale = "count_model_with_log_link_and_exposure_offset",
        evaluation_scale = "final_rate_or_ratio_scale_after_offset_removal",
        outcome_scale = output_scale
      ))
    }

    if (model_name == "linear_log") {
      return(list(
        fit_scale = "log_final_outcome_scale",
        evaluation_scale = "final_outcome_scale_after_exponentiation",
        outcome_scale = output_scale
      ))
    }

    if (model_name == "linear_identity") {
      return(list(
        fit_scale = "final_outcome_scale_identity",
        evaluation_scale = "final_outcome_scale",
        outcome_scale = output_scale
      ))
    }

    list(
      fit_scale = "unknown",
      evaluation_scale = "unknown",
      outcome_scale = output_scale
    )
  }

  check_prediction_validity <- function(pred, health_indicator_type) {
    if (length(pred) != nrow(d) || any(!is.finite(pred))) {
      return(list(ok = FALSE, reason = "Predictions are missing, non-finite, or have the wrong length."))
    }

    tol <- sqrt(.Machine$double.eps)

    if (health_indicator_type == "proportion") {
      if (any(pred < -tol | pred > 1 + tol)) {
        return(list(ok = FALSE, reason = "Predictions are outside the valid proportion range [0, 1]."))
      }
    }

    if (health_indicator_type == "percentage") {
      if (any(pred < -tol | pred > 100 + tol)) {
        return(list(ok = FALSE, reason = "Predictions are outside the valid percentage range [0, 100]."))
      }
    }

    if (health_indicator_type %in% c("rate", "ratio")) {
      if (any(pred < -tol)) {
        return(list(ok = FALSE, reason = "Predictions are negative on the final rate/ratio scale."))
      }
    }

    list(ok = TRUE, reason = NA_character_)
  }

  with_captured_warnings <- function(expr) {
    warnings <- character(0)
    value <- withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(value = value, warnings = unique(warnings))
  }

  has_required_data_for_model <- function(model_name) {
    if (model_name == "linear_log" && !all(is.finite(d$log_response))) {
      return(list(
        ok = FALSE,
        reason = "`linear_log` requires strictly positive observed outcome values on the final scale. Action: use another model or remove zero-valued outcomes if methodologically justified."
      ))
    }

    if (model_name %in% c("binomial", "quasibinomial")) {
      if (input_source != "counts" && model_name == "binomial") {
        return(list(ok = FALSE, reason = "`binomial` requires count data and is not available for indicator-population input."))
      }
      if (health_indicator_type %in% c("rate", "ratio")) {
        return(list(ok = FALSE, reason = "Binomial-family models are only valid for proportion or percentage indicators."))
      }
    }

    if (model_name %in% c("poisson", "quasipoisson", "negative_binomial")) {
      if (input_source != "counts") {
        return(list(ok = FALSE, reason = "Count models require numerator and denominator/exposure variables."))
      }
      if (!health_indicator_type %in% c("rate", "ratio")) {
        return(list(ok = FALSE, reason = "Poisson-type count models are only valid for rate or ratio indicators in this SII implementation."))
      }
      if (any(!is.finite(d$model_offset))) {
        return(list(ok = FALSE, reason = "Count-rate models require a finite exposure offset for every unit."))
      }
    }

    list(ok = TRUE, reason = NA_character_)
  }

  # ---------------------------------------------------------------------------
  # 3. Fit one model and return a standardized fitted-model object.
  # ---------------------------------------------------------------------------

  fit_one_model <- function(model_name) {
    scales <- model_scales(model_name)
    precheck <- has_required_data_for_model(model_name)

    if (!isTRUE(precheck$ok)) {
      return(list(
        model_name = model_name,
        object = NULL,
        predictions = rep(NA_real_, nrow(d)),
        record = empty_model_record(
          model_name = model_name,
          fit_scale = scales$fit_scale,
          evaluation_scale = scales$evaluation_scale,
          outcome_scale = scales$outcome_scale,
          status = "excluded_before_fit",
          fitted = FALSE,
          converged = FALSE,
          prediction_valid = FALSE,
          reason_excluded = precheck$reason,
          message = precheck$reason
        ),
        warnings = character(0)
      ))
    }

    fit_result <- tryCatch(
      with_captured_warnings({
        if (model_name == "binomial") {
          stats::glm(
            cbind(numerator, count_failures) ~ social_position,
            family = stats::binomial(),
            data = d
          )
        } else if (model_name == "quasibinomial") {
          if (input_source == "counts") {
            stats::glm(
              cbind(numerator, count_failures) ~ social_position,
              family = stats::quasibinomial(),
              data = d
            )
          } else {
            stats::glm(
              model_response ~ social_position,
              family = stats::quasibinomial(),
              weights = population_weight,
              data = d
            )
          }
        } else if (model_name == "poisson") {
          stats::glm(
            model_response ~ social_position + offset(model_offset),
            family = stats::poisson(),
            data = d
          )
        } else if (model_name == "quasipoisson") {
          stats::glm(
            model_response ~ social_position + offset(model_offset),
            family = stats::quasipoisson(),
            data = d
          )
        } else if (model_name == "negative_binomial") {
          if (!requireNamespace("MASS", quietly = TRUE)) {
            stop("Package `MASS` is required for `negative_binomial`. Action: install MASS or remove `negative_binomial` from `models`.")
          }
          MASS::glm.nb(
            model_response ~ social_position + offset(model_offset),
            data = d
          )
        } else if (model_name == "linear_log") {
          stats::lm(
            log_response ~ social_position,
            weights = population_weight,
            data = d
          )
        } else if (model_name == "linear_identity") {
          stats::lm(
            health_indicator_output ~ social_position,
            weights = population_weight,
            data = d
          )
        } else {
          stop("Unknown model name: ", model_name)
        }
      }),
      error = function(e) {
        list(error = conditionMessage(e))
      }
    )

    if (!is.null(fit_result$error)) {
      return(list(
        model_name = model_name,
        object = NULL,
        predictions = rep(NA_real_, nrow(d)),
        record = empty_model_record(
          model_name = model_name,
          fit_scale = scales$fit_scale,
          evaluation_scale = scales$evaluation_scale,
          outcome_scale = scales$outcome_scale,
          status = "fit_error",
          fitted = FALSE,
          converged = FALSE,
          prediction_valid = FALSE,
          reason_excluded = "Model fitting failed.",
          message = paste0("Model fitting failed. Action: inspect the data or remove this model. Error: ", fit_result$error)
        ),
        warnings = character(0)
      ))
    }

    fit <- fit_result$value
    fit_warnings <- fit_result$warnings

    converged <- TRUE
    if (!is.null(fit$converged)) {
      converged <- isTRUE(fit$converged)
    }

    predictions <- tryCatch({
      if (model_name %in% c("binomial", "quasibinomial")) {
        p <- as.numeric(stats::predict(fit, newdata = d, type = "response"))
        if (health_indicator_type == "percentage") 100 * p else p
      } else if (model_name %in% c("poisson", "quasipoisson", "negative_binomial")) {
        predicted_counts <- as.numeric(stats::predict(fit, newdata = d, type = "response"))
        predicted_counts / exp(d$model_offset)
      } else if (model_name == "linear_log") {
        as.numeric(exp(stats::predict(fit, newdata = d)))
      } else if (model_name == "linear_identity") {
        as.numeric(stats::predict(fit, newdata = d))
      } else {
        rep(NA_real_, nrow(d))
      }
    }, error = function(e) {
      res <- rep(NA_real_, nrow(d))
      attr(res, "prediction_error") <- conditionMessage(e)
      res
    })

    if (is.character(predictions) && length(predictions) == 1L) {
      prediction_error <- predictions
      predictions <- rep(NA_real_, nrow(d))
    } else {
      prediction_error <- NA_character_
    }

    pred_check <- check_prediction_validity(predictions, health_indicator_type)
    prediction_valid <- isTRUE(pred_check$ok)

    status <- if (converged && prediction_valid) {
      "fitted"
    } else if (!converged) {
      "not_converged"
    } else {
      "invalid_predictions"
    }

    reason_excluded <- NA_character_
    message <- NA_character_
    if (!converged) {
      reason_excluded <- "Model did not report convergence."
      message <- "Model was fitted but did not report convergence. Action: inspect model diagnostics or choose another model."
    }
    if (!prediction_valid) {
      reason_excluded <- pred_check$reason
      message <- paste0(pred_check$reason, " Action: inspect this model or allow selection among other candidate models.")
    }
    if (!is.na(prediction_error)) {
      reason_excluded <- "Prediction failed."
      message <- paste0("Prediction failed. Action: inspect this model or remove it. Error: ", prediction_error)
    }
    if (length(fit_warnings) > 0L && is.na(message)) {
      message <- paste(fit_warnings, collapse = " | ")
    }

    list(
      model_name = model_name,
      object = fit,
      predictions = predictions,
      record = empty_model_record(
        model_name = model_name,
        fit_scale = scales$fit_scale,
        evaluation_scale = scales$evaluation_scale,
        outcome_scale = scales$outcome_scale,
        status = status,
        fitted = TRUE,
        converged = converged,
        prediction_valid = prediction_valid,
        reason_excluded = reason_excluded,
        message = message
      ),
      warnings = fit_warnings
    )
  }

  # ---------------------------------------------------------------------------
  # 4. Fit all candidate models independently.
  # ---------------------------------------------------------------------------

  fitted_list <- stats::setNames(lapply(candidate_models, fit_one_model), candidate_models)

  model_records <- do.call(rbind, lapply(fitted_list, function(x) x$record))
  rownames(model_records) <- NULL

  prediction_columns <- lapply(fitted_list, function(x) x$predictions)
  predictions <- data.frame(
    Unit_ID = d$unit_id,
    Social_Position = d$social_position,
    Observed = d$health_indicator_output,
    Weight = d$population_weight,
    stringsAsFactors = FALSE
  )
  for (nm in names(prediction_columns)) {
    predictions[[paste0("Pred_", nm)]] <- prediction_columns[[nm]]
  }

  n_fitted <- sum(model_records$Fitted)
  n_converged <- sum(model_records$Converged)
  n_viable <- sum(model_records$Fitted & model_records$Converged & model_records$Prediction_Valid)

  if (n_viable == 0L) {
    reasons <- paste(
      paste0(model_records$Model, ": ", model_records$Message %||% model_records$Reason_Excluded),
      collapse = " | "
    )
    stop(
      "None of the candidate SII models produced valid fitted predictions. ",
      "Action: inspect `models`, indicator type, and outcome values; try a simpler compatible model such as `linear_identity` when methodologically appropriate. ",
      "Details: ", reasons,
      call. = FALSE
    )
  }

  # Keep only successful model objects in a named list for downstream modules.
  model_objects <- lapply(fitted_list, function(x) x$object)
  names(model_objects) <- names(fitted_list)

  model_warnings <- lapply(fitted_list, function(x) x$warnings)
  names(model_warnings) <- names(fitted_list)

  # ---------------------------------------------------------------------------
  # 5. Return auditable fitted-model object.
  # ---------------------------------------------------------------------------

  list(
    models = model_objects,
    model_table = model_records,
    predictions = predictions,
    prepared_data = d,
    metadata = list(
      health_indicator_type = health_indicator_type,
      input_source = input_source,
      candidate_models = candidate_models,
      model_set_source = validation$model_set_source,
      model_selection_metric = validation$model_selection_metric,
      measure = if (!is.null(validation$measure)) validation$measure else "sii",
      estimand = if (!is.null(validation$estimand)) validation$estimand else "SII = Y_advantaged - Y_disadvantaged",
      estimand_type = if (!is.null(validation$estimand_type)) validation$estimand_type else "absolute_difference",
      extreme_contrast_rule = if (!is.null(validation$extreme_contrast_rule)) validation$extreme_contrast_rule else "f(advantaged_position) - f(disadvantaged_position)",
      prediction_rule = "All fitted-model predictions are returned on the final outcome scale used for model-selection metrics.",
      transformed_model_rule = "Models fitted on log scales are back-transformed before prediction comparison. Logit-transformed linear models are not part of this version.",
      count_rate_prediction_rule = "Poisson-type count models with exposure offsets are converted back to the final rate/ratio scale by dividing predicted counts by exp(offset).",
      n_candidate_models = length(candidate_models),
      n_fitted_models = n_fitted,
      n_converged_models = n_converged,
      n_viable_models = n_viable
    ),
    diagnostics = list(
      warnings_by_model = model_warnings,
      excluded_models = model_records$Model[!model_records$Fitted | !model_records$Converged | !model_records$Prediction_Valid],
      fitted_models = model_records$Model[model_records$Fitted],
      viable_models = model_records$Model[model_records$Fitted & model_records$Converged & model_records$Prediction_Valid]
    )
  )
}
