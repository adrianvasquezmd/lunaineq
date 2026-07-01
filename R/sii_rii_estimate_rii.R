#' The confidence interval is then exponentiated:
#'   CI = exp(log(RII) +/- z * SE_log_RII)
#'
#' This avoids invalid negative confidence limits for a ratio.
#'
#' Important model-specific rule:
#'   For binomial/quasibinomial models, RII is a ratio of predicted probabilities:
#'     p_advantaged / p_disadvantaged
#'   It is NOT an odds ratio. The function therefore does not exponentiate the
#'   logit-scale coefficient or logit-scale contrast.
#'
#' @keywords internal
.sii_rii_estimate_rii <- function(selected, validation, conf_level = NULL, vcov_type = "HC1") {

  # ---------------------------------------------------------------------------
  # 1. Validate upstream contracts.
  # ---------------------------------------------------------------------------

  required_selected_fields <- c(
    "selected_model_name",
    "selected_model",
    "model_selection_table",
    "predictions",
    "prepared_data",
    "metadata",
    "diagnostics"
  )
  missing_selected_fields <- setdiff(required_selected_fields, names(selected))
  if (length(missing_selected_fields) > 0L) {
    stop(
      "Invalid `selected` object: missing field(s): ",
      paste(missing_selected_fields, collapse = ", "), ". ",
      "Action: call the model-selection helper before `.sii_rii_estimate_rii()`.",
      call. = FALSE
    )
  }

  required_validation_fields <- c(
    "health_indicator_type",
    "input_source",
    "social_position_method",
    "conf_level"
  )
  missing_validation_fields <- setdiff(required_validation_fields, names(validation))
  if (length(missing_validation_fields) > 0L) {
    stop(
      "Invalid `validation` object: missing field(s): ",
      paste(missing_validation_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_validate_inputs(..., measure = 'rii')` before `.sii_rii_estimate_rii()`.",
      call. = FALSE
    )
  }

  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop(
      "Package `emmeans` is required to compute the Wald RII contrast in this version. ",
      "Action: install `emmeans` or use a resampling interval once `rii_resample.R` is available.",
      call. = FALSE
    )
  }

  model_name <- as.character(selected$selected_model_name)
  model <- selected$selected_model
  d <- selected$prepared_data

  if (length(model_name) != 1L || is.na(model_name) || !nzchar(model_name)) {
    stop(
      "Invalid selected model name. ",
      "Action: inspect model-selection output.",
      call. = FALSE
    )
  }
  if (is.null(model)) {
    stop(
      "Invalid selected model: selected model object is NULL. ",
      "Action: inspect model fitting and model-selection outputs.",
      call. = FALSE
    )
  }
  if (!is.data.frame(d) || nrow(d) < 2L) {
    stop(
      "Invalid prepared data in selected object: at least two ecological units are required. ",
      "Action: inspect data-preparation output.",
      call. = FALSE
    )
  }

  required_data_columns <- c(
    "social_position",
    "health_indicator_output",
    "population_weight"
  )
  missing_data_columns <- setdiff(required_data_columns, names(d))
  if (length(missing_data_columns) > 0L) {
    stop(
      "Invalid prepared data in selected object: missing column(s): ",
      paste(missing_data_columns, collapse = ", "), ". ",
      "Action: rebuild the workflow from the data-preparation helper.",
      call. = FALSE
    )
  }

  health_indicator_type <- as.character(validation$health_indicator_type)
  input_source <- as.character(validation$input_source)
  social_position_method <- as.character(validation$social_position_method)

  if (!health_indicator_type %in% c("proportion", "percentage", "rate", "ratio")) {
    stop("Internal error: unsupported `health_indicator_type` in validation.", call. = FALSE)
  }
  if (!input_source %in% c("counts", "indicator_population")) {
    stop("Internal error: unsupported `input_source` in validation.", call. = FALSE)
  }
  if (!social_position_method %in% c("classic", "bounded")) {
    stop("Internal error: unsupported `social_position_method` in validation.", call. = FALSE)
  }

  if (is.null(conf_level)) {
    conf_level <- validation$conf_level
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      is.na(conf_level) || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop(
      "`conf_level` must be a single numeric value between 0 and 1. ",
      "Action: use a value such as 0.95.",
      call. = FALSE
    )
  }

  if (!is.character(vcov_type) || length(vcov_type) != 1L ||
      is.na(vcov_type) || !nzchar(vcov_type)) {
    stop(
      "`vcov_type` must be a single non-empty character value. ",
      "Action: use `vcov_type = 'HC1'` unless you have a specific reason to change it.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Internal helpers.
  # ---------------------------------------------------------------------------

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  get_vcov <- function(object, type) {
    robust_used <- FALSE
    warning_message <- NA_character_

    V <- NULL
    if (requireNamespace("sandwich", quietly = TRUE)) {
      V <- tryCatch(
        sandwich::vcovHC(object, type = type),
        error = function(e) {
          warning_message <<- paste0(
            "Robust sandwich covariance could not be computed (", conditionMessage(e), "). ",
            "Falling back to `stats::vcov()`."
          )
          NULL
        },
        warning = function(w) {
          warning_message <<- paste0(
            "Robust sandwich covariance produced a warning (", conditionMessage(w), "). ",
            "Falling back to `stats::vcov()`."
          )
          NULL
        }
      )
      if (!is.null(V)) {
        robust_used <- TRUE
      }
    } else {
      warning_message <- paste0(
        "Package `sandwich` is not installed. ",
        "Using `stats::vcov()` for the Wald interval."
      )
    }

    if (is.null(V)) {
      V <- tryCatch(
        stats::vcov(object),
        error = function(e) {
          stop(
            "Could not compute a covariance matrix for the selected model (", conditionMessage(e), "). ",
            "Action: inspect the selected model or use bootstrap/jackknife once available.",
            call. = FALSE
          )
        }
      )
    }

    list(
      vcov = V,
      robust_used = robust_used,
      warning = warning_message
    )
  }

  get_positions <- function(prepared_data, method) {
    if (method == "classic") {
      return(list(disadvantaged = 0, advantaged = 1))
    }

    observed <- prepared_data$social_position
    observed <- observed[is.finite(observed)]
    if (length(observed) < 2L) {
      stop(
        "Cannot compute bounded RII positions because fewer than two finite social positions are available. ",
        "Action: inspect data-preparation output.",
        call. = FALSE
      )
    }

    list(
      disadvantaged = min(observed),
      advantaged = max(observed)
    )
  }

  make_extreme_predictions <- function(values, positions) {
    data.frame(
      Position_Role = c("Disadvantaged", "Advantaged"),
      Social_Position = c(positions$disadvantaged, positions$advantaged),
      Predicted_Value = as.numeric(values),
      stringsAsFactors = FALSE
    )
  }

  require_positive_extremes <- function(values, model_name) {
    if (length(values) != 2L || any(!is.finite(values)) || any(values <= 0)) {
      stop(
        "RII cannot be estimated because at least one extreme predicted value is non-positive or non-finite. ",
        "RII requires Y_disadvantaged > 0 and Y_advantaged > 0. ",
        "Selected model: ", model_name, ". ",
        "Action: inspect `results_extremes`, use `social_position_method = 'bounded'`, ",
        "exclude models that can predict non-positive values such as `linear_identity`, ",
        "or use a model with a positive response scale when appropriate.",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  extract_emmeans_prediction <- function(emm_object, candidate_names = c("response", "prob", "rate", "emmean")) {
    out <- as.data.frame(summary(emm_object, infer = c(FALSE, FALSE)))
    pred_col <- intersect(candidate_names, names(out))[1L]
    if (is.na(pred_col)) {
      stop(
        "Could not extract response-scale predictions from `emmeans`. ",
        "Action: inspect the selected model and `emmeans` output.",
        call. = FALSE
      )
    }
    as.numeric(out[[pred_col]])
  }

  get_prediction_vcov <- function(emm_object) {
    V_theta <- tryCatch(
      as.matrix(stats::vcov(emm_object)),
      error = function(e) {
        stop(
          "Could not extract the covariance matrix of the two extreme predictions from `emmeans` (",
          conditionMessage(e), "). ",
          "Action: inspect selected model diagnostics.",
          call. = FALSE
        )
      }
    )

    if (!all(dim(V_theta) == c(2L, 2L)) || any(!is.finite(V_theta))) {
      stop(
        "Invalid covariance matrix for the two extreme predictions. ",
        "Action: inspect the selected model or choose a different model.",
        call. = FALSE
      )
    }

    V_theta
  }

  delta_log_ratio <- function(theta, V_theta, label) {
    if (!requireNamespace("msm", quietly = TRUE)) {
      stop(
        "Package `msm` is required to compute the Wald interval for RII when ",
        "the model does not provide a direct log-ratio contrast. ",
        "This applies to ", label, ". ",
        "Action: install `msm`, choose a log-link model when appropriate, or use bootstrap/jackknife once available.",
        call. = FALSE
      )
    }

    if (length(theta) != 2L || any(!is.finite(theta)) || any(theta <= 0)) {
      stop(
        "RII log-ratio delta method requires two strictly positive response-scale predictions. ",
        "Action: inspect extreme predictions or use a different model.",
        call. = FALSE
      )
    }

    as.numeric(msm::deltamethod(~ log(x2 / x1), theta, V_theta))
  }

  extract_link_contrast <- function(contrast_object, level) {
    out <- as.data.frame(summary(contrast_object, infer = c(TRUE, TRUE), level = level))

    estimate_col <- intersect(c("estimate", "emmean"), names(out))[1L]
    se_col <- intersect(c("SE", "se"), names(out))[1L]
    lower_col <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(out))[1L]
    upper_col <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(out))[1L]

    if (any(is.na(c(estimate_col, se_col, lower_col, upper_col)))) {
      stop(
        "Could not extract log-RII estimate, SE, and confidence limits from the `emmeans` contrast object. ",
        "Action: inspect the selected model and `emmeans` output.",
        call. = FALSE
      )
    }

    list(
      log_estimate = as.numeric(out[[estimate_col]][1L]),
      se_log = as.numeric(out[[se_col]][1L]),
      log_lower = as.numeric(out[[lower_col]][1L]),
      log_upper = as.numeric(out[[upper_col]][1L]),
      raw = out
    )
  }

  # ---------------------------------------------------------------------------
  # 3. Determine disadvantaged and advantaged positions.
  # ---------------------------------------------------------------------------

  positions <- get_positions(d, social_position_method)

  if (!is.finite(positions$disadvantaged) || !is.finite(positions$advantaged) ||
      positions$disadvantaged >= positions$advantaged) {
    stop(
      "Invalid RII positions: disadvantaged position must be less than advantaged position. ",
      "Action: inspect social position construction in data-preparation output.",
      call. = FALSE
    )
  }

  at_positions <- list(
    social_position = c(positions$disadvantaged, positions$advantaged)
  )

  # ---------------------------------------------------------------------------
  # 4. Compute covariance matrix.
  # ---------------------------------------------------------------------------

  vc <- get_vcov(model, vcov_type)
  V <- vc$vcov
  warnings <- character(0)
  if (!is.na(vc$warning)) {
    warnings <- c(warnings, vc$warning)
  }

  # ---------------------------------------------------------------------------
  # 5. Estimate extremes and Wald SE/CI.
  # ---------------------------------------------------------------------------

  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha / 2)

  estimation_engine <- NA_character_
  contrast_raw <- NULL
  extreme_predictions <- NULL
  estimate <- NA_real_
  log_estimate <- NA_real_
  se_log <- NA_real_
  ci_lower <- NA_real_
  ci_upper <- NA_real_
  log_ci_lower <- NA_real_
  log_ci_upper <- NA_real_

  if (model_name %in% c("poisson", "quasipoisson", "negative_binomial")) {

    # For log-link count models with offset = 0:
    #   log(Y_R) = alpha + beta * R
    # Therefore:
    #   log(RII) = log(Y_A) - log(Y_D)
    # This is a valid log-rate/log-ratio contrast. It is not an odds ratio.

    emm_link <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      type = "link",
      offset = 0,
      vcov. = V
    )
    ctr <- emmeans::contrast(emm_link, method = list(log_RII = c(-1, 1)))
    cs <- extract_link_contrast(ctr, conf_level)

    emm_resp <- emmeans::regrid(emm_link, transform = "response")
    pred_values <- extract_emmeans_prediction(
      emm_resp,
      candidate_names = c("rate", "response", "emmean")
    )

    require_positive_extremes(pred_values, model_name)

    log_estimate <- cs$log_estimate
    se_log <- cs$se_log
    log_ci_lower <- cs$log_lower
    log_ci_upper <- cs$log_upper
    estimate <- exp(log_estimate)
    ci_lower <- exp(log_ci_lower)
    ci_upper <- exp(log_ci_upper)
    contrast_raw <- cs$raw
    extreme_predictions <- make_extreme_predictions(pred_values, positions)
    estimation_engine <- "emmeans log-link contrast with offset = 0; RII = exp(log-rate contrast)"

  } else if (model_name == "linear_log") {

    # For transformed linear-log models:
    #   log(Y_R) = alpha + beta * R
    # Therefore:
    #   log(RII) = log(Y_A) - log(Y_D)
    # This can be estimated directly on the fitted transformed scale.

    emm_log <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      vcov. = V
    )
    ctr <- emmeans::contrast(emm_log, method = list(log_RII = c(-1, 1)))
    cs <- extract_link_contrast(ctr, conf_level)

    emm_sum <- as.data.frame(summary(emm_log, infer = c(FALSE, FALSE)))
    em_col <- intersect(c("emmean", "response"), names(emm_sum))[1L]
    if (is.na(em_col)) {
      stop("Could not extract log-scale predictions from `emmeans`.", call. = FALSE)
    }

    theta_log <- as.numeric(emm_sum[[em_col]])
    if (length(theta_log) != 2L || any(!is.finite(theta_log))) {
      stop(
        "Invalid log-scale predictions for the selected `linear_log` model. ",
        "Action: inspect selected model diagnostics.",
        call. = FALSE
      )
    }

    pred_values <- exp(theta_log)
    require_positive_extremes(pred_values, model_name)

    log_estimate <- cs$log_estimate
    se_log <- cs$se_log
    log_ci_lower <- cs$log_lower
    log_ci_upper <- cs$log_upper
    estimate <- exp(log_estimate)
    ci_lower <- exp(log_ci_lower)
    ci_upper <- exp(log_ci_upper)
    contrast_raw <- cs$raw
    extreme_predictions <- make_extreme_predictions(pred_values, positions)
    estimation_engine <- "emmeans contrast on fitted log scale; RII = exp(log-response contrast)"

  } else if (model_name %in% c("binomial", "quasibinomial")) {

    # Critical RII rule:
    # For logit-link models, exp(link-scale contrast) is an odds ratio.
    # RII is NOT the odds ratio. RII is:
    #   p_advantaged / p_disadvantaged
    # Therefore we obtain response-scale predicted probabilities and compute
    # the log ratio of probabilities using a package-based delta method.

    emm_link <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      type = "link",
      vcov. = V
    )
    emm_resp <- emmeans::regrid(emm_link, transform = "response")

    pred_prob <- extract_emmeans_prediction(
      emm_resp,
      candidate_names = c("prob", "response", "emmean")
    )
    require_positive_extremes(pred_prob, model_name)

    V_prob <- get_prediction_vcov(emm_resp)
    se_log <- delta_log_ratio(pred_prob, V_prob, "binomial/quasibinomial probability ratios")

    log_estimate <- log(pred_prob[2L] / pred_prob[1L])
    log_ci_lower <- log_estimate - z_value * se_log
    log_ci_upper <- log_estimate + z_value * se_log
    estimate <- exp(log_estimate)
    ci_lower <- exp(log_ci_lower)
    ci_upper <- exp(log_ci_upper)

    scale_factor <- if (health_indicator_type == "percentage") 100 else 1
    pred_values <- pred_prob * scale_factor

    contrast_raw <- data.frame(
      contrast = "log_RII",
      estimate = log_estimate,
      SE = se_log,
      lower.CL = log_ci_lower,
      upper.CL = log_ci_upper,
      response_ratio = estimate,
      stringsAsFactors = FALSE
    )
    extreme_predictions <- make_extreme_predictions(pred_values, positions)
    estimation_engine <- paste(
      "emmeans response-scale probabilities + msm::deltamethod for log(p_advantaged / p_disadvantaged);",
      "not an odds-ratio contrast"
    )

  } else if (model_name == "linear_identity") {

    # Identity-scale models can predict non-positive values.
    # That may be tolerable for an absolute difference, but it is not tolerable
    # for RII. We therefore enforce strict positivity at both extremes.

    emm <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      vcov. = V
    )

    pred_values <- extract_emmeans_prediction(
      emm,
      candidate_names = c("emmean", "response")
    )
    require_positive_extremes(pred_values, model_name)

    V_pred <- get_prediction_vcov(emm)
    se_log <- delta_log_ratio(pred_values, V_pred, "identity-scale predicted-value ratios")

    log_estimate <- log(pred_values[2L] / pred_values[1L])
    log_ci_lower <- log_estimate - z_value * se_log
    log_ci_upper <- log_estimate + z_value * se_log
    estimate <- exp(log_estimate)
    ci_lower <- exp(log_ci_lower)
    ci_upper <- exp(log_ci_upper)

    contrast_raw <- data.frame(
      contrast = "log_RII",
      estimate = log_estimate,
      SE = se_log,
      lower.CL = log_ci_lower,
      upper.CL = log_ci_upper,
      response_ratio = estimate,
      stringsAsFactors = FALSE
    )
    extreme_predictions <- make_extreme_predictions(pred_values, positions)
    estimation_engine <- "emmeans identity-scale predictions + msm::deltamethod for log(Y_advantaged / Y_disadvantaged)"

  } else {
    stop(
      "Unsupported selected model for RII estimation: ", model_name, ". ",
      "Action: inspect candidate model validation and fitting outputs.",
      call. = FALSE
    )
  }

  if (!is.finite(estimate) || estimate <= 0 ||
      !is.finite(log_estimate) ||
      !is.finite(se_log) || se_log < 0 ||
      !is.finite(ci_lower) || ci_lower <= 0 ||
      !is.finite(ci_upper) || ci_upper <= 0 ||
      !is.finite(log_ci_lower) || !is.finite(log_ci_upper)) {
    stop(
      "The Wald RII estimate or interval contains non-finite or non-positive values. ",
      "Action: inspect the selected model, consider a different candidate model, ",
      "or use bootstrap/jackknife once available.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Assemble output.
  # ---------------------------------------------------------------------------

  results_overall <- data.frame(
    Metric = "Relative Index of Inequality",
    Estimate = estimate,
    Log_Estimate = log_estimate,
    Standard_Error = se_log,
    Standard_Error_Scale = "log",
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    Log_CI_Lower = log_ci_lower,
    Log_CI_Upper = log_ci_upper,
    Confidence_Level = conf_level,
    CI_Method = "wald",
    Interval_Scale = "log",
    Selected_Model = model_name,
    Health_Indicator_Type = health_indicator_type,
    Input_Source = input_source,
    Social_Position_Method = social_position_method,
    Disadvantaged_Position = positions$disadvantaged,
    Advantaged_Position = positions$advantaged,
    Disadvantaged_Predicted_Value = extreme_predictions$Predicted_Value[1L],
    Advantaged_Predicted_Value = extreme_predictions$Predicted_Value[2L],
    stringsAsFactors = FALSE
  )

  metadata <- c(
    selected$metadata,
    list(
      measure = "rii",
      estimand = "RII = Y_advantaged / Y_disadvantaged",
      estimand_type = "relative_ratio",
      estimation_engine = estimation_engine,
      contrast_definition = "RII = Y_advantaged / Y_disadvantaged; log_RII = log(Y_advantaged) - log(Y_disadvantaged)",
      covariance_type_requested = vcov_type,
      robust_covariance_used = vc$robust_used,
      covariance_source = if (vc$robust_used) {
        paste0("sandwich::vcovHC(type = '", vcov_type, "')")
      } else {
        "stats::vcov()"
      },
      wald_interval_rule = "Wald confidence intervals for RII are computed on the log-ratio scale and exponentiated.",
      binomial_model_rule = paste(
        "For binomial/quasibinomial models, RII is computed as a ratio of response-scale predicted probabilities;",
        "the logit-scale exponentiated contrast is an odds ratio and is not used as RII."
      ),
      positivity_rule = "RII requires strictly positive predicted values at both disadvantaged and advantaged social-position extremes."
    )
  )

  diagnostics <- c(
    selected$diagnostics,
    list(
      estimate_warnings = unique(warnings),
      package_based_wald = TRUE,
      manual_gradient_used = FALSE,
      emmeans_required = TRUE,
      msm_required_for_response_scale_log_ratio = model_name %in% c("binomial", "quasibinomial", "linear_identity"),
      extreme_prediction_order = "row 1 = disadvantaged, row 2 = advantaged",
      positivity_check = "Both extreme predicted values are finite and strictly positive.",
      sign_and_ratio_check = "Estimate equals Advantaged_Predicted_Value / Disadvantaged_Predicted_Value.",
      odds_ratio_warning = if (model_name %in% c("binomial", "quasibinomial")) {
        "RII is a probability ratio, not an odds ratio; logit-scale exponentiated contrasts are not used."
      } else {
        NA_character_
      }
    )
  )

  list(
    results_overall = results_overall,
    extreme_predictions = extreme_predictions,
    contrast = contrast_raw,
    covariance_matrix = V,
    selected_model_name = model_name,
    selected_model = model,
    model_selection_table = selected$model_selection_table,
    predictions = selected$predictions,
    prepared_data = d,
    metadata = metadata,
    diagnostics = diagnostics
  )
}
