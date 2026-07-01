#'     intended back-transformation.
#'
#' The contrast is always explicit: c(-1, 1), meaning:
#'   -1 * Y_disadvantaged + 1 * Y_advantaged.
#'
#' @keywords internal
.sii_rii_estimate_sii <- function(selected, validation, conf_level = NULL, vcov_type = "HC1") {

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
      "Action: call `.sii_rii_select_model()` before `.sii_rii_estimate_sii()`.",
      call. = FALSE
    )
  }

  required_validation_fields <- c(
    "health_indicator_type",
    "input_source",
    "social_position_method",
    "conf_level",
    "estimand"
  )
  missing_validation_fields <- setdiff(required_validation_fields, names(validation))
  if (length(missing_validation_fields) > 0L) {
    stop(
      "Invalid `validation` object: missing field(s): ",
      paste(missing_validation_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_validate_inputs()` before `.sii_rii_estimate_sii()`.",
      call. = FALSE
    )
  }

  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop(
      "Package `emmeans` is required to compute the Wald SII contrast in this version. ",
      "Action: install `emmeans` or use a resampling interval once `sii_resample.R` is available.",
      call. = FALSE
    )
  }

  model_name <- as.character(selected$selected_model_name)
  model <- selected$selected_model
  d <- selected$prepared_data

  if (length(model_name) != 1L || is.na(model_name) || !nzchar(model_name)) {
    stop(
      "Invalid selected model name. ",
      "Action: inspect `.sii_rii_select_model()` output.",
      call. = FALSE
    )
  }
  if (is.null(model)) {
    stop(
      "Invalid selected model: selected model object is NULL. ",
      "Action: inspect `.sii_rii_fit_models()` and `.sii_rii_select_model()` outputs.",
      call. = FALSE
    )
  }
  if (!is.data.frame(d) || nrow(d) < 2L) {
    stop(
      "Invalid prepared data in selected object: at least two ecological units are required. ",
      "Action: inspect `.sii_rii_prepare_data()` output.",
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
      "Action: rebuild the workflow from `.sii_rii_prepare_data()`.",
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
        "Cannot compute bounded SII positions because fewer than two finite social positions are available. ",
        "Action: inspect `.sii_rii_prepare_data()` output.",
        call. = FALSE
      )
    }

    list(
      disadvantaged = min(observed),
      advantaged = max(observed)
    )
  }

  extract_contrast_summary <- function(contrast_object, level) {
    out <- as.data.frame(summary(contrast_object, infer = c(TRUE, TRUE), level = level))

    estimate_col <- intersect(c("estimate", "response", "emmean"), names(out))[1]
    se_col <- intersect(c("SE", "se"), names(out))[1]
    lower_col <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(out))[1]
    upper_col <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(out))[1]

    if (any(is.na(c(estimate_col, se_col, lower_col, upper_col)))) {
      stop(
        "Could not extract estimate, SE, and confidence limits from the `emmeans` contrast object. ",
        "Action: inspect the selected model and `emmeans` output.",
        call. = FALSE
      )
    }

    list(
      estimate = as.numeric(out[[estimate_col]][1L]),
      se = as.numeric(out[[se_col]][1L]),
      lower = as.numeric(out[[lower_col]][1L]),
      upper = as.numeric(out[[upper_col]][1L]),
      raw = out
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

  # ---------------------------------------------------------------------------
  # 3. Determine disadvantaged and advantaged positions.
  # ---------------------------------------------------------------------------

  positions <- get_positions(d, social_position_method)

  if (!is.finite(positions$disadvantaged) || !is.finite(positions$advantaged) ||
      positions$disadvantaged >= positions$advantaged) {
    stop(
      "Invalid SII positions: disadvantaged position must be less than advantaged position. ",
      "Action: inspect social position construction in `.sii_rii_prepare_data()`.",
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
  # 5. Estimate extremes and Wald SE/CI using package-based contrasts.
  # ---------------------------------------------------------------------------

  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha / 2)

  estimation_engine <- NA_character_
  contrast_raw <- NULL
  extreme_predictions <- NULL
  estimate <- NA_real_
  se <- NA_real_
  ci_lower <- NA_real_
  ci_upper <- NA_real_

  if (model_name %in% c("binomial", "quasibinomial")) {

    emm <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      type = "link",
      vcov. = V
    )
    emm_resp <- emmeans::regrid(emm, transform = "response")
    ctr <- emmeans::contrast(emm_resp, method = list(SII = c(-1, 1)))
    cs <- extract_contrast_summary(ctr, conf_level)

    scale_factor <- if (health_indicator_type == "percentage") 100 else 1
    extreme_df <- as.data.frame(summary(emm_resp, infer = c(FALSE, FALSE)))
    pred_col <- intersect(c("prob", "response", "emmean"), names(extreme_df))[1]
    if (is.na(pred_col)) {
      stop("Could not extract response-scale predictions from `emmeans`.", call. = FALSE)
    }
    pred_values <- as.numeric(extreme_df[[pred_col]]) * scale_factor

    estimate <- cs$estimate * scale_factor
    se <- cs$se * scale_factor
    ci_lower <- cs$lower * scale_factor
    ci_upper <- cs$upper * scale_factor
    contrast_raw <- cs$raw
    extreme_predictions <- make_extreme_predictions(pred_values, positions)
    estimation_engine <- "emmeans::emmeans + emmeans::regrid + explicit contrast c(-1, 1)"

  } else if (model_name %in% c("poisson", "quasipoisson", "negative_binomial")) {

    emm <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      type = "link",
      offset = 0,
      vcov. = V
    )
    emm_resp <- emmeans::regrid(emm, transform = "response")
    ctr <- emmeans::contrast(emm_resp, method = list(SII = c(-1, 1)))
    cs <- extract_contrast_summary(ctr, conf_level)

    extreme_df <- as.data.frame(summary(emm_resp, infer = c(FALSE, FALSE)))
    pred_col <- intersect(c("rate", "response", "emmean"), names(extreme_df))[1]
    if (is.na(pred_col)) {
      stop("Could not extract response-scale rate predictions from `emmeans`.", call. = FALSE)
    }
    pred_values <- as.numeric(extreme_df[[pred_col]])

    estimate <- cs$estimate
    se <- cs$se
    ci_lower <- cs$lower
    ci_upper <- cs$upper
    contrast_raw <- cs$raw
    extreme_predictions <- make_extreme_predictions(pred_values, positions)
    estimation_engine <- "emmeans::emmeans with offset = 0 + emmeans::regrid + explicit contrast c(-1, 1)"

  } else if (model_name == "linear_identity") {

    emm <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      vcov. = V
    )
    ctr <- emmeans::contrast(emm, method = list(SII = c(-1, 1)))
    cs <- extract_contrast_summary(ctr, conf_level)

    extreme_df <- as.data.frame(summary(emm, infer = c(FALSE, FALSE)))
    pred_col <- intersect(c("emmean", "response"), names(extreme_df))[1]
    if (is.na(pred_col)) {
      stop("Could not extract identity-scale predictions from `emmeans`.", call. = FALSE)
    }
    pred_values <- as.numeric(extreme_df[[pred_col]])

    estimate <- cs$estimate
    se <- cs$se
    ci_lower <- cs$lower
    ci_upper <- cs$upper
    contrast_raw <- cs$raw
    extreme_predictions <- make_extreme_predictions(pred_values, positions)
    estimation_engine <- "emmeans::emmeans + explicit contrast c(-1, 1)"

  } else if (model_name == "linear_log") {

    if (!requireNamespace("msm", quietly = TRUE)) {
      stop(
        "Package `msm` is required to compute the Wald interval for transformed linear models ",
        "(`linear_log`) using `msm::deltamethod()`. ",
        "Action: install `msm`, choose a non-transformed model, or use bootstrap/jackknife once available.",
        call. = FALSE
      )
    }

    emm <- emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = at_positions,
      vcov. = V
    )
    emm_sum <- as.data.frame(summary(emm, infer = c(FALSE, FALSE)))
    em_col <- intersect(c("emmean", "response"), names(emm_sum))[1]
    if (is.na(em_col)) {
      stop("Could not extract transformed-scale predictions from `emmeans`.", call. = FALSE)
    }

    theta <- as.numeric(emm_sum[[em_col]])
    if (length(theta) != 2L || any(!is.finite(theta))) {
      stop(
        "Invalid transformed-scale predictions for the selected transformed model. ",
        "Action: inspect selected model diagnostics.",
        call. = FALSE
      )
    }

    V_theta <- tryCatch(
      as.matrix(stats::vcov(emm)),
      error = function(e) {
        stop(
          "Could not extract the covariance matrix of transformed-scale predictions from `emmeans` (",
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

    pred_values <- exp(theta)
    estimate <- pred_values[2L] - pred_values[1L]
    se <- as.numeric(msm::deltamethod(~ exp(x2) - exp(x1), theta, V_theta))
    estimation_engine <- "emmeans::emmeans on log scale + msm::deltamethod for exp(x2) - exp(x1)"

    ci_lower <- estimate - z_value * se
    ci_upper <- estimate + z_value * se
    contrast_raw <- data.frame(
      contrast = "SII",
      estimate = estimate,
      SE = se,
      lower.CL = ci_lower,
      upper.CL = ci_upper,
      stringsAsFactors = FALSE
    )
    extreme_predictions <- make_extreme_predictions(pred_values, positions)

  } else {
    stop(
      "Unsupported selected model for SII estimation: ", model_name, ". ",
      "Action: inspect candidate model validation and fitting outputs.",
      call. = FALSE
    )
  }

  if (!is.finite(estimate) || !is.finite(se) || se < 0 ||
      !is.finite(ci_lower) || !is.finite(ci_upper)) {
    stop(
      "The Wald SII estimate or interval contains non-finite values. ",
      "Action: inspect the selected model, consider a different candidate model, or use bootstrap/jackknife once available.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Assemble output.
  # ---------------------------------------------------------------------------

  results_overall <- data.frame(
    Metric = "Slope Index of Inequality",
    Estimate = estimate,
    Standard_Error = se,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    Confidence_Level = conf_level,
    CI_Method = "wald",
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
      estimand = "SII = Y_advantaged - Y_disadvantaged",
      estimation_engine = estimation_engine,
      contrast_definition = "Explicit contrast c(-1, 1): -Y_disadvantaged + Y_advantaged",
      covariance_type_requested = vcov_type,
      robust_covariance_used = vc$robust_used,
      covariance_source = if (vc$robust_used) {
        paste0("sandwich::vcovHC(type = '", vcov_type, "')")
      } else {
        "stats::vcov()"
      },
      transformed_linear_model_rule = paste(
        "For linear_log, emmeans estimates the two extremes on the fitted transformed scale;",
        "msm::deltamethod() is used for the nonlinear back-transformed SII contrast."
      )
    )
  )

  diagnostics <- c(
    selected$diagnostics,
    list(
      estimate_warnings = unique(warnings),
      package_based_wald = TRUE,
      manual_gradient_used = FALSE,
      emmeans_required = TRUE,
      msm_required_for_transformed_linear_models = model_name == "linear_log",
      extreme_prediction_order = "row 1 = disadvantaged, row 2 = advantaged",
      sign_check = "Estimate equals Advantaged_Predicted_Value - Disadvantaged_Predicted_Value"
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
