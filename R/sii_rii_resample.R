#' Iterative Variance Estimation Engine for SII/RII Families
#'
#' @description
#' `sii_rii_resample` governs the stochastic non-parametric replication topology
#' required for bootstrapped and jackknifed standard errors.
#'
#' The engine operates recursively across the primary ecological unit vector.
#' For each iteration, it reconstructs the fractional rank (ridit) vector and refits
#' the optimal selected regression distribution (`glm`, `lm`, or `glm.nb`), yielding
#' non-parametric bounds immune to complex error-structure misspecifications.
#'
#' Shared convention:
#'   social_position = 0 -> most disadvantaged
#'   social_position = 1 -> most advantaged
#'
#' SII estimand:
#'   SII = Y_advantaged - Y_disadvantaged
#'
#' RII estimand:
#'   RII = Y_advantaged / Y_disadvantaged
#'
#' For RII, endpoint predictions must be strictly positive. Bootstrap intervals
#' are percentile intervals on the RII scale. Jackknife intervals are calculated
#' on the log-RII scale and exponentiated, to avoid impossible negative ratio
#' intervals.
#'
#' @importFrom stats qnorm aggregate glm as.formula binomial quasibinomial poisson quasipoisson lm predict sd quantile
#' @keywords internal
#' @noRd
.sii_rii_resample <- function(selected,
                              validation,
                              point_estimate = NULL,
                              ci_method = NULL,
                              n_boot = NULL,
                              conf_level = NULL,
                              seed = NULL) {

  # ---------------------------------------------------------------------------
  # 1. Validate upstream contracts.
  # ---------------------------------------------------------------------------

  required_selected_fields <- c(
    "selected_model_name",
    "selected_model",
    "prepared_data",
    "metadata",
    "diagnostics"
  )
  missing_selected_fields <- setdiff(required_selected_fields, names(selected))
  if (length(missing_selected_fields) > 0L) {
    stop(
      "Invalid `selected` object: missing field(s): ",
      paste(missing_selected_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_select_model()` before `.sii_rii_resample()`.",
      call. = FALSE
    )
  }

  required_validation_fields <- c(
    "health_indicator_type",
    "input_source",
    "social_position_method",
    "ci_method",
    "n_boot",
    "conf_level",
    "seed",
    "estimand"
  )
  missing_validation_fields <- setdiff(required_validation_fields, names(validation))
  if (length(missing_validation_fields) > 0L) {
    stop(
      "Invalid `validation` object: missing field(s): ",
      paste(missing_validation_fields, collapse = ", "), ". ",
      "Action: call `.sii_rii_validate_inputs()` before `.sii_rii_resample()`.",
      call. = FALSE
    )
  }

  measure <- validation$measure
  if (is.null(measure)) {
    if (!is.null(validation$estimand) &&
        grepl("^RII", validation$estimand, ignore.case = TRUE)) {
      measure <- "rii"
    } else {
      measure <- "sii"
    }
  }
  measure <- match.arg(as.character(measure), choices = c("sii", "rii"))

  measure_label <- if (measure == "sii") {
    "Slope Index of Inequality"
  } else {
    "Relative Index of Inequality"
  }

  estimand_text <- if (measure == "sii") {
    "SII = Y_advantaged - Y_disadvantaged"
  } else {
    "RII = Y_advantaged / Y_disadvantaged"
  }

  if (is.null(ci_method)) {
    ci_method <- validation$ci_method
  }
  ci_method <- as.character(ci_method)
  if (length(ci_method) != 1L || !ci_method %in% c("bootstrap", "jackknife")) {
    stop(
      "`.sii_rii_resample()` only handles `ci_method = 'bootstrap'` or `ci_method = 'jackknife'`. ",
      "Action: use the measure-specific estimate helper for Wald intervals, or call `.sii_rii_resample()` only for resampling methods.",
      call. = FALSE
    )
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

  if (is.null(n_boot)) {
    n_boot <- validation$n_boot
  }
  if (ci_method == "bootstrap") {
    if (!is.numeric(n_boot) || length(n_boot) != 1L ||
        is.na(n_boot) || !is.finite(n_boot) ||
        n_boot < 100 || abs(n_boot - round(n_boot)) > sqrt(.Machine$double.eps)) {
      stop(
        "`n_boot` must be an integer-like numeric value >= 100 when `ci_method = 'bootstrap'`. ",
        "Action: use a value such as 2000 or 10000.",
        call. = FALSE
      )
    }
    n_boot <- as.integer(round(n_boot))
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop(
        "Package `boot` is required for `ci_method = 'bootstrap'`. ",
        "Action: install it with `install.packages('boot')`, or use `ci_method = 'wald'` or `ci_method = 'jackknife'`.",
        call. = FALSE
      )
    }
  }

  if (ci_method == "jackknife" &&
      !requireNamespace("bootstrap", quietly = TRUE)) {
    stop(
      "Package `bootstrap` is required for `ci_method = 'jackknife'` in this version. ",
      "Action: install it with `install.packages('bootstrap')`, or use `ci_method = 'wald'` or `ci_method = 'bootstrap'`.",
      call. = FALSE
    )
  }

  if (is.null(seed)) {
    seed <- validation$seed
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L ||
        is.na(seed) || !is.finite(seed) || seed < 0 ||
        abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
      stop(
        "`seed` must be NULL or a single non-negative integer-like numeric value. ",
        "Action: use a value such as 123, or leave it as NULL.",
        call. = FALSE
      )
    }
    seed <- as.integer(round(seed))
  }

  alpha <- 1 - conf_level
  z_value <- stats::qnorm(1 - alpha / 2)

  model_name <- selected$selected_model_name
  d0 <- selected$prepared_data

  required_data_columns <- c(
    "unit_id",
    "equity_order",
    "population_weight",
    "health_indicator_internal",
    "health_indicator_output",
    "model_response",
    "social_position"
  )
  missing_data_columns <- setdiff(required_data_columns, names(d0))
  if (length(missing_data_columns) > 0L) {
    stop(
      "Invalid `selected$prepared_data`: missing column(s): ",
      paste(missing_data_columns, collapse = ", "), ". ",
      "Action: call `.sii_rii_prepare_data()` and `.sii_rii_fit_models()` before resampling.",
      call. = FALSE
    )
  }

  health_indicator_type <- validation$health_indicator_type
  input_source <- validation$input_source
  social_position_method <- validation$social_position_method

  # ---------------------------------------------------------------------------
  # 2. Local helpers.
  # ---------------------------------------------------------------------------

  add_warning <- function(warnings, message) {
    unique(c(warnings, message))
  }

  refit_ridit <- function(d) {
    # Reconstruct the relative social position in the resampled data.
    # Duplicated units in bootstrap samples remain duplicated rows and jointly
    # contribute their resampled population/exposure weights.
    total_weight <- sum(d$population_weight)
    if (!is.finite(total_weight) || total_weight <= 0) {
      stop("Resampled total population/exposure weight is not positive.", call. = FALSE)
    }

    ord <- order(d$equity_order, d$unit_id, seq_len(nrow(d)))
    d <- d[ord, , drop = FALSE]
    d$original_resample_row <- seq_len(nrow(d))

    block <- stats::aggregate(
      population_weight ~ equity_order,
      data = d,
      FUN = sum
    )
    block <- block[order(block$equity_order), , drop = FALSE]
    block$cumulative_weight_before_block <- c(0, head(cumsum(block$population_weight), -1L))
    block$cumulative_weight_after_block <- cumsum(block$population_weight)
    block$social_position <- (
      block$cumulative_weight_before_block + 0.5 * block$population_weight
    ) / total_weight

    d$social_position <- block$social_position[
      match(d$equity_order, block$equity_order)
    ]

    d <- d[order(d$original_resample_row), , drop = FALSE]
    d$original_resample_row <- NULL

    d
  }

  endpoint_positions <- function(d) {
    if (identical(social_position_method, "classic")) {
      c(disadvantaged = 0, advantaged = 1)
    } else {
      c(
        disadvantaged = min(d$social_position),
        advantaged = max(d$social_position)
      )
    }
  }

  fit_one_model <- function(d, model_name) {
    model <- NULL

    if (model_name == "binomial") {
      if (input_source != "counts" ||
          !health_indicator_type %in% c("proportion", "percentage")) {
        stop("`binomial` is only valid for proportion/percentage with counts.", call. = FALSE)
      }
      model <- stats::glm(
        stats::as.formula("cbind(numerator, count_failures) ~ social_position"),
        data = d,
        family = stats::binomial()
      )
    } else if (model_name == "quasibinomial") {
      if (!health_indicator_type %in% c("proportion", "percentage")) {
        stop("`quasibinomial` is only valid for proportion/percentage outcomes.", call. = FALSE)
      }
      if (input_source == "counts") {
        model <- stats::glm(
          stats::as.formula("cbind(numerator, count_failures) ~ social_position"),
          data = d,
          family = stats::quasibinomial()
        )
      } else {
        model <- stats::glm(
          stats::as.formula("health_indicator_internal ~ social_position"),
          data = d,
          family = stats::quasibinomial(),
          weights = population_weight
        )
      }
    } else if (model_name == "poisson") {
      if (input_source != "counts" ||
          !health_indicator_type %in% c("rate", "ratio")) {
        stop("`poisson` is only valid for rate/ratio with counts.", call. = FALSE)
      }
      model <- stats::glm(
        stats::as.formula("model_response ~ social_position + offset(model_offset)"),
        data = d,
        family = stats::poisson()
      )
    } else if (model_name == "quasipoisson") {
      if (input_source != "counts" ||
          !health_indicator_type %in% c("rate", "ratio")) {
        stop("`quasipoisson` is only valid for rate/ratio with counts.", call. = FALSE)
      }
      model <- stats::glm(
        stats::as.formula("model_response ~ social_position + offset(model_offset)"),
        data = d,
        family = stats::quasipoisson()
      )
    } else if (model_name == "negative_binomial") {
      if (input_source != "counts" ||
          !health_indicator_type %in% c("rate", "ratio")) {
        stop("`negative_binomial` is only valid for rate/ratio with counts.", call. = FALSE)
      }
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop(
          "Package `MASS` is required to refit `negative_binomial` during resampling. ",
          "Action: install `MASS`, or select a different model.",
          call. = FALSE
        )
      }
      model <- MASS::glm.nb(
        stats::as.formula("model_response ~ social_position + offset(model_offset)"),
        data = d
      )
    } else if (model_name == "linear_log") {
      if (any(!is.finite(d$log_response))) {
        stop("`linear_log` cannot be refitted because at least one resampled outcome is not strictly positive.", call. = FALSE)
      }
      model <- stats::lm(
        stats::as.formula("log_response ~ social_position"),
        data = d,
        weights = population_weight
      )
    } else if (model_name == "linear_identity") {
      model <- stats::lm(
        stats::as.formula("health_indicator_output ~ social_position"),
        data = d,
        weights = population_weight
      )
    } else {
      stop(
        "Unsupported selected model for resampling: `", model_name, "`. ",
        "Action: inspect candidate-model validation and model-selection output.",
        call. = FALSE
      )
    }

    model
  }

  predict_endpoint <- function(model, model_name, position) {
    nd <- data.frame(social_position = position)

    if (model_name %in% c("binomial", "quasibinomial")) {
      p <- stats::predict(model, newdata = nd, type = "response")
      value <- if (health_indicator_type == "percentage") 100 * p else p
    } else if (model_name %in% c("poisson", "quasipoisson", "negative_binomial")) {
      # With offset = 0, type = "response" returns the modeled scaled
      # rate/ratio, not a count for an arbitrary denominator.
      nd$model_offset <- 0
      value <- stats::predict(model, newdata = nd, type = "response")
    } else if (model_name == "linear_log") {
      value <- exp(stats::predict(model, newdata = nd))
    } else if (model_name == "linear_identity") {
      value <- stats::predict(model, newdata = nd)
    } else {
      stop("Unsupported model in endpoint prediction.", call. = FALSE)
    }

    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) {
      stop("Endpoint prediction is not a single finite numeric value.", call. = FALSE)
    }
    value
  }

  validate_endpoint_prediction <- function(value) {
    if (!is.finite(value)) {
      return(FALSE)
    }
    if (health_indicator_type == "proportion") {
      return(value >= 0 && value <= 1)
    }
    if (health_indicator_type == "percentage") {
      return(value >= 0 && value <= 100)
    }
    if (health_indicator_type %in% c("rate", "ratio")) {
      return(value >= 0)
    }
    FALSE
  }

  estimand_from_prepared_resample <- function(d_resampled, return_log = FALSE) {
    d_resampled <- refit_ridit(d_resampled)
    positions <- endpoint_positions(d_resampled)

    model <- fit_one_model(d_resampled, model_name)

    y_disadvantaged <- predict_endpoint(
      model = model,
      model_name = model_name,
      position = unname(positions["disadvantaged"])
    )
    y_advantaged <- predict_endpoint(
      model = model,
      model_name = model_name,
      position = unname(positions["advantaged"])
    )

    if (!validate_endpoint_prediction(y_disadvantaged) ||
        !validate_endpoint_prediction(y_advantaged)) {
      stop("At least one endpoint prediction is outside the valid outcome range.", call. = FALSE)
    }

    if (measure == "sii") {
      value <- y_advantaged - y_disadvantaged
      if (return_log) {
        return(NA_real_)
      }
      return(value)
    }

    # RII requires strict positivity at both endpoints. This is checked here,
    # not in input validation, because endpoint predictions depend on the
    # selected and refitted model.
    if (y_disadvantaged <= 0 || y_advantaged <= 0) {
      stop(
        "RII cannot be computed because at least one endpoint prediction is not strictly positive.",
        call. = FALSE
      )
    }

    log_value <- log(y_advantaged) - log(y_disadvantaged)
    if (return_log) {
      return(log_value)
    }
    exp(log_value)
  }

  warnings <- character(0)


  extract_point_estimate <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }

    if (is.list(x) && !is.data.frame(x) && !is.null(x$results_overall)) {
      x <- x$results_overall
    }

    if (is.data.frame(x)) {
      if (!"Estimate" %in% names(x)) {
        stop(
          "`point_estimate` was supplied as a data frame but does not contain an `Estimate` column. ",
          "Action: pass a single numeric estimate or the full estimate object returned by the measure-specific estimate helper.",
          call. = FALSE
        )
      }
      x <- x$Estimate[1L]
    }

    if (is.list(x)) {
      stop(
        "`point_estimate` was supplied as a list but the numeric estimate could not be extracted. ",
        "Action: pass `estimate$results_overall$Estimate[1L]`, the full estimate object returned by the measure-specific estimate helper, or leave `point_estimate = NULL`.",
        call. = FALSE
      )
    }

    suppressWarnings(as.numeric(x))
  }

  if (is.null(point_estimate)) {
    point_estimate <- tryCatch(
      estimand_from_prepared_resample(d0),
      error = function(e) NA_real_
    )
    if (!is.finite(point_estimate)) {
      stop(
        "Could not compute the observed ", toupper(measure),
        " estimate through the resampling estimator. ",
        "Action: inspect the selected model, transformed outcomes, and endpoint predictions.",
        call. = FALSE
      )
    }
    warnings <- add_warning(
      warnings,
      "`point_estimate` was not supplied; it was recomputed using the resampling estimator."
    )
  } else {
    point_estimate <- extract_point_estimate(point_estimate)
    if (length(point_estimate) != 1L || !is.finite(point_estimate)) {
      stop(
        "`point_estimate` must be NULL, a single finite numeric value, ",
        "or an estimate object containing `results_overall$Estimate`. ",
        "Action: pass `estimate$results_overall$Estimate[1L]`, the full estimate object, or leave `point_estimate = NULL`.",
        call. = FALSE
      )
    }
  }

  point_log_estimate <- NA_real_
  if (measure == "rii") {
    if (point_estimate <= 0) {
      stop(
        "`point_estimate` must be strictly positive for RII. ",
        "Action: inspect endpoint predictions or leave `point_estimate = NULL`.",
        call. = FALSE
      )
    }
    point_log_estimate <- log(point_estimate)
  }

  # ---------------------------------------------------------------------------
  # 3. Bootstrap ecological-unit resampling.
  # ---------------------------------------------------------------------------

  if (ci_method == "bootstrap") {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    statistic <- function(data, indices) {
      d_b <- data[indices, , drop = FALSE]
      value <- tryCatch(
        estimand_from_prepared_resample(d_b),
        error = function(e) NA_real_,
        warning = function(w) NA_real_
      )
      value <- suppressWarnings(as.numeric(value))
      if (length(value) != 1L || !is.finite(value)) NA_real_ else value
    }

    boot_object <- boot::boot(
      data = d0,
      statistic = statistic,
      R = n_boot
    )

    replicates <- suppressWarnings(as.numeric(boot_object$t[, 1]))
    valid <- is.finite(replicates)
    if (measure == "rii") {
      valid <- valid & replicates > 0
    }
    n_valid <- sum(valid)

    if (n_valid < max(30L, ceiling(0.50 * n_boot))) {
      stop(
        "Too few valid bootstrap replicates were obtained. ",
        "Valid replicates: ", n_valid, " of ", n_boot, ". ",
        "Action: inspect model convergence and endpoint positivity; consider a simpler model or use `ci_method = 'wald'`.",
        call. = FALSE
      )
    }

    if (n_valid < n_boot) {
      warnings <- add_warning(
        warnings,
        paste0(
          "Some bootstrap replicates failed and were excluded from the interval calculation: ",
          n_boot - n_valid, " of ", n_boot, "."
        )
      )
    }

    valid_values <- replicates[valid]
    se <- stats::sd(valid_values)
    ci <- stats::quantile(
      valid_values,
      probs = c(alpha / 2, 1 - alpha / 2),
      names = FALSE,
      na.rm = TRUE,
      type = 6
    )
    ci_lower <- ci[1]
    ci_upper <- ci[2]

    log_estimate <- NA_real_
    se_log <- NA_real_
    log_ci_lower <- NA_real_
    log_ci_upper <- NA_real_
    interval_scale <- "estimand_scale"
    interval_type <- "percentile_ecological_unit_bootstrap"

    if (measure == "rii") {
      log_estimate <- point_log_estimate
      log_values <- log(valid_values)
      se_log <- stats::sd(log_values)
      log_ci <- stats::quantile(
        log_values,
        probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE,
        na.rm = TRUE,
        type = 6
      )
      log_ci_lower <- log_ci[1]
      log_ci_upper <- log_ci[2]
      interval_scale <- "ratio_scale_percentile_with_log_scale_diagnostics"
    }

    results_replicates <- data.frame(
      Replicate = seq_along(replicates),
      Replicate_Type = "bootstrap_ecological_unit",
      Estimate = replicates,
      Log_Estimate = if (measure == "rii") ifelse(replicates > 0, log(replicates), NA_real_) else NA_real_,
      Status = ifelse(valid, "valid", "failed"),
      Message = ifelse(valid, NA_character_, "Bootstrap sample failed to refit model or produce valid endpoint predictions."),
      stringsAsFactors = FALSE
    )

    results_overall <- data.frame(
      Metric = measure_label,
      Estimate = point_estimate,
      Log_Estimate = log_estimate,
      Standard_Error = se,
      Standard_Error_Scale = if (measure == "rii") "ratio_scale" else "difference_scale",
      Standard_Error_Log = se_log,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      Log_CI_Lower = log_ci_lower,
      Log_CI_Upper = log_ci_upper,
      Confidence_Level = conf_level,
      CI_Method = "bootstrap",
      Interval_Type = interval_type,
      Interval_Scale = interval_scale,
      Selected_Model = model_name,
      Health_Indicator_Type = health_indicator_type,
      Input_Source = input_source,
      Social_Position_Method = social_position_method,
      stringsAsFactors = FALSE
    )

    return(list(
      results_overall = results_overall,
      results_replicates = results_replicates,
      resampling_object = boot_object,
      metadata = list(
        measure = measure,
        estimand = estimand_text,
        estimand_type = if (measure == "sii") "absolute_difference" else "relative_ratio",
        resampling_method = "ecological_unit_bootstrap",
        resampling_unit = "ecological_unit",
        model_refitting_rule = "The selected model class is refitted in every bootstrap sample; model selection is not repeated.",
        ridit_rule = "The relative social position is reconstructed in every bootstrap sample using resampled population/exposure weights and the original social ordering.",
        interval_rule = if (measure == "sii") {
          "Percentile interval from valid bootstrap replicate estimates on the SII scale."
        } else {
          "Percentile interval from valid bootstrap replicate estimates on the RII scale; log-scale diagnostics are also reported."
        },
        n_requested_replicates = n_boot,
        n_valid_replicates = n_valid,
        seed = seed
      ),
      diagnostics = list(
        warnings = warnings,
        n_failed_replicates = n_boot - n_valid,
        failure_rate = (n_boot - n_valid) / n_boot,
        selected_model = model_name,
        ratio_positivity_check = if (measure == "rii") {
          "All valid RII bootstrap replicates require strictly positive endpoint predictions."
        } else {
          NA_character_
        }
      )
    ))
  }

  # ---------------------------------------------------------------------------
  # 4. Jackknife ecological-unit delete-one resampling.
  # ---------------------------------------------------------------------------

  theta <- function(retained_indices) {
    retained_indices <- as.integer(retained_indices)
    d_j <- d0[retained_indices, , drop = FALSE]
    value <- tryCatch(
      estimand_from_prepared_resample(d_j),
      error = function(e) NA_real_,
      warning = function(w) NA_real_
    )
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) NA_real_ else value
  }

  jack_object <- tryCatch(
    bootstrap::jackknife(x = seq_len(nrow(d0)), theta = theta),
    error = function(e) {
      warnings <<- add_warning(
        warnings,
        paste0(
          "`bootstrap::jackknife()` failed: ", conditionMessage(e),
          " The delete-one values were computed explicitly."
        )
      )
      NULL
    }
  )

  jack_values <- rep(NA_real_, nrow(d0))
  for (i in seq_len(nrow(d0))) {
    keep <- setdiff(seq_len(nrow(d0)), i)
    jack_values[i] <- theta(keep)
  }

  valid <- is.finite(jack_values)
  if (measure == "rii") {
    valid <- valid & jack_values > 0
  }
  n_valid <- sum(valid)
  m <- length(jack_values)

  if (n_valid < max(3L, ceiling(0.75 * m))) {
    stop(
      "Too few valid jackknife delete-one estimates were obtained. ",
      "Valid estimates: ", n_valid, " of ", m, ". ",
      "Action: inspect influential units, model convergence, and endpoint positivity; consider `ci_method = 'wald'` or a simpler model.",
      call. = FALSE
    )
  }

  if (n_valid < m) {
    warnings <- add_warning(
      warnings,
      paste0(
        "Some jackknife delete-one estimates failed and were excluded from the SE calculation: ",
        m - n_valid, " of ", m, "."
      )
    )
  }

  valid_values <- jack_values[valid]

  if (measure == "sii") {
    jack_mean <- mean(valid_values)
    se <- sqrt((n_valid - 1) / n_valid * sum((valid_values - jack_mean)^2))
    ci_lower <- point_estimate - z_value * se
    ci_upper <- point_estimate + z_value * se

    log_estimate <- NA_real_
    se_log <- NA_real_
    log_ci_lower <- NA_real_
    log_ci_upper <- NA_real_
    standard_error_scale <- "difference_scale"
    interval_scale <- "difference_scale"
    interval_rule <- "Normal approximation using delete-one jackknife standard error on the SII difference scale."
  } else {
    log_values <- log(valid_values)
    jack_mean <- mean(log_values)
    se_log <- sqrt((n_valid - 1) / n_valid * sum((log_values - jack_mean)^2))

    log_estimate <- point_log_estimate
    log_ci_lower <- log_estimate - z_value * se_log
    log_ci_upper <- log_estimate + z_value * se_log
    ci_lower <- exp(log_ci_lower)
    ci_upper <- exp(log_ci_upper)

    se <- NA_real_
    standard_error_scale <- "log_ratio_scale"
    interval_scale <- "log_ratio_scale_exponentiated"
    interval_rule <- "Normal approximation using delete-one jackknife standard error on the log-RII scale, exponentiated to the RII scale."
  }

  results_replicates <- data.frame(
    Replicate = seq_along(jack_values),
    Replicate_Type = "jackknife_delete_one",
    Omitted_Unit = as.character(d0$unit_id),
    Estimate = jack_values,
    Log_Estimate = if (measure == "rii") ifelse(jack_values > 0, log(jack_values), NA_real_) else NA_real_,
    Status = ifelse(valid, "valid", "failed"),
    Message = ifelse(valid, NA_character_, "Delete-one sample failed to refit model or produce valid endpoint predictions."),
    stringsAsFactors = FALSE
  )

  results_overall <- data.frame(
    Metric = measure_label,
    Estimate = point_estimate,
    Log_Estimate = log_estimate,
    Standard_Error = se,
    Standard_Error_Scale = standard_error_scale,
    Standard_Error_Log = se_log,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    Log_CI_Lower = log_ci_lower,
    Log_CI_Upper = log_ci_upper,
    Confidence_Level = conf_level,
    CI_Method = "jackknife",
    Interval_Type = "wald_from_delete_one_ecological_jackknife_se",
    Interval_Scale = interval_scale,
    Selected_Model = model_name,
    Health_Indicator_Type = health_indicator_type,
    Input_Source = input_source,
    Social_Position_Method = social_position_method,
    stringsAsFactors = FALSE
  )

  list(
    results_overall = results_overall,
    results_replicates = results_replicates,
    resampling_object = jack_object,
    metadata = list(
      measure = measure,
      estimand = estimand_text,
      estimand_type = if (measure == "sii") "absolute_difference" else "relative_ratio",
      resampling_method = "ecological_unit_jackknife_delete_one",
      resampling_unit = "ecological_unit",
      model_refitting_rule = "The selected model class is refitted in every delete-one sample; model selection is not repeated.",
      ridit_rule = "The relative social position is reconstructed in every delete-one sample using remaining population/exposure weights and the original social ordering.",
      interval_rule = interval_rule,
      n_delete_one_estimates = m,
      n_valid_delete_one_estimates = n_valid,
      jackknife_mean = jack_mean
    ),
    diagnostics = list(
      warnings = warnings,
      n_failed_replicates = m - n_valid,
      failure_rate = (m - n_valid) / m,
      selected_model = model_name,
      ratio_positivity_check = if (measure == "rii") {
        "All valid RII jackknife estimates require strictly positive endpoint predictions."
      } else {
        NA_character_
      }
    )
  )
}
