#' Validation and Configuration Engine for Regression-Based Inequality Indices
#'
#' @description
#' `sii_rii_validate_inputs` is the strict upstream entry point for `sii_ineqeco()` and
#' `rii_ineqeco()`. It performs exhaustive structural and statistical typing validations
#' on ecological-unit vectors, mapping analytical pathways to model specifications based
#' on the user's measure estimand (`sii` or `rii`).
#'
#' This module enforces constraints on population distributions, verifies variance methodologies,
#' identifies outcome-bound violations, and locks the analytical configuration object passed
#' to downstream regression and prediction steps.
#'
#' @noRd
.sii_rii_validate_inputs <- function(data,
                                     measure = c("sii", "rii"),
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
                                 verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # 1. Small local validators used only inside this input-validation module.
  # ---------------------------------------------------------------------------

  is_scalar_character <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }

  is_integerish <- function(x, tol = sqrt(.Machine$double.eps)) {
    is.numeric(x) && all(is.finite(x)) && all(abs(x - round(x)) < tol)
  }

  normalize_model_names <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x <- gsub("\\s+", " ", x)

    aliases <- c(
      "linear" = "linear_identity",
      "identity" = "linear_identity",
      "linear identity" = "linear_identity",
      "linear_identity" = "linear_identity",
      "linear ln(y)" = "linear_log",
      "linear log(y)" = "linear_log",
      "linear log" = "linear_log",
      "linear_log" = "linear_log",
      "negative binomial" = "negative_binomial",
      "negative_binomial" = "negative_binomial",
      "negbin" = "negative_binomial",
      "binomial" = "binomial",
      "quasibinomial" = "quasibinomial",
      "poisson" = "poisson",
      "quasipoisson" = "quasipoisson"
    )

    out <- unname(aliases[x])
    out[is.na(out)] <- x[is.na(out)]
    out
  }

  allowed_models_for_case <- function(type, source) {
    if (source == "counts" && type %in% c("proportion", "percentage")) {
      return(c("binomial", "quasibinomial", "linear_identity"))
    }
    if (source == "counts" && type %in% c("rate", "ratio")) {
      return(c("poisson", "quasipoisson", "negative_binomial", "linear_log", "linear_identity"))
    }
    if (source == "indicator_population" && type %in% c("proportion", "percentage")) {
      return(c("quasibinomial", "linear_identity"))
    }
    if (source == "indicator_population" && type %in% c("rate", "ratio")) {
      return(c("linear_log", "linear_identity"))
    }
    stop("Internal error: unsupported indicator type/input source combination.", call. = FALSE)
  }

  validate_aic_bic_comparability <- function(metric, candidate_models, type, source) {
    if (!(metric %in% c("AIC", "BIC"))) {
      return(list(ok = TRUE, message = NA_character_))
    }

    if (source == "counts" && type %in% c("rate", "ratio") &&
        all(candidate_models %in% c("poisson", "negative_binomial"))) {
      return(list(ok = TRUE,
                  message = "AIC/BIC selection is restricted to likelihood-comparable count models."))
    }

    if (source == "counts" && type %in% c("proportion", "percentage") &&
        all(candidate_models %in% c("binomial"))) {
      return(list(ok = TRUE,
                  message = "AIC/BIC selection is available for the binomial model; no cross-family comparison is performed."))
    }

    list(
      ok = FALSE,
      message = paste0(
        "AIC/BIC cannot be used to compare the requested models because they are not likelihood-comparable. ",
        "Use MAE, MSE, RMSE, or MAPE, or provide a likelihood-comparable subset of models. ",
        "For rates/ratios with counts, comparable AIC/BIC candidates are: poisson, negative_binomial. ",
        "For proportions/percentages with counts, AIC/BIC is only valid for binomial-only selection."
      )
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Validate global scalar arguments.
  # ---------------------------------------------------------------------------

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or tibble.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`data` has zero rows.", call. = FALSE)
  }

  health_indicator_type <- match.arg(tolower(trimws(as.character(health_indicator_type))),
                                     choices = c("proportion", "percentage", "rate", "ratio"))
  social_position_method <- match.arg(tolower(trimws(as.character(social_position_method))),
                                      choices = c("classic", "bounded"))
  ci_method <- match.arg(tolower(trimws(as.character(ci_method))),
                         choices = c("wald", "bootstrap", "jackknife"))

  model_selection_metric <- toupper(trimws(as.character(model_selection_metric)[1L]))
  if (!(model_selection_metric %in% c("MSE", "MAE", "RMSE", "MAPE", "AIC", "BIC"))) {
    stop(
      "`model_selection_metric` must be one of: MSE, MAE, RMSE, MAPE, AIC, BIC. ",
      "Action: use MSE as the default robust choice for comparing different model families on the final outcome scale, ",
      "or use AIC/BIC only with likelihood-comparable model subsets.",
      call. = FALSE
    )
  }

  if (!is.logical(higher_ineq_is_favorable) || length(higher_ineq_is_favorable) != 1L ||
      is.na(higher_ineq_is_favorable)) {
    stop("`higher_ineq_is_favorable` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single numeric value between 0 and 1.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0 ||
        abs(seed - round(seed)) > sqrt(.Machine$double.eps)) {
      stop("`seed` must be NULL or a single non-negative integer-like numeric value.", call. = FALSE)
    }
    seed <- as.integer(round(seed))
  }

  if (ci_method == "bootstrap") {
    if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) ||
        n_boot < 100 || abs(n_boot - round(n_boot)) > sqrt(.Machine$double.eps)) {
      stop("`n_boot` must be a single integer-like numeric value >= 100 when `ci_method = 'bootstrap'`.",
           call. = FALSE)
    }
    n_boot <- as.integer(round(n_boot))
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop("Package `boot` is required when `ci_method = 'bootstrap'`. Please install it.",
           call. = FALSE)
    }
  } else {
    n_boot <- as.integer(round(n_boot))
  }

  if (ci_method == "wald") {
    if (!requireNamespace("emmeans", quietly = TRUE)) {
      stop(
        "Package `emmeans` is required when `ci_method = 'wald'`. ",
        "Action: install `emmeans`, or use `ci_method = 'bootstrap'` or `ci_method = 'jackknife'`.",
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Validate variable-name arguments and decide the input source.
  # ---------------------------------------------------------------------------

  if (!is_scalar_character(equity_stratifier_var)) {
    stop("`equity_stratifier_var` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is_scalar_character(analysis_unit_var)) {
    stop("`analysis_unit_var` must be a single non-empty character string.", call. = FALSE)
  }

  validation_warnings <- character(0)
  ignored_arguments <- character(0)

  has_complete_counts <- !is.null(health_numerator_var) && !is.null(health_denominator_var)
  has_partial_counts <- xor(!is.null(health_numerator_var), !is.null(health_denominator_var))
  has_complete_indicator_population <- !is.null(health_indicator_var) && !is.null(population_weights_var)
  has_partial_indicator_population <- xor(!is.null(health_indicator_var), !is.null(population_weights_var))

  if (has_partial_counts) {
    stop(
      "Incomplete count-based input: `health_numerator_var` and `health_denominator_var` must be supplied together. ",
      "Action: either provide both variables, or omit both and use `health_indicator_var` + `population_weights_var`.",
      call. = FALSE
    )
  }

  if (has_partial_indicator_population) {
    stop(
      "Incomplete indicator-population input: `health_indicator_var` and `population_weights_var` must be supplied together. ",
      "Action: either provide both variables, or omit both and use `health_numerator_var` + `health_denominator_var`.",
      call. = FALSE
    )
  }

  if (has_complete_counts && (!is_scalar_character(health_numerator_var) ||
                              !is_scalar_character(health_denominator_var))) {
    stop(
      "For count-based SII, both `health_numerator_var` and `health_denominator_var` must be single non-empty character strings. ",
      "Action: pass column names such as `health_numerator_var = \"deaths\"` and `health_denominator_var = \"population\"`.",
      call. = FALSE
    )
  }

  if (has_complete_indicator_population && (!is_scalar_character(health_indicator_var) ||
                                            !is_scalar_character(population_weights_var))) {
    stop(
      "For indicator-population SII, both `health_indicator_var` and `population_weights_var` must be single non-empty character strings. ",
      "Action: pass column names such as `health_indicator_var = \"rate\"` and `population_weights_var = \"population\"`.",
      call. = FALSE
    )
  }

  if (!has_complete_counts && !has_complete_indicator_population) {
    stop(
      "No complete input source was supplied. ",
      "Action: provide either `health_numerator_var` + `health_denominator_var`, or `health_indicator_var` + `population_weights_var`.",
      call. = FALSE
    )
  }

  # Priority rule: if numerator and denominator are supplied, counts are the
  # source of the point estimate and model family. Indicator/population columns,
  # if also supplied, are not used by this version.
  input_source <- if (has_complete_counts) "counts" else "indicator_population"

  if (has_complete_counts && has_complete_indicator_population) {
    validation_warnings <- c(
      validation_warnings,
      paste0(
        "Both count variables and indicator-population variables were supplied. ",
        "This version uses `health_numerator_var` + `health_denominator_var` as the primary source and ignores ",
        "`health_indicator_var` + `population_weights_var` for point estimation, model fitting, weighting, and uncertainty."
      )
    )
    ignored_arguments <- c(ignored_arguments, "health_indicator_var", "population_weights_var")
  }

  required_vars <- c(equity_stratifier_var, analysis_unit_var)
  if (input_source == "counts") {
    required_vars <- c(required_vars, health_numerator_var, health_denominator_var)
  } else {
    required_vars <- c(required_vars, health_indicator_var, population_weights_var)
  }

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop("The following variables were not found in `data`: ",
         paste(missing_vars, collapse = ", "), ".", call. = FALSE)
  }

  duplicated_roles <- names(which(table(required_vars) > 1L))
  if (length(duplicated_roles) > 0L) {
    stop("Each variable role must refer to a distinct column. Duplicated role columns: ",
         paste(duplicated_roles, collapse = ", "), ".", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 4. Extract and validate required columns without modifying the data.
  # ---------------------------------------------------------------------------

  unit_id <- data[[analysis_unit_var]]
  equity <- data[[equity_stratifier_var]]

  if (any(is.na(unit_id))) {
    stop("Missing values detected in `analysis_unit_var`.", call. = FALSE)
  }
  if (anyDuplicated(as.character(unit_id)) > 0L) {
    stop("`analysis_unit_var` must uniquely identify ecological units; duplicated unit IDs were found.",
         call. = FALSE)
  }
  if (!is.numeric(equity) || any(!is.finite(equity))) {
    stop("`equity_stratifier_var` must be numeric, finite, and non-missing.", call. = FALSE)
  }
  if (length(unique(equity)) < 2L) {
    stop("`equity_stratifier_var` must contain at least two distinct values to estimate a social gradient.",
         call. = FALSE)
  }

  if (input_source == "counts") {
    numerator <- data[[health_numerator_var]]
    denominator <- data[[health_denominator_var]]

    if (!is.numeric(numerator) || any(!is.finite(numerator))) {
      stop("`health_numerator_var` must be numeric, finite, and non-missing.", call. = FALSE)
    }
    if (!is.numeric(denominator) || any(!is.finite(denominator))) {
      stop("`health_denominator_var` must be numeric, finite, and non-missing.", call. = FALSE)
    }
    if (any(numerator < 0)) {
      stop("`health_numerator_var` must contain values >= 0.", call. = FALSE)
    }
    if (any(denominator <= 0)) {
      stop("`health_denominator_var` must contain values > 0.", call. = FALSE)
    }
    if (!is_integerish(numerator)) {
      stop("`health_numerator_var` must contain integer counts. This function does not silently round counts.",
           call. = FALSE)
    }
    if (health_indicator_type %in% c("proportion", "percentage") && !is_integerish(denominator)) {
      stop("`health_denominator_var` must contain integer denominators for binomial proportion/percentage models. This function does not silently round denominators.",
           call. = FALSE)
    }

    if (health_indicator_type %in% c("proportion", "percentage") && any(numerator > denominator)) {
      stop("For `proportion` or `percentage` with counts, numerator values must be <= denominator values.",
           call. = FALSE)
    }

    if (health_indicator_type %in% c("rate", "ratio")) {
      if (!is.numeric(rate_scaling_factor) || length(rate_scaling_factor) != 1L ||
          is.na(rate_scaling_factor) || !is.finite(rate_scaling_factor) || rate_scaling_factor <= 0) {
        stop("`rate_scaling_factor` must be a single numeric value > 0 for `rate` or `ratio` with counts.",
             call. = FALSE)
      }
    } else {
      rate_scaling_factor <- 1
    }

    observed_outcome <- if (health_indicator_type %in% c("proportion", "percentage")) {
      numerator / denominator * ifelse(health_indicator_type == "percentage", 100, 1)
    } else {
      numerator / denominator * rate_scaling_factor
    }
    model_weight <- denominator
  }

  if (input_source == "indicator_population") {
    indicator <- data[[health_indicator_var]]
    population_weight <- data[[population_weights_var]]

    if (!is.numeric(indicator) || any(!is.finite(indicator))) {
      stop("`health_indicator_var` must be numeric, finite, and non-missing.", call. = FALSE)
    }
    if (!is.numeric(population_weight) || any(!is.finite(population_weight))) {
      stop("`population_weights_var` must be numeric, finite, and non-missing.", call. = FALSE)
    }
    if (any(population_weight <= 0)) {
      stop("`population_weights_var` must contain values > 0.", call. = FALSE)
    }

    if (health_indicator_type == "proportion" && any(indicator < 0 | indicator > 1)) {
      stop("For `proportion`, `health_indicator_var` must be on the 0-1 scale.", call. = FALSE)
    }
    if (health_indicator_type == "percentage" && any(indicator < 0 | indicator > 100)) {
      stop("For `percentage`, `health_indicator_var` must be on the 0-100 scale.", call. = FALSE)
    }
    if (health_indicator_type %in% c("rate", "ratio") && any(indicator < 0)) {
      stop("For `rate` or `ratio`, `health_indicator_var` must contain values >= 0.", call. = FALSE)
    }

    if (health_indicator_type %in% c("rate", "ratio")) {
      if (is.null(rate_scaling_factor)) {
        rate_scaling_factor <- 1
      }
      if (!is.numeric(rate_scaling_factor) || length(rate_scaling_factor) != 1L ||
          is.na(rate_scaling_factor) || !is.finite(rate_scaling_factor) || rate_scaling_factor <= 0) {
        stop("`rate_scaling_factor`, when supplied for `rate` or `ratio`, must be a single numeric value > 0.",
             call. = FALSE)
      }
    } else {
      rate_scaling_factor <- 1
    }

    observed_outcome <- indicator
    model_weight <- population_weight
  }

  # ---------------------------------------------------------------------------
  # 5. Validate model-selection metric against the observed outcome.
  # ---------------------------------------------------------------------------

  if (model_selection_metric == "MAPE" && any(observed_outcome == 0)) {
    stop(
      "`model_selection_metric = 'MAPE'` is not available because the observed outcome contains zero values. ",
      "Action: use MSE, RMSE, or MAE instead. MSE is the default and is comparable across candidate models when evaluated on the final outcome scale.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Resolve allowed/candidate models and validate model comparability.
  # ---------------------------------------------------------------------------

  allowed_models <- allowed_models_for_case(health_indicator_type, input_source)

  if (is.null(models)) {
    candidate_models <- allowed_models
    model_set_source <- "all_allowed_models"
  } else {
    if (!is.character(models) || length(models) < 1L || any(is.na(models)) ||
        any(!nzchar(trimws(models)))) {
      stop("`models` must be NULL or a non-empty character vector of model names.", call. = FALSE)
    }
    candidate_models <- unique(normalize_model_names(models))
    model_set_source <- "user_requested_models"
  }

  invalid_models <- setdiff(candidate_models, allowed_models)
  if (length(invalid_models) > 0L) {
    stop(
      "The following models are not allowed for `health_indicator_type = '", health_indicator_type,
      "'` with input source `", input_source, "`: ", paste(invalid_models, collapse = ", "),
      ". Allowed models are: ", paste(allowed_models, collapse = ", "), ". ",
      "Action: use `models = NULL` to evaluate all allowed models automatically, or choose only from the allowed set above.",
      call. = FALSE
    )
  }

  if ("negative_binomial" %in% candidate_models && !requireNamespace("MASS", quietly = TRUE)) {
    stop("Package `MASS` is required for `negative_binomial` models. Please install it.",
         call. = FALSE)
  }

  if (any(candidate_models %in% c("quasibinomial", "quasipoisson")) &&
      model_selection_metric %in% c("AIC", "BIC")) {
    stop(
      "AIC/BIC cannot be used with quasi-likelihood models (`quasibinomial` or `quasipoisson`) because they do not provide a comparable full likelihood. ",
      "Action: use MSE, RMSE, MAE, or MAPE; or remove quasi models and provide a likelihood-comparable subset with `models = ...`.",
      call. = FALSE
    )
  }

  comparability <- validate_aic_bic_comparability(
    metric = model_selection_metric,
    candidate_models = candidate_models,
    type = health_indicator_type,
    source = input_source
  )
  if (!isTRUE(comparability$ok)) {
    stop(comparability$message, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 7. Validate model-specific data requirements.
  # ---------------------------------------------------------------------------

  linear_log_data_available <- !("linear_log" %in% candidate_models && any(observed_outcome <= 0))


  model_data_flags <- list(
    linear_log_data_available = linear_log_data_available,
    linear_log_message = if (linear_log_data_available) NA_character_ else
      "`linear_log` requires strictly positive observed outcome values and should be excluded during model fitting."
  )

  if (!is.null(models)) {
    requested_unavailable <- character(0)
    if ("linear_log" %in% candidate_models && !linear_log_data_available) {
      requested_unavailable <- c(requested_unavailable, "linear_log")
    }
    if (length(requested_unavailable) == length(candidate_models)) {
      stop(
        "All requested models are incompatible with the observed outcome values: ",
        paste(requested_unavailable, collapse = ", "), ". ",
        "Action: use `models = NULL` to let the function evaluate all allowed models and exclude incompatible ones, ",
        "or include `linear_identity` when appropriate.",
        call. = FALSE
      )
    }
  }

  if (length(unique(observed_outcome)) < 2L) {
    stop("The observed health outcome must contain at least two distinct values to estimate a gradient.",
         call. = FALSE)
  }

  if (length(unique(model_weight)) < 1L || sum(model_weight) <= 0) {
    stop("The population/exposure weights must sum to a positive value.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 8. Return a compact, auditable validation object for downstream modules.
  # ---------------------------------------------------------------------------

  estimand_label <- if (measure == "sii") {
    "SII = Y_advantaged - Y_disadvantaged"
  } else {
    "RII = Y_advantaged / Y_disadvantaged"
  }

  estimand_type <- if (measure == "sii") "absolute_difference" else "relative_ratio"
  extreme_contrast_rule <- if (measure == "sii") {
    "f(advantaged_position) - f(disadvantaged_position)"
  } else {
    "f(advantaged_position) / f(disadvantaged_position)"
  }
  wald_interval_scale <- if (measure == "sii") "outcome_difference_scale" else "log_ratio_scale"

  if (measure == "rii") {
    validation_warnings <- c(
      validation_warnings,
      paste0(
        "RII requires strictly positive predicted values at the disadvantaged and advantaged social-position extremes. ",
        "This cannot be guaranteed during input validation and must be checked after model fitting/estimation."
      )
    )
  }

  list(
    measure = measure,
    health_indicator_type = health_indicator_type,
    input_source = input_source,
    variable_names = list(
      analysis_unit_var = analysis_unit_var,
      equity_stratifier_var = equity_stratifier_var,
      health_indicator_var = health_indicator_var,
      population_weights_var = population_weights_var,
      health_numerator_var = health_numerator_var,
      health_denominator_var = health_denominator_var
    ),
    higher_ineq_is_favorable = higher_ineq_is_favorable,
    rate_scaling_factor = rate_scaling_factor,
    allowed_models = allowed_models,
    candidate_models = candidate_models,
    model_set_source = model_set_source,
    model_selection_metric = model_selection_metric,
    model_selection_metric_definition = paste0(
      "Predictive metrics MAE/MSE/RMSE/MAPE are computed on the final outcome scale and weighted by the ",
      "population/exposure weights. AIC/BIC are only allowed for likelihood-comparable model sets."
    ),
    aic_bic_comparability_note = comparability$message,
    social_position_method = social_position_method,
    ci_method = ci_method,
    n_boot = n_boot,
    conf_level = conf_level,
    seed = seed,
    verbose = verbose,
    estimand = estimand_label,
    estimand_type = estimand_type,
    extreme_contrast_rule = extreme_contrast_rule,
    wald_interval_scale = wald_interval_scale,
    outcome_scale = if (health_indicator_type == "percentage") "0_to_100" else if (health_indicator_type == "proportion") "0_to_1" else "original_scaled_indicator",
    diagnostics = list(
      n_units = nrow(data),
      n_candidate_models = length(candidate_models),
      metric_is_weighted = model_selection_metric %in% c("MAE", "MSE", "RMSE", "MAPE"),
      metric_weighting_rule = "Predictive model-selection metrics are population/exposure-weighted.",
      rii_positivity_rule = if (measure == "rii") "RII requires predicted disadvantaged and advantaged values to be strictly positive; checked downstream." else NA_character_,
      external_standard_errors = "Not part of this version: `health_se_var` is intentionally absent from the interface.",
      warnings = validation_warnings,
      ignored_arguments = unique(ignored_arguments),
      model_data_flags = model_data_flags
    )
  )
}
