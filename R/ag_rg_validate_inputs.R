#' Input Validation and Standardization Engine for AG/RG Functions
#'
#' @description
#' `ag_rg_validate_inputs` enforces strict epidemiological and statistical typing constraints
#' on all user-supplied data prior to inequality estimation.
#'
#' It performs systematic checks on input dimensions, numeric bounds, logical flags, and missing data.
#' It resolves complex estimation pathways ("scenarios") based strictly on the specified arguments,
#' differentiating between parametric delta-method scenarios (raw counts or provided standard errors)
#' and deterministic scenarios.
#'
#' If numerator and denominator are supplied, count-based estimation is
#' prioritized. Health indicator and population-weight arguments, if also
#' supplied, are validated as column names but ignored for point estimation.
#'
#' Grouping validation is conditional:
#'
#' * `territorial_method` is required and validated only when
#'   `grouping_approach = "territorial"`.
#' * For `weighted_cut`, `weighted_midpoint`, and `fractional`,
#'   `territorial_method` is ignored and standardized to `NA_character_`.
#'
#' Resampling-package checks are conditional:
#'
#' * `boot` is required only for `ci_method = "bootstrap"`.
#' * `bootstrap` is required only for `ci_method = "jackknife"`.
#'
#' @importFrom stats complete.cases
#' @keywords internal
#' @noRd
.ag_rg_validate_inputs <- function(data,
                                     health_indicator_type,
                                     health_indicator_var,
                                     population_weights_var,
                                     health_numerator_var,
                                     health_denominator_var,
                                     health_se_var,
                                     analysis_unit_var,
                                     equity_stratifier_var,
                                     higher_ineq_is_favorable,
                                     rate_scaling_factor,
                                     grouping_approach,
                                     territorial_method,
                                     n_groups,
                                     ci_method,
                                     n_boot,
                                     conf_level,
                                     verbose) {

  stop_if <- function(condition, message) {
    if (isTRUE(condition)) stop(message, call. = FALSE)
  }

  log_warning <- character(0)
  log_message <- character(0)

  add_warning <- function(message) {
    log_warning <<- c(log_warning, message)
    if (isTRUE(verbose)) warning(message, call. = FALSE, immediate. = TRUE)
    invisible(NULL)
  }

  add_message <- function(message) {
    log_message <<- c(log_message, message)
    if (isTRUE(verbose)) message(message)
    invisible(NULL)
  }

  # ---------------------------------------------------------------------------
  # 1. Basic data checks
  # ---------------------------------------------------------------------------
  stop_if(missing(data) || is.null(data), "`data` must be supplied as a data frame.")
  stop_if(!is.data.frame(data), "`data` must be a data frame with one row per ecological unit.")
  stop_if(nrow(data) < 2, "`data` must contain at least two ecological units.")

  stop_if(is.null(equity_stratifier_var), "`equity_stratifier_var` must be supplied.")

  get_required_column <- function(var_name, arg_name) {
    stop_if(is.null(var_name) || !is.character(var_name) || length(var_name) != 1,
            sprintf("`%s` must be a single character string naming a column in `data`.", arg_name))
    stop_if(!(var_name %in% names(data)),
            sprintf("Column `%s` supplied through `%s` was not found in `data`.", var_name, arg_name))
    data[[var_name]]
  }

  get_optional_column <- function(var_name, arg_name) {
    if (is.null(var_name)) return(NULL)
    stop_if(!is.character(var_name) || length(var_name) != 1,
            sprintf("`%s` must be NULL or a single character string naming a column in `data`.", arg_name))
    stop_if(!(var_name %in% names(data)),
            sprintf("Column `%s` supplied through `%s` was not found in `data`.", var_name, arg_name))
    data[[var_name]]
  }

  # Validate all non-NULL user-supplied column names. Some columns may later be
  # ignored when counts are prioritized, but invalid supplied names are still
  # treated as user errors.
  equity_stratifier <- get_required_column(equity_stratifier_var, "equity_stratifier_var")
  analysis_unit <- get_optional_column(analysis_unit_var, "analysis_unit_var")
  health_indicator <- get_optional_column(health_indicator_var, "health_indicator_var")
  population_weights <- get_optional_column(population_weights_var, "population_weights_var")
  health_numerator <- get_optional_column(health_numerator_var, "health_numerator_var")
  health_denominator <- get_optional_column(health_denominator_var, "health_denominator_var")
  health_se <- get_optional_column(health_se_var, "health_se_var")

  # ---------------------------------------------------------------------------
  # 2. Argument choices
  # ---------------------------------------------------------------------------
  health_indicator_type <- match.arg(
    health_indicator_type,
    choices = c("rate", "ratio", "proportion", "percentage")
  )

  grouping_approach <- match.arg(
    grouping_approach,
    choices = c("territorial", "weighted_cut", "weighted_midpoint", "fractional")
  )

  if (identical(grouping_approach, "territorial")) {
    if (is.null(territorial_method) || (length(territorial_method) == 1 && is.na(territorial_method))) {
      stop(
        "`territorial_method` must be supplied when `grouping_approach = 'territorial'`. ",
        "Use one of: 'quantile', 'kmeans', 'fisher', 'manual', or 'equal'.",
        call. = FALSE
      )
    }
    territorial_method <- match.arg(
      territorial_method,
      choices = c("quantile", "kmeans", "fisher", "manual", "equal")
    )
  } else {
    # Important: territorial_method does not apply to non-territorial grouping.
    # It must not be validated or required here.
    territorial_method <- NA_character_
  }

  ci_method_requested <- match.arg(
    ci_method,
    choices = c("delta", "bootstrap", "jackknife")
  )

  stop_if(is.null(higher_ineq_is_favorable) || !is.logical(higher_ineq_is_favorable) ||
            length(higher_ineq_is_favorable) != 1 || is.na(higher_ineq_is_favorable),
          "`higher_ineq_is_favorable` must be a single TRUE/FALSE value.")

  stop_if(!is.logical(verbose) || length(verbose) != 1 || is.na(verbose),
          "`verbose` must be a single TRUE/FALSE value.")

  n_groups <- suppressWarnings(as.integer(n_groups))
  stop_if(is.na(n_groups) || n_groups < 2,
          "`n_groups` must be an integer >= 2.")

  stop_if(!is.numeric(conf_level) || length(conf_level) != 1 || !is.finite(conf_level) ||
            conf_level <= 0 || conf_level >= 1,
          "`conf_level` must be a numeric value between 0 and 1, e.g., 0.95.")

  if (identical(ci_method_requested, "bootstrap")) {
    n_boot <- suppressWarnings(as.integer(n_boot))
    stop_if(is.na(n_boot) || n_boot < 100,
            "`n_boot` must be an integer >= 100 when `ci_method = 'bootstrap'`.")
    if (!requireNamespace("boot", quietly = TRUE)) {
      stop(
        "Package 'boot' is required for `ci_method = 'bootstrap'`. ",
        "Please install it with install.packages('boot').",
        call. = FALSE
      )
    }
  } else {
    # Keep a numeric value for downstream outputs, but do not require a minimum
    # when bootstrap is not requested.
    n_boot <- suppressWarnings(as.integer(n_boot))
    if (is.na(n_boot)) n_boot <- NA_integer_
  }

  if (identical(ci_method_requested, "jackknife") &&
      !requireNamespace("bootstrap", quietly = TRUE)) {
    stop(
      "Package 'bootstrap' is required for `ci_method = 'jackknife'`. ",
      "Please install it with install.packages('bootstrap').",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 3. Scale and variance family
  # ---------------------------------------------------------------------------
  if (health_indicator_type %in% c("rate", "ratio")) {
    stop_if(is.null(rate_scaling_factor) || !is.numeric(rate_scaling_factor) ||
              length(rate_scaling_factor) != 1 || !is.finite(rate_scaling_factor) ||
              rate_scaling_factor <= 0,
            "For rates/ratios, `rate_scaling_factor` must be a numeric scalar > 0, e.g., 1000 or 100000.")
    final_multiplier <- rate_scaling_factor
    variance_family <- "poisson"
  } else if (identical(health_indicator_type, "percentage")) {
    final_multiplier <- 100
    variance_family <- "binomial"
  } else {
    final_multiplier <- 1
    variance_family <- "binomial"
  }

  # ---------------------------------------------------------------------------
  # 4. Scenario detection from supplied arguments
  # ---------------------------------------------------------------------------
  has_numerator_arg <- !is.null(health_numerator_var)
  has_denominator_arg <- !is.null(health_denominator_var)
  has_indicator_arg <- !is.null(health_indicator_var)
  has_population_arg <- !is.null(population_weights_var)
  has_se_arg <- !is.null(health_se_var)

  stop_if(xor(has_numerator_arg, has_denominator_arg),
          "`health_numerator_var` and `health_denominator_var` must be supplied together or both be NULL.")

  has_counts_args <- has_numerator_arg && has_denominator_arg

  if (!has_counts_args) {
    stop_if(xor(has_indicator_arg, has_population_arg),
            "When count variables are not supplied, `health_indicator_var` and `population_weights_var` must be supplied together or both be NULL.")
    stop_if(!(has_indicator_arg && has_population_arg),
            "Provide either (`health_numerator_var` and `health_denominator_var`) or (`health_indicator_var` and `population_weights_var`).")
  }

  # Count-based estimation takes priority whenever numerator and denominator are
  # supplied. Additional indicator/population arguments are allowed but ignored
  # for point estimation, after their column names have been validated above.
  use_counts <- has_counts_args
  use_indicator_population <- !use_counts && has_indicator_arg && has_population_arg
  use_provided_se <- has_se_arg

  if (use_counts && (has_indicator_arg || has_population_arg)) {
    add_warning(paste0(
      "Count variables were supplied, so count-based point estimation is used. ",
      "Any supplied health indicator or population-weight variables are ignored for point estimation."
    ))
  }

  scenario <- if (use_counts && !use_provided_se) {
    "counts_no_se"
  } else if (use_counts && use_provided_se) {
    "counts_with_se"
  } else if (use_indicator_population && !use_provided_se) {
    "indicator_population_no_se"
  } else {
    "indicator_population_with_se"
  }

  # ---------------------------------------------------------------------------
  # 5. Length and numeric checks
  # ---------------------------------------------------------------------------
  stop_if(!is.numeric(equity_stratifier),
          "The column supplied through `equity_stratifier_var` must be numeric.")

  n <- length(equity_stratifier)
  if (!is.null(analysis_unit)) {
    stop_if(length(analysis_unit) != n,
            "The column supplied through `analysis_unit_var` must have the same length as `equity_stratifier_var`.")
  }

  if (use_counts) {
    stop_if(length(health_numerator) != n || length(health_denominator) != n,
            "The columns supplied through `health_numerator_var`, `health_denominator_var`, and `equity_stratifier_var` must have the same length.")
    stop_if(!is.numeric(health_numerator) || !is.numeric(health_denominator),
            "The columns supplied through `health_numerator_var` and `health_denominator_var` must be numeric.")
    stop_if(any(health_numerator < 0, na.rm = TRUE),
            "The column supplied through `health_numerator_var` cannot contain negative values.")
    stop_if(any(health_denominator <= 0, na.rm = TRUE),
            "The column supplied through `health_denominator_var` must contain values > 0.")
  } else {
    stop_if(length(health_indicator) != n || length(population_weights) != n,
            "The columns supplied through `health_indicator_var`, `population_weights_var`, and `equity_stratifier_var` must have the same length.")
    stop_if(!is.numeric(health_indicator) || !is.numeric(population_weights),
            "The columns supplied through `health_indicator_var` and `population_weights_var` must be numeric.")
    stop_if(any(population_weights <= 0, na.rm = TRUE),
            "The column supplied through `population_weights_var` must contain values > 0.")
  }

  if (use_provided_se) {
    stop_if(length(health_se) != n,
            "The column supplied through `health_se_var` must have the same length as `equity_stratifier_var`.")
    stop_if(!is.numeric(health_se),
            "The column supplied through `health_se_var` must be numeric.")
    stop_if(any(health_se < 0, na.rm = TRUE),
            "The column supplied through `health_se_var` cannot contain negative values.")
    stop_if(all(is.na(health_se)),
            "`health_se_var` was supplied, but all SE values are missing. Use `health_se_var = NULL` if no standard errors are available.")
  }

  # ---------------------------------------------------------------------------
  # 6. Standardized internal data frame
  # ---------------------------------------------------------------------------
  if (use_counts) {
    df <- data.frame(
      unit_id = if (is.null(analysis_unit)) as.character(seq_len(n)) else as.character(analysis_unit),
      numerator = as.numeric(health_numerator),
      denominator = as.numeric(health_denominator),
      population_weight = as.numeric(health_denominator),
      equity_stratifier = as.numeric(equity_stratifier),
      stringsAsFactors = FALSE
    )
    df$indicator_internal <- df$numerator / df$denominator
    df$health_indicator <- df$indicator_internal * final_multiplier
  } else {
    df <- data.frame(
      unit_id = if (is.null(analysis_unit)) as.character(seq_len(n)) else as.character(analysis_unit),
      numerator = NA_real_,
      denominator = NA_real_,
      population_weight = as.numeric(population_weights),
      equity_stratifier = as.numeric(equity_stratifier),
      stringsAsFactors = FALSE
    )
    df$health_indicator <- as.numeric(health_indicator)
    df$indicator_internal <- df$health_indicator / final_multiplier
  }

  df$health_se <- if (use_provided_se) as.numeric(health_se) else NA_real_
  df$se_internal <- if (use_provided_se) df$health_se / final_multiplier else NA_real_

  # Standardized ordering variable.
  #
  # Contract used by .ag_rg_group_data():
  #   smaller strat_order values correspond to greater disadvantage;
  #   larger strat_order values correspond to greater advantage.
  #
  # Therefore, after sorting ascending by strat_order:
  #   group 1 = most disadvantaged / least advantaged;
  #   group K = most advantaged / least disadvantaged.
  #
  # If higher values of the equity stratifier are favorable (e.g., wealth),
  # the original stratifier already runs from disadvantage to advantage.
  # If higher values are unfavorable (e.g., unmet basic needs), the sign is
  # reversed so that the most disadvantaged units come first.
  df$strat_order <- if (isTRUE(higher_ineq_is_favorable)) {
    df$equity_stratifier
  } else {
    -df$equity_stratifier
  }

  # ---------------------------------------------------------------------------
  # 7. Missingness handling
  # ---------------------------------------------------------------------------
  required_vars <- c("indicator_internal", "population_weight", "equity_stratifier", "strat_order")
  if (use_counts) required_vars <- c(required_vars, "numerator", "denominator")

  n_before <- nrow(df)
  df <- df[stats::complete.cases(df[, required_vars, drop = FALSE]), , drop = FALSE]
  n_dropped <- n_before - nrow(df)

  if (n_dropped > 0) {
    add_warning(sprintf(
      "%s ecological unit(s) were removed because of missing values required for point estimation.",
      n_dropped
    ))
  }

  if (use_provided_se && any(is.na(df$se_internal))) {
    stop(
      "`health_se_var` contains missing values among included units. Missing SE values cannot be treated as zero.",
      call. = FALSE
    )
  }

  stop_if(nrow(df) < 2, "At least two complete ecological units are required.")

  # ---------------------------------------------------------------------------
  # 8. Return standardized object
  # ---------------------------------------------------------------------------
  list(
    df = df,
    scenario = scenario,
    final_multiplier = final_multiplier,
    variance_family = variance_family,
    use_counts = use_counts,
    grouping_approach = grouping_approach,
    territorial_method = territorial_method,
    ci_method_requested = ci_method_requested,
    n_boot = n_boot,
    n_original = n,
    log_warning = log_warning,
    log_message = log_message
  )
}
