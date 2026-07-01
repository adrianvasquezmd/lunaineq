#' Relative Ecological Inequality Gap for Aggregated Data
#'
#' @description
#' `rg_ineqeco()` computes the Relative Gap (RG) in a health indicator between the most
#' disadvantaged and most advantaged ecological groups. This metric quantifies the
#' relative scale (ratio) of inequality across population strata.
#'
#' The estimand is formally defined as:
#'
#' \deqn{RG = Y_D / Y_A}
#'
#' where \eqn{Y_D} denotes the estimate of the health indicator in the most
#' socioeconomically disadvantaged group, and \eqn{Y_A} denotes the estimate in
#' the most advantaged group.
#'
#' @param data A data frame containing the grouped ecological data.
#' @param health_indicator_type Character. The scale of the health indicator: `"rate"`, `"ratio"`, `"proportion"`, or `"percentage"`.
#' @param health_indicator_var Character. The name of the column containing the point estimates of the health indicator.
#' @param population_weights_var Character. The name of the column containing the population sizes or analytic weights for each ecological unit.
#' @param health_numerator_var Character. The name of the column containing the event counts (numerator). Required if raw counts are used.
#' @param health_denominator_var Character. The name of the column containing the population sizes (denominator). Required if raw counts are used.
#' @param health_se_var Character. The name of the column containing standard errors for the health indicator. Used for delta-method analytical variance.
#' @param analysis_unit_var Character. The name of the column identifying each ecological unit (e.g., municipality, province).
#' @param equity_stratifier_var Character. The name of the column containing the socioeconomic stratifier (e.g., income, education). Must be numeric or ordered factor.
#' @param higher_ineq_is_favorable Logical. If `TRUE`, implies that a higher indicator value is favorable. Defaults to `FALSE` (e.g., mortality rates).
#' @param rate_scaling_factor Numeric. A scaling factor to adjust rates (e.g., `1000` for per 1,000 population). If `NULL`, it defaults based on `health_indicator_type`.
#' @param grouping_approach Character. Method to group the ecological units: `"territorial"`, `"weighted_cut"`, `"weighted_midpoint"`, or `"fractional"`.
#' @param territorial_method Character. If `grouping_approach` is `"territorial"`, the algorithm to build quantiles: `"quantile"`, `"kmeans"`, `"fisher"`, `"manual"`, or `"equal"`.
#' @param manual_breaks Numeric vector. If `territorial_method` is `"manual"`, the precise cut-points for grouping.
#' @param n_groups Integer. The number of quantiles to form (e.g., 5 for quintiles). Defaults to `5`.
#' @param ci_method Character. The method for confidence interval estimation: `"delta"`, `"bootstrap"`, or `"jackknife"`.
#' @param n_boot Integer. The number of bootstrap replicates if `ci_method = "bootstrap"`. Defaults to `2000`.
#' @param conf_level Numeric. The confidence level for the intervals. Defaults to `0.95`.
#' @param seed Integer. An optional seed for reproducible resampling.
#' @param verbose Logical. Whether to print methodological notes and warnings to the console. Defaults to `TRUE`.
#'
#' @return An S3 object of class `c("rg_ineqeco", "ag_rg_ineqeco", "list")` containing:
#' \item{point_estimate}{The point estimate of the relative gap.}
#' \item{se}{The standard error of the estimate (often on the log-ratio scale).}
#' \item{ci_lower}{The lower bound of the confidence interval.}
#' \item{ci_upper}{The upper bound of the confidence interval.}
#' \item{methodological_notes}{Detailed explanations of the analytical approach and variance estimator used.}
#' \item{results_groups}{A data frame summarizing the grouped data, including group-specific estimates and variance.}
#'
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' \dontrun{
#' # Use the built-in paho_data
#' data(paho_data)
#' df_2020 <- subset(paho_data, year == 2020)
#' 
#' # Calculate relative gap using jackknife resampling
#' rg <- rg_ineqeco(
#'   data = df_2020,
#'   health_indicator_type = "rate",
#'   health_numerator_var = "maternal_death",
#'   health_denominator_var = "live_births",
#'   rate_scaling_factor = 100000,
#'   analysis_unit_var = "state",
#'   equity_stratifier_var = "ubn",
#'   grouping_approach = "fractional",
#'   n_groups = 5,
#'   ci_method = "jackknife"
#' )
#' print(rg)
#' }
rg_ineqeco <- function(data,
                       health_indicator_type = c("rate", "ratio", "proportion", "percentage"),
                       health_indicator_var = NULL,
                       population_weights_var = NULL,
                       health_numerator_var = NULL,
                       health_denominator_var = NULL,
                       health_se_var = NULL,
                       analysis_unit_var = NULL,
                       equity_stratifier_var = NULL,
                       higher_ineq_is_favorable = FALSE,
                       rate_scaling_factor = NULL,
                       grouping_approach = c("territorial", "weighted_cut", "weighted_midpoint", "fractional"),
                       territorial_method = c("quantile", "kmeans", "fisher", "manual", "equal"),
                       manual_breaks = NULL,
                       n_groups = 5,
                       ci_method = c("delta", "bootstrap", "jackknife"),
                       n_boot = 2000,
                       conf_level = 0.95,
                       seed = NULL,
                       verbose = TRUE) {

  health_indicator_type <- match.arg(health_indicator_type)
  grouping_approach <- match.arg(grouping_approach)
  territorial_method <- if (identical(grouping_approach, "territorial")) {
    match.arg(territorial_method)
  } else {
    NA_character_
  }
  ci_method <- match.arg(ci_method)



  val <- .ag_rg_validate_inputs(
    data = data,
    health_indicator_type = health_indicator_type,
    health_indicator_var = health_indicator_var,
    population_weights_var = population_weights_var,
    health_numerator_var = health_numerator_var,
    health_denominator_var = health_denominator_var,
    health_se_var = health_se_var,
    analysis_unit_var = analysis_unit_var,
    equity_stratifier_var = equity_stratifier_var,
    higher_ineq_is_favorable = higher_ineq_is_favorable,
    rate_scaling_factor = rate_scaling_factor,
    grouping_approach = grouping_approach,
    territorial_method = territorial_method,
    n_groups = n_groups,
    ci_method = ci_method,
    n_boot = n_boot,
    conf_level = conf_level,
    verbose = verbose
  )

  df <- val$df
  scenario <- val$scenario
  final_multiplier <- val$final_multiplier
  variance_family <- val$variance_family
  use_counts <- val$use_counts
  ci_method_used <- val$ci_method_requested

  grouper <- function(d) {
    .ag_rg_group_data(
      df = d,
      use_counts = use_counts,
      grouping_approach = grouping_approach,
      territorial_method = territorial_method,
      manual_breaks = manual_breaks,
      n_groups = n_groups
    )
  }

  grouped_df <- grouper(df)

  lower_group <- attr(grouped_df, "lower_group")
  upper_group <- attr(grouped_df, "upper_group")
  n_groups_observed <- attr(grouped_df, "n_groups_observed")

  est <- .rg_estimate_gap(
    grouped_df = grouped_df,
    scenario = scenario,
    use_counts = use_counts,
    variance_family = variance_family,
    final_multiplier = final_multiplier,
    conf_level = conf_level
  )

  point_gap <- est$point_gap
  se_out <- est$se_out
  ci_low <- est$ci_low
  ci_high <- est$ci_high
  disadvantaged_value <- est$disadvantaged_value
  advantaged_value <- est$advantaged_value

  interval_interpretation <- if (identical(ci_method_used, "delta")) {
    if (identical(scenario, "counts_no_se")) {
      "inferential_delta_from_counts_log_ratio_scale"
    } else if (scenario %in% c("counts_with_se", "indicator_population_with_se")) {
      "inferential_delta_from_provided_se_log_ratio_scale"
    } else {
      "descriptive_point_estimate_only"
    }
  } else {
    NA_character_
  }

  res <- .ag_rg_resample(
    grouped_df = grouped_df,
    base_units = df,
    scenario = scenario,
    variance_family = variance_family,
    ci_method_used = ci_method_used,
    n_boot = n_boot,
    conf_level = conf_level,
    metric_multiplier = est$metric_multiplier,
    metric_calculator = est$metric_calculator,
    grouper = grouper,
    seed = seed
  )

  if (ci_method_used %in% c("bootstrap", "jackknife")) {
    se_out <- res$se_out
    ci_low <- res$ci_low
    ci_high <- res$ci_high

    if (identical(ci_method_used, "jackknife") && is.finite(point_gap)) {
      jk <- res$resampling_results$estimate
      jk <- jk[is.finite(jk) & jk > 0]
      if (length(jk) >= 3) {
        log_jk <- log(jk)
        se_log <- sqrt((length(log_jk) - 1) / length(log_jk) * sum((log_jk - mean(log_jk))^2))
        z <- stats::qnorm(1 - (1 - conf_level) / 2)
        ci_low <- exp(log(point_gap) - z * se_log)
        ci_high <- exp(log(point_gap) + z * se_log)
      }
    }

    interval_interpretation <- if (
      identical(ci_method_used, "bootstrap") &&
        identical(scenario, "counts_no_se")
    ) {
      "inferential_parametric_bootstrap_from_counts"
    } else {
      "ecological_sensitivity"
    }
  }

  results_groups <- .ag_rg_group_summary(
    grouped_df = grouped_df,
    scenario = scenario,
    use_counts = use_counts,
    final_multiplier = final_multiplier,
    variance_family = variance_family,
    ci_method_used = ci_method_used,
    interval_interpretation = interval_interpretation
  )

  out <- .ag_rg_format_outputs(
    metric_name = "relative_ecological_gap",
    point_gap = point_gap,
    se_out = se_out,
    ci_low = ci_low,
    ci_high = ci_high,
    conf_level = conf_level,
    ci_method_used = ci_method_used,
    interval_interpretation = interval_interpretation,
    lower_group = lower_group,
    upper_group = upper_group,
    disadvantaged_value = disadvantaged_value,
    advantaged_value = advantaged_value,
    grouped_df = grouped_df,
    resampling_results = res$resampling_results,
    n_groups = n_groups_observed,
    n_original = val$n_original,
    grouping_approach = grouping_approach,
    territorial_method = territorial_method,
    scenario = scenario,
    use_counts = use_counts,
    variance_family = variance_family,
    higher_ineq_is_favorable = higher_ineq_is_favorable,
    final_multiplier = final_multiplier,
    log_message = val$log_message,
    log_warning = c(val$log_warning, res$log_warning),
    results_groups = results_groups,
    methodological_notes_specific = "The relative gap is computed as disadvantaged group value divided by advantaged group value. Delta and jackknife intervals for relative gaps are computed on the log-ratio scale when possible."
  )

  class(out) <- c("rg_ineqeco", "ag_rg_ineqeco", class(out))
  out
}
