#' Internal Point and Delta-Interval Estimator for Absolute Ecological Gaps
#'
#' @description
#' Core estimation engine for the Absolute Gap (AG), invoked post-validation and
#' ecological grouping. It computes the absolute difference between the most
#' disadvantaged and most advantaged groups: \eqn{AG = Y_D - Y_A}.
#'
#' Variance estimation strictly follows the `health_se_var` pathway:
#' - If user-supplied standard errors are provided, group-level variance is pooled
#'   via weighted independent-error propagation.
#' - If raw event counts are provided, binomial or Poisson delta-method standard
#'   errors are computed based on the specified variance family.
#' - If only descriptive rates/populations are provided without variance priors,
#'   analytical SE propagation is bypassed, deferring to resampling inference if requested.
#'
#' @param grouped_df A structured data frame returned by `.ag_rg_group_data()`.
#' @param scenario Character denoting the data availability scenario.
#' @param use_counts Logical indicating whether raw events and denominators are used.
#' @param variance_family Character denoting the assumed underlying distribution (`"binomial"` or `"poisson"`).
#' @param final_multiplier Numeric scaling factor applied strictly to the point estimate vector.
#' @param conf_level Numeric denoting the confidence bound (default: 0.95).
#'
#' @return A complex list of parameter estimates, bounds, and closures for resampling.
#'
#' @importFrom stats weighted.mean qnorm
#' @noRd

.ag_estimate_gap <- function(grouped_df,
                             scenario,
                             use_counts,
                             variance_family = c("binomial", "poisson"),
                             final_multiplier = 1,
                             conf_level = 0.95) {

  variance_family <- match.arg(variance_family)

  if (missing(grouped_df) || !is.data.frame(grouped_df)) {
    stop("`grouped_df` must be a data frame returned by .ag_rg_group_data().", call. = FALSE)
  }

  valid_groups <- sort(unique(grouped_df$group[!is.na(grouped_df$group)]))
  if (length(valid_groups) < 2) {
    stop("At least two valid ecological groups are required.", call. = FALSE)
  }

  lower_group <- attr(grouped_df, "lower_group")
  upper_group <- attr(grouped_df, "upper_group")
  if (is.null(lower_group) || is.null(upper_group)) {
    lower_group <- min(valid_groups)
    upper_group <- max(valid_groups)
  }

  group_estimate_internal <- function(gdata, group_val) {
    g <- gdata[gdata$group == group_val, , drop = FALSE]
    if (nrow(g) == 0) return(NA_real_)

    if (isTRUE(use_counts)) {
      den <- sum(g$allocated_denominator, na.rm = TRUE)
      if (!is.finite(den) || den <= 0) return(NA_real_)
      num <- sum(g$allocated_numerator, na.rm = TRUE)
      return(num / den)
    }

    w <- g$allocated_population
    y <- g$indicator_internal
    if (length(w) == 0 || all(is.na(w)) || sum(w, na.rm = TRUE) <= 0) return(NA_real_)
    stats::weighted.mean(y, w = w, na.rm = TRUE)
  }

  group_se_internal_delta <- function(gdata, group_val) {
    g <- gdata[gdata$group == group_val, , drop = FALSE]
    if (nrow(g) == 0) return(NA_real_)

    if (identical(scenario, "counts_no_se")) {
      den <- sum(g$allocated_denominator, na.rm = TRUE)
      num <- sum(g$allocated_numerator, na.rm = TRUE)
      if (!is.finite(den) || den <= 0) return(NA_real_)
      p_hat <- num / den

      if (identical(variance_family, "binomial")) {
        return(sqrt(p_hat * (1 - p_hat) / den))
      }

      return(sqrt(p_hat / den))
    }

    if (scenario %in% c("counts_with_se", "indicator_population_with_se")) {
      w <- g$allocated_population
      v <- (g$se_internal)^2

      if (length(w) == 0 || all(is.na(w)) || sum(w, na.rm = TRUE) <= 0) return(NA_real_)
      if (any(is.na(v))) return(NA_real_)

      ww <- w / sum(w, na.rm = TRUE)
      return(sqrt(sum((ww^2) * v, na.rm = TRUE)))
    }

    NA_real_
  }

  y_d_internal <- group_estimate_internal(grouped_df, lower_group)
  y_a_internal <- group_estimate_internal(grouped_df, upper_group)

  y_d <- y_d_internal * final_multiplier
  y_a <- y_a_internal * final_multiplier

  se_d_internal <- group_se_internal_delta(grouped_df, lower_group)
  se_a_internal <- group_se_internal_delta(grouped_df, upper_group)

  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  point <- y_d - y_a
  metric_multiplier <- final_multiplier

  if (is.finite(se_d_internal) && is.finite(se_a_internal)) {
    se_out <- sqrt(se_d_internal^2 + se_a_internal^2) * final_multiplier
    ci_low <- point - z * se_out
    ci_high <- point + z * se_out
  } else {
    se_out <- NA_real_
    ci_low <- NA_real_
    ci_high <- NA_real_
  }

  metric_calculator <- function(gdata) {
    group_estimate_internal(gdata, min(gdata$group, na.rm = TRUE)) -
      group_estimate_internal(gdata, max(gdata$group, na.rm = TRUE))
  }

  interval_scale <- "difference_scale"
  standard_error_log <- NA_real_
  log_estimate <- NA_real_
  log_ci_low <- NA_real_
  log_ci_high <- NA_real_

  list(
    point_gap = as.numeric(point),
    se_out = as.numeric(se_out),
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    disadvantaged_value = as.numeric(y_d),
    advantaged_value = as.numeric(y_a),
    disadvantaged_value_internal = as.numeric(y_d_internal),
    advantaged_value_internal = as.numeric(y_a_internal),
    disadvantaged_se_internal = as.numeric(se_d_internal),
    advantaged_se_internal = as.numeric(se_a_internal),
    log_estimate = as.numeric(log_estimate),
    standard_error_log = as.numeric(standard_error_log),
    log_ci_low = as.numeric(log_ci_low),
    log_ci_high = as.numeric(log_ci_high),
    metric_calculator = metric_calculator,
    metric_multiplier = metric_multiplier,
    interval_scale = interval_scale
  )
}
