#' Internal Ecological Grouping and Summary Constructor
#'
#' @description
#' `ag_rg_group_data` implements the core aggregation and grouping algorithms for ecological data.
#' It maps continuous ecological units onto a discrete ordinal scale based on their
#' socioeconomic stratifier. Supported algorithms include unweighted territorial methods
#' (`quantile`, `kmeans`, `fisher`), user-defined deterministic bounds (`manual`),
#' population-weighted empirical quantiles (`weighted_cut`), and fractional population allocation (`fractional`).
#'
#' `ag_rg_group_summary` subsequently compiles these assignments, applying the specified
#' variance family and generating group-level point estimates and propagated delta-method standard errors.
#'
#' @details
#' These functions assume pre-validated and standardized inputs from `.ag_rg_validate_inputs`.
#'
#' @importFrom Hmisc wtd.quantile
#' @importFrom classInt classIntervals
#' @noRd

.ag_rg_group_data <- function(df,
                                use_counts,
                                grouping_approach = c("territorial", "weighted_cut", "weighted_midpoint", "fractional"),
                                territorial_method = c("quantile", "kmeans", "fisher", "manual", "equal"),
                                manual_breaks = NULL,
                                n_groups = 5) {

  # ---- Basic validation -------------------------------------------------------
  if (missing(df) || !is.data.frame(df)) {
    stop("`df` must be a data frame with standardized ecological-unit columns.", call. = FALSE)
  }

  if (missing(use_counts) || !is.logical(use_counts) || length(use_counts) != 1 || is.na(use_counts)) {
    stop("`use_counts` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  grouping_approach <- match.arg(grouping_approach)

  if (grouping_approach == "territorial") {
    if (is.null(territorial_method)) {
      stop("`territorial_method` cannot be NULL when `grouping_approach = 'territorial'`.", call. = FALSE)
    }
    territorial_method <- match.arg(territorial_method)
  } else {
    territorial_method <- NA_character_
  }

  n_groups <- as.integer(n_groups)
  if (is.na(n_groups) || n_groups < 2) {
    stop("`n_groups` must be an integer >= 2.", call. = FALSE)
  }

  required_cols <- c(
    "unit_id", "population_weight", "equity_stratifier", "strat_order",
    "indicator_internal", "health_indicator", "health_se", "se_internal"
  )
  if (isTRUE(use_counts)) {
    required_cols <- c(required_cols, "numerator", "denominator")
  }

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "`df` is missing required standardized column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(df$population_weight) || any(df$population_weight <= 0, na.rm = TRUE)) {
    stop("`df$population_weight` must be numeric and strictly positive for all included units.", call. = FALSE)
  }
  if (!is.numeric(df$strat_order)) {
    stop("`df$strat_order` must be numeric.", call. = FALSE)
  }

  # ---- Internal weighted quantile -------------------------------------------
  # Kept inside the grouping constructor because it is only used for group
  # construction and should not become a separate package-level helper.
  weighted_quantile <- function(x, w, probs) {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      stop(
        "Package `Hmisc` is required for `grouping_approach = 'weighted_cut'`. ",
        "Install it with install.packages('Hmisc').",
        call. = FALSE
      )
    }
    as.numeric(Hmisc::wtd.quantile(x = x, weights = w, probs = probs, na.rm = TRUE))
  }

  # ---- Order units and compute population positions --------------------------
  d <- df[order(df$strat_order, df$unit_id), , drop = FALSE]
  rownames(d) <- NULL

  pop_total <- sum(d$population_weight, na.rm = TRUE)
  if (!is.finite(pop_total) || pop_total <= 0) {
    stop("Total population/reference weight must be positive.", call. = FALSE)
  }

  d$population_share <- d$population_weight / pop_total
  d$cum_population_share <- cumsum(d$population_share)
  d$prev_cum_population_share <- c(0, head(d$cum_population_share, -1))
  d$ridit <- d$prev_cum_population_share + d$population_share / 2

  # ---- Territorial grouping --------------------------------------------------
  if (grouping_approach == "territorial") {
    if (!requireNamespace("classInt", quietly = TRUE)) {
      stop(
        "Package `classInt` is required for territorial grouping. ",
        "Install it with install.packages('classInt').",
        call. = FALSE
      )
    }

    style <- switch(
      territorial_method,
      quantile = "quantile",
      kmeans   = "kmeans",
      fisher   = "fisher",
      equal    = "equal",
      manual   = "fixed"
    )

    if (territorial_method == "manual") {
      if (is.null(manual_breaks)) {
        stop("`manual_breaks` is required when `territorial_method = 'manual'`.", call. = FALSE)
      }
      breaks_in <- sort(unique(as.numeric(manual_breaks)))
      if (length(breaks_in) < 1 || any(!is.finite(breaks_in))) {
        stop("`manual_breaks` must contain finite numeric cut points.", call. = FALSE)
      }

      if (length(breaks_in) == n_groups + 1) {
        breaks <- breaks_in
      } else {
        breaks <- c(-Inf, breaks_in, Inf)
      }
      actual_groups <- length(breaks) - 1
    } else {
      ci <- tryCatch(
        classInt::classIntervals(d$strat_order, n = n_groups, style = style),
        error = function(e) {
          stop(sprintf("Unable to create territorial groups: %s", e$message), call. = FALSE)
        }
      )
      breaks <- ci$brks
      actual_groups <- length(breaks) - 1
    }

    if (length(unique(breaks)) < length(breaks)) {
      stop(
        "Grouping failed because duplicated break values were generated. ",
        "Try fewer groups or a different grouping method.",
        call. = FALSE
      )
    }

    d$group <- as.integer(cut(
      d$strat_order,
      breaks = breaks,
      include.lowest = TRUE,
      labels = seq_len(actual_groups)
    ))

    d$allocation_fraction <- 1
    d$allocated_population <- d$population_weight
    d$allocated_denominator <- if (isTRUE(use_counts)) d$denominator else NA_real_
    d$allocated_numerator <- if (isTRUE(use_counts)) d$numerator else NA_real_
    d$grouping_break_low <- breaks[pmax(d$group, 1)]
    d$grouping_break_high <- breaks[pmax(d$group, 1) + 1]
    d$grouping_method_detail <- territorial_method
  }

  # ---- Weighted cut grouping -------------------------------------------------
  if (grouping_approach == "weighted_cut") {
    probs <- seq(0, 1, length.out = n_groups + 1)
    breaks <- unique(weighted_quantile(d$strat_order, d$population_weight, probs))

    if (length(breaks) < n_groups + 1) {
      stop(
        "Weighted cut failed: not enough unique weighted break values to form ",
        "the requested number of groups.",
        call. = FALSE
      )
    }

    d$group <- as.integer(cut(
      d$strat_order,
      breaks = breaks,
      include.lowest = TRUE,
      labels = seq_len(n_groups)
    ))

    d$allocation_fraction <- 1
    d$allocated_population <- d$population_weight
    d$allocated_denominator <- if (isTRUE(use_counts)) d$denominator else NA_real_
    d$allocated_numerator <- if (isTRUE(use_counts)) d$numerator else NA_real_
    d$grouping_break_low <- breaks[pmax(d$group, 1)]
    d$grouping_break_high <- breaks[pmax(d$group, 1) + 1]
    d$grouping_method_detail <- "weighted_quantile_cut"
  }

  # ---- Weighted midpoint grouping -------------------------------------------
  if (grouping_approach == "weighted_midpoint") {
    d$group <- floor(d$ridit * n_groups) + 1
    d$group <- pmin(pmax(d$group, 1), n_groups)

    d$allocation_fraction <- 1
    d$allocated_population <- d$population_weight
    d$allocated_denominator <- if (isTRUE(use_counts)) d$denominator else NA_real_
    d$allocated_numerator <- if (isTRUE(use_counts)) d$numerator else NA_real_
    d$grouping_break_low <- (d$group - 1) / n_groups
    d$grouping_break_high <- d$group / n_groups
    d$grouping_method_detail <- "population_weighted_midpoint_ridit"
  }

  # ---- Fractional grouping ---------------------------------------------------
  if (grouping_approach == "fractional") {
    group_limits <- data.frame(
      group = seq_len(n_groups),
      group_lower = (seq_len(n_groups) - 1) / n_groups,
      group_upper = seq_len(n_groups) / n_groups
    )

    idx <- expand.grid(row_id = seq_len(nrow(d)), group = seq_len(n_groups))
    out <- d[idx$row_id, , drop = FALSE]
    gl <- group_limits[idx$group, , drop = FALSE]

    overlap <- pmax(
      0,
      pmin(out$cum_population_share, gl$group_upper) -
        pmax(out$prev_cum_population_share, gl$group_lower)
    )

    fraction <- ifelse(out$population_share > 0, overlap / out$population_share, 0)

    out$group <- gl$group
    out$allocation_fraction <- fraction
    out$allocated_population <- out$population_weight * fraction
    out$allocated_denominator <- if (isTRUE(use_counts)) out$denominator * fraction else NA_real_
    out$allocated_numerator <- if (isTRUE(use_counts)) out$numerator * fraction else NA_real_
    out$grouping_break_low <- gl$group_lower
    out$grouping_break_high <- gl$group_upper
    out$grouping_method_detail <- "fractional_population_allocation"

    d <- out[out$allocation_fraction > 0, , drop = FALSE]
    rownames(d) <- NULL
  }

  # ---- Harmonize returned columns -------------------------------------------
  d$grouping_approach <- grouping_approach

  required_output_cols <- c(
    "unit_id", "group", "grouping_approach", "grouping_method_detail",
    "allocation_fraction", "allocated_population", "allocated_numerator",
    "allocated_denominator", "population_weight", "population_share", "ridit",
    "equity_stratifier", "strat_order", "indicator_internal", "health_indicator",
    "health_se", "se_internal", "numerator", "denominator", "grouping_break_low",
    "grouping_break_high", "cum_population_share", "prev_cum_population_share"
  )

  for (nm in setdiff(required_output_cols, names(d))) {
    d[[nm]] <- NA_real_
  }

  valid_groups <- sort(unique(d$group[!is.na(d$group)]))
  if (length(valid_groups) < 2) {
    stop("Grouping resulted in fewer than two valid groups.", call. = FALSE)
  }

  attr(d, "valid_groups") <- valid_groups
  attr(d, "lower_group") <- min(valid_groups)
  attr(d, "upper_group") <- max(valid_groups)
  attr(d, "grouping_approach") <- grouping_approach
  attr(d, "territorial_method") <- territorial_method
  attr(d, "n_groups_requested") <- n_groups
  attr(d, "n_groups_observed") <- length(valid_groups)

  d
}


# Internal ecological group summary constructor
#
# This function is intended for use by ag_ineqeco() and rg_ineqeco() after
# .ag_rg_group_data(). It estimates group-level health indicators and, when
# possible, group-level standard errors. It does not compute the final absolute
# or relative inequality metric; those remain visible inside the public functions.

.ag_rg_group_summary <- function(grouped_df,
                                   scenario,
                                   use_counts,
                                   final_multiplier = 1,
                                   variance_family = c("binomial", "poisson"),
                                   ci_method_used = NA_character_,
                                   interval_interpretation = NA_character_) {

  # ---- Basic validation -------------------------------------------------------
  if (missing(grouped_df) || !is.data.frame(grouped_df)) {
    stop("`grouped_df` must be a data frame returned by .ag_rg_group_data().", call. = FALSE)
  }

  if (missing(use_counts) || !is.logical(use_counts) || length(use_counts) != 1 || is.na(use_counts)) {
    stop("`use_counts` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  allowed_scenarios <- c(
    "counts_no_se", "counts_with_se",
    "indicator_population_no_se", "indicator_population_with_se"
  )
  if (missing(scenario) || !(scenario %in% allowed_scenarios)) {
    stop("`scenario` must be one of: counts_no_se, counts_with_se, indicator_population_no_se, indicator_population_with_se.", call. = FALSE)
  }

  variance_family <- match.arg(variance_family)

  if (!is.numeric(final_multiplier) || length(final_multiplier) != 1 || !is.finite(final_multiplier) || final_multiplier <= 0) {
    stop("`final_multiplier` must be a positive numeric scalar.", call. = FALSE)
  }

  required_cols <- c(
    "unit_id", "group", "allocated_population", "allocation_fraction",
    "indicator_internal", "se_internal"
  )
  if (isTRUE(use_counts)) {
    required_cols <- c(required_cols, "allocated_numerator", "allocated_denominator")
  }

  missing_cols <- setdiff(required_cols, names(grouped_df))
  if (length(missing_cols) > 0) {
    stop(
      "`grouped_df` is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  groups <- sort(unique(grouped_df$group[!is.na(grouped_df$group)]))
  if (length(groups) < 2) {
    stop("At least two valid groups are required to build group summaries.", call. = FALSE)
  }

  lower_group <- min(groups)
  upper_group <- max(groups)
  total_allocated_population <- sum(grouped_df$allocated_population, na.rm = TRUE)

  if (!is.finite(total_allocated_population) || total_allocated_population <= 0) {
    stop("Total allocated population/reference weight must be positive.", call. = FALSE)
  }

  # ---- Local group estimators -------------------------------------------------
  # Kept local because they are part of the group summary calculation and do not
  # need to exist as separate package-level helpers.
  group_estimate_internal <- function(gdata, group_value) {
    sub <- gdata[gdata$group == group_value, , drop = FALSE]
    if (nrow(sub) == 0) return(NA_real_)

    if (isTRUE(use_counts)) {
      num <- sum(sub$allocated_numerator, na.rm = TRUE)
      den <- sum(sub$allocated_denominator, na.rm = TRUE)
      if (!is.finite(den) || den <= 0) return(NA_real_)
      return(num / den)
    }

    w <- sub$allocated_population
    y <- sub$indicator_internal
    ok <- is.finite(w) & is.finite(y) & w > 0
    if (!any(ok)) return(NA_real_)
    sum(w[ok] * y[ok], na.rm = TRUE) / sum(w[ok], na.rm = TRUE)
  }

  group_se_internal_delta <- function(gdata, group_value) {
    sub <- gdata[gdata$group == group_value, , drop = FALSE]
    if (nrow(sub) == 0) return(NA_real_)

    if (scenario == "indicator_population_no_se") {
      return(NA_real_)
    }

    if (scenario %in% c("counts_with_se", "indicator_population_with_se")) {
      W <- sum(sub$allocated_population, na.rm = TRUE)
      if (!is.finite(W) || W <= 0) return(NA_real_)
      a <- sub$allocated_population / W
      return(sqrt(sum((a^2) * (sub$se_internal^2), na.rm = TRUE)))
    }

    # counts_no_se: analytical approximation from aggregated counts.
    num <- sum(sub$allocated_numerator, na.rm = TRUE)
    den <- sum(sub$allocated_denominator, na.rm = TRUE)
    if (!is.finite(den) || den <= 0) return(NA_real_)

    y <- num / den
    if (!is.finite(y)) return(NA_real_)

    if (variance_family == "binomial") {
      if (y < 0 || y > 1) return(NA_real_)
      return(sqrt(y * (1 - y) / den))
    }

    sqrt(num) / den
  }

  # ---- Build table ------------------------------------------------------------
  rows <- lapply(groups, function(g) {
    sub <- grouped_df[grouped_df$group == g, , drop = FALSE]

    numerator <- if (isTRUE(use_counts)) sum(sub$allocated_numerator, na.rm = TRUE) else NA_real_
    denominator <- if (isTRUE(use_counts)) sum(sub$allocated_denominator, na.rm = TRUE) else NA_real_
    pop <- sum(sub$allocated_population, na.rm = TRUE)

    health_internal <- group_estimate_internal(grouped_df, g)
    se_internal <- group_se_internal_delta(grouped_df, g)

    data.frame(
      group = g,
      group_role = ifelse(
        g == lower_group,
        "most_disadvantaged",
        ifelse(g == upper_group, "most_advantaged", "intermediate")
      ),
      grouping_approach = if ("grouping_approach" %in% names(grouped_df)) {
        unique(as.character(sub$grouping_approach))[1]
      } else {
        NA_character_
      },
      grouping_method_detail = if ("grouping_method_detail" %in% names(grouped_df)) {
        unique(as.character(sub$grouping_method_detail))[1]
      } else {
        NA_character_
      },
      n_units_contributing = length(unique(sub$unit_id)),
      fractional_units = sum(sub$allocation_fraction, na.rm = TRUE),
      population_weight = pop,
      population_share = pop / total_allocated_population,
      numerator = numerator,
      denominator = denominator,
      health_estimate = health_internal * final_multiplier,
      standard_error = se_internal * final_multiplier,
      ci_method_used = ci_method_used,
      interval_interpretation = interval_interpretation,
      units = paste(unique(sub$unit_id), collapse = ", "),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  attr(out, "lower_group") <- lower_group
  attr(out, "upper_group") <- upper_group
  attr(out, "scenario") <- scenario
  attr(out, "use_counts") <- use_counts
  attr(out, "final_multiplier") <- final_multiplier
  attr(out, "variance_family") <- variance_family

  out
}
