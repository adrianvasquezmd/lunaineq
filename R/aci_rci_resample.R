#' Iterative Variance Estimation Engine for ACI/RCI Families
#'
#' @description
#' `aci_rci_resample` governs the stochastic non-parametric replication topology
#' required for bootstrapped and jackknifed standard errors of the Concentration Indices.
#'
#' The engine operates recursively across the ecological unit vector, strictly maintaining
#' grouped-level integrity while reconstructing ridits and empirical covariances across each
#' stochastically perturbed set.
#'
#' Resampling is performed over ecological units/groups, not individuals.
#' Smoothed concentration curves are not recomputed and are not used by this
#' function.
#'
#' Important behavior:
#' \itemize{
#'   \item Bootstrap resamples ecological units/groups with replacement.
#'   \item Jackknife removes one ecological unit/group at a time.
#'   \item `health_se_var`, if supplied, is not used directly by non-parametric
#'     bootstrap or jackknife. It is retained in resampled data only for
#'     traceability and possible downstream diagnostics.
#'   \item Each replicate rebuilds population shares, fractional ranks,
#'     cumulative shares, RCI, ACI, and bounded-outcome corrections.
#' }
#'
#' @param prepared A list returned by `.aci_rci_prepare_data()`.
#' @param config A validated configuration list returned by
#'   `.aci_rci_validate_inputs()`.
#' @param method Character. Either `"bootstrap"` or `"jackknife"`.
#' @param n_boot Integer. Number of bootstrap replicates when
#'   `method = "bootstrap"`.
#' @param seed Integer or NULL. Optional random seed for bootstrap
#'   reproducibility. Ignored for jackknife.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list with:
#' \describe{
#'   \item{replicates}{A data frame with one row per valid replicate and one
#'     column per metric.}
#'   \item{warnings}{Warnings generated during resampling.}
#'   \item{diagnostics}{Counts of requested, valid, failed, and skipped
#'     replicates.}
#' }
#'
#' @importFrom stats setNames
#' @importFrom utils head
#' @keywords internal
#' @noRd
#'
.aci_rci_resample <- function(
  prepared,
  config,
  method = c("bootstrap", "jackknife"),
  n_boot = NULL,
  seed = NULL,
  verbose = FALSE
) {
  warnings_out <- character()
  diagnostics <- list()

  add_warning <- function(x) {
    warnings_out <<- c(warnings_out, x)
    if (isTRUE(config$verbose)) warning(x, call. = FALSE, immediate. = TRUE)
    invisible(NULL)
  }

  stop_resample <- function(x) {
    stop(x, call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Basic checks
  # -------------------------------------------------------------------------

  method <- match.arg(method)

  if (!is.list(prepared)) {
    stop_resample("`prepared` must be the list returned by `.aci_rci_prepare_data()`.")
  }

  if (!is.list(config)) {
    stop_resample("`config` must be the list returned by `.aci_rci_validate_inputs()`.")
  }

  if (is.null(prepared$data_units) || !is.data.frame(prepared$data_units)) {
    stop_resample("`prepared$data_units` must be a data frame.")
  }

  data_units <- prepared$data_units

  required_columns <- c(
    ".row_id",
    ".unit",
    ".health",
    ".health_se",
    ".population_weight",
    ".equity",
    ".equity_order"
  )

  missing_columns <- setdiff(required_columns, names(data_units))

  if (length(missing_columns) > 0L) {
    stop_resample(
      paste0(
        "`prepared$data_units` is missing required internal column(s): ",
        paste(missing_columns, collapse = ", "),
        ". Please run `.aci_rci_prepare_data()` before `.aci_rci_resample()`."
      )
    )
  }

  n_units <- nrow(data_units)

  if (n_units < 2L) {
    stop_resample("At least two ecological units/groups are required for resampling.")
  }

  if (method == "jackknife" && n_units < 3L) {
    stop_resample("Jackknife resampling requires at least three ecological units/groups.")
  }

  if (method == "bootstrap") {
    if (is.null(n_boot)) {
      n_boot <- config$n_boot
    }

    if (is.null(n_boot) || !is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot)) {
      stop_resample("`n_boot` must be a single numeric value when `method = \"bootstrap\"`.")
    }

    n_boot <- as.integer(n_boot)

    if (n_boot < 1L) {
      stop_resample("`n_boot` must be at least 1 when `method = \"bootstrap\"`.")
    }

    if (n_boot < 100L) {
      add_warning(
        "`n_boot` is less than 100. Bootstrap confidence intervals may be unstable."
      )
    }
  } else {
    n_boot <- NA_integer_
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop_resample("`verbose` must be a single TRUE/FALSE value.")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop_resample("`seed` must be NULL or a single numeric value.")
    }

    if (method == "bootstrap") {
      set.seed(as.integer(seed))
    }
  }

  if (!is.null(config$health_se_var)) {
    add_warning(
      paste0(
        "`health_se_var` is present, but non-parametric ",
        method,
        " resampling does not use ecological-unit standard errors directly. ",
        "Uncertainty is quantified by resampling ecological units/groups."
      )
    )
  }

  if (method == "bootstrap") {
    add_warning(
      "Bootstrap resampling will sample ecological units/groups with replacement and reconstruct ranks and cumulative shares in each replicate."
    )
  } else {
    add_warning(
      "Jackknife resampling will remove one ecological unit/group at a time and reconstruct ranks and cumulative shares in each replicate."
    )
  }

  # -------------------------------------------------------------------------
  # Helper to rebuild prepared object from already-standardized internal data
  # -------------------------------------------------------------------------

  rebuild_prepared_from_units <- function(unit_data) {
    if (!is.data.frame(unit_data)) {
      stop("Internal resample data must be a data frame.", call. = FALSE)
    }

    if (nrow(unit_data) < 2L) {
      stop("Fewer than two units in replicate.", call. = FALSE)
    }

    if (any(!is.finite(unit_data$.health))) {
      stop("Non-finite health values in replicate.", call. = FALSE)
    }

    if (any(!is.finite(unit_data$.population_weight)) ||
        any(unit_data$.population_weight <= 0)) {
      stop("Invalid population weights in replicate.", call. = FALSE)
    }

    if (any(!is.finite(unit_data$.equity_order))) {
      stop("Non-finite equity order values in replicate.", call. = FALSE)
    }

    # Stable ordering after resampling. For bootstrap duplicates, `.resample_row`
    # preserves the replicate row order so repeated units remain distinct rows.
    if (!".resample_row" %in% names(unit_data)) {
      unit_data$.resample_row <- seq_len(nrow(unit_data))
    }

    order_index <- order(unit_data$.equity_order, unit_data$.resample_row, na.last = NA)
    unit_data <- unit_data[order_index, , drop = FALSE]
    rownames(unit_data) <- NULL

    total_weight <- sum(unit_data$.population_weight)

    if (!is.finite(total_weight) || total_weight <= 0) {
      stop("Total population/exposure weight is not positive in replicate.", call. = FALSE)
    }

    weighted_health_total <- sum(unit_data$.population_weight * unit_data$.health)
    mean_health <- weighted_health_total / total_weight

    if (!is.finite(weighted_health_total) || !is.finite(mean_health) || mean_health <= 0) {
      stop("Replicate produced non-positive or non-finite mean health.", call. = FALSE)
    }

    unit_data$.population_share <- unit_data$.population_weight / total_weight
    unit_data$.cum_population <- cumsum(unit_data$.population_share)
    unit_data$.cum_population_lag <- c(0, utils::head(unit_data$.cum_population, -1L))

    # Fractional rank with tie handling harmonized with
    # `.aci_rci_prepare_data()` and the original luna implementation.
    unit_data$.fractional_rank <- NA_real_
    tie_groups <- split(seq_len(nrow(unit_data)), unit_data$.equity_order)

    for (idx in tie_groups) {
      tied_pop_start <- min(unit_data$.cum_population_lag[idx])
      tied_pop_end <- max(unit_data$.cum_population[idx])
      unit_data$.fractional_rank[idx] <- tied_pop_start + 0.5 * (tied_pop_end - tied_pop_start)
    }

    unit_data$.weighted_health <- unit_data$.population_weight * unit_data$.health
    unit_data$.health_share <- unit_data$.weighted_health / weighted_health_total
    unit_data$.cum_health <- cumsum(unit_data$.health_share)
    unit_data$.cum_health_lag <- c(0, utils::head(unit_data$.cum_health, -1L))

    curve_data <- data.frame(
      .curve_point = seq_len(nrow(unit_data) + 1L),
      .cum_population = c(0, unit_data$.cum_population),
      .cum_health_empirical = c(0, unit_data$.cum_health),
      stringsAsFactors = FALSE
    )

    summary <- list(
      n_units = nrow(unit_data),
      total_population_weight = total_weight,
      weighted_health_total = weighted_health_total,
      mean_health = mean_health,
      health_indicator_internal_scale = if (!is.null(prepared$summary$health_indicator_internal_scale)) {
        prepared$summary$health_indicator_internal_scale
      } else {
        NA_character_
      },
      health_lower_bound = config$health_lower_bound,
      health_upper_bound = config$health_upper_bound,
      bounded_outcome = isTRUE(config$bounded_outcome)
    )

    list(
      data_units = unit_data,
      curve_data = curve_data,
      summary = summary,
      warnings = character(),
      diagnostics = list(
        n_rows_prepared = nrow(unit_data),
        total_population_weight = total_weight,
        weighted_health_total = weighted_health_total,
        weighted_mean_health = mean_health
      )
    )
  }

  # -------------------------------------------------------------------------
  # Helper to compute all metrics in one replicate
  # -------------------------------------------------------------------------

  compute_replicate_metrics <- function(unit_data) {
    replicate_prepared <- rebuild_prepared_from_units(unit_data)
    u <- replicate_prepared$data_units
    m_h <- replicate_prepared$summary$mean_health
    
    rep_rci <- (2 / m_h) * sum(u$.population_share * u$.health * u$.fractional_rank) - 1
    rep_aci <- m_h * rep_rci
    
    rep_wagstaff <- NA_real_
    rep_erreygers <- NA_real_
    
    bounded_outcome <- isTRUE(config$bounded_outcome)
    lower_bound <- config$health_lower_bound
    upper_bound <- config$health_upper_bound
    correction_methods <- config$correction_methods
    if (is.null(correction_methods)) correction_methods <- character()
    correction_methods <- unique(tolower(correction_methods))
    
    if (bounded_outcome && length(correction_methods) > 0L) {
      if (!is.null(lower_bound) && !is.null(upper_bound)) {
        bounds_range <- as.numeric(upper_bound) - as.numeric(lower_bound)
        if ("wagstaff" %in% correction_methods) {
          denom <- (upper_bound - m_h) * (m_h - lower_bound)
          if (is.finite(denom) && abs(denom) > 1e-10) {
            rep_wagstaff <- (m_h * bounds_range / denom) * rep_rci
          }
        }
        if ("erreygers" %in% correction_methods) {
          if (is.finite(bounds_range) && bounds_range > 1e-10) {
            rep_erreygers <- (4 * m_h / bounds_range) * rep_rci
          }
        }
      }
    }

    out <- c(
      RCI = rep_rci,
      ACI = rep_aci,
      Wagstaff_Corrected_RCI = rep_wagstaff,
      Erreygers_Corrected_RCI = rep_erreygers,
      Mean_Health = m_h,
      N_Units = nrow(u),
      Total_Weight = replicate_prepared$summary$total_population_weight
    )

    as.numeric(out) |>
      stats::setNames(names(out))
  }

  # -------------------------------------------------------------------------
  # Generate replicate index sets
  # -------------------------------------------------------------------------

  if (method == "bootstrap") {
    n_requested <- n_boot
    replicate_indices <- vector("list", n_requested)

    for (b in seq_len(n_requested)) {
      replicate_indices[[b]] <- sample.int(n_units, size = n_units, replace = TRUE)
    }

    replicate_type <- "bootstrap"
  } else {
    n_requested <- n_units
    replicate_indices <- vector("list", n_requested)

    for (j in seq_len(n_requested)) {
      replicate_indices[[j]] <- setdiff(seq_len(n_units), j)
    }

    replicate_type <- "jackknife"
  }

  # -------------------------------------------------------------------------
  # Run replicates
  # -------------------------------------------------------------------------

  replicate_rows <- vector("list", n_requested)
  failure_messages <- character()
  failed_replicates <- integer()
  valid_count <- 0L

  for (r in seq_len(n_requested)) {
    if (verbose && (r == 1L || r %% 100L == 0L || r == n_requested)) {
      message(
        "ACI/RCI ", replicate_type, " replicate ",
        r, " of ", n_requested
      )
    }

    idx <- replicate_indices[[r]]
    unit_data <- data_units[idx, , drop = FALSE]
    unit_data$.resample_row <- seq_len(nrow(unit_data))

    metric_values <- tryCatch(
      compute_replicate_metrics(unit_data),
      error = function(e) {
        failure_messages <<- c(failure_messages, conditionMessage(e))
        failed_replicates <<- c(failed_replicates, r)
        NULL
      }
    )

    if (is.null(metric_values)) {
      next
    }

    valid_count <- valid_count + 1L

    replicate_rows[[valid_count]] <- data.frame(
      Replicate_ID = r,
      Replicate_Type = replicate_type,
      RCI = metric_values[["RCI"]],
      ACI = metric_values[["ACI"]],
      Wagstaff_Corrected_RCI = metric_values[["Wagstaff_Corrected_RCI"]],
      Erreygers_Corrected_RCI = metric_values[["Erreygers_Corrected_RCI"]],
      Mean_Health = metric_values[["Mean_Health"]],
      N_Units = metric_values[["N_Units"]],
      Total_Weight = metric_values[["Total_Weight"]],
      stringsAsFactors = FALSE
    )
  }

  if (valid_count == 0L) {
    stop_resample(
      paste0(
        "All ", replicate_type,
        " replicates failed. Resampling-based confidence intervals cannot be computed."
      )
    )
  }

  replicates <- do.call(rbind, replicate_rows[seq_len(valid_count)])
  rownames(replicates) <- NULL

  n_failed <- length(failed_replicates)

  if (n_failed > 0L) {
    add_warning(
      paste0(
        n_failed,
        " ",
        replicate_type,
        " replicate(s) failed and were excluded from confidence interval calculations."
      )
    )
  }

  # Warning threshold: if more than 20% fail, resampling intervals may be
  # unreliable. This is intentionally a warning rather than a stop because some
  # bounded correction metrics may fail while base RCI/ACI remain valid.
  failure_rate <- n_failed / n_requested

  if (failure_rate > 0.20) {
    add_warning(
      paste0(
        "More than 20% of ",
        replicate_type,
        " replicates failed. Confidence intervals may be unstable. ",
        "Review data sparsity, zero/negative mean health, and bounded correction factors."
      )
    )
  }

  # Metric-specific valid counts. Some bounded corrections may be NA even when
  # the base replicate was valid.
  metric_columns <- c(
    "RCI",
    "ACI",
    "Wagstaff_Corrected_RCI",
    "Erreygers_Corrected_RCI"
  )

  metric_valid_counts <- stats::setNames(
    vapply(
      metric_columns,
      function(col) sum(is.finite(replicates[[col]])),
      numeric(1)
    ),
    metric_columns
  )

  for (metric_name in names(metric_valid_counts)) {
    if (metric_valid_counts[[metric_name]] == 0L) {
      add_warning(
        paste0(
          "No valid resampled estimates were available for metric `",
          metric_name,
          "`."
        )
      )
    } else if (metric_valid_counts[[metric_name]] < max(2L, floor(0.50 * valid_count))) {
      add_warning(
        paste0(
          "Fewer than half of valid replicates produced finite estimates for metric `",
          metric_name,
          "`. Confidence intervals for this metric may be unavailable or unstable."
        )
      )
    }
  }

  # -------------------------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------------------------

  diagnostics$method <- method
  diagnostics$replicate_type <- replicate_type
  diagnostics$n_units_original <- n_units
  diagnostics$n_replicates_requested <- n_requested
  diagnostics$n_replicates_valid <- valid_count
  diagnostics$n_replicates_failed <- n_failed
  diagnostics$failure_rate <- failure_rate
  diagnostics$failed_replicates <- failed_replicates
  diagnostics$unique_failure_messages <- unique(failure_messages)
  diagnostics$metric_valid_counts <- metric_valid_counts
  diagnostics$health_se_supplied <- !is.null(config$health_se_var)
  diagnostics$health_se_used_directly <- FALSE
  diagnostics$smoothed_curves_recomputed <- FALSE

  list(
    replicates = replicates,
    warnings = unique(warnings_out),
    diagnostics = diagnostics
  )
}
