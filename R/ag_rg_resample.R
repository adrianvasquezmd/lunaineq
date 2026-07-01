#' Resampling and Sensitivity Engine for AG/RG Ecological Gap Functions
#'
#' @description
#' Generates robust, non-parametric or parametric variance approximations for the Absolute and Relative Gaps.
#'
#' Methodological pathways:
#' - **Parametric Bootstrap (`counts_no_se`)**: Simulates grouped binomial/Poisson events using the provided population parameters and empirical probabilities.
#' - **Ecological Bootstrap**: Performs standard resampling with replacement across primary ecological units (rows), automatically re-running the grouping clustering algorithm per draw.
#' - **Delete-One Jackknife**: A leave-one-out variance approximation over ecological units, suited for smaller deterministic datasets.
#'
#' @importFrom stats quantile rbinom rpois sd qnorm
#' @importFrom boot boot boot.ci
#' @importFrom bootstrap jackknife
#' @noRd

.ag_rg_resample <- function(grouped_df,
                            base_units,
                            scenario,
                            variance_family,
                            ci_method_used,
                            n_boot,
                            conf_level,
                            metric_multiplier,
                            metric_calculator,
                            grouper,
                            seed = NULL) {

  se_out <- NA_real_
  ci_low <- NA_real_
  ci_high <- NA_real_
  boot_method_label <- NA_character_
  log_warning <- character(0)

  ci_percentile <- function(x, conf_level) {
    alpha <- 1 - conf_level
    as.numeric(stats::quantile(
      x,
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE,
      names = FALSE,
      type = 6
    ))
  }

  resampling_results <- data.frame(
    replicate = integer(0),
    method = character(0),
    omitted_unit = character(0),
    estimate = numeric(0),
    stringsAsFactors = FALSE
  )

  if (!(ci_method_used %in% c("bootstrap", "jackknife"))) {
    return(list(
      se_out = se_out,
      ci_low = ci_low,
      ci_high = ci_high,
      resampling_results = resampling_results,
      log_warning = log_warning
    ))
  }

  if (ci_method_used == "bootstrap") {
    if (!is.null(seed)) set.seed(seed)

    if (identical(scenario, "counts_no_se")) {
      original_unit_rows <- !duplicated(grouped_df$unit_id)
      units_for_sim <- grouped_df[original_unit_rows, , drop = FALSE]

      if (identical(variance_family, "binomial")) {
        non_integer_den <- abs(units_for_sim$denominator - floor(units_for_sim$denominator)) >
          sqrt(.Machine$double.eps)
        if (any(non_integer_den, na.rm = TRUE)) {
          stop("Binomial bootstrap requires integer-valued denominators.", call. = FALSE)
        }
      }

      boot_stat_count <- function(data, indices = NULL) {
        metric_calculator(data) * metric_multiplier
      }

      boot_ran_count <- function(data, mle) {
        bg <- data
        original_unit_rows <- !duplicated(bg$unit_id)
        units <- bg[original_unit_rows, , drop = FALSE]

        if (identical(mle$variance_family, "binomial")) {
          p <- pmin(pmax(units$numerator / units$denominator, 0), 1)
          sim_num_unit <- stats::rbinom(
            n = nrow(units),
            size = units$denominator,
            prob = p
          )
        } else {
          sim_num_unit <- stats::rpois(
            n = nrow(units),
            lambda = pmax(units$numerator, 0)
          )
        }

        sim_map <- data.frame(
          unit_id = units$unit_id,
          sim_numerator = sim_num_unit,
          stringsAsFactors = FALSE
        )

        idx <- match(bg$unit_id, sim_map$unit_id)
        bg$numerator <- sim_map$sim_numerator[idx]
        bg$allocated_numerator <- bg$numerator * bg$allocation_fraction
        bg
      }

      boot_out <- boot::boot(
        data = grouped_df,
        statistic = boot_stat_count,
        R = n_boot,
        sim = "parametric",
        ran.gen = boot_ran_count,
        mle = list(variance_family = variance_family)
      )

      boot_estimates <- as.numeric(boot_out$t[, 1])
      boot_method_label <- "parametric_bootstrap_counts_boot_package"
    } else {
      boot_stat_ecological <- function(data, indices) {
        bd <- data[indices, , drop = FALSE]
        bd$unit_id <- paste0(bd$unit_id, "__draw", seq_along(indices))

        bg <- tryCatch(grouper(bd), error = function(e) NULL)
        if (is.null(bg)) return(NA_real_)

        metric_calculator(bg) * metric_multiplier
      }

      boot_out <- boot::boot(
        data = base_units,
        statistic = boot_stat_ecological,
        R = n_boot,
        sim = "ordinary"
      )

      boot_estimates <- as.numeric(boot_out$t[, 1])
      boot_method_label <- "ecological_bootstrap_units_boot_package"
    }

    valid_estimates <- boot_estimates[is.finite(boot_estimates)]
    if (length(valid_estimates) < max(30, ceiling(0.50 * n_boot))) {
      stop(
        "Too few valid bootstrap replicates were produced. ",
        "Check grouping, extreme group estimates, and whether relative gaps have non-positive denominators.",
        call. = FALSE
      )
    }

    boot_ci <- tryCatch(
      boot::boot.ci(boot_out, conf = conf_level, type = "perc"),
      error = function(e) NULL
    )

    if (!is.null(boot_ci) && !is.null(boot_ci$percent) && ncol(boot_ci$percent) >= 5) {
      ci_low <- as.numeric(boot_ci$percent[1, 4])
      ci_high <- as.numeric(boot_ci$percent[1, 5])
    } else {
      log_warning <- c(
        log_warning,
        "boot::boot.ci() could not compute percentile confidence limits; percentile limits were computed from the bootstrap replicate distribution using stats::quantile()."
      )
      qs <- ci_percentile(valid_estimates, conf_level)
      ci_low <- qs[1]
      ci_high <- qs[2]
    }

    se_out <- stats::sd(valid_estimates, na.rm = TRUE)

    resampling_results <- data.frame(
      replicate = seq_len(n_boot),
      method = boot_method_label,
      omitted_unit = NA_character_,
      estimate = boot_estimates,
      stringsAsFactors = FALSE
    )
  }

  if (ci_method_used == "jackknife") {
    unit_ids <- base_units$unit_id

    jack_statistic <- function(keep) {
      keep <- as.integer(keep)
      jd <- base_units[keep, , drop = FALSE]
      if (nrow(jd) < 2) return(NA_real_)

      jg <- tryCatch(grouper(jd), error = function(e) NULL)
      if (is.null(jg)) return(NA_real_)

      metric_calculator(jg) * metric_multiplier
    }

    jack_out <- bootstrap::jackknife(
      x = seq_len(nrow(base_units)),
      theta = jack_statistic
    )

    jk_estimates <- as.numeric(jack_out$jack.values)
    se_out <- as.numeric(jack_out$jack.se)

    resampling_results <- data.frame(
      replicate = seq_along(unit_ids),
      method = "ecological_jackknife_units_bootstrap_package",
      omitted_unit = unit_ids,
      estimate = jk_estimates,
      stringsAsFactors = FALSE
    )
  }

  list(
    se_out = se_out,
    ci_low = ci_low,
    ci_high = ci_high,
    resampling_results = resampling_results,
    log_warning = log_warning
  )
}
