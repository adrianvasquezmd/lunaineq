#' Absolute inequality gap (two input modes)
#'
#' Computes the absolute inequality gap between the extreme social groups defined based on a
#' numeric equity stratifier. The function operates in two modes: \strong{(Mode A)} using a numerator and
#' denominator (\code{health_numerator} and \code{health_denominator}), which is the preferred
#' and most accurate method, as it relies directly on event counts;
#' and \strong{(Mode B)} using an aggregated health indicator (\code{health_indicator})
#' and population weights (\code{population_weights}), which approximates the gap based on weighted group means.
#' The function also provides automatic interpretations of the results in
#' four languages: Spanish, English, French, and Portuguese.
#'
#' In both cases, the social stratifier is grouped into \code{n_groups} using the selected method (e.g., quantiles),
#' ensuring that \emph{group 1} corresponds to the most disadvantaged group and the last group to the most advantaged
#' (as defined by \code{higher_ineq_is_favorable}). The absolute gap is calculated as the difference in the health indicator \eqn{y}
#' —defined as the mean value of the indicator within each group— between the most disadvantaged and the most advantaged groups: \eqn{\mathrm{AG} = y_{\text{disadvantaged}} - y_{\text{advantaged}}}.
#'
#' @param health_indicator_type Character string: either \code{"proportion"} for bounded indicators in \eqn{[0,1]}, or \code{"rate"} if already scaled (e.g., deaths per 100,000); required in all cases.
#' @param health_indicator Numeric vector with the health indicator values (e.g., \code{c(0.12, 0.34)} for proportions or \code{c(120, 345)} for rates); aggregated value per unit; used only in Mode B.
#' @param population_weights Numeric vector > 0 representing population size or weight per unit (e.g., \code{c(1000, 5000, 1200)}); used only in Mode B.
#' @param health_numerator Numeric vector with event counts (e.g., \code{c(5, 12, 8)} for deaths); used only in Mode A and preferred for robust estimation.
#' @param health_denominator Numeric vector > 0 associated with \code{health_numerator} (e.g., \code{c(500, 1000, 800)} representing population at risk); used only in Mode A.
#' @param equity_stratifier Numeric vector representing the social stratifier (e.g., poverty rate, education index); required in all cases.
#' @param higher_ineq_is_favorable Logical; \code{TRUE} if higher values of \code{equity_stratifier} indicate better social conditions (e.g., Gross Domestic Product per capita), or \code{FALSE} if they indicate worse conditions (e.g., Unsatisfied Basic Needs).
#' @param rate_scaling_factor Numeric scalar > 0 applied only when \code{health_indicator_type = "rate"}; commonly \code{1000} or \code{100000} (e.g., maternal mortality ratio per 100000 live births).
#' @param grouping_method String defining how to group the social stratifier: \code{"quantiles"}, \code{"equal intervals"}, \code{"manual intervals"}, \code{"k-means"}, \code{"hierarchical clustering"}, \code{"bagging clustering"}, \code{"fisher-jenks"}, \code{"jenks"}, or \code{"maximum differences"}.
#' @param n_groups Integer ≥ 2 specifying the number of social groups to create; ignored if \code{grouping_method = "manual intervals"}.
#' @param manual_breaks Numeric vector defining custom cut points (e.g., \code{c(20, 40)} defines intervals \code{0–20}, \code{20–40}, and \code{40+}); used only when \code{grouping_method = "manual intervals"}.
#' @param conf_level Numeric scalar between 0 and 1 (e.g., \code{0.95}) specifying the confidence level for inequality estimates; default is \code{0.95}.
#' @param language_interpretation Character string indicating the desired output language for interpretation and note; allowed values are \code{"en"}, \code{"es"}, \code{"fr"}, \code{"pt"} or full names (\code{"english"}, \code{"spanish"}, \code{"french"}, \code{"portuguese"}); default is \code{"en"}.
#'
#' @return A \code{list} with the following components:
#' \item{summary_table}{A \code{tibble} with the inequality metric name (\code{inequality_metric}), the estimated absolute gap (\code{value}), and its confidence interval bounds (\code{ci_lower}, \code{ci_upper}).}
#' \item{interpretation}{Character string with the automatic interpretation of the result in the selected language.}
#' \item{note}{Character string with the explanatory note in the selected language, based on the type of health indicator.}
#' \item{group_summary}{\code{tibble} summarizing the indicator for each social group. Contains either weighted means (Mode B) or raw counts and denominators (Mode A).}
#' \item{global_health_mean}{Overall mean of the health indicator: weighted average (Mode B) or global rate (Mode A).}
#' \item{data}{Final \code{tibble} used for the calculation, including processed variables and groupings.}
#'
#' @details
#' \subsection{A) Using numerator and denominator (preferred)}{
#' When \code{health_numerator} and \code{health_denominator} are provided (Mode A), the function uses actual counts to compute the inequality gap and its confidence interval.
#'
#' \itemize{
#'   \item \strong{Proportions}: The absolute gap is calculated as the difference between the proportion of events in the most disadvantaged and the most advantaged groups. Confidence intervals are computed using \code{PropCIs::diffscoreci()}, which implements score-based intervals for the difference of two independent proportions.
#'
#'   \item \strong{Rates}: The absolute gap is computed as the rate difference between the two extreme groups. Counts are assumed to follow a Poisson distribution, and confidence intervals are estimated using \code{ratesci::scoreci(distrib = "poi", contrast = "RD")}, which uses score-based methods for the rate difference. Final results are scaled by \code{rate_scaling_factor} (e.g., 1,000 or 100,000).
#'
#'   \item \emph{This mode does not rely on population weights}: the inference is based on actual event counts, making it more statistically robust.
#' }
#' }
#'
#' \subsection{B) With aggregated indicator weighted by population (approximate)}{
#' \itemize{
#'   \item Based on \code{health_indicator} and \code{population_weights}, weighted group means
#'         of the health indicator are calculated.
#'   \item \strong{Proportions}: Wald-type CI for difference in proportions, with
#'         \code{n} approximated using \code{population_weights}.
#'   \item \strong{Rates}: Wald-type CI for difference in “scaled” rates;
#'         \code{rate_scaling_factor} is used to express the result per 1,000/100,000.
#'   \item \emph{Warning}: these CIs may under- or overestimate the
#'         variance when sample sizes are small, heterogeneity is high,
#'         or prior aggregation has occurred; prefer mode A when feasible.
#' }}
#'
#' @note Mode B was developed as an alternative for situations where only aggregated
#' health indicators (e.g., proportions or rates) are available, and the disaggregated
#' data (numerator and denominator) is not accessible. Although Mode A is more precise,
#' Mode B still allows for valid inequality comparisons when only summary data is available.
#'
#' @examples
#' # ——— Mode A (preferred: numerator and denominator) ———
#'
#' # Health indicator: Rate. Suppose data1 with:
#' #   maternal_deaths  = number of maternal deaths.
#' #   live_births      = number of live births.
#' #   ubn_percent      = % of population with at least one unmet basic need (higher = worse social status).
#' data(data1)
#' ag_luna(
#'   health_indicator_type    = "rate",
#'   health_numerator         = data1$maternal_deaths,
#'   health_denominator       = data1$live_births,
#'   equity_stratifier        = data1$ubn_percent,
#'   higher_ineq_is_favorable = FALSE,
#'   rate_scaling_factor      = 100000,
#'   grouping_method          = "quantiles",
#'   n_groups                 = 4,
#'   conf_level               = 0.95,
#'   language_interpretation = "en"
#' )
#'
#' # Health indicator: Proportion. Suppose data2 with:
#' #   skilled_births   = number of births attended by skilled personnel.
#' #   total_births     = total number of births.
#' #   ubn_index        = % of population with at least one unmet basic need (higher = worse social status).
#' data(data2)
#' ag_luna(
#'   health_indicator_type    = "proportion",
#'   health_numerator         = data2$skilled_births,
#'   health_denominator       = data2$total_births,
#'   equity_stratifier        = data2$ubn_percent,
#'   higher_ineq_is_favorable = FALSE,
#'   grouping_method          = "quantiles",
#'   n_groups                 = 4,
#'   conf_level               = 0.95,
#'   language_interpretation = "en"
#' )
#'
#' # ——— Mode B (aggregated indicator + population as weight) ———
#'
#' # Suppose data1 with:
#' #   mmr       = maternal mortality ratio per 100,000 live births.
#' #   population = population weight of the unit of analysis.
#' #   ubn_index = % of population with at least one unmet basic need.
#' data(data1)
#' ag_luna(
#'   health_indicator_type    = "rate",
#'   health_indicator         = data1$mmr,
#'   equity_stratifier        = data1$ubn_percent,
#'   population_weights       = data1$population_zone,
#'   higher_ineq_is_favorable = FALSE,
#'   rate_scaling_factor      = 100000,
#'   grouping_method          = "quantiles",
#'   n_groups                 = 4,
#'   conf_level               = 0.95,
#'   language_interpretation = "en"
#' )
#'
#' # Suppose data2 with:
#' #   skilled_births_prop = proportion of births attended by skilled personnel (in [0,1]).
#' #   population           = population weight of the unit of analysis.
#' #   ubn_index            = % of population with at least one unmet basic need.
#' data(data2)
#' ag_luna(
#'   health_indicator_type    = "proportion",
#'   health_indicator         = data2$skilled_births_prop,
#'   equity_stratifier        = data2$ubn_percent,
#'   population_weights       = data2$population,
#'   higher_ineq_is_favorable = FALSE,
#'   grouping_method          = "quantiles",
#'   n_groups                 = 4,
#'   conf_level               = 0.95,
#'   language_interpretation = "en"
#' )
#'
#' @seealso \code{PropCIs::diffscoreci}, \code{ratesci::scoreci}, \code{classInt::classIntervals}
#'
#' @references
#' World Health Organization. Handbook on health inequality monitoring with a special focus on low- and middle-income countries. 2013.
#' Available from: https://www.who.int/publications/i/item/9789241548632
#'
#' @author
#' Adrián Vázquez-Mejía, MD MSc(c) (\email{apvasquez.md@gmail.com})
#' \cr Independent Consultant
#'
#' Oscar J. Mújica, MD MPH (\email{mujicaos@paho.org})
#' \cr Pan American Health Organization (PAHO), Washington D.C.
#'
#' Antonio Sanhueza, PhD (\email{sanhueza@paho.org})
#' \cr Pan American Health Organization (PAHO), Washington D.C.
#'
#' @export

ag_luna <- function(health_indicator_type = NULL,
               health_indicator = NULL,
               population_weights = NULL,
               health_numerator = NULL,
               health_denominator = NULL,
               equity_stratifier = NULL,
               higher_ineq_is_favorable = NULL,
               rate_scaling_factor = NULL,
               grouping_method = NULL,
               n_groups = NULL,
               manual_breaks = NULL,
               conf_level = NULL,
               language_interpretation = NULL) {

  # ───────────────────────────────────────────────
  # Argument validation (input consistency checks)
  # ───────────────────────────────────────────────

  # Check health indicator type
  if (is.null(health_indicator_type) || !(health_indicator_type %in% c("proportion", "rate"))) {
    stop("`health_indicator_type` is required and must be either 'proportion' or 'rate'.")
  }

  # Check presence health indicator & (health numerator | health denominator)
  if (is.null(health_indicator) && (is.null(health_numerator) || is.null(health_denominator))) {
    stop("You must provide either `health_indicator`, or both `health_numerator` and `health_denominator`.")
  }

  # Check presence health numerator & health denominator
  use_num_den <- !is.null(health_numerator) && !is.null(health_denominator)

  # Mode A: numerator/denominator provided
  if (use_num_den) {

    # Check health_numerator
    if (!is.null(health_numerator)) {
      if (!is.numeric(health_numerator)) {
        stop("`health_numerator` must be numeric.")
      }
      if (any(health_numerator < 0, na.rm = TRUE)) {
        stop("`health_numerator` cannot contain negative values. All values must be ≥ 0.")
      }
      if (any(health_numerator %% 1 != 0, na.rm = TRUE)) {
        warning("`health_numerator` contains non-integer values. These will be rounded to the nearest whole number.")
        health_numerator <- round(health_numerator)
      }
    }

    # Check health denominator
    if (!is.null(health_denominator)) {
      if (!is.numeric(health_denominator)) {
        stop("`health_denominator` must be numeric.")
      }
      if (any(health_denominator <= 0, na.rm = TRUE)) {
        stop("`health_denominator` must contain values > 0 only.")
      }
      if (any(health_denominator %% 1 != 0, na.rm = TRUE)) {
        warning("`health_denominator` contains non-integer values. These will be rounded to the nearest whole number.")
        health_denominator <- round(health_denominator)
      }
    }

    #Check health indicator calculated from numerator/denominator
    if (health_indicator_type == "proportion") {
      ratio_calc <- health_numerator / health_denominator
      if (any(ratio_calc < 0 | ratio_calc > 1, na.rm = TRUE)) {
        stop("Health indicator calculated from numerator/denominator fall outside the range [0,1]. Check inputs.")
      }
    }
    if (health_indicator_type == "rate") {
      ratio_calc <- health_numerator / health_denominator * rate_scaling_factor
      if (any(ratio_calc < 0, na.rm = TRUE)) {
        stop("Health indicator values calculated from numerator/denominator fall outside the range [0, Inf]. Check inputs.")
      }
    }
    if (!is.null(health_indicator)) {
      warning("`health_numerator` and `health_denominator` detected. health_indicator will be ignored because health_numerator and health_denominator were provided. The indicator will be derived from those inputs instead.")
    }

    #Check length of inputs
    if (length(health_numerator) != length(health_denominator)) {
      stop("`health_numerator` must have the same length as `health_denominator`.")
    }

    if (length(health_numerator) != length(equity_stratifier)) {
      stop("`health_numerator` must have the same length as `equity_stratifier`.")
    }
  }

  # Mode B: numerator/denominator not provided
  if (!use_num_den) {

    if (!is.null(health_numerator)) {
      warning(
        "`health_numerator` was provided without `health_denominator` and will be ignored."
      )
      health_numerator <- NULL
    }

    if (!is.null(health_denominator)) {
      warning(
        "`health_denominator` was provided without `health_numerator` and will be ignored."
      )
      health_denominator <- NULL
    }

    # Check health indicator
    if (is.null(health_indicator)) stop("Missing `health_indicator`. Provide it when raw numerator/denominator are not available.")
    if (!is.numeric(health_indicator)) stop("`health_indicator` must be numeric.")
    if (any(is.na(health_indicator))) stop("Missing values detected in `health_indicator`.")
    if (health_indicator_type == "proportion") {
      if (any(health_indicator < 0 | health_indicator > 1, na.rm = TRUE)) {
        stop("`health_indicator` values must be between 0 and 1 for type 'proportion'.")
      }
    } else if (health_indicator_type == "rate") {
      if (is.null(rate_scaling_factor)) {
        stop("For rates, a numeric `rate_scaling_factor` > 0 (e.g., 1000 or 100000) must be specified.")
      }
      if (any(health_indicator < 0, na.rm = TRUE)) {
        stop("`health_indicator` values must be ≥ 0 for type 'rate'.")
      }
    }

    # Check population weights
    if (is.null(population_weights)) {
      stop("`population_weights` must be provided when `health_numerator` and `health_denominator` are not.")
    }
    if (!is.numeric(population_weights)) stop("`population_weights` must be numeric.")
    if (any(population_weights <= 0, na.rm = TRUE)) stop("`population_weights` must be strictly greater than 0.")
    if (any(is.na(population_weights))) stop("Missing values detected in `population_weights`.")
    if (!use_num_den && any(population_weights %% 1 != 0)) {
      warning("population_weights contain non-integer values. CIs computed with simple Wald approximation may be unreliable; prefer numerator/denominator when possible.")
    }

    # Check length of inputs (health indicator, population weights and equity stratifier)
    if (length(population_weights) != length(health_indicator)) {
      stop("`population_weights` must have the same length as `health_indicator`.")
    }
    if (length(population_weights) != length(equity_stratifier)) {
      stop("`population_weights` must have the same length as `equity_stratifier`.")
    }
  }

  # Check equity stratifier
  if (is.null(equity_stratifier)) {
    stop("You must provide `equity_stratifier`.")
  }
  if (!is.numeric(equity_stratifier)) stop("`equity_stratifier` must be numeric.")
  if (any(is.na(equity_stratifier))) stop("Missing values detected in `equity_stratifier`.")

  # Check higher equity stratifier is favorable
  if (is.null(higher_ineq_is_favorable)) {
    stop("Please specify `higher_ineq_is_favorable` (TRUE if higher values = more favorable).")
  }
  if (!is.logical(higher_ineq_is_favorable)) {
    stop("`higher_ineq_is_favorable` must be logical (TRUE/FALSE).")
  }

  # Check health indicator type
  if (health_indicator_type == "rate") {
    if (is.null(rate_scaling_factor) || !is.numeric(rate_scaling_factor) ||
        length(rate_scaling_factor) != 1 || rate_scaling_factor <= 0) {
      stop("For `health_indicator_type = 'rate'`, you must provide a numeric `rate_scaling_factor` > 0 (e.g., 1000 or 100000).")
    }
  }

  # Warn if scaling factor is ignored for proportions
  if (health_indicator_type == "proportion") {
    if (!is.null(rate_scaling_factor)) {
      warning("`rate_scaling_factor` is ignored because `health_indicator_type` is set to `proportion`. It only applies to rates.")
    }
    rate_scaling_factor <- 1
  }

  # Check grouping method
  if (is.null(grouping_method)) {
    warning("No `grouping_method` specified. Defaulting to 'quantiles'.")
    grouping_method <- "quantiles"
  }

  allowed_methods <- c(
    "manual intervals", "equal intervals", "quantiles", "k-means",
    "hierarchical clustering", "bagging clustering", "fisher-jenks", "jenks",
    "maximum differences"
  )

  if (!(tolower(grouping_method) %in% tolower(allowed_methods))) {
    stop("Unrecognized `grouping_method`. Allowed values are: ",
         paste(allowed_methods, collapse = ", "), ".")
  }

  if (tolower(grouping_method) == "manual intervals" && is.null(manual_breaks)) {
    stop("When using 'manual intervals', you must provide `manual_breaks`.")
  }

  # Check number of groups
  if (is.null(n_groups)) {
    warning("No `n_groups` specified. Defaulting to 4.")
    n_groups <- 4L
  } else {
    if (!is.numeric(n_groups) || length(n_groups) != 1) {
      stop("`n_groups` must be a single numeric integer greater than or equal to 2.")
    }
    if (n_groups < 2 || n_groups != as.integer(n_groups)) {
      stop("`n_groups` must be an integer greater than or equal to 2.")
    }
    n_groups <- as.integer(n_groups)
  }

  # Check manual breakpoints
  if (!is.null(manual_breaks)) {
    if (is.character(manual_breaks) && length(manual_breaks) == 1) {
      manual_breaks <- as.numeric(unlist(strsplit(gsub("\\s+", " ", manual_breaks), "[, ]+")))
    } else {
      manual_breaks <- as.numeric(manual_breaks)
    }

    if (any(is.na(manual_breaks))) stop("`manual_breaks` contains non-numeric values.")
    if (length(manual_breaks) == 0) stop("`manual_breaks` cannot be empty.")
    if (any(duplicated(manual_breaks))) stop("`manual_breaks` must not contain duplicate values.")

    if (!is.null(equity_stratifier)) {
      if (any(manual_breaks <= min(equity_stratifier, na.rm = TRUE) |
              manual_breaks >= max(equity_stratifier, na.rm = TRUE))) {
        stop("Some `manual_breaks` are outside the range of `equity_stratifier`. Check your breakpoints.")
      }
    }
  }

  # Check confidence level
  if (is.null(conf_level)) {
    warning("No `conf_level` specified. Defaulting to 0.95.")
    conf_level <- 0.95
  } else {
    conf_level <- as.numeric(conf_level)
    if (is.na(conf_level) || length(conf_level) != 1 || !(conf_level > 0 && conf_level < 1)) {
      stop("`conf_level` must be a numeric value between 0 and 1 (e.g., 0.95).")
    }
  }

  # Check language interpretation
  alias_map <- c(en = "english", es = "spanish", pt = "portuguese", fr = "french")
  full_to_alias <- setNames(names(alias_map), alias_map)

  if (missing(language_interpretation) || is.null(language_interpretation) ||
      !nzchar(trimws(as.character(language_interpretation)))) {
    warning("`language_interpretation` not provided: defaulting to 'en'.")
    language_interpretation <- "en"
  } else {
    language_interpretation <- tolower(trimws(as.character(language_interpretation)))
  }

  if (language_interpretation %in% names(full_to_alias)) {
    language_interpretation <- unname(full_to_alias[language_interpretation])
  } else if (language_interpretation %in% names(alias_map)) {
  } else {
    stop("Unrecognized `language_interpretation`. Allowed aliases: ",
         paste(names(alias_map), collapse = ", "),
         ". Allowed full names: ",
         paste(unname(alias_map), collapse = ", "), ".")
  }

  # ────────────────────────────────
  # Base data frame preparation
  # ────────────────────────────────
  # If both numerator and denominator are provided, construct health_indicator from them;
  # otherwise, use the health_indicator provided directly.
  health_indicator_final <- if (!is.null(health_numerator) && !is.null(health_denominator)) {
    if (health_indicator_type == "proportion") {
      as.numeric(health_numerator / health_denominator)
    } else { # "rate"
      as.numeric((health_numerator / health_denominator) * rate_scaling_factor)
    }
  } else {
    as.numeric(health_indicator)
  }

  if (!use_num_den) {
    df <- tibble::tibble(
      health_indicator = health_indicator_final,
      equity_stratifier = as.numeric(equity_stratifier),
      population_weights = as.numeric(population_weights)
    )
  } else {
    df <- tibble::tibble(
      health_indicator   = health_indicator_final,
      health_numerator   = health_numerator,
      health_denominator = health_denominator,
      equity_stratifier  = as.numeric(equity_stratifier)
    )
  }

  # Create social position variable for ranking and grouping purposes.
  # equity_stratifier_order is used to rank and group units. If higher values of the stratifier
  # are considered more favorable (higher_ineq_is_favorable = TRUE), the order is preserved.
  # Otherwise, the scale is reversed to ensure that group 1 always corresponds to the most disadvantaged.
  df <- df %>%
    mutate(equity_stratifier_order = if (isTRUE(higher_ineq_is_favorable)) equity_stratifier else -equity_stratifier) %>%
    arrange(equity_stratifier_order)

  # ──────────────────────────────────────────
  # Safe creation of breaks for grouping
  # ──────────────────────────────────────────
  safe_class_intervals <- function(x, n, style, fixedBreaks = NULL) {
    tryCatch({
      if (!is.null(fixedBreaks)) {
        ci <- classInt::classIntervals(x, n = n, style = style, fixedBreaks = fixedBreaks)
      } else {
        ci <- classInt::classIntervals(x, n = n, style = style)
      }
      ci$brks
    }, error = function(e) {
      stop(sprintf("Unable to create breaks using style = '%s'. Check the method or the data. Original error: %s", style, e$message))
    })
  }

  # Map human-readable method names to classInt styles (if needed)
  style <- switch(tolower(grouping_method),
                  "manual intervals"     = "fixed",
                  "equal intervals"      = "equal",
                  "quantiles"            = "quantile",
                  "k-means"              = "kmeans",
                  "hierarchical clustering" = "hclust",
                  "bagging clustering"    = "bclust",
                  "fisher-jenks"         = "fisher",
                  "jenks"                = "jenks",
                  "maximum differences"       = "maximum",
                  stop("Unrecognized grouping_method.")
  )

  # Compute breaks for social stratification
  if (style == "fixed") {
    fixedBreaks <- c(-Inf, sort(manual_breaks), Inf)
    breaks_vec <- safe_class_intervals(df$equity_stratifier_order, n = length(fixedBreaks) - 1, style = "fixed", fixedBreaks = fixedBreaks)
  } else {
    breaks_vec <- safe_class_intervals(df$equity_stratifier_order, n = n_groups, style = style)
  }

  # Assign social groups based on equity_stratifier_order and computed breaks
  df <- df %>%
    mutate(group = cut(equity_stratifier_order, breaks = breaks_vec, include.lowest = TRUE)) %>%
    mutate(group = as.integer(group))

  # Check that grouping resulted in at least two non-empty groups
  if (n_distinct(df$group) < 2) stop("Grouping resulted in fewer than two valid groups. Consider changing `grouping_method` or `n_groups`.")

  # ──────────────────────────────────────────────────────────────────────────────
  # Global weighted means for the health indicator and the equity stratifier
  # ──────────────────────────────────────────────────────────────────────────────
  if (!use_num_den) {
    global_health_mean <- weighted.mean(df$health_indicator, w = df$population_weights, na.rm = TRUE)
  } else {
    total_num <- sum(df$health_numerator, na.rm = TRUE)
    total_den <- sum(df$health_denominator, na.rm = TRUE)

    global_health_mean <- if (health_indicator_type == "proportion") {
      total_num / total_den
    } else if (health_indicator_type == "rate") {
      (total_num / total_den) * rate_scaling_factor
    }
  }

  # ──────────────────────────────────────────────────────────────────────────────
  # Group-level summary
  # - pop_group: total population within each group
  # - health_indicator_group: weighted intra-group mean of the health indicator
  # ──────────────────────────────────────────────────────────────────────────────

  if (!use_num_den) {
    group_summary <- df %>%
      group_by(group) %>%
      summarise(
        health_indicator_group = weighted.mean(health_indicator, w = population_weights, na.rm = TRUE),
        pop_group              = sum(population_weights, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(
        population_weights = pop_group,
        health_indicator   = health_indicator_group
      ) %>%
      arrange(group)
  } else {
    group_summary <- df %>%
      group_by(group) %>%
      summarise(
        health_numerator   = sum(health_numerator, na.rm = TRUE),   # total numerator in group
        health_denominator = sum(health_denominator, na.rm = TRUE), # total denominator in group
        .groups = "drop"
      ) %>%
      mutate(
        health_indicator = (health_numerator / health_denominator) * rate_scaling_factor
      ) %>%
      arrange(group)
  }

  # ──────────────────────────────
  # Absolute Gap Calculation
  # ──────────────────────────────

  min_g <- min(group_summary$group)
  max_g <- max(group_summary$group)

  if (!use_num_den) {

    z <- qnorm(1 - (1 - conf_level) / 2)

    if (health_indicator_type == "proportion") {
      # Approximation for independent proportions
      p1 <- group_summary$health_indicator[group_summary$group == min_g]
      n1 <- group_summary$population_weights[group_summary$group == min_g]
      p2 <- group_summary$health_indicator[group_summary$group == max_g]
      n2 <- group_summary$population_weights[group_summary$group == max_g]

      # Standard error of the difference between two independent proportions
      se <- sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2))
      # Confidence interval for the absolute gap
      value   <- p1 - p2
      ci_low  <- value - z * se
      ci_high <- value + z * se

    } else if (health_indicator_type == "rate") {
      # Approximation for rates
      x1 <- group_summary$health_indicator[group_summary$group == min_g] *
        group_summary$population_weights[group_summary$group == min_g] / rate_scaling_factor
      n1 <- group_summary$population_weights[group_summary$group == min_g]
      x2 <- group_summary$health_indicator[group_summary$group == max_g] *
        group_summary$population_weights[group_summary$group == max_g] / rate_scaling_factor
      n2 <- group_summary$population_weights[group_summary$group == max_g]

      # Combined rate across both groups
      p_comb <- (x1 + x2) / (n1 + n2)

      # Standard error of the scaled rate difference (Wald method)
      se <- sqrt(p_comb * (1 - p_comb) * (1 / n1 + 1 / n2))
      value   <- (x1 / n1) - (x2 / n2)
      ci_low  <- value - z * se
      ci_high <- value + z * se

      value   <- value   * rate_scaling_factor
      ci_low  <- ci_low  * rate_scaling_factor
      ci_high <- ci_high * rate_scaling_factor
    }

  } else {

    if (health_indicator_type == "proportion") {
      x1 <- sum(df$health_numerator[df$group == min_g], na.rm = TRUE)
      n1 <- sum(df$health_denominator[df$group == min_g], na.rm = TRUE)
      x2 <- sum(df$health_numerator[df$group == max_g], na.rm = TRUE)
      n2 <- sum(df$health_denominator[df$group == max_g], na.rm = TRUE)

      ci     <- PropCIs::diffscoreci(x1, n1, x2, n2, conf.level = conf_level)
      value  <- (x1 / n1) - (x2 / n2)
      ci_low <- as.numeric(ci$conf.int[1])
      ci_high <- as.numeric(ci$conf.int[2])

    } else if (health_indicator_type == "rate") {
      x1 <- sum(df$health_numerator[df$group == min_g], na.rm = TRUE)
      n1 <- sum(df$health_denominator[df$group == min_g], na.rm = TRUE)
      x2 <- sum(df$health_numerator[df$group == max_g], na.rm = TRUE)
      n2 <- sum(df$health_denominator[df$group == max_g], na.rm = TRUE)

      res <- ratesci::scoreci(
        x1 = as.numeric(x1),
        n1 = as.numeric(n1),
        x2 = as.numeric(x2),
        n2 = as.numeric(n2),
        distrib  = "poi",      # Poisson distribution for count data
        contrast = "RD",       # Rate difference
        level    = conf_level,
        skew     = TRUE
      )

      value   <- as.numeric(res$estimates[1, "est"])   * rate_scaling_factor
      ci_low  <- as.numeric(res$estimates[1, "lower"]) * rate_scaling_factor
      ci_high <- as.numeric(res$estimates[1, "upper"]) * rate_scaling_factor
    }
  }

  # Confidence interval check
  if (!is.finite(ci_low) || !is.finite(ci_high)) stop("Invalid confidence intervals (NA/Inf).")

  summary_table <- tibble::tibble(
    inequality_metric = "Absolute gap",
    value      = value,
    ci_lower   = ci_low,
    ci_upper   = ci_high
  )

  # ──────────────────────────────
  # Automatic interpretation
  # ──────────────────────────────
  factor_legible <- format(rate_scaling_factor, scientific = FALSE, big.mark = ",")

  # Interpretation in English
  interpretation_english <- if (health_indicator_type == "rate") {
    if (value > 0) {
      paste0("The absolute gap is ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%CI [", round(ci_low, 2), "; ", round(ci_high, 2), "]. This indicates that the most disadvantaged group has ",
             round(value, 2), " more events per ", factor_legible, " individuals in the base population than the most advantaged group.")
    } else if (value < 0) {
      paste0("The absolute gap is ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%CI [", round(ci_low, 2), "; ", round(ci_high, 2), "]. This indicates that the most advantaged group has ",
             abs(round(value, 2)), " more events per ", factor_legible, " individuals in the base population than the most disadvantaged group.")
    } else {
      paste0("The absolute gap is 0. There are no absolute differences per ", factor_legible, " individuals between the extreme social groups.")
    }
  } else if (health_indicator_type == "proportion") {
    if (value > 0) {
      paste0("The absolute gap is ", round(value * 100, 2), " percentage points (pp), ", round(conf_level * 100, 1),
             "%CI [", round(ci_low * 100, 2), "pp; ", round(ci_high * 100, 2), "pp]. This indicates that the health indicator in the most disadvantaged group is higher by ", round(value * 100, 2), " percentage points (pp) compared to the most advantaged group.")
    } else if (value < 0) {
      paste0("The absolute gap is ", round(value * 100, 2), " percentage points (pp), ", round(conf_level * 100, 1),
             "%CI [", round(ci_low * 100, 2), "pp; ", round(ci_high * 100, 2), "pp]. This indicates that the health indicator in the most advantaged group is higher by ", abs(round(value * 100, 2)), " percentage points (pp) compared to the most disadvantaged group.")
    } else {
      paste0("The absolute gap is 0 pp. There are no absolute differences between the extreme groups.")
    }
  }

  interpretation_english <- paste0(
    interpretation_english,
    if (ci_low <= 0 && ci_high >= 0) {
      " However, this difference is not statistically significant."
    } else {
      " This difference is statistically significant."
    }
  )

  # Interpretation in Spanish
  interpretation_spanish <- if (health_indicator_type == "rate") {
    if (value > 0) {
      paste0("La brecha absoluta es ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%IC [", round(ci_low, 2), "; ", round(ci_high, 2), "]. Esto indica que el grupo más desfavorecido tiene ",
             round(value, 2), " más eventos por cada ", factor_legible, " individuos de la población base que el grupo social más favorecido.")
    } else if (value < 0) {
      paste0("La brecha absoluta es ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%IC [", round(ci_low, 2), "; ", round(ci_high, 2), "]. Esto indica que el grupo más favorecido tiene ",
             abs(round(value, 2)), " más eventos por cada ", factor_legible, " individuos de la población base que el grupo social más desfavorecido.")
    } else {
      paste0("La brecha absoluta es 0. No hay diferencias absolutas por cada ", factor_legible, " individuos de la población entre los grupos sociales extremos.")
    }
  } else if (health_indicator_type == "proportion") {
    if (value > 0) {
      paste0("La brecha absoluta es ", round(value * 100, 2), " puntos porcentuales (pp), ", round(conf_level * 100, 1),
             "%IC [", round(ci_low * 100, 2), "pp; ", round(ci_high * 100, 2), "pp]. Esto indica que el indicador de salud en el grupo más desfavorecido es mayor por ", round(value * 100, 2), " puntos porcentuales (pp) en comparación con el grupo más favorecido.")
    } else if (value < 0) {
      paste0("La brecha absoluta es ", round(value * 100, 2), " puntos porcentuales (pp), ", round(conf_level * 100, 1),
             "%IC [", round(ci_low * 100, 2), "pp; ", round(ci_high * 100, 2), "pp]. Esto indica que el indicador de salud en el grupo más favorecido es mayor por ", abs(round(value * 100, 2)), " puntos porcentuales (pp) en comparación con el grupo más desfavorecido.")
    } else {
      paste0("La brecha absoluta es 0 pp. No hay diferencias absolutas entre los grupos extremos.")
    }
  }

  interpretation_spanish <- paste0(
    interpretation_spanish,
    if (ci_low <= 0 && ci_high >= 0) {
      " Sin embargo, esta diferencia no es estadísticamente significativa."
    } else {
      " Esta diferencia es estadísticamente significativa."
    }
  )

  # Interprétation en French
  interpretation_french <- if (health_indicator_type == "rate") {
    if (value > 0) {
      paste0("L'écart absolu est de ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%IC [", round(ci_low, 2), " ; ", round(ci_high, 2), "]. Cela indique que le groupe le plus défavorisé présente ",
             round(value, 2), " événements de plus pour chaque ", factor_legible, " individus de la population de base que le groupe social le plus favorisé.")
    } else if (value < 0) {
      paste0("L'écart absolu est de ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%IC [", round(ci_low, 2), " ; ", round(ci_high, 2), "]. Cela indique que le groupe le plus favorisé présente ",
             abs(round(value, 2)), " événements de plus pour chaque ", factor_legible, " individus de la population de base que le groupe social le plus défavorisé.")
    } else {
      paste0("L'écart absolu est de 0. Il n'y a pas de différences absolues pour chaque ", factor_legible, " individus de la population entre les groupes sociaux extrêmes.")
    }
  } else if (health_indicator_type == "proportion") {
    if (value > 0) {
      paste0("L'écart absolu est de ", round(value * 100, 2), " points de pourcentage (pp), ", round(conf_level * 100, 1),
             "%IC [", round(ci_low * 100, 2), "pp ; ", round(ci_high * 100, 2), "pp]. Cela indique que l'indicateur de santé dans le groupe le plus défavorisé est supérieur de ",
             round(value * 100, 2), " points de pourcentage (pp) par rapport au groupe le plus favorisé.")
    } else if (value < 0) {
      paste0("L'écart absolu est de ", round(value * 100, 2), " points de pourcentage (pp), ", round(conf_level * 100, 1),
             "%IC [", round(ci_low * 100, 2), "pp ; ", round(ci_high * 100, 2), "pp]. Cela indique que l'indicateur de santé dans le groupe le plus favorisé est supérieur de ",
             abs(round(value * 100, 2)), " points de pourcentage (pp) par rapport au groupe le plus défavorisé.")
    } else {
      paste0("L'écart absolu est de 0 pp. Il n'y a pas de différences absolues entre les groupes extrêmes.")
    }
  }

  interpretation_french <- paste0(
    interpretation_french,
    if (ci_low <= 0 && ci_high >= 0) {
      " Cependant, cette différence n’est pas statistiquement significative."
    } else {
      " Cette différence est statistiquement significative."
    }
  )

  # Interpretation in Portuguese
  interpretation_portuguese <- if (health_indicator_type == "rate") {
    if (value > 0) {
      paste0("A diferença absoluta é de ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%IC [", round(ci_low, 2), "; ", round(ci_high, 2), "]. Isso indica que o grupo mais desfavorecido apresenta ",
             round(value, 2), " eventos a mais por cada ", factor_legible, " indivíduos da população de base do que o grupo social mais favorecido.")
    } else if (value < 0) {
      paste0("A diferença absoluta é de ", round(value, 2), ", ", round(conf_level * 100, 1),
             "%IC [", round(ci_low, 2), "; ", round(ci_high, 2), "]. Isso indica que o grupo mais favorecido apresenta ",
             abs(round(value, 2)), " eventos a mais por cada ", factor_legible, " indivíduos da população de base do que o grupo social mais desfavorecido.")
    } else {
      paste0("A diferença absoluta é 0. Não há diferenças absolutas por cada ", factor_legible, " indivíduos da população entre os grupos sociais extremos.")
    }
  } else if (health_indicator_type == "proportion") {
    if (value > 0) {
      paste0("A diferença absoluta é de ", round(value * 100, 2), " pontos percentuais (pp), ", round(conf_level * 100, 1),
             "%IC [", round(ci_low * 100, 2), "pp; ", round(ci_high * 100, 2), "pp]. Isso indica que o indicador de saúde no grupo mais desfavorecido é maior em ", round(value * 100, 2), " pontos percentuais (pp) em comparação com o grupo mais favorecido.")
    } else if (value < 0) {
      paste0("A diferença absoluta é de ", round(value * 100, 2), " pontos percentuais (pp), ", round(conf_level * 100, 1),
             "%IC [", round(ci_low * 100, 2), "pp; ", round(ci_high * 100, 2), "pp]. Isso indica que o indicador de saúde no grupo mais favorecido é maior em ", abs(round(value * 100, 2)), " pontos percentuais (pp) em comparação com o grupo mais desfavorecido.")
    } else {
      paste0("A diferença absoluta é 0 pp. Não há diferenças absolutas entre os grupos extremos.")
    }
  }

  interpretation_portuguese <- paste0(
    interpretation_portuguese,
    if (ci_low <= 0 && ci_high >= 0) {
      " No entanto, essa diferença não é estatisticamente significativa."
    } else {
      " Essa diferença é estatisticamente significativa."
    }
  )

  # Note based on indicator type
  note_spanish <- if (health_indicator_type == "rate") {
    "Nota: Los eventos corresponden al numerador y la población base al denominador del indicador de salud. Ejemplo: RMM."
  } else if (health_indicator_type == "proportion") {
    "Nota: El valor representa una proporción (por ejemplo, partos atendidos por personal calificado)."
  }

  note_english <- if (health_indicator_type == "rate") {
    "Note: The events correspond to the numerator and the base population to the denominator of the health indicator. Example: MMR."
  } else if (health_indicator_type == "proportion") {
    "Note: The value represents a proportion (e.g., deliveries attended by skilled health personnel)."
  }

  note_portuguese <- if (health_indicator_type == "rate") {
    "Nota: Os eventos correspondem ao numerador e a população de base ao denominador do indicador de saúde. Exemplo: RMM."
  } else if (health_indicator_type == "proportion") {
    "Nota: O valor representa uma proporção (por exemplo, partos assistidos por pessoal qualificado)."
  }

  note_french <- if (health_indicator_type == "rate") {
    "Remarque : Les événements correspondent au numérateur et la population de base au dénominateur de l’indicateur de santé. Exemple : TMM."
  } else if (health_indicator_type == "proportion") {
    "Remarque : La valeur représente une proportion (par exemple, accouchements assistés par du personnel qualifié)."
  }

  if (language_interpretation == "en") {
    interpretation <- interpretation_english
    note <- note_english
  } else if  (language_interpretation == "es") {
    interpretation <- interpretation_spanish
    note <- note_spanish
  } else if  (language_interpretation == "pt") {
    interpretation <- interpretation_portuguese
    note <- note_portuguese
  } else if  (language_interpretation == "fr") {
    interpretation <- interpretation_french
    note <- note_french
  }

  # ────────────
  # Final output
  # ────────────
  return(list(
    summary_table = summary_table,
    interpretation = interpretation,
    note = note,
    group_summary = group_summary,
    global_health_mean = global_health_mean,
    data = df
  ))
}

