#' Relative Concentration Index (two input modes and optional corrections)
#'
#' Computes the Relative Concentration Index to quantify socioeconomic-related inequalities in a health indicator.
#' The function operates in two modes: \strong{(Mode A)} using a numerator and denominator (\code{health_numerator} and \code{health_denominator}),
#' which is the preferred and most accurate method as it relies directly on event counts; and \strong{(Mode B)} using an aggregated health indicator
#' (\code{health_indicator}) and population weights (\code{population_weights}), which is used when disaggregated data are not available.
#' The RCI is estimated using either the Fuller method (area under a fitted symmetric Lorenz curve) or the Wagstaff & Kakwani method
#' (regression approach). Optionally, the function computes correction variants (Wagstaff adjustment and Erreygers correction) when the health indicator is bounded (e.g., a proportion)
#' . The function provides automatic interpretations in four languages: Spanish, English, French, and Portuguese.
#'
#' @param health_indicator_type Character string: either \code{"proportion"} for bounded indicators in \eqn{[0,1]}, or \code{"rate"} if already scaled (e.g., deaths per 100,000); required in all cases.
#' @param health_indicator Numeric vector with the health indicator values (e.g., \code{c(0.12, 0.34)} for proportions or \code{c(120, 345)} for rates); aggregated value per unit; used only in Mode B.
#' @param population_weights Numeric vector > 0 representing population size or weight per unit (e.g., \code{c(1000, 5000, 1200)}); used only in Mode B.
#' @param health_numerator Numeric vector with event counts (e.g., \code{c(5, 12, 8)} for deaths); used only in Mode A and preferred for robust estimation.
#' @param health_denominator Numeric vector > 0 associated with \code{health_numerator} (e.g., \code{c(500, 1000, 800)} representing population at risk); used only in Mode A.
#' @param equity_stratifier Numeric vector representing the social stratifier (e.g., poverty rate, education index); required in all cases.
#' @param higher_ineq_is_favorable Logical; \code{TRUE} if higher values of \code{equity_stratifier} indicate better social conditions (e.g., Gross Domestic Product per capita), or \code{FALSE} if they indicate worse conditions (e.g., Unsatisfied Basic Needs).
#' @param rate_scaling_factor Numeric scalar > 0 applied only when \code{health_indicator_type = "rate"}; commonly \code{1000} or \code{100000} (e.g., maternal mortality ratio per 100000 live births).
#' @param calculation_method Character string specifying the method used to compute the RCI. Allowed values are \code{"Fuller"} (smoothed Lorenz curve approach) and \code{"WK"} (Wagstaff & Kakwani regression/smoothing approach). Default is \code{"WK"}.
#' @param apply_corrections Logical; \code{TRUE} to compute additional corrected indices (Wagstaff adjustment and Erreygers correction) in addition to the standard RCI; default \code{FALSE}.
#' @param smoothing_factor Numeric scalar (e.g., \code{0.5}) controlling the smoothness of the concentration curve. Used only when \code{calculation_method = "WK"}; ignored otherwise.
#' @param integration_method Character string specifying the optimization algorithm for parameter estimation when \code{calculation_method = "Fuller"}. Allowed values include \code{"BFGS"}, \code{"CG"}, \code{"Nelder-Mead"}, \code{"L-BFGS-B"}, \code{"nlm"}, \code{"nlminb"}, \code{"spg"}, \code{"ucminf"}, \code{"Rcgmin"}, \code{"Rvmmin"}, \code{"newuoa"}, \code{"bobyqa"}, \code{"nmkb"}, and \code{"hjkb"}. Default is \code{"BFGS"}.
#' @param conf_level Numeric scalar between 0 and 1 (e.g., \code{0.95}) specifying the confidence level for the RCI estimate. Default is \code{0.95}.
#' @param language_interpretation Character string indicating the desired output language for interpretation and note; allowed values are \code{"en"}, \code{"es"}, \code{"fr"}, \code{"pt"} or full names (\code{"english"}, \code{"spanish"}, \code{"french"}, \code{"portuguese"}); default is \code{"en"}.
#'
#' @return A \code{list} with the following components:
#' \item{summary}{A \code{tibble} with one or more rows containing: the RCI estimate(s) (\code{estimate}), the confidence interval bounds (\code{ci_lower}, \code{ci_upper}), a short \code{metric} label and the \code{method} used. When \code{apply_corrections = TRUE} additional rows for the Wagstaff and Erreygers adjustments are included.}
#' \item{interpretation}{Character string (or character vector when corrections are applied) with the automatic interpretation of the result(s) in the selected language.}
#' \item{note}{Character string with the explanatory note in the selected language, based on the type of health indicator.}
#' \item{global_health_mean}{Overall mean of the health indicator: weighted average (Mode B) or global rate (Mode A).}
#' \item{data}{Final \code{tibble} used for the calculation, including processed variables and the estimated concentration/Lorenz curve.}
#'
#' @details
#'
#' \subsection{Estimation methods}{
#' The Relative Concentration Index (RCI) can be estimated using two methods:
#'
#' \describe{
#'   \item{\strong{Wagstaff & Kakwani (WK) method}}{
#'   The method estimates the index from the observed distribution by smoothing the empirical concentration curve and fitting a linear regression to obtain the RCI. The discrete expression used to summarise concentration is
#'
#'   \deqn{\mathrm{RCI} \;=\; \frac{\sum_j p_j \cdot (2X_j - 1) \cdot y_j}{\mu}}
#'
#'   where \eqn{p_j} is the population share of subgroup \eqn{j}; \eqn{y_j} is the value of the health indicator in subgroup \eqn{j}; \eqn{X_j} is the relative rank of subgroup \eqn{j}, calculated as \eqn{X_j = \sum_{i=1}^{j} p_i - 0.5\,p_j}; and \eqn{\mu} denotes the global mean of the health indicator. The fitted regression provides the RCI estimate and its confidence interval.
#'   }
#'
#'   \item{\strong{Fuller method}}{
#'   The method fits a symmetric Lorenz-type curve to the cumulative distribution and derives the RCI as the relative area between the fitted curve and the line of equality. The curve used is
#'
#'   \deqn{y = \frac{1}{\left( \frac{1}{e^{k - 1}} - 1 \right)} \left( e^{\frac{x}{k - x}} - 1 \right)}
#'
#'   where \eqn{x} denotes cumulative population share, \eqn{y} denotes cumulative health indicator share, and \eqn{k} is a shape parameter estimated by minimizing squared deviations between observed and fitted values.
#'
#'   The RCI corresponds to the area between this fitted curve and the line of equality.
#'   }
#' }
#' }
#'
#' \subsection{Correction methods}{
#' When requested, the function returns two corrected indices derived from the standard RCI: the Wagstaff adjustment and the Erreygers correction.
#' The Wagstaff adjustment is implemented here as \eqn{RCI_{Wagstaff} = RCI / (1 - \mu)}, where \eqn{\mu} denotes the global mean of the health indicator; this rescales the RCI to account for upper bounds imposed by the mean when the indicator is bounded. The Erreygers correction is implemented here as \eqn{RCI_{Erreygers} = RCI \times (4\mu)}, producing an absolute-scale index (in percentage points when \eqn{\mu} is a proportion) suitable for comparisons of bounded indicators. Users should consult the original references for methodological details.
#' }
#'
#' @note Mode B is an approximation used when only aggregated indicators and population weights are available; variance estimates may be less precise than those from Mode A (numerator/denominator).
#'
#' @examples
#' # ——— Mode A (preferred: numerator and denominator) ———
#'
#' # Health indicator: Rate. Suppose data1 with:
#' #   maternal_deaths  = number of maternal deaths.
#' #   live_births      = number of live births.
#' #   ubn_percent      = % of population with at least one unmet basic need (higher = worse social status).
#' data(data1)
#' rci_luna(
#'   health_indicator_type    = "rate",
#'   health_numerator         = data1$maternal_deaths,
#'   health_denominator       = data1$live_births,
#'   equity_stratifier        = data1$ubn_percent,
#'   higher_ineq_is_favorable = FALSE,
#'   rate_scaling_factor      = 100000,
#'   calculation_method       = "WK",
#'   smoothing_factor         = 0.6,
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' # Health indicator: Proportion. Suppose data2 with:
#' #   skilled_births   = number of births attended by skilled personnel.
#' #   total_births     = total number of births.
#' #   ubn_index        = % of population with at least one unmet basic need (higher = worse social status).
#' data(data2)
#' rci_luna(
#'   health_indicator_type    = "proportion",
#'   health_numerator         = data2$skilled_births,
#'   health_denominator       = data2$total_births,
#'   equity_stratifier        = data2$ubn_percent,
#'   higher_ineq_is_favorable = FALSE,
#'   calculation_method       = "Fuller",
#'   apply_corrections        = TRUE,
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' # ——— Mode B (aggregated indicator + population as weight) ———
#'
#' # Suppose data1 with:
#' #   mmr       = maternal mortality ratio per 100,000 live births.
#' #   population = population weight of the unit of analysis.
#' #   ubn_index = % of population with at least one unmet basic need.
#' data(data1)
#' rci_luna(
#'   health_indicator_type    = "rate",
#'   health_indicator         = data1$mmr,
#'   equity_stratifier        = data1$ubn_percent,
#'   population_weights       = data1$population_zone,
#'   higher_ineq_is_favorable = FALSE,
#'   rate_scaling_factor      = 100000,
#'   calculation_method       = "WK",
#'   smoothing_factor         = 0.6,
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' # Suppose data2 with:
#' #   skilled_births_prop = proportion of births attended by skilled personnel (in [0,1]).
#' #   population           = population weight of the unit of analysis.
#' #   ubn_index            = % of population with at least one unmet basic need.
#' data(data2)
#' rci_luna(
#'   health_indicator_type    = "proportion",
#'   health_indicator         = data2$skilled_births_prop,
#'   equity_stratifier        = data2$ubn_percent,
#'   population_weights       = data2$population,
#'   higher_ineq_is_favorable = FALSE,
#'   calculation_method       = "Fuller",
#'   apply_corrections        = TRUE,
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' @seealso \code{stats::smooth.spline}, \code{stats::glm}, \code{optimx::optimx}
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

rci_luna <- function(
    health_indicator_type = NULL,
    health_indicator = NULL,
    population_weights = NULL,
    health_numerator = NULL,
    health_denominator = NULL,
    equity_stratifier = NULL,
    higher_ineq_is_favorable = NULL,
    rate_scaling_factor = NULL,
    calculation_method = NULL,
    apply_corrections = NULL,
    smoothing_factor = NULL,
    integration_method = NULL,
    conf_level = NULL,
    language_interpretation = NULL
) {

  # ───────────────────────────────────────────────
  # Argument validation (input consistency checks)
  # ───────────────────────────────────────────────
  # Check health indicator
  if (is.null(health_indicator_type) || !(health_indicator_type %in% c("proportion", "rate"))) {
    stop("`health_indicator_type` is required and must be either 'proportion' or 'rate'.")
  }

  if (is.null(conf_level)) {
    warning("`conf_level` not provided: defaulting to 0.95.")
    conf_level <- 0.95
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single numeric value between 0 and 1 (e.g. 0.95).")
  }

  if (is.null(calculation_method)) {
    warning("`calculation_method` not provided: defaulting to 'Fuller' (numerical integration).")
    calculation_method <- "WK"
  } else {
    if (!is.character(calculation_method) || length(calculation_method) != 1) {
      stop("`calculation_method` must be a single string: 'Fuller' or 'WK'.")
    }
    allowed_methods <- c("Fuller", "WK")
    if (!(calculation_method %in% allowed_methods)) {
      stop("`calculation_method` must be one of: ", paste(allowed_methods, collapse = ", "), ".")
    }
  }

  if (is.null(apply_corrections)) {
    warning("`apply_corrections` not provided: defaulting to FALSE.")
    apply_corrections <- FALSE
  } else if (!is.logical(apply_corrections) || length(apply_corrections) != 1) {
    stop("`apply_corrections` must be TRUE or FALSE.")
  }

  if (calculation_method == "WK") {
    if (is.null(smoothing_factor)) {
      warning("`smoothing_factor` not provided: defaulting to 0.5 (used only for calculation_method = 'WK').")
      smoothing_factor <- 0.5
    }
    if (!is.numeric(smoothing_factor) || length(smoothing_factor) != 1) {
      stop("`smoothing_factor` must be a single numeric value (e.g. 0.5) when calculation_method = 'WK'.")
    }
  } else {
    if (!is.null(smoothing_factor)) {
      warning("`smoothing_factor` provided but ignored because calculation_method != 'WK'.")
    }
    smoothing_factor <- NULL
  }

  # If user requests the 'best' integration method, enable search across available methods
  if (identical(integration_method, "best")) {
    warning("Searching across available optimization methods (this may increase runtime).")
    search_best_integration <- TRUE
  } else {
    search_best_integration <- FALSE
  }

  if (is.null(integration_method)) {
    warning("`integration_method` not provided: defaulting to 'BFGS'.")
    integration_method <- "BFGS"
  }

  allowed_integration <- c("BFGS","CG","Nelder-Mead","L-BFGS-B","nlm","nlminb","spg","ucminf",
                           "Rcgmin","Rvmmin","newuoa","bobyqa","nmkb","hjkb")
  if (!(integration_method %in% allowed_integration)) {
    stop("`integration_method` not allowed. Choose one of: ", paste(allowed_integration, collapse = ", "))
  }

  if (is.null(higher_ineq_is_favorable)) {
    stop("`higher_ineq_is_favorable` is required (TRUE if higher stratifier values = more favorable).")
  }
  if (!is.logical(higher_ineq_is_favorable) || length(higher_ineq_is_favorable) != 1) {
    stop("`higher_ineq_is_favorable` must be logical (TRUE/FALSE).")
  }

  # scaling factor rules
  if (health_indicator_type == "proportion") {
    if (!is.null(rate_scaling_factor)) warning("`rate_scaling_factor` is ignored because `health_indicator_type` is 'proportion'. It only applies to rates.")
    rate_scaling_factor <- 1
  } else {
    if (is.null(rate_scaling_factor)) stop("For `health_indicator_type = 'rate'` supply `rate_scaling_factor` (e.g., 1000 or 100000).")
    if (!is.numeric(rate_scaling_factor) || length(rate_scaling_factor) != 1 || rate_scaling_factor <= 0)
      stop("`rate_scaling_factor` must be a single positive numeric value.")
  }

  # equity stratifier validation
  if (is.null(equity_stratifier)) stop("`equity_stratifier` is required.")
  if (!is.numeric(equity_stratifier) || any(is.na(equity_stratifier))) stop("`equity_stratifier` must be numeric and complete.")

  ## Mode detection (must provide one of the two valid input modes)
  use_num_den <- !is.null(health_numerator) && !is.null(health_denominator)
  use_indicator_pop <- !use_num_den && !is.null(health_indicator) && !is.null(population_weights)

  if (!(use_num_den || use_indicator_pop)) {
    stop("You must provide either health_numerator & health_denominator (preferred) OR health_indicator & population_weights (approximate).")
  }

  ## Mode A checks
  if (use_num_den) {
    if (!is.numeric(health_numerator) || !is.numeric(health_denominator)) stop("health_numerator and health_denominator must be numeric.")
    if (any(is.na(health_numerator)) || any(is.na(health_denominator))) stop("Missing values in health_numerator or health_denominator.")
    if (any(health_denominator <= 0, na.rm = TRUE)) stop("All health_denominator values must be > 0.")
    if (any(health_numerator < 0, na.rm = TRUE)) stop("health_numerator must be >= 0.")
    if (!is.null(health_indicator)) {
      warning("Both counts and indicator provided: the indicator will be ignored and derived from numerator/denominator.")
      health_indicator <- NULL
    }
    if (length(health_numerator) != length(health_denominator)) stop("numerator and denominator must have same length.")
    if (length(health_numerator) != length(equity_stratifier)) stop("numerator must have same length as equity_stratifier.")
  }

  ## Mode B checks
  if (use_indicator_pop) {
    if (!is.numeric(health_indicator) || any(is.na(health_indicator))) stop("health_indicator must be numeric and complete.")
    if (!is.numeric(population_weights) || any(is.na(population_weights)) || any(population_weights <= 0)) stop("population_weights must be numeric, >0, and complete.")
    if (health_indicator_type == "proportion" && any(health_indicator < 0 | health_indicator > 1, na.rm = TRUE))
      stop("For proportions, health_indicator must be in [0,1].")
    if (length(health_indicator) != length(population_weights)) stop("health_indicator and population_weights must have same length.")
    if (length(population_weights) != length(equity_stratifier)) stop("population_weights and equity_stratifier must have same length.")
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
  } else if (!(language_interpretation %in% names(alias_map))) {
    stop("Unrecognized `language_interpretation`. Allowed aliases: ",
         paste(names(alias_map), collapse = ", "),
         ". Allowed full names: ",
         paste(unname(alias_map), collapse = ", "), ".")
  }

  # ──────────────────────────────
  # Base data frame preparation
  # ──────────────────────────────
  if (use_num_den) {
    # indicator scaled according to type
    if (health_indicator_type == "proportion") {
      indicator_vec <- health_numerator / health_denominator
    } else {
      indicator_vec <- (health_numerator / health_denominator) * rate_scaling_factor
    }
    df_master <- tibble::tibble(
      indicator    = as.numeric(indicator_vec),
      equity       = as.numeric(equity_stratifier),
      numerator    = as.numeric(round(health_numerator)),
      denominator  = as.numeric(round(health_denominator)),
      pop_use      = as.numeric(round(health_denominator))   # population to use = denominator when available
    )
  } else {
    df_master <- tibble::tibble(
      indicator = as.numeric(health_indicator),
      equity    = as.numeric(equity_stratifier),
      pop_use   = as.numeric(population_weights)
    )
  }

  # ──────────────────────────────────────
  # Ordering by the inequality variable
  # ──────────────────────────────────────
  df_master <- df_master %>%
    mutate(equity_order = if (isTRUE(higher_ineq_is_favorable)) equity else -equity) %>%
    arrange(equity_order)

  # ──────────────────────────────────────────────────────
  # Compute population used and social position (ridit)
  # ──────────────────────────────────────────────────────
  df_master <- df_master %>%
    mutate(wpop  = pop_use / sum(pop_use, na.rm = TRUE),
             cwpop = cumsum(wpop),
             cumw  = cumsum(pop_use),
             cumw1 = lag(cumw, default = 0),
             sumw  = sum(pop_use, na.rm = TRUE)) %>%
      group_by(equity_order) %>%
      mutate(cumwr         = max(cumw, na.rm = TRUE),
             cumwr1        = min(cumw1, na.rm = TRUE),
             social_position = (cumwr1 + 0.5 * (cumwr - cumwr1)) / sumw) %>%
      ungroup() %>%
      select(-cumw, -cumw1, -sumw, -cumwr, -cumwr1)

  # ────────────────────────────────────────────────
  # Global weighted means for the health indicator
  # ────────────────────────────────────────────────
  if (use_num_den) {
    total_num <- sum(df_master$numerator, na.rm = TRUE)
    total_den <- sum(df_master$denominator, na.rm = TRUE)
    global_health_mean <- if (health_indicator_type == "proportion") total_num / total_den else (total_num / total_den) * rate_scaling_factor
  } else {
    global_health_mean <- weighted.mean(df_master$indicator, w = df_master$pop_use, na.rm = TRUE)
  }

  # ────────────────────────────────
  # Branch by calculation_method
  # ────────────────────────────────

  if (calculation_method == "Fuller") {

    df <- df_master %>%
      arrange(equity_order)

    if ("numerator" %in% names(df)) {
      df <- df %>% mutate(f_indicator = numerator)
    } else {
      df <- df %>% mutate(f_indicator = indicator * pop_use / rate_scaling_factor)
    }

    df <- df %>% mutate(w_indicator = f_indicator / sum(f_indicator, na.rm = TRUE),
                        cwi = cumsum(w_indicator))

    func_lorenz <- function(x, k) (exp(x / (k - x)) - 1) / (exp(1 / (k - 1)) - 1)
    loss_lorenz <- function(k, data) sum((data$y - func_lorenz(data$x, k))^2)

    p0 <- c(10, -10)
    results_opt <- purrr::map_dfr(p0, function(p) {
      optimx::optimx(par = p,
                     fn = loss_lorenz,
                     data = rbind(c(0,0), df %>% dplyr::select(y = cwi, x = cwpop)),
                     method = if (search_best_integration) "all.methods" else integration_method,
                     control = list(all.methods = search_best_integration, save.failures = TRUE, maxit = 2500)) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "method") %>%
        mutate(par0 = p)
    })

    best_fit <- results_opt %>% filter(convcode == 0) %>% arrange(value) %>% slice_head(n = 1)
    if (nrow(best_fit) == 0) best_fit <- results_opt %>% arrange(value) %>% slice_head(n = 1)
    k_hat <- as.numeric(best_fit$p1[1])
    method_used <- as.character(best_fit$method[1])

    grid_x <- seq(0, 1, by = 0.01)
    grid_y <- func_lorenz(grid_x, k_hat)
    rci_est <- 2 * sum((grid_x - grid_y) * 0.01)

    df <- df %>% mutate(func_lorenz = func_lorenz(cwpop, k_hat),
                        q = cwi / global_health_mean,
                        a = (indicator / global_health_mean) * (2 * social_position - 1 - rci_est) + 2 - lag(q, 1, default = 0) - q,
                        fa2 = wpop * a^2)

    n <- nrow(df)
    rci_se <- sqrt(abs((1 / n) * (sum(df$fa2, na.rm = TRUE) - (1 + rci_est)^2))) / sqrt(n)
    z <- qnorm(1 - (1 - conf_level) / 2)
    rci_low <- rci_est - z * rci_se
    rci_up  <- rci_est + z * rci_se

    results_tbl <- tibble::tibble(metric = "Relative concentration index",
                                  value = rci_est,
                                  ci_lower = rci_low,
                                  ci_upper = rci_up,
                                  method = "Fuller")




    if (isTRUE(apply_corrections)) {
      mu <- global_health_mean
      den <- 1 - mu; if (den <= .Machine$double.eps) den <- .Machine$double.eps

      est_wag <- rci_est / den
      ci_wag_low <- rci_low / den
      ci_wag_up <- rci_up / den

      est_err <- rci_est * (4 * mu)
      ci_err_low <- rci_est * (4 * mu)
      ci_err_up <- rci_up * (4 * mu)

      results_tbl <- bind_rows(results_tbl,
                               tibble::tibble(metric = "Relative concentration index (Wagstaff adj.)",
                                              value = est_wag,
                                              ci_lower = ci_wag_low,
                                              ci_upper = ci_wag_up,
                                              method = "Fuller"),
                               tibble::tibble(metric = "Relative concentration index (Erreygers adj.)",
                                              value = est_err,
                                              ci_lower = ci_err_low,
                                              ci_upper = ci_err_up,
                                              method = "Fuller"))
    }

  } else if (calculation_method == "WK") {

    df <- df_master %>%
      arrange(equity_order)

    if ("numerator" %in% names(df)) {
      df <- df %>% mutate(f_indicator = numerator)
    } else {
      df <- df %>% mutate(f_indicator = indicator * pop_use / rate_scaling_factor)
    }

    df <- df %>% mutate(w_indicator = f_indicator / sum(f_indicator, na.rm = TRUE),
                        cwi = cumsum(w_indicator))


    #Smoothing
    nodos <- tibble::tibble(cwpop = c(0, df$cwpop, 1),
                            cwi = c(0, df$cwi, 1)) %>%
             arrange(cwpop) %>%
             distinct(cwpop, .keep_all = TRUE)

    fit_smooth <- stats::smooth.spline(x = nodos$cwpop, y = nodos$cwi, spar = smoothing_factor)
    df$func_lorenz <- stats::predict(fit_smooth, x = df$cwpop)$y
    df$func_lorenz <- pmin(pmax(df$func_lorenz, 0), 1)
    iso <- stats::isoreg(df$cwpop, df$func_lorenz)
    df$func_lorenz <- pmin(pmax(approx(iso$x, iso$yf, xout = df$cwpop)$y, 0), 1)

    df <- df %>%
      mutate(
        sumw     = sum(pop_use, na.rm = TRUE),
        sigma1   = sum((pop_use / sumw) * ((social_position - 0.5)^2), na.rm = TRUE),
        lhs      = 2 * sigma1 * (indicator / global_health_mean),
        rhs      = social_position
      )

    mod <- tryCatch(
      glm(lhs ~ rhs, family = gaussian(), data = df, weights = pop_use),
      error = function(e) stop("Regression for WK failed: ", e$message)
    )

    rci_est <- as.numeric(coef(mod)[["rhs"]])
    ci_mod <- tryCatch(confint.default(mod, parm = "rhs", level = conf_level), error = function(e) c(NA_real_, NA_real_))

    rci_low <- ci_mod[1]
    rci_up <- ci_mod[2]

    results_tbl <- tibble::tibble(metric = "Relative concentration index",
                                  value = rci_est,
                                  ci_lower = rci_low,
                                  ci_upper = rci_up,
                                  method = "Wagstaff & Kakwani")

    if (apply_corrections) {
      den <- 1 - global_health_mean; if (den <= .Machine$double.eps) den <- .Machine$double.eps

      # Wagstaff
      df <- df %>% mutate(lhs_wag = lhs / den)
      mod_wag <- glm(lhs_wag ~ rhs, family = gaussian(), data = df, weights = pop_use)
      est_wag <- as.numeric(coef(mod_wag)[["rhs"]])
      ci_wag  <- tryCatch(confint.default(mod_wag, parm = "rhs", level = conf_level),
                          error = function(e) c(NA_real_, NA_real_))

      ci_wag_low <- ci_wag[1]
      ci_wag_up <- ci_wag[2]

      # Erreygers
      k_val <- 4 * global_health_mean
      df <- df %>% mutate(lhs_err = lhs * k_val)
      mod_err <- glm(lhs_err ~ rhs, family = gaussian(), data = df, weights = pop_use)
      est_err <- as.numeric(coef(mod_err)[["rhs"]])
      ci_err  <- tryCatch(confint.default(mod_err, parm = "rhs", level = conf_level),
                          error = function(e) c(NA_real_, NA_real_))

      ci_err_low <- ci_err[1]
      ci_err_up <- ci_err[2]

      results_tbl <- bind_rows(results_tbl,
                               tibble::tibble(metric = "Relative concentration index (Wagstaff adj.)",
                                              value = est_wag,
                                              ci_lower = ci_wag_low,
                                              ci_upper = ci_wag_up,
                                              method = "Wagstaff & Kakwani"),
                               tibble::tibble(metric = "Relative concentration index (Erreygers adj.)",
                                              value = est_err,
                                              ci_lower = ci_err_low,
                                              ci_upper = ci_err_up,
                                              method = "Wagstaff & Kakwani"))
    }
  }

  df <- df %>%
    select(health_indicator = indicator,
           equity_stratifier = equity,
           cumulative_population = cwpop,
           social_position,
           cumulative_event_health  = cwi,
           func_lorenz
           )

  # ────────────────
  # Interpretations
  # ────────────────

  if (!apply_corrections) {

    part_main_english <- if (rci_est < 0) {
      paste0(
        "The Relative Concentration Index is ", round(rci_est * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(rci_low * 100, 2), "%; ", round(rci_up * 100, 2), "%]. ",
        "This means that the Health Indicator is concentrated among the most disadvantaged. ",
        "In relative terms, the concentration among the most disadvantaged corresponds to ", abs(round(rci_est * 100, 2)),
        "% of the maximum possible concentration."
      )
    } else if (rci_est > 0) {
      paste0(
        "The Relative Concentration Index is ", round(rci_est * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(rci_low * 100, 2), "%; ", round(rci_up * 100, 2), "%]. ",
        "This means that the Health Indicator is concentrated among the most advantaged. ",
        "In relative terms, the concentration among the most advantaged corresponds to ", abs(round(rci_est * 100, 2)),
        "% of the maximum possible concentration."
      )
    } else {
      paste0(
        "The Relative Concentration Index is ", round(rci_est * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(rci_low * 100, 2), "%; ", round(rci_up * 100, 2), "%]. ",
        "This indicates that, in relative terms, the Health Indicator is distributed equitably across socioeconomic groups, with no concentration in either the most or least disadvantaged."
      )
    }

    part_sig_english <- if (rci_low <= 0 & rci_up >= 0) {
      " However, this concentration is not statistically significant."
    } else {
      " Moreover, this concentration is statistically significant."
    }

    interpretation_english <- paste0(part_main_english, part_sig_english)

    note_english <- if (rci_est < 0) {
      "Note: The maximum possible concentration would occur if the entire Health Indicator were concentrated among the most disadvantaged group."
    } else {
      "Note: The maximum possible concentration would occur if the entire Health Indicator were concentrated among the most advantaged group."
    }

  } else {

    part_main_english_1 <- if (rci_est < 0) {
      paste0(
        "The Relative Concentration Index (Standard) is ", round(rci_est * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(rci_low * 100, 2), "%; ", round(rci_up * 100, 2), "%]. ",
        "This means that the Health Indicator is concentrated among the most disadvantaged. ",
        "In relative terms, the concentration among the most disadvantaged corresponds to ", abs(round(rci_est * 100, 2)),
        "% of the maximum possible concentration."
      )
    } else if (rci_est > 0) {
      paste0(
        "The Relative Concentration Index (Standard) is ", round(rci_est * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(rci_low * 100, 2), "%; ", round(rci_up * 100, 2), "%]. ",
        "This means that the Health Indicator is concentrated among the most advantaged. ",
        "In relative terms, the concentration among the most advantaged corresponds to ", abs(round(rci_est * 100, 2)),
        "% of the maximum possible concentration."
      )
    } else {
      paste0(
        "The Relative Concentration Index (Standard) is ", round(rci_est * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(rci_low * 100, 2), "%; ", round(rci_up * 100, 2), "%]. ",
        "This indicates that, in relative terms, the Health Indicator is equitably distributed across socioeconomic groups, without concentration among either the most or least disadvantaged."
      )
    }

    part_sig_english_1 <- if (rci_low <= 0 & rci_up >= 0) {
      " However, this concentration is not statistically significant."
    } else {
      " Moreover, this concentration is statistically significant."
    }

    interpretation_english_1 <- paste0(part_main_english_1, part_sig_english_1)

    part_main_english_2 <- if (est_wag < 0) {
      paste0(
        "The Wagstaff Concentration Index is ", round(est_wag * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(ci_wag_low * 100, 2), "%; ", round(ci_wag_up * 100, 2), "%]. ",
        "This negative value indicates concentration of the indicator among the most disadvantaged groups. ",
        "In relative terms, the magnitude corresponds to ", abs(round(est_wag * 100, 2)),
        "% of the maximum feasible inequality given the mean of the indicator."
      )
    } else if (est_wag > 0) {
      paste0(
        "The Wagstaff Concentration Index is ", round(est_wag * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(ci_wag_low * 100, 2), "%; ", round(ci_wag_up * 100, 2), "%]. ",
        "This positive value indicates concentration of the indicator among the most advantaged groups. ",
        "In relative terms, the magnitude corresponds to ", abs(round(est_wag * 100, 2)),
        "% of the maximum feasible inequality given the mean of the indicator."
      )
    } else {
      paste0(
        "The Wagstaff Concentration Index is ", round(est_wag * 100, 2),
        "%, ", conf_level * 100, "%CI [", round(ci_wag_low * 100, 2), "%; ", round(ci_wag_up * 100, 2), "%]. ",
        "The null value is consistent with the absence of a socioeconomic gradient. ",
        "In relative terms, the magnitude corresponds to 0% of the maximum feasible inequality given the mean of the indicator."
      )
    }

    part_sig_english_2 <- if (ci_wag_low <= 0 & ci_wag_up >= 0) {
      " However, the estimated inequality is not statistically significant."
    } else {
      " Moreover, the estimated inequality is statistically significant."
    }

    interpretation_english_2 <- paste0(part_main_english_2, part_sig_english_2)

    part_main_english_3 <- if (est_err < 0) {
      paste0(
        "The Erreygers Concentration Index is ", round(est_err * 100, 2),
        " pp, ", conf_level * 100, "%CI [", round(ci_err_low * 100, 2), " pp; ", round(ci_err_up * 100, 2), " pp]. ",
        "This negative value implies that, across the socioeconomic gradient, the health indicator is on average ",
        abs(round(est_err * 100, 2)), " percentage points higher among the lower socioeconomic strata than among the upper strata, ",
        "reflecting a pro-poor absolute concentration of the indicator."
      )
    } else if (est_err > 0) {
      paste0(
        "The Erreygers Concentration Index is ", round(est_err * 100, 2),
        " pp, ", conf_level * 100, "%CI [", round(ci_err_low * 100, 2), " pp; ", round(ci_err_up * 100, 2), " pp]. ",
        "This positive value implies that, across the socioeconomic gradient, the health indicator is on average ",
        abs(round(est_err * 100, 2)), " percentage points higher among the upper socioeconomic strata than among the lower strata, ",
        "reflecting a pro-rich absolute concentration of the indicator."
      )
    } else {
      paste0(
        "The Erreygers Concentration Index is ", round(est_err * 100, 2),
        " pp, ", conf_level * 100, "%CI [", round(ci_err_low * 100, 2), " pp; ", round(ci_err_up * 100, 2), " pp]. ",
        "The value close to zero indicates that there are no relevant absolute differences in the indicator level between the extremes of the socioeconomic gradient."
      )
    }

    part_sig_english_3 <- if (ci_err_low <= 0 & ci_err_up >= 0) {
      " However, the estimated inequality is not statistically significant."
    } else {
      " Moreover, the estimated inequality is statistically significant."
    }

    interpretation_english_3 <- paste0(part_main_english_3, part_sig_english_3)

    note_english <- if (rci_est < 0) {
      "Note: The maximum possible concentration would occur if the entire Health Indicator were concentrated among the most disadvantaged group."
    } else {
      "Note: The maximum possible concentration would occur if the entire Health Indicator were concentrated among the most advantaged group."
    }

    interpretation_english <- c(interpretation_english_1, interpretation_english_2, interpretation_english_3)
  }

  if (!apply_corrections) {

    part_main_spanish <- if (rci_est < 0) {
      paste0(
        "El Índice de Concentración Relativo es ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Esto significa que el Indicador de Salud se concentra en los más desfavorecidos. ",
        "En términos relativos, la concentración en los más desfavorecidos es ", abs(round(rci_est*100, 2)),
        "% de la concentración máxima posible."
      )
    } else if (rci_est > 0) {
      paste0(
        "El Índice de Concentración Relativo es ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Esto significa que el Indicador de Salud se concentra en los más favorecidos. ",
        "En términos relativos, la concentración en los más favorecidos es ", abs(round(rci_est*100, 2)),
        "% de la concentración máxima posible."
      )
    } else if (rci_est == 0) {
      paste0(
        "El Índice de Concentración Relativo es de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Esto indica que, en términos relativos, el Indicador de Salud está distribuido de manera equitativa entre los grupos socioeconómicos, sin concentrarse ni en los más desfavorecidos ni en los menos desfavorecidos."
      )
    }

    part_sig_spanish <- if (rci_low <= 0 & rci_up >= 0) {
      paste0(" Sin embargo, esta concentración no es estadísticamente significativa.")
    } else {
      paste0(" Además, esta concentración es estadísticamente significativa.")
    }

    interpretation_spanish <- paste0(part_main_spanish, part_sig_spanish)

    # Nota de interpretación
    note_spanish <- if (rci_est < 0) {
      paste0(
        "Nota: La concentración máxima posible es aquella que ocurriría si todo el Indicador de Salud se concentraría en el extremo de los más desfavorecidos."
      )
    } else {
      paste0(
        "Nota: La concentración máxima posible es aquella que ocurriría si todo el Indicador de Salud se concentraría en el extremo de los más favorecidos."
      )
    }

  } else {
    part_main_spanish_1 <- if (rci_est < 0) {
      paste0(
        "El Índice de Concentración Relativo (Estándar) es ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Esto significa que el Indicador de Salud se concentra en los más desfavorecidos. ",
        "En términos relativos, la concentración en los más desfavorecidos es ", abs(round(rci_est*100, 2)),
        "% de la concentración máxima posible."
      )
    } else if (rci_est > 0) {
      paste0(
        "El Índice de Concentración Relativo (Estándar) es ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Esto significa que el Indicador de Salud se concentra en los más favorecidos. ",
        "En términos relativos, la concentración en los más favorecidos es ", abs(round(rci_est*100, 2)),
        "% de la concentración máxima posible."
      )
    } else if (rci_est == 0) {
      paste0(
        "El Índice de Concentración Relativo (Estándar) es de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Esto indica que, en términos relativos, el Indicador de Salud está distribuido de manera equitativa entre los grupos socioeconómicos, sin concentrarse ni en los más desfavorecidos ni en los menos desfavorecidos."
      )
    }

    part_sig_spanish_1 <- if (rci_low <= 0 & rci_up >= 0) {
      paste0(" Sin embargo, esta concentración no es estadísticamente significativa.")
    } else {
      paste0(" Además, esta concentración es estadísticamente significativa.")
    }

    interpretation_spanish_1 <- paste0(part_main_spanish_1, part_sig_spanish_1)

    part_main_spanish_2 <- if (est_wag < 0) {
      paste0(
        "El Índice de Concentración según Wagstaff es ", round(est_wag*100, 2),
        "%, ", conf_level*100, "%IC [", round(ci_wag_low*100, 2), "%; ", round(ci_wag_up*100, 2), "%]. ",
        "Este valor negativo indica concentración del indicador en los grupos más desfavorecidos. ",
        "En términos relativos, la magnitud equivale al ", abs(round(est_wag*100, 2)),
        "% de la desigualdad máxima factible dada la media del indicador."
      )
    } else if (est_wag > 0) {
      paste0(
        "El Índice de Concentración según Wagstaff es ", round(est_wag*100, 2),
        "%, ", conf_level*100, "%IC [", round(ci_wag_low*100, 2), "%; ", round(ci_wag_up*100, 2), "%]. ",
        "Este valor positivo indica concentración del indicador en los grupos más favorecidos. ",
        "En términos relativos, la magnitud equivale al ", abs(round(est_wag*100, 2)),
        "% de la desigualdad máxima factible dada la media del indicador."
      )
    } else {
      paste0(
        "El Índice de Concentración según Wagstaff es ", round(est_wag*100, 2),
        "%, ", conf_level*100, "%IC [", round(ci_wag_low*100, 2), "%; ", round(ci_wag_up*100, 2), "%]. ",
        "El valor nulo es consistente con ausencia de gradiente socioeconómico. ",
        "En términos relativos, la magnitud equivale al 0% de la desigualdad máxima factible dada la media del indicador."
      )
    }

    part_sig_spanish_2 <- if (ci_wag_low <= 0 & ci_wag_up >= 0) {
      " No obstante, la desigualdad estimada no es estadísticamente significativa."
    } else {
      " Además, la desigualdad estimada es estadísticamente significativa."
    }

    interpretation_spanish_2 <- paste0(part_main_spanish_2, part_sig_spanish_2)

    part_main_spanish_3 <- if (est_err < 0) {
      paste0(
        "El índice de Concentración de Erreygers es ", round(est_err*100, 2),
        "pp, ", conf_level*100, "% IC [", round(ci_err_low*100, 2), "pp; ", round(ci_err_up*100, 2), "pp]. ",
        "Este valor negativo implica que, a lo largo del gradiente socioeconómico, el indicador de salud es en promedio ",
        abs(round(est_err*100, 2)), " puntos porcentuales más alto en los estratos socioeconómicos inferiores que en los superiores, ",
        "lo que refleja una concentración absoluta pro-pobre del indicador."
      )
    } else if (est_err > 0) {
      paste0(
        "El índice de Concentración de Erreygers es ", round(est_err*100, 2),
        "pp, ", conf_level*100, "% IC [", round(ci_err_low*100, 2), "pp; ", round(ci_err_up*100, 2), "pp]. ",
        "Este valor positivo implica que, a lo largo del gradiente socioeconómico, el indicador de salud es en promedio ",
        abs(round(est_err*100, 2)), " puntos porcentuales más alto en los estratos socioeconómicos superiores que en los inferiores, ",
        "lo que refleja una concentración absoluta pro-rico del indicador."
      )
    } else {
      paste0(
        "El índice de Concentración de Erreygers es ", round(est_err*100, 2),
        "pp, ", conf_level*100, "% IC [", round(ci_err_low*100, 2), "pp; ", round(ci_err_up*100, 2), "pp]. ",
        "El valor cercano a cero indica que no se observan diferencias absolutas relevantes en el nivel del indicador entre los extremos del gradiente socioeconómico."
      )
    }

    part_sig_spanish_3 <- if (ci_err_low <= 0 & ci_err_up >= 0) {
      " No obstante, la desigualdad estimada no alcanza significación estadística."
    } else {
      " Además, la desigualdad estimada es estadísticamente significativa."
    }

    interpretation_spanish_3 <- paste0(part_main_spanish_3, part_sig_spanish_3)


    note_spanish <- if (rci_est < 0) {
      paste0(
        "Nota: La concentración máxima posible es aquella que ocurriría si todo el Indicador de Salud se concentraría en el extremo de los más desfavorecidos."
      )
    } else {
      paste0(
        "Nota: La concentración máxima posible es aquella que ocurriría si todo el Indicador de Salud se concentraría en el extremo de los más favorecidos."
      )
    }

    interpretation_spanish = c(interpretation_spanish_1, interpretation_spanish_2, interpretation_spanish_3)
  }

  if (!apply_corrections) {

    part_main_french <- if (rci_est < 0) {
      paste0(
        "L’Indice de Concentration Relatif est de ", round(rci_est*100, 2),
        "%, IC à ", conf_level*100, "% [", round(rci_low*100, 2), "% ; ", round(rci_up*100, 2), "%]. ",
        "Cela signifie que l’indicateur de santé est concentré parmi les plus défavorisés. ",
        "En termes relatifs, la concentration chez les plus défavorisés correspond à ", abs(round(rci_est*100, 2)),
        "% de la concentration maximale possible."
      )
    } else if (rci_est > 0) {
      paste0(
        "L’Indice de Concentration Relatif est de ", round(rci_est*100, 2),
        "%, IC à ", conf_level*100, "% [", round(rci_low*100, 2), "% ; ", round(rci_up*100, 2), "%]. ",
        "Cela signifie que l’indicateur de santé est concentré parmi les plus favorisés. ",
        "En termes relatifs, la concentration chez les plus favorisés correspond à ", abs(round(rci_est*100, 2)),
        "% de la concentration maximale possible."
      )
    } else {
      paste0(
        "L’Indice de Concentration Relatif est de ", round(rci_est*100, 2),
        "%, IC à ", conf_level*100, "% [", round(rci_low*100, 2), "% ; ", round(rci_up*100, 2), "%]. ",
        "Cela indique que, en termes relatifs, l’indicateur de santé est distribué de manière équitable entre les groupes socioéconomiques, sans concentration vers les plus favorisés ni les plus défavorisés."
      )
    }

    part_sig_french <- if (rci_low <= 0 & rci_up >= 0) {
      " Cependant, cette concentration n’est pas statistiquement significative."
    } else {
      " De plus, cette concentration est statistiquement significative."
    }

    interpretation_french <- paste0(part_main_french, part_sig_french)

    note_french <- if (rci_est < 0) {
      "Note : La concentration maximale possible correspondrait à une situation où tout l’indicateur de santé serait concentré chez les plus défavorisés."
    } else {
      "Note : La concentration maximale possible correspondrait à une situation où tout l’indicateur de santé serait concentré chez les plus favorisés."
    }

  } else {

    part_main_french_1 <- if (rci_est < 0) {
      paste0(
        "L’Indice de Concentration Relatif (Standardisé) est de ", round(rci_est*100, 2),
        "%, IC à ", conf_level*100, "% [", round(rci_low*100, 2), "% ; ", round(rci_up*100, 2), "%]. ",
        "Cela signifie que l’indicateur de santé est concentré parmi les plus défavorisés. ",
        "En termes relatifs, cette concentration chez les plus défavorisés représente ", abs(round(rci_est*100, 2)),
        "% de la concentration maximale possible."
      )
    } else if (rci_est > 0) {
      paste0(
        "L’Indice de Concentration Relatif (Standardisé) est de ", round(rci_est*100, 2),
        "%, IC à ", conf_level*100, "% [", round(rci_low*100, 2), "% ; ", round(rci_up*100, 2), "%]. ",
        "Cela signifie que l’indicateur de santé est concentré parmi les plus favorisés. ",
        "En termes relatifs, cette concentration chez les plus favorisés représente ", abs(round(rci_est*100, 2)),
        "% de la concentration maximale possible."
      )
    } else {
      paste0(
        "L’Indice de Concentration Relatif (Standardisé) est de ", round(rci_est*100, 2),
        "%, IC à ", conf_level*100, "% [", round(rci_low*100, 2), "% ; ", round(rci_up*100, 2), "%]. ",
        "Cela indique une répartition équitable de l’indicateur de santé entre les groupes socioéconomiques, sans concentration marquée."
      )
    }

    part_sig_french_1 <- if (rci_low <= 0 & rci_up >= 0) {
      " Cependant, cette concentration n’est pas statistiquement significative."
    } else {
      " De plus, cette concentration est statistiquement significative."
    }

    interpretation_french_1 <- paste0(part_main_french_1, part_sig_french_1)

    part_main_french_2 <- if (est_wag < 0) {
      paste0(
        "L’Indice de Concentration de Wagstaff est de ", round(est_wag*100, 2),
        "%, IC à ", conf_level*100, "% [", round(ci_wag_low*100, 2), "% ; ", round(ci_wag_up*100, 2), "%]. ",
        "Cette valeur négative indique une concentration de l’indicateur parmi les plus défavorisés. ",
        "En termes relatifs, cela correspond à ", abs(round(est_wag*100, 2)),
        "% de l’inégalité maximale réalisable selon la moyenne de l’indicateur."
      )
    } else if (est_wag > 0) {
      paste0(
        "L’Indice de Concentration de Wagstaff est de ", round(est_wag*100, 2),
        "%, IC à ", conf_level*100, "% [", round(ci_wag_low*100, 2), "% ; ", round(ci_wag_up*100, 2), "%]. ",
        "Cette valeur positive indique une concentration de l’indicateur parmi les plus favorisés. ",
        "En termes relatifs, cela correspond à ", abs(round(est_wag*100, 2)),
        "% de l’inégalité maximale réalisable selon la moyenne de l’indicateur."
      )
    } else {
      paste0(
        "L’Indice de Concentration de Wagstaff est de ", round(est_wag*100, 2),
        "%, IC à ", conf_level*100, "% [", round(ci_wag_low*100, 2), "% ; ", round(ci_wag_up*100, 2), "%]. ",
        "Cette valeur nulle est cohérente avec l’absence de gradient socioéconomique."
      )
    }

    part_sig_french_2 <- if (ci_wag_low <= 0 & ci_wag_up >= 0) {
      " Cependant, l’inégalité estimée n’est pas statistiquement significative."
    } else {
      " De plus, l’inégalité estimée est statistiquement significative."
    }

    interpretation_french_2 <- paste0(part_main_french_2, part_sig_french_2)

    part_main_french_3 <- if (est_err < 0) {
      paste0(
        "L’Indice de Concentration d’Erreygers est de ", round(est_err*100, 2),
        " points de pourcentage, IC à ", conf_level*100, "% [", round(ci_err_low*100, 2), " pp ; ", round(ci_err_up*100, 2), " pp]. ",
        "Cette valeur négative indique qu’en moyenne, l’indicateur est ", abs(round(est_err*100, 2)),
        " points de pourcentage plus élevé dans les groupes socioéconomiques les plus défavorisés, ",
        "ce qui reflète une concentration absolue pro-pauvres."
      )
    } else if (est_err > 0) {
      paste0(
        "L’Indice de Concentration d’Erreygers est de ", round(est_err*100, 2),
        " points de pourcentage, IC à ", conf_level*100, "% [", round(ci_err_low*100, 2), " pp ; ", round(ci_err_up*100, 2), " pp]. ",
        "Cette valeur positive indique qu’en moyenne, l’indicateur est ", abs(round(est_err*100, 2)),
        " points de pourcentage plus élevé dans les groupes socioéconomiques les plus favorisés, ",
        "ce qui reflète une concentration absolue pro-riches."
      )
    } else {
      paste0(
        "L’Indice de Concentration d’Erreygers est de ", round(est_err*100, 2),
        " points de pourcentage, IC à ", conf_level*100, "% [", round(ci_err_low*100, 2), " pp ; ", round(ci_err_up*100, 2), " pp]. ",
        "Une valeur proche de zéro suggère l’absence de disparités absolues marquées entre les extrémités du gradient socioéconomique."
      )
    }

    part_sig_french_3 <- if (ci_err_low <= 0 & ci_err_up >= 0) {
      " Cependant, l’inégalité estimée n’est pas statistiquement significative."
    } else {
      " De plus, l’inégalité estimée est statistiquement significative."
    }

    interpretation_french_3 <- paste0(part_main_french_3, part_sig_french_3)

    note_french <- if (rci_est < 0) {
      "Note : La concentration maximale possible correspondrait à une situation où tout l’indicateur de santé serait concentré chez les plus défavorisés."
    } else {
      "Note : La concentration maximale possible correspondrait à une situation où tout l’indicateur de santé serait concentré chez les plus favorisés."
    }

    interpretation_french <- c(interpretation_french_1, interpretation_french_2, interpretation_french_3)
  }


  if (!apply_corrections) {

    part_main_portuguese <- if (rci_est < 0) {
      paste0(
        "O Índice de Concentração Relativo é de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Isso significa que o Indicador de Saúde está concentrado nos grupos mais desfavorecidos. ",
        "Em termos relativos, a concentração nos mais desfavorecidos equivale a ", abs(round(rci_est*100, 2)),
        "% da concentração máxima possível."
      )
    } else if (rci_est > 0) {
      paste0(
        "O Índice de Concentração Relativo é de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Isso significa que o Indicador de Saúde está concentrado nos grupos mais favorecidos. ",
        "Em termos relativos, a concentração nos mais favorecidos equivale a ", abs(round(rci_est*100, 2)),
        "% da concentração máxima possível."
      )
    } else {
      paste0(
        "O Índice de Concentração Relativo é de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Isso indica que, em termos relativos, o Indicador de Saúde está distribuído de forma equitativa entre os grupos socioeconômicos, sem concentração nos mais nem nos menos favorecidos."
      )
    }

    part_sig_portuguese <- if (rci_low <= 0 & rci_up >= 0) {
      " No entanto, essa concentração não é estatisticamente significativa."
    } else {
      " Além disso, essa concentração é estatisticamente significativa."
    }

    interpretation_portuguese <- paste0(part_main_portuguese, part_sig_portuguese)

    note_portuguese <- if (rci_est < 0) {
      "Nota: A concentração máxima possível ocorreria se todo o Indicador de Saúde estivesse concentrado no extremo dos mais desfavorecidos."
    } else {
      "Nota: A concentração máxima possível ocorreria se todo o Indicador de Saúde estivesse concentrado no extremo dos mais favorecidos."
    }

  } else {

    part_main_portuguese_1 <- if (rci_est < 0) {
      paste0(
        "O Índice de Concentração Relativo (Padrão) é de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Isso significa que o Indicador de Saúde está concentrado nos grupos mais desfavorecidos. ",
        "Em termos relativos, a concentração nos mais desfavorecidos equivale a ", abs(round(rci_est*100, 2)),
        "% da concentração máxima possível."
      )
    } else if (rci_est > 0) {
      paste0(
        "O Índice de Concentração Relativo (Padrão) é de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Isso significa que o Indicador de Saúde está concentrado nos grupos mais favorecidos. ",
        "Em termos relativos, a concentração nos mais favorecidos equivale a ", abs(round(rci_est*100, 2)),
        "% da concentração máxima possível."
      )
    } else {
      paste0(
        "O Índice de Concentração Relativo (Padrão) é de ", round(rci_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
        "Isso indica que, em termos relativos, o Indicador de Saúde está distribuído de forma equitativa entre os grupos socioeconômicos, sem concentração nos mais nem nos menos favorecidos."
      )
    }

    part_sig_portuguese_1 <- if (rci_low <= 0 & rci_up >= 0) {
      " No entanto, essa concentração não é estatisticamente significativa."
    } else {
      " Além disso, essa concentração é estatisticamente significativa."
    }

    interpretation_portuguese_1 <- paste0(part_main_portuguese_1, part_sig_portuguese_1)

    part_main_portuguese_2 <- if (est_wag < 0) {
      paste0(
        "O Índice de Concentração segundo Wagstaff é de ", round(est_wag*100, 2),
        "%, ", conf_level*100, "%IC [", round(ci_wag_low*100, 2), "%; ", round(ci_wag_up*100, 2), "%]. ",
        "Este valor negativo indica concentração do indicador nos grupos mais desfavorecidos. ",
        "Em termos relativos, a magnitude equivale a ", abs(round(est_wag*100, 2)),
        "% da desigualdade máxima possível, dada a média do indicador."
      )
    } else if (est_wag > 0) {
      paste0(
        "O Índice de Concentração segundo Wagstaff é de ", round(est_wag*100, 2),
        "%, ", conf_level*100, "%IC [", round(ci_wag_low*100, 2), "%; ", round(ci_wag_up*100, 2), "%]. ",
        "Este valor positivo indica concentração do indicador nos grupos mais favorecidos. ",
        "Em termos relativos, a magnitude equivale a ", abs(round(est_wag*100, 2)),
        "% da desigualdade máxima possível, dada a média do indicador."
      )
    } else {
      paste0(
        "O Índice de Concentração segundo Wagstaff é de ", round(est_wag*100, 2),
        "%, ", conf_level*100, "%IC [", round(ci_wag_low*100, 2), "%; ", round(ci_wag_up*100, 2), "%]. ",
        "O valor nulo é consistente com ausência de gradiente socioeconômico. ",
        "Em termos relativos, a magnitude equivale a 0% da desigualdade máxima possível, dada a média do indicador."
      )
    }

    part_sig_portuguese_2 <- if (ci_wag_low <= 0 & ci_wag_up >= 0) {
      " No entanto, a desigualdade estimada não é estatisticamente significativa."
    } else {
      " Além disso, a desigualdade estimada é estatisticamente significativa."
    }

    interpretation_portuguese_2 <- paste0(part_main_portuguese_2, part_sig_portuguese_2)

    part_main_portuguese_3 <- if (est_err < 0) {
      paste0(
        "O Índice de Concentração de Erreygers é de ", round(est_err*100, 2),
        "pp, ", conf_level*100, "%IC [", round(ci_err_low*100, 2), "pp; ", round(ci_err_up*100, 2), "pp]. ",
        "Este valor negativo implica que, ao longo do gradiente socioeconômico, o indicador de saúde é em média ",
        abs(round(est_err*100, 2)), " pontos percentuais mais alto nos estratos socioeconômicos inferiores do que nos superiores, ",
        "refletindo uma concentração absoluta pró-pobre do indicador."
      )
    } else if (est_err > 0) {
      paste0(
        "O Índice de Concentração de Erreygers é de ", round(est_err*100, 2),
        "pp, ", conf_level*100, "%IC [", round(ci_err_low*100, 2), "pp; ", round(ci_err_up*100, 2), "pp]. ",
        "Este valor positivo implica que, ao longo do gradiente socioeconômico, o indicador de saúde é em média ",
        abs(round(est_err*100, 2)), " pontos percentuais mais alto nos estratos socioeconômicos superiores do que nos inferiores, ",
        "refletindo uma concentração absoluta pró-rico do indicador."
      )
    } else {
      paste0(
        "O Índice de Concentração de Erreygers é de ", round(est_err*100, 2),
        "pp, ", conf_level*100, "%IC [", round(ci_err_low*100, 2), "pp; ", round(ci_err_up*100, 2), "pp]. ",
        "O valor próximo de zero indica ausência de desigualdades absolutas relevantes no nível do indicador entre os extremos do gradiente socioeconômico."
      )
    }

    part_sig_portuguese_3 <- if (ci_err_low <= 0 & ci_err_up >= 0) {
      " No entanto, a desigualdade estimada não é estatisticamente significativa."
    } else {
      " Além disso, a desigualdade estimada é estatisticamente significativa."
    }

    interpretation_portuguese_3 <- paste0(part_main_portuguese_3, part_sig_portuguese_3)

    note_portuguese <- if (rci_est < 0) {
      "Nota: A concentração máxima possível ocorreria se todo o Indicador de Saúde estivesse concentrado no extremo dos mais desfavorecidos."
    } else {
      "Nota: A concentração máxima possível ocorreria se todo o Indicador de Saúde estivesse concentrado no extremo dos mais favorecidos."
    }

    interpretation_portuguese <- c(
      interpretation_portuguese_1,
      interpretation_portuguese_2,
      interpretation_portuguese_3
    )
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
    summary = results_tbl,
    interpretation = interpretation,
    note = note,
    global_health_mean = global_health_mean,
    data = df
  ))
}

