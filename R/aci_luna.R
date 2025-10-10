#' Absolute Concentration Index (two input modes and two methods)
#'
#' Computes the Absolute Concentration Index to quantify socioeconomic-related inequalities in a health indicator.
#' The function operates in two modes: \strong{(Mode A)} using a numerator and denominator (\code{health_numerator} and \code{health_denominator}),
#' which is the preferred and most accurate method as it relies directly on event counts; and \strong{(Mode B)} using an aggregated health indicator
#' (\code{health_indicator}) and population weights (\code{population_weights}), which is used when disaggregated data are not available.
#' The ACI is estimated using either the Fuller method (area under the symmetric Lorenz curve)
#' or the Wagstaff & Kakwani method (regression approach). The function also provides automatic interpretations of the results in
#' four languages: Spanish, English, French, and Portuguese.
#'
#' @param health_indicator_type Character string: either \code{"proportion"} for bounded indicators in \eqn{[0,1]}, or \code{"rate"} if already scaled (e.g., deaths per 100,000); required in all cases.
#' @param health_indicator Numeric vector with the health indicator values (e.g., \code{c(0.12, 0.34)} for proportions or \code{c(120, 345)} for rates); aggregated value per unit; used only in Mode B.
#' @param population_weights Numeric vector > 0 representing population size or weight per unit (e.g., \code{c(1000, 5000, 1200)}); used only in Mode B.
#' @param health_numerator Numeric vector with event counts (e.g., \code{c(5, 12, 8)} for deaths); used only in Mode A and preferred for robust estimation.
#' @param health_denominator Numeric vector > 0 associated with \code{health_numerator} (e.g., \code{c(500, 1000, 800)} representing population at risk); used only in Mode A.
#' @param equity_stratifier Numeric vector representing the social stratifier (e.g., poverty rate, education index); required in all cases.
#' @param higher_ineq_is_favorable Logical; \code{TRUE} if higher values of \code{equity_stratifier} indicate better social conditions (e.g., Gross Domestic Product per capita), or \code{FALSE} if they indicate worse conditions (e.g., Unsatisfied Basic Needs).
#' @param rate_scaling_factor Numeric scalar > 0 applied only when \code{health_indicator_type = "rate"}; commonly \code{1000} or \code{100000} (e.g., maternal mortality ratio per 100000 live births).
#' @param calculation_method Character string specifying the method used to compute the ACI. Allowed values are \code{"Fuller"} (smoothed Lorenz curve approach) and \code{"WK"} (Wagstaff & Kakwani regression approach). Default is \code{"WK"}.
#' @param smoothing_factor Numeric scalar (e.g., \code{0.5}) controlling the smoothness of the curve. Used only when \code{calculation_method = "WK"}; ignored otherwise.
#' @param integration_method Character string specifying the optimization algorithm for parameter estimation when \code{calculation_method = "Fuller"}. Allowed values include \code{"BFGS"}, \code{"CG"}, \code{"Nelder-Mead"}, \code{"L-BFGS-B"}, \code{"nlm"}, \code{"nlminb"}, \code{"spg"}, \code{"ucminf"}, \code{"Rcgmin"}, \code{"Rvmmin"}, \code{"newuoa"}, \code{"bobyqa"}, \code{"nmkb"}, and \code{"hjkb"}. Default is \code{"BFGS"}.
#' @param conf_level Numeric scalar between 0 and 1 (e.g., \code{0.95}) specifying the confidence level for the ACI estimate. Default is \code{0.95}.
#' @param language_interpretation Character string indicating the desired output language for interpretation and note; allowed values are \code{"en"}, \code{"es"}, \code{"fr"}, \code{"pt"} or full names (\code{"english"}, \code{"spanish"}, \code{"french"}, \code{"portuguese"}); default is \code{"en"}.
#'
#' @return A \code{list} with the following components:
#' \item{summary_table}{A \code{tibble} with the Absolute Concentration Index estimate (\code{estimate}), its confidence interval bounds (\code{ci_lower}, \code{ci_upper}), and the method used (\code{method}).}
#' \item{interpretation}{Character string with the automatic interpretation of the result in the selected language.}
#' \item{note}{Character string with the explanatory note in the selected language, based on the type of health indicator.}
#' \item{global_health_mean}{Overall mean of the health indicator: weighted average (Mode B) or global rate (Mode A).}
#' \item{data}{Final \code{tibble} used for the calculation, including processed variables and the estimated Lorenz curve.}
#'
#' @details
#'
#' \subsection{Estimation methods}{
#' The Absolute Concentration Index (ACI) can be estimated using two methods:
#'
#' \describe{
#'   \item{\strong{Wagstaff & Kakwani (WK) method}}{
#'   This method estimates the ACI directly from the observed data without assuming any specific functional form for the concentration curve. The ACI is computed using the formula:
#'
#'   \deqn{\mathrm{ACI} \;=\; \sum_j p_j \cdot (2X_j - 1) \cdot y_j}
#'
#'   where \eqn{p_j} is the population share of subgroup \eqn{j}; \eqn{y_j} is the value of the health indicator in subgroup \eqn{j}; and \eqn{X_j} is the relative rank of subgroup \eqn{j}, calculated as \eqn{X_j = \sum_{i=1}^{j} p_i - 0.5\,p_j}.
#'   }
#'
#'   \item{\strong{Fuller method}}{
#'   This method estimates the ACI by fitting a symmetric Lorenz curve (Murray, 1996). The Lorenz function is:
#'
#'   \deqn{y = \frac{1}{\left( \frac{1}{e^{k - 1}} - 1 \right)} \left( e^{\frac{x}{k - x}} - 1 \right)}
#'
#'   where \eqn{x} is the cumulative population share; \eqn{y} is the cumulative health indicator share; and \eqn{k} is the shape parameter estimated by minimizing the squared error between observed and fitted values.
#'
#'   The ACI corresponds to the area between the fitted curve and the line of equality, scaled by the mean of the health indicator.
#'   }
#' }
#' }
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
#' aci_luna(
#'   health_indicator_type    = "rate",
#'   unit_analysis            = data1$zone,
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
#' aci_luna(
#'   health_indicator_type    = "proportion",
#'   unit_analysis            = data2$zone,
#'   health_numerator         = data2$skilled_births,
#'   health_denominator       = data2$total_births,
#'   equity_stratifier        = data2$ubn_percent,
#'   higher_ineq_is_favorable = FALSE,
#'   calculation_method       = "Fuller",
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
#' aci_luna(
#'   health_indicator_type    = "rate",
#'   unit_analysis            = data1$zone,
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
#' aci_luna(
#'   health_indicator_type    = "proportion",
#'   unit_analysis            = data2$zone,
#'   health_indicator         = data2$skilled_births_prop,
#'   equity_stratifier        = data2$ubn_percent,
#'   population_weights       = data2$population,
#'   higher_ineq_is_favorable = FALSE,
#'   calculation_method       = "Fuller",
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
#' @export

aci_luna <- function(
    health_indicator_type = NULL,
    unit_analysis = NULL,
    health_indicator = NULL,
    population_weights = NULL,
    health_numerator = NULL,
    health_denominator = NULL,
    equity_stratifier = NULL,
    higher_ineq_is_favorable = NULL,
    rate_scaling_factor = NULL,
    calculation_method = NULL,
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

  # Check Unit analysis
  if (is.null(unit_analysis)) {
    stop("You must provide `unit_analysis`.")
  }
  if (any(is.na(unit_analysis))) stop("Missing values detected in `unit_analysis`.")

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
      unit_analysis = as.character(unit_analysis),
      indicator    = as.numeric(indicator_vec),
      equity       = as.numeric(equity_stratifier),
      numerator    = as.numeric(round(health_numerator)),
      denominator  = as.numeric(round(health_denominator)),
      pop_use      = as.numeric(round(health_denominator))   # population to use = denominator when available
    )
  } else {
    df_master <- tibble::tibble(
      unit_analysis = as.character(unit_analysis),
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
    rci_est <- 2 * sum((grid_x - grid_y) * 0.01) * global_health_mean

    df <- df %>% mutate(func_lorenz = func_lorenz(cwpop, k_hat),
                        q = cwi / global_health_mean,
                        a = (indicator / global_health_mean) * (2 * social_position - 1 - rci_est) + 2 - lag(q, 1, default = 0) - q,
                        fa2 = wpop * a^2)

    n <- nrow(df)
    rci_se <- sqrt(abs((1 / n) * (sum(df$fa2, na.rm = TRUE) - (1 + rci_est)^2))) / sqrt(n)
    z <- qnorm(1 - (1 - conf_level) / 2)
    rci_low <- rci_est - z * rci_se *global_health_mean
    rci_up  <- rci_est + z * rci_se *global_health_mean

    results_tbl <- tibble::tibble(metric = "Absolute concentration index",
                                  value = rci_est,
                                  ci_lower = rci_low,
                                  ci_upper = rci_up,
                                  method = "Fuller")

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


    fit_cobs <- cobs(nodos$cwpop, nodos$cwi, constraint = "increase", nknots = smoothing_factor)
    pred <- predict(fit_cobs, z = df$cwpop)$fit
    df$func_lorenz <- pmin(pmax(pred, 0), 1)

    # fit_smooth <- stats::smooth.spline(x = nodos$cwpop, y = nodos$cwi, spar = smoothing_factor)
    # df$func_lorenz <- stats::predict(fit_smooth, x = df$cwpop)$y
    # df$func_lorenz <- pmin(pmax(df$func_lorenz, 0), 1)
    # iso <- stats::isoreg(df$cwpop, df$func_lorenz)
    # df$func_lorenz <- pmin(pmax(approx(iso$x, iso$yf, xout = df$cwpop)$y, 0), 1)

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

    rci_est <- as.numeric(coef(mod)[["rhs"]]) *global_health_mean
    ci_mod <- tryCatch(confint.default(mod, parm = "rhs", level = conf_level), error = function(e) c(NA_real_, NA_real_))

    rci_low <- ci_mod[1] *global_health_mean
    rci_up <- ci_mod[2] *global_health_mean

    results_tbl <- tibble::tibble(metric = "Absolute concentration index",
                                  value = rci_est,
                                  ci_lower = rci_low,
                                  ci_upper = rci_up,
                                  method = "Wagstaff & Kakwani")
  }

  df <- df %>%
    select(unit_analysis = unit_analysis,
           health_indicator = indicator,
           equity_stratifier = equity,
           cumulative_population = cwpop,
           social_position,
           cumulative_event_health  = cwi,
           func_lorenz
    )

  # ────────────────
  # Interpretations
  # ────────────────

  factor_legible <- format(rate_scaling_factor, scientific = FALSE, big.mark = ",")


  part_main_english <- if (rci_est < 0 & health_indicator_type == "rate") {
    paste0(
      "The Absolute Concentration Index is ", round(rci_est, 2),
      ", ", conf_level*100, "%CI [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "This means that the health indicator is concentrated among the most disadvantaged. ",
      "In absolute terms, ", abs(round(rci_est, 2)), " events per ", factor_legible,
      " people in the base population are concentrated in the most disadvantaged groups, compared to a scenario with no socioeconomic differences."
    )
  } else if (rci_est > 0 & health_indicator_type == "rate") {
    paste0(
      "The Absolute Concentration Index is ", round(rci_est, 2),
      ", ", conf_level*100, "%CI [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "This means that the health indicator is concentrated among the most advantaged. ",
      "In absolute terms, ", abs(round(rci_est, 2)), " events per ", factor_legible,
      " people in the base population are concentrated in the most advantaged groups, compared to a scenario with no socioeconomic differences."
    )
  } else if (rci_est == 0 & health_indicator_type == "rate") {
    paste0(
      "The Absolute Concentration Index is ", round(rci_est*100, 2),
      "%, ", conf_level*100, "%CI [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
      "This indicates that, in absolute terms, the health indicator is equally distributed across socioeconomic groups, with no concentration in either the most disadvantaged or the most advantaged."
    )
  } else if (rci_est < 0 & health_indicator_type == "proportion") {
    paste0(
      "The Absolute Concentration Index is ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%CI [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "This means that the health indicator is concentrated among the most disadvantaged. ",
      "In absolute terms, ", abs(round(rci_est*100, 2)), " percentage points of the health indicator are concentrated in the most disadvantaged groups, compared to a scenario with no socioeconomic differences."
    )
  } else if (rci_est > 0 & health_indicator_type == "proportion") {
    paste0(
      "The Absolute Concentration Index is ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%CI [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "This means that the health indicator is concentrated among the most advantaged. ",
      "In absolute terms, ", abs(round(rci_est*100, 2)), " percentage points of the health indicator are concentrated in the most advantaged groups, compared to a scenario with no socioeconomic differences."
    )
  } else if (rci_est == 0 & health_indicator_type == "proportion") {
    paste0(
      "The Absolute Concentration Index is ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%CI [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "This indicates that, in absolute terms, the health indicator is equally distributed across socioeconomic groups, with no concentration in either the most disadvantaged or the most advantaged."
    )
  }

  part_sig_english <- if (rci_low <= 0 & rci_up >= 0) {
    paste0(" However, this concentration is not statistically significant.")
  } else {
    paste0(" Moreover, this concentration is statistically significant.")
  }

  interpretation_english <- paste0(part_main_english, part_sig_english)

  note_english <- if (health_indicator_type == "rate") {
    paste0(
      "Note: Events correspond to the numerator of the health indicator, and the base population to its denominator. ",
      "For example, in the case of the maternal mortality ratio (MMR), events refer to maternal deaths and the base population refers to live births."
    )
  } else {
    paste0(
      "Note: The health indicator value corresponds to the evaluated proportion. For example: Proportion of births attended by skilled health personnel."
    )
  }


  part_main_spanish <- if (rci_est < 0 & health_indicator_type == "rate") {
    paste0(
      "El Índice de Concentración Absoluto es ", round(rci_est, 2),
      ", ", conf_level*100, "%IC [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "Esto significa que el indicador de salud se concentra en los más desfavorecidos. ",
      "En términos absolutos, ", abs(round(rci_est, 2)), " eventos por cada ", factor_legible,
      " personas de la población base se concentran en los grupos más desfavorecidos, en comparación con un escenario sin diferencias socioeconómicas."
    )
  } else if (rci_est > 0 & health_indicator_type == "rate") {
    paste0(
      "El Índice de Concentración Absoluto es ", round(rci_est, 2),
      ", ", conf_level*100, "%IC [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "Esto significa que el indicador de salud se concentra en los más favorecidos. ",
      "En términos absolutos, ", abs(round(rci_est, 2)), " eventos por cada ", factor_legible,
      " personas de la población base se concentran en los grupos más favorecidos, en comparación con un escenario sin diferencias socioeconómicas."
    )
  } else if (rci_est == 0 & health_indicator_type == "rate") {
    paste0(
      "El índice de Concentración Absoluto es de ", round(rci_est*100, 2),
      "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
      "Esto indica que, en términos absolutos, el indicador de salud está distribuido de manera equitativa entre los grupos socioeconómicos, sin concentrarse ni en los más desfavorecidos ni en los menos desfavorecidos."
    )
  } else if (rci_est < 0 & health_indicator_type == "proportion") {
    paste0(
      "El índice de Concentración Absoluto es ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "Esto significa que el indicador de salud se concentra en los más desfavorecidos. ",
      "En términos absolutos, ", abs(round(rci_est*100, 2)), " puntos porcentuales del indicador de salud están concentrados en los grupos más desfavorecidos, en comparación con un escenario sin diferencias socioeconómicas."
    )
  } else if (rci_est > 0 & health_indicator_type == "proportion") {
    paste0(
      "El índice de Concentración Absoluto es ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "Esto significa que el indicador de salud se concentra en los más favorecidos. ",
      "En términos absolutos, ", abs(round(rci_est*100, 2)), " puntos porcentuales del indicador de salud están concentrados en los grupos más favorecidos, en comparación con un escenario sin diferencias socioeconómicas."
    )
  } else if (rci_est == 0 & health_indicator_type == "proportion") {
    paste0(
      "El índice de Concentración Absoluto es de ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "Esto indica que, en términos absolutos, el indicador de salud está distribuido de manera equitativa entre los grupos socioeconómicos, sin concentrarse ni en los más desfavorecidos ni en los menos desfavorecidos."
    )
  }

  part_sig_spanish <- if (rci_low <= 0 & rci_up >= 0) {
    paste0(" Sin embargo, esta concentración no es estadísticamente significativa.")
  } else {
    paste0(" Además, esta concentración es estadísticamente significativa.")
  }

  interpretation_spanish <- paste0(part_main_spanish, part_sig_spanish)

  note_spanish <- if (health_indicator_type == "rate") {
    paste0(
      "Nota: Los eventos corresponden al numerador del Indicador de Salud y la población base al denominador del Indicador de Salud. ",
      "Por ejemplo, en el caso de la RMM (Razón de Mortalidad Materna), los eventos serían las muertes maternas y la población base serían los nacidos vivos."
    )
  } else {
    paste0(
      "Nota: El valor del Indicador de Salud corresponde a la proporción evaluada. Por ejemplo: Proporción de partos atendidos por personal calificado."
    )
  }


  part_main_french <- if (rci_est < 0 & health_indicator_type == "rate") {
    paste0(
      "L'indice de concentration absolue est de ", round(rci_est, 2),
      ", ", conf_level*100, "%IC [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "Cela signifie que l'indicateur de santé est concentré chez les plus défavorisés. ",
      "En termes absolus, ", abs(round(rci_est, 2)), " événements pour ", factor_legible,
      " personnes de la population de base sont concentrés dans les groupes les plus défavorisés, par rapport à un scénario sans différences socioéconomiques."
    )
  } else if (rci_est > 0 & health_indicator_type == "rate") {
    paste0(
      "L'indice de concentration absolue est de ", round(rci_est, 2),
      ", ", conf_level*100, "%IC [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "Cela signifie que l'indicateur de santé est concentré chez les plus favorisés. ",
      "En termes absolus, ", abs(round(rci_est, 2)), " événements pour ", factor_legible,
      " personnes de la population de base sont concentrés dans les groupes les plus favorisés, par rapport à un scénario sans différences socioéconomiques."
    )
  } else if (rci_est == 0 & health_indicator_type == "rate") {
    paste0(
      "L'indice de concentration absolue est de ", round(rci_est*100, 2),
      "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
      "Cela indique qu'en termes absolus, l'indicateur de santé est réparti équitablement entre les groupes socioéconomiques, sans concentration chez les plus ou les moins défavorisés."
    )
  } else if (rci_est < 0 & health_indicator_type == "proportion") {
    paste0(
      "L'indice de concentration absolue est de ", round(rci_est*100, 2),
      " pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), " pp; ", round(rci_up*100, 2), " pp]. ",
      "Cela signifie que l'indicateur de santé est concentré chez les plus défavorisés. ",
      "En termes absolus, ", abs(round(rci_est*100, 2)), " points de pourcentage de l'indicateur de santé sont concentrés dans les groupes les plus défavorisés, par rapport à un scénario sans différences socioéconomiques."
    )
  } else if (rci_est > 0 & health_indicator_type == "proportion") {
    paste0(
      "L'indice de concentration absolue est de ", round(rci_est*100, 2),
      " pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), " pp; ", round(rci_up*100, 2), " pp]. ",
      "Cela signifie que l'indicateur de santé est concentré chez les plus favorisés. ",
      "En termes absolus, ", abs(round(rci_est*100, 2)), " points de pourcentage de l'indicateur de santé sont concentrés dans les groupes les plus favorisés, par rapport à un scénario sans différences socioéconomiques."
    )
  } else if (rci_est == 0 & health_indicator_type == "proportion") {
    paste0(
      "L'indice de concentration absolue est de ", round(rci_est*100, 2),
      " pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), " pp; ", round(rci_up*100, 2), " pp]. ",
      "Cela indique qu'en termes absolus, l'indicateur de santé est réparti équitablement entre les groupes socioéconomiques, sans concentration chez les plus ou les moins défavorisés."
    )
  }

  part_sig_french <- if (rci_low <= 0 & rci_up >= 0) {
    paste0(" Cependant, cette concentration n'est pas statistiquement significative.")
  } else {
    paste0(" De plus, cette concentration est statistiquement significative.")
  }

  interpretation_french <- paste0(part_main_french, part_sig_french)

  note_french <- if (health_indicator_type == "rate") {
    paste0(
      "Remarque : Les événements correspondent au numérateur de l’indicateur de santé et la population de base au dénominateur. ",
      "Par exemple, dans le cas du Taux de mortalité maternelle, les événements seraient les décès maternels et la population de base serait les naissances vivantes."
    )
  } else {
    paste0(
      "Remarque : La valeur de l’indicateur de santé correspond à la proportion évaluée. Par exemple : proportion d’accouchements assistés par du personnel qualifié."
    )
  }


  part_main_portuguese <- if (rci_est < 0 & health_indicator_type == "rate") {
    paste0(
      "O Índice de Concentração Absoluta é ", round(rci_est, 2),
      ", ", conf_level*100, "%IC [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "Isso significa que o indicador de saúde está concentrado entre os mais desfavorecidos. ",
      "Em termos absolutos, ", abs(round(rci_est, 2)), " eventos por cada ", factor_legible,
      " pessoas da população base estão concentrados entre os grupos mais desfavorecidos, em comparação com um cenário sem diferenças socioeconômicas."
    )
  } else if (rci_est > 0 & health_indicator_type == "rate") {
    paste0(
      "O Índice de Concentração Absoluta é ", round(rci_est, 2),
      ", ", conf_level*100, "%IC [", round(rci_low, 2), "; ", round(rci_up, 2), "]. ",
      "Isso significa que o indicador de saúde está concentrado entre os mais favorecidos. ",
      "Em termos absolutos, ", abs(round(rci_est, 2)), " eventos por cada ", factor_legible,
      " pessoas da população base estão concentrados entre os grupos mais favorecidos, em comparação com um cenário sem diferenças socioeconômicas."
    )
  } else if (rci_est == 0 & health_indicator_type == "rate") {
    paste0(
      "O Índice de Concentração Absoluta é ", round(rci_est*100, 2),
      "%, ", conf_level*100, "%IC [", round(rci_low*100, 2), "%; ", round(rci_up*100, 2), "%]. ",
      "Isso indica que, em termos absolutos, o indicador de saúde está distribuído de maneira equitativa entre os grupos socioeconômicos, sem se concentrar nem nos mais desfavorecidos nem nos mais favorecidos."
    )
  } else if (rci_est < 0 & health_indicator_type == "proportion") {
    paste0(
      "O Índice de Concentração Absoluta é ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "Isso significa que o indicador de saúde está concentrado entre os mais desfavorecidos. ",
      "Em termos absolutos, ", abs(round(rci_est*100, 2)), " pontos percentuais do indicador de saúde estão concentrados entre os grupos mais desfavorecidos, em comparação com um cenário sem diferenças socioeconômicas."
    )
  } else if (rci_est > 0 & health_indicator_type == "proportion") {
    paste0(
      "O Índice de Concentração Absoluta é ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "Isso significa que o indicador de saúde está concentrado entre os mais favorecidos. ",
      "Em termos absolutos, ", abs(round(rci_est*100, 2)), " pontos percentuais do indicador de saúde estão concentrados entre os grupos mais favorecidos, em comparação com um cenário sem diferenças socioeconômicas."
    )
  } else if (rci_est == 0 & health_indicator_type == "proportion") {
    paste0(
      "O Índice de Concentração Absoluta é ", round(rci_est*100, 2),
      "pp, ", conf_level*100, "%IC [", round(rci_low*100, 2), "pp; ", round(rci_up*100, 2), "pp]. ",
      "Isso indica que, em termos absolutos, o indicador de saúde está distribuído de maneira equitativa entre os grupos socioeconômicos, sem se concentrar nem nos mais desfavorecidos nem nos mais favorecidos."
    )
  }

  part_sig_portuguese <- if (rci_low <= 0 & rci_up >= 0) {
    paste0(" No entanto, essa concentração não é estatisticamente significativa.")
  } else {
    paste0(" Além disso, essa concentração é estatisticamente significativa.")
  }

  interpretation_portuguese <- paste0(part_main_portuguese, part_sig_portuguese)

  note_portuguese <- if (health_indicator_type == "rate") {
    paste0(
      "Nota: Os eventos correspondem ao numerador do Indicador de Saúde e a população base ao denominador do Indicador de Saúde. ",
      "Por exemplo, no caso da Razão de Mortalidade Materna (RMM), os eventos seriam as mortes maternas e a população base seriam os nascidos vivos."
    )
  } else {
    paste0(
      "Nota: O valor do Indicador de Saúde corresponde à proporção avaliada. Por exemplo: Proporção de partos assistidos por pessoal qualificado."
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
    summary_table = results_tbl,
    interpretation = interpretation,
    note = note,
    global_health_mean = global_health_mean,
    data = df
  ))
}

