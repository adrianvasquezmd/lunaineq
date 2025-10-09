#' Relative index of inequality (two input modes)
#'
#' Computes the Relative Index of Inequality (RII) between the extreme social groups defined based on a
#' numeric equity stratifier. The function operates in two modes: \strong{(Mode A)} using a numerator and
#' denominator (\code{health_numerator} and \code{health_denominator}), which is the preferred
#' and most accurate method, as it relies directly on event counts;
#' and \strong{(Mode B)} using an aggregated health indicator (\code{health_indicator})
#' and population weights (\code{population_weights}), which approximates the index based on weighted group means.
#' The function also provides automatic interpretations of the results in
#' four languages: Spanish, English, French, and Portuguese.
#'
#' In both cases, the social stratifier is transformed into a relative social position (ridit)
#' ranging from 0 (most disadvantaged) to 1 (most advantaged), ensuring that \emph{position 0}
#' corresponds to the most disadvantaged and \emph{position 1} to the most advantaged
#' (as defined by \code{higher_ineq_is_favorable}). The RII is calculated as the relative difference in the health indicator \eqn{y}
#' —predicted from the regression model at the two extremes— between the most disadvantaged and the most advantaged positions: \eqn{\mathrm{RII} = \hat{y}_{\text{advantaged}} / \hat{y}_{\text{disadvantaged}}}.
#'
#' @param health_indicator_type Character string: either \code{"proportion"} for bounded indicators in \eqn{[0,1]}, or \code{"rate"} if already scaled (e.g., deaths per 100,000); required in all cases.
#' @param health_indicator Numeric vector with the health indicator values (e.g., \code{c(0.12, 0.34)} for proportions or \code{c(120, 345)} for rates); aggregated value per unit; used only in Mode B.
#' @param population_weights Numeric vector > 0 representing population size or weight per unit (e.g., \code{c(1000, 5000, 1200)}); used only in Mode B.
#' @param health_numerator Numeric vector with event counts (e.g., \code{c(5, 12, 8)} for deaths); used only in Mode A and preferred for robust estimation.
#' @param health_denominator Numeric vector > 0 associated with \code{health_numerator} (e.g., \code{c(500, 1000, 800)} representing population at risk); used only in Mode A.
#' @param equity_stratifier Numeric vector representing the social stratifier (e.g., poverty rate, education index); required in all cases.
#' @param higher_ineq_is_favorable Logical; \code{TRUE} if higher values of \code{equity_stratifier} indicate better social conditions (e.g., Gross Domestic Product per capita), or \code{FALSE} if they indicate worse conditions (e.g., Unsatisfied Basic Needs).
#' @param rate_scaling_factor Numeric scalar > 0 applied only when \code{health_indicator_type = "rate"}; commonly \code{1000} or \code{100000} (e.g., maternal mortality ratio per 100000 live births).
#' @param models Character vector specifying the regression models to evaluate; for proportions: default is \code{c("binomial", "linear logit(y)")} (Mode A) and \code{c("quasibinomial", "linear logit(y)")} (Mode B), allowed models include \code{"linear"}, \code{"linear ln(y)"}, \code{"linear logit(y)"}, \code{"binomial"}, \code{"quasibinomial"}; for rates: default is \code{c("poisson", "negative binomial")} (Mode A) and \code{c("linear ln(y)")} (Mode B), allowed models include \code{"linear"}, \code{"linear ln(y)"}, \code{"poisson"}, \code{"negative binomial"}.
#' @param model_evaluation_metric Character string defining the metric for best model selection; allowed values: \code{"MAE"}, \code{"MSE"}, \code{"RMSE"}, \code{"MAPE"}, \code{"AIC"}; default is \code{"MSE"}.
#' @param conf_level Numeric scalar between 0 and 1 (e.g., \code{0.95}) specifying the confidence level for inequality estimates; default is \code{0.95}.
#' @param language_interpretation Character string indicating the desired output language for interpretation and note; allowed values are \code{"en"}, \code{"es"}, \code{"fr"}, \code{"pt"} or full names (\code{"english"}, \code{"spanish"}, \code{"french"}, \code{"portuguese"}); default is \code{"en"}.
#'
#' @return A \code{list} with components:
#' \item{summary_table}{A \code{tibble} with \code{inequality_metric}, \code{selected_model}, \code{estimate}, \code{ci_lower}, \code{ci_upper}.}
#' \item{interpretation}{Character: automatic interpretation in the selected language.}
#' \item{note}{Character: explanatory note/alerts in the selected language.}
#' \item{global_health_mean}{Numeric: overall mean of the health indicator (weighted or derived from num/den).}
#' \item{data}{\code{tibble} used for calculations. Contains at least: \code{health_indicator}, \code{equity_stratifier}, \code{pop_use}, \code{social_position} (ridit), fitted prediction columns for candidate models, and \code{final_model_pred}.}
#' \item{data_with_extremes}{A \code{data.frame} / \code{tibble} that appends the extreme social positions (0 and 1) and their predicted values (from emmeans). Useful for plotting extremes.}
#'
#' @details
#' \subsection{A) Using numerator and denominator (preferred)}{
#' This method calculates the health indicator directly from event counts and denominators, ensuring statistical robustness.
#'
#' \itemize{
#' \item \strong{Proportions}: The RII is computed using binomial, linear or logit-linear regression models.
#' \item \strong{Rates}: The RII is estimated from Poisson, negative binomial, or log-linear regression models, scaled by \code{rate_scaling_factor}.
#' \item \emph{This mode does not require population weights}: estimates are based on observed event counts.
#' }
#' }
#'
#' \subsection{B) Using aggregated indicator and population weights (approximate)}{
#' This method approximates the RII based on weighted means.
#'
#' \itemize{
#' \item \strong{Proportions}: Quasibinomial or logit-linear models are used.
#' \item \strong{Rates}: Linear or log-linear regression models with weights are applied.
#' \item \emph{Warning}: Results may be biased when using small samples or aggregated data. Prefer Mode A when possible.
#' }
#' }
#'
#' @note Mode B is designed for cases where only summary indicators and population weights are available. However, Mode A should be used whenever individual counts are accessible, as it yields more accurate and robust estimates.
#'
#' @examples
#' # ——— Mode A (preferred: numerator and denominator) ———
#'
#' # Example with a rate indicator (e.g., maternal mortality ratio)
#' data(data1)
#' rii_luna(
#'   health_indicator_type    = "rate",
#'   unit_analysis            = data1$zone,
#'   health_numerator         = data1$maternal_deaths,
#'   health_denominator       = data1$live_births,
#'   equity_stratifier        = data1$ubn_percent,
#'   higher_ineq_is_favorable = FALSE,
#'   rate_scaling_factor      = 100000,
#'   models = c('poisson', 'negative binomial'),
#'   model_evaluation_metric = 'MSE',
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' # Example with a proportion indicator (e.g., skilled birth attendance)
#' data(data2)
#' rii_luna(
#'   unit_analysis            = data2$zone,
#'   health_numerator         = data2$skilled_births,
#'   health_denominator       = data2$total_births,
#'   equity_stratifier        = data2$ubn_percent,
#'   higher_ineq_is_favorable = FALSE,
#'   health_indicator_type    = "proportion",
#'   models = c('binomial', 'linear logit(y)'),
#'   model_evaluation_metric = 'MSE',
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' # ——— Mode B (aggregated indicator + population as weight) ———
#'
#' # Example with rate indicator
#' data(data1)
#' rii_luna(
#'   unit_analysis            = data1$zone,
#'   health_indicator         = data1$mmr,
#'   equity_stratifier        = data1$ubn_percent,
#'   population_weights       = data1$population_zone,
#'   higher_ineq_is_favorable = FALSE,
#'   health_indicator_type    = "rate",
#'   models = c('linear ln(y)'),
#'   model_evaluation_metric = 'MSE',
#'   rate_scaling_factor      = 100000,
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' # Example with proportion indicator
#' data(data2)
#' rii_luna(
#'   unit_analysis            = data2$zone,
#'   health_indicator         = data2$skilled_births_prop,
#'   equity_stratifier        = data2$ubn_percent,
#'   population_weights       = data2$population,
#'   higher_ineq_is_favorable = FALSE,
#'   health_indicator_type    = "proportion",
#'   models = c('quasibinomial', 'linear logit(y)'),
#'   model_evaluation_metric = 'MSE',
#'   conf_level               = 0.95,
#'   language_interpretation  = "en"
#' )
#'
#' @seealso \code{emmeans::emmeans}, \code{sandwich::vcovHC}, \code{Metrics::mae}
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

rii_luna <- function(health_indicator_type = NULL,
                unit_analysis = NULL,
                health_indicator = NULL,
                population_weights = NULL,
                health_numerator = NULL,
                health_denominator = NULL,
                equity_stratifier = NULL,
                higher_ineq_is_favorable = NULL,
                rate_scaling_factor = NULL,
                models = NULL,
                model_evaluation_metric = NULL,
                conf_level = NULL,
                language_interpretation = NULL) {

  # ───────────────────────────────────────────────
  # Argument validation (input consistency checks)
  # ───────────────────────────────────────────────

  # Check presence health numerator & health denominator
  use_num_den <- !is.null(health_numerator) && !is.null(health_denominator)

  # Check health indicator
  if (is.null(health_indicator_type) || !(health_indicator_type %in% c("proportion", "rate"))) {
    stop("`health_indicator_type` is required and must be either 'proportion' or 'rate'.")
  }

  # Check model evaluation metric
  if (is.null(model_evaluation_metric)) {
    model_evaluation_metric <- "MSE"
    warning("`model_evaluation_metric` not provided: defaulting to 'MSE'.")
  }
  allowed_metrics <- c("MAE", "MSE", "RMSE", "MAPE", "AIC")
  if (!(model_evaluation_metric %in% allowed_metrics)) {
    stop("`model_evaluation_metric` must be one of: ", paste(allowed_metrics, collapse = ", "), ".")
  }

  # Check confidence level
  if (is.null(conf_level)) {
    conf_level <- 0.95
    warning("`conf_level` not provided: defaulting to 0.95.")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single numeric value between 0 and 1 (e.g. 0.95).")
  }

  # Check inequality stratifier
  if (is.null(equity_stratifier)) {
    stop("You must provide `equity_stratifier`.")
  }
  if (!is.numeric(equity_stratifier)) stop("`equity_stratifier` must be numeric.")
  if (any(is.na(equity_stratifier))) stop("Missing values detected in `equity_stratifier`.")

  # Check Unit analysis
  if (is.null(unit_analysis)) {
    stop("You must provide `unit_analysis`.")
  }
  if (any(is.na(unit_analysis))) stop("Missing values detected in `unit_analysis`.")

  # Check higher equity stratifier is favorable
  if (is.null(higher_ineq_is_favorable)) {
    stop("`higher_ineq_is_favorable` is required (TRUE if higher stratifier values = more favorable).")
  }
  if (!is.logical(higher_ineq_is_favorable)) {
    stop("`higher_ineq_is_favorable` must be logical (TRUE/FALSE).")
  }

  # Warn if scaling factor is ignored for proportions
  if (health_indicator_type == "proportion") {
    if (!is.null(rate_scaling_factor)) {
      warning("`rate_scaling_factor` is ignored because `health_indicator_type` is 'proportion'. It only applies to rates.")
    }
    rate_scaling_factor <- 1
  }

  # Check health indicator type
  if (is.null(rate_scaling_factor) || !is.numeric(rate_scaling_factor) ||
      length(rate_scaling_factor) != 1 || rate_scaling_factor <= 0) {
    stop("For `health_indicator_type = 'rate'` you must provide a numeric `rate_scaling_factor` > 0 (e.g., 1000 or 100000).")
  }

  # Mode A: numerator/denominator provided
  if (use_num_den) {

    # Numerator/denominator are numeric and complete
    if (!is.numeric(health_numerator) || !is.numeric(health_denominator)) {
      stop("`health_numerator` and `health_denominator` must be numeric when provided.")
    }
    if (any(is.na(health_numerator)) || any(is.na(health_denominator))) {
      stop("Missing values detected in `health_numerator` or `health_denominator`.")
    }

    # Check health numerator and health denominador
    if (any(health_numerator < 0, na.rm = TRUE)) {
      stop("`health_numerator` contains negative values. All values must be >= 0.")
    }
    if (any(health_denominator <= 0, na.rm = TRUE)) {
      stop("All values in `health_denominator` must be > 0.")
    }
    if (any(health_numerator %% 1 != 0, na.rm = TRUE)) {
      warning("`health_numerator` contains non-integer values. These will be rounded to the nearest integer.")
      health_numerator <- round(health_numerator)
    }
    if (any(health_denominator %% 1 != 0, na.rm = TRUE)) {
      warning("`health_denominator` contains non-integer values. These will be rounded to the nearest integer.")
      health_denominator <- round(health_denominator)
    }

    #Check health indicator calculated from numerator/denominator
    if (health_indicator_type == "proportion") {
      ratio_calc <- health_numerator / health_denominator
      if (any(is.na(ratio_calc)) || any(ratio_calc < 0 | ratio_calc > 1, na.rm = TRUE)) {
        stop("Calculated proportions from numerator/denominator fall outside [0,1]. Check `health_numerator` and `health_denominator`.")
      }
    } else {
      ratio_calc_rate <- (health_numerator / health_denominator) * rate_scaling_factor
      if (any(is.na(ratio_calc_rate)) || any(ratio_calc_rate < 0, na.rm = TRUE)) {
        stop("Calculated rates from numerator/denominator are invalid. Check `health_numerator`, `health_denominator`, and `rate_scaling_factor`.")
      }
    }
    if (!is.null(health_indicator)) {
      warning("`health_numerator` and `health_denominator` provided: `health_indicator` will be ignored and derived from num/den.")
      health_indicator <- NULL
    }

    #Check length of inputs
    if (length(health_numerator) != length(health_denominator)) {
      stop("`health_numerator` and `health_denominator` must have the same length.")
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
    if (is.null(health_indicator)) stop("`health_indicator` is required when numerator/denominator are not provided.")
    if (!is.numeric(health_indicator)) stop("`health_indicator` must be numeric.")
    if (any(is.na(health_indicator))) stop("Missing values detected in `health_indicator`.")

    if (health_indicator_type == "proportion") {
      if (any(health_indicator < 0 | health_indicator > 1, na.rm = TRUE)) {
        stop("For `health_indicator_type = 'proportion'`, `health_indicator` values must be between 0 and 1.")
      }
    } else {
      if (any(health_indicator < 0, na.rm = TRUE)) {
        stop("`health_indicator` (rates) must be >= 0.")
      }
    }

    # Check population weights
    if (is.null(population_weights)) stop("`population_weights` must be provided when numerator/denominator are not.")
    if (!is.numeric(population_weights)) stop("`population_weights` must be numeric.")
    if (any(is.na(population_weights))) stop("Missing values detected in `population_weights`.")
    if (any(population_weights <= 0, na.rm = TRUE)) stop("`population_weights` must contain values strictly greater than 0.")
    if (!all(population_weights %% 1 == 0)) {
      warning("`population_weights` contain non-integer values. CI approximations may be unreliable; prefer numerator/denominator when possible.")
    }

    #Check length of inputs
    if (length(population_weights) != length(health_indicator)) {
      stop("`population_weights` must have the same length as `health_indicator`.")
    }
    if (length(population_weights) != length(equity_stratifier)) {
      stop("`population_weights` must have the same length as `equity_stratifier`.")
    }
  }

  # Model defaults and allowed models depending on scenario
  if (!use_num_den && health_indicator_type == "proportion") {
    default_models <- c("quasibinomial", "linear logit(y)")
    allowed_models <- c("quasibinomial", "linear logit(y)", "linear ln(y)", "linear")
    scenario_text <- "proportion without numerator/denominator"
  } else if (use_num_den && health_indicator_type == "proportion") {
    default_models <- c("binomial", "linear logit(y)")
    allowed_models <- c("binomial", "linear logit(y)", "linear ln(y)", "linear")
    scenario_text <- "proportion with numerator/denominator"
  } else if (!use_num_den && health_indicator_type == "rate") {
    default_models <- c("linear ln(y)")
    allowed_models <- c("linear ln(y)", "linear")
    scenario_text <- "rate without numerator/denominator"
  } else if (use_num_den && health_indicator_type == "rate") {
    default_models <- c("poisson", "negative binomial")
    allowed_models <- c("poisson", "negative binomial", "linear ln(y)", "linear")
    scenario_text <- "rate with numerator/denominator"
  }

  if (is.null(models) || length(models) == 0) {
    models <- default_models
    warning(sprintf("`models` not provided: defaulting to %s for scenario: %s.",
                    paste0("c('", paste(default_models, collapse = "', '"), "')"),
                    scenario_text))
  }

  models <- as.character(models)
  models <- trimws(models)
  models <- unique(models)

  invalid_models <- setdiff(models, allowed_models)
  if (length(invalid_models) > 0) {
    stop(sprintf("Models not allowed in this scenario (%s): %s. Allowed models: %s.",
                 scenario_text,
                 paste(invalid_models, collapse = ", "),
                 paste(allowed_models, collapse = ", ")))
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

  # ────────────────────────────────
  # Base data frame preparation
  # ────────────────────────────────
  if (!use_num_den) {
    df <- tibble::tibble(
      unit_analysis = as.character(unit_analysis),
      health_indicator = as.numeric(health_indicator),
      equity_stratifier = as.numeric(equity_stratifier),
      population_weights = as.numeric(population_weights)
    )
  } else {
    df <- tibble::tibble(
      unit_analysis = as.character(unit_analysis),
      health_indicator = if (health_indicator_type == "proportion") {
        as.numeric(health_numerator) / as.numeric(health_denominator)
      } else {
        (as.numeric(health_numerator) / as.numeric(health_denominator)) * rate_scaling_factor
      },
      equity_stratifier   = as.numeric(equity_stratifier),
      health_numerator    = round(as.numeric(health_numerator)),
      health_denominator  = round(as.numeric(health_denominator))
    )
  }

  # ──────────────────────────────────────
  # Ordering by the inequality variable
  # ──────────────────────────────────────
  df <- df %>%
    mutate(equity_stratifier_order = if (isTRUE(higher_ineq_is_favorable)) equity_stratifier else -equity_stratifier) %>%
    arrange(equity_stratifier_order)

  # ──────────────────────────────────────
  # Compute population used and social position (ridit midpoint)
  # ──────────────────────────────────────
  df <- df %>%
    mutate(
      pop_use = if (!use_num_den) population_weights else health_denominator
    ) %>%
    arrange(equity_stratifier_order) %>%
    mutate(
      cumw  = cumsum(pop_use),
      cumw1 = lag(cumw, default = 0),
      sumw  = sum(pop_use)
    ) %>%
    group_by(equity_stratifier_order) %>%
    mutate(
      cumwr  = max(cumw),
      cumwr1 = min(cumw1),
      social_position = (cumwr1 + 0.5 * (cumwr - cumwr1)) / sumw
    ) %>%
    ungroup() %>%
    select(-cumw, -cumw1, -sumw, -cumwr, -cumwr1)

  # ──────────────────────────────────────────────────────────────────────────────
  # Global weighted means for the health indicator and the equity stratifier
  # ──────────────────────────────────────────────────────────────────────────────
  if (!use_num_den) {
    global_health_mean <- weighted.mean(df$health_indicator, w = df$pop_use, na.rm = TRUE)
  } else {
    total_num <- sum(df$health_numerator, na.rm = TRUE)
    total_den <- sum(df$health_denominator, na.rm = TRUE)

    global_health_mean <- if (health_indicator_type == "proportion") {
      total_num / total_den
    } else if (health_indicator_type == "rate") {
      (total_num / total_den) * rate_scaling_factor
    }
  }

  # ───────────────────────────────────────────────────────────
  # Fit candidate models and compute predictions / fit metrics
  # ───────────────────────────────────────────────────────────
  model_objects <- list()

  if (health_indicator_type == "rate") {
    eps <- 1e-3
  } else { # "proportion"
    eps <- 1e-6
  }

  # Rate models
  if (health_indicator_type == "rate") {

    # Count models (require numerator/denominator)
    if (use_num_den) {
      model_objects$model_poisson <- tryCatch(
        glm(
          health_numerator ~ social_position + offset(log(health_denominator / rate_scaling_factor)),
          data = df,
          family = poisson()
        ),
        error = function(e) stop("Failed to fit Poisson model: ", e$message)
      )

      model_objects$model_negbin <- tryCatch(
        MASS::glm.nb(
          health_numerator ~ social_position + offset(log(health_denominator / rate_scaling_factor)),
          data = df
        ),
        error = function(e) stop("Failed to fit negative binomial model: ", e$message)
      )
    }

    # Linear-type models (weighted by pop_use)
    model_objects$model_linear <- tryCatch(
      glm(
        health_indicator ~ social_position,
        data = df,
        weights = pop_use,
        family = gaussian()
      ),
      error = function(e) stop("Failed to fit linear model: ", e$message)
    )

    model_objects$model_linear_lny <- tryCatch(
      glm(
        log(health_indicator) ~ social_position,
        data = df %>% mutate(health_indicator = ifelse(health_indicator == 0, eps, health_indicator)),
        weights = pop_use,
        family = gaussian()
      ),
      error = function(e) stop("Failed to fit linear (ln(y)) model: ", e$message)
    )

    # initialize prediction columns
    df$pred_model_poisson <- NA_real_
    df$pred_model_negbin  <- NA_real_
    df$pred_model_linear  <- NA_real_
    df$pred_model_linear_lny <- NA_real_

    # rowwise predictions
    if (!is.null(model_objects$model_poisson)) {
      # predict returns expected counts; divide by (denominator / factor) to get rate
      df$pred_model_poisson <- predict(model_objects$model_poisson, type = "response") / (df$health_denominator / rate_scaling_factor)
    }
    if (!is.null(model_objects$model_negbin)) {
      df$pred_model_negbin <- predict(model_objects$model_negbin, type = "response") / (df$health_denominator / rate_scaling_factor)
    }
    df$pred_model_linear <- predict(model_objects$model_linear, type = "response")
    df$pred_model_linear_lny <- exp(predict(model_objects$model_linear_lny, type = "link"))

    # predictions at extremes (social_position = 0 and 1)
    nd0 <- data.frame(social_position = 0, health_denominator = 1, rate_scaling_factor = 1)
    nd1 <- data.frame(social_position = 1, health_denominator = 1, rate_scaling_factor = 1)

    p0_lin <- as.numeric(predict(model_objects$model_linear, newdata = nd0, type = "response"))
    p1_lin <- as.numeric(predict(model_objects$model_linear, newdata = nd1, type = "response"))

    p0_lny <- exp(as.numeric(predict(model_objects$model_linear_lny, newdata = nd0, type = "link")))
    p1_lny <- exp(as.numeric(predict(model_objects$model_linear_lny, newdata = nd1, type = "link")))

    p0_pois <- if (!is.null(model_objects$model_poisson)) as.numeric(predict(model_objects$model_poisson, newdata = nd0, type = "response")) else NA_real_
    p1_pois <- if (!is.null(model_objects$model_poisson)) as.numeric(predict(model_objects$model_poisson, newdata = nd1, type = "response")) else NA_real_

    p0_negb <- if (!is.null(model_objects$model_negbin)) as.numeric(predict(model_objects$model_negbin, newdata = nd0, type = "response")) else NA_real_
    p1_negb <- if (!is.null(model_objects$model_negbin)) as.numeric(predict(model_objects$model_negbin, newdata = nd1, type = "response")) else NA_real_

    # --- Fit metrics table ---
    fit_metrics <- data.frame(
      model = c("linear", "linear ln(y)", if (use_num_den) "poisson" else NA, if (use_num_den) "negative binomial" else NA),
      MAE = c(
        mae(df$health_indicator, df$pred_model_linear),
        mae(df$health_indicator, df$pred_model_linear_lny),
        if (!is.null(model_objects$model_poisson)) mae(df$health_indicator, df$pred_model_poisson) else NA_real_,
        if (!is.null(model_objects$model_negbin)) mae(df$health_indicator, df$pred_model_negbin) else NA_real_
      ),
      MSE = c(
        mse(df$health_indicator, df$pred_model_linear),
        mse(df$health_indicator, df$pred_model_linear_lny),
        if (!is.null(model_objects$model_poisson)) mse(df$health_indicator, df$pred_model_poisson) else NA_real_,
        if (!is.null(model_objects$model_negbin)) mse(df$health_indicator, df$pred_model_negbin) else NA_real_
      ),
      RMSE = c(
        rmse(df$health_indicator, df$pred_model_linear),
        rmse(df$health_indicator, df$pred_model_linear_lny),
        if (!is.null(model_objects$model_poisson)) rmse(df$health_indicator, df$pred_model_poisson) else NA_real_,
        if (!is.null(model_objects$model_negbin)) rmse(df$health_indicator, df$pred_model_negbin) else NA_real_
      ),
      MAPE = c(
        mape(df$health_indicator, df$pred_model_linear),
        mape(df$health_indicator, df$pred_model_linear_lny),
        if (!is.null(model_objects$model_poisson)) mape(df$health_indicator, df$pred_model_poisson) else NA_real_,
        if (!is.null(model_objects$model_negbin)) mape(df$health_indicator, df$pred_model_negbin) else NA_real_
      ),
      AIC = c(
        model_objects$model_linear$aic,
        model_objects$model_linear_lny$aic,
        if (!is.null(model_objects$model_poisson)) model_objects$model_poisson$aic else NA_real_,
        if (!is.null(model_objects$model_negbin)) model_objects$model_negbin$aic else NA_real_
      ),
      infeasible_estimates = c(
        any(df$pred_model_linear < 0, na.rm = TRUE),
        any(df$pred_model_linear_lny < 0, na.rm = TRUE),
        if (!is.null(model_objects$model_poisson)) any(df$pred_model_poisson < 0, na.rm = TRUE) else FALSE,
        if (!is.null(model_objects$model_negbin)) any(df$pred_model_negbin < 0, na.rm = TRUE) else FALSE
      ),
      pred_pos0 = c(p0_lin, p0_lny, p0_pois, p0_negb),
      pred_pos1 = c(p1_lin, p1_lny, p1_pois, p1_negb),
      indicator_has_zero = rep(any(df$health_indicator == 0, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )

    fit_metrics1 <- fit_metrics %>% dplyr::filter(!is.na(model))

    fit_metrics2 <- fit_metrics1 %>%
      mutate(
        pos01_infeasible = (pred_pos0 < 0) | (pred_pos1 < 0),
        model_viable = ifelse(pos01_infeasible | infeasible_estimates, FALSE, TRUE)
      )

  } else if (health_indicator_type == "proportion") {
    # Proportion models ----------------------------------------------------

    model_objects <- list()

    # Binomial (requires counts)
    if (use_num_den) {
      model_objects$model_binomial <- tryCatch(
        glm(
          cbind(health_numerator, health_denominator - health_numerator) ~ social_position,
          data = df,
          family = binomial(link = "logit")
        ),
        error = function(e) stop("Failed to fit binomial model: ", e$message)
      )
    }

    # Linear models (weighted)
    model_objects$model_linear <- tryCatch(
      glm(
        health_indicator ~ social_position,
        data = df,
        weights = pop_use,
        family = gaussian()
      ),
      error = function(e) stop("Failed to fit linear model: ", e$message)
    )

    model_objects$model_linear_lny <- tryCatch(
      glm(
        log(health_indicator) ~ social_position,
        data = df %>% mutate(health_indicator = ifelse(health_indicator == 0, eps, health_indicator)),
        weights = pop_use,
        family = gaussian()
      ),
      error = function(e) stop("Failed to fit linear ln(y) model: ", e$message)
    )

    # linear model on logit(y): transform y safely into (eps, 1-eps)
    model_objects$model_linear_logit <- tryCatch(
      glm(
        logit(health_indicator) ~ social_position,
        data = df %>% mutate(health_indicator = ifelse(health_indicator == 0, eps, health_indicator),
                             health_indicator = ifelse(health_indicator == 1, 1 - eps, health_indicator)),
        weights = pop_use,
        family = gaussian()
      ),
      error = function(e) stop("Failed to fit linear logit(y) model: ", e$message)
    )

    # Quasibinomial (alternative when counts not available or to allow overdispersion)
    model_objects$model_quasibinomial <- tryCatch(
      glm(
        health_indicator ~ social_position,
        data = df,
        weights = pop_use,
        family = quasibinomial(link = "logit")
      ),
      error = function(e) stop("Failed to fit quasibinomial model: ", e$message)
    )

    # init prediction cols
    df$pred_model_binomial <- NA_real_
    df$pred_model_quasi   <- NA_real_
    df$pred_model_linear  <- NA_real_
    df$pred_model_linear_lny <- NA_real_
    df$pred_model_linear_logit <- NA_real_

    # rowwise predictions
    df$pred_model_linear <- predict(model_objects$model_linear, type = "response")
    df$pred_model_linear_lny <- exp(predict(model_objects$model_linear_lny, type = "link"))
    df$pred_model_linear_logit <- plogis(predict(model_objects$model_linear_logit, type = "response"))
    df$pred_model_quasi <- predict(model_objects$model_quasibinomial, type = "response")

    if (!is.null(model_objects$model_binomial)) {
      df$pred_model_binomial <- predict(model_objects$model_binomial, type = "response")
    }

    # extreme predictions
    nd0 <- data.frame(social_position = 0)
    nd1 <- data.frame(social_position = 1)

    p0_lin  <- as.numeric(predict(model_objects$model_linear, newdata = nd0, type = "response"))
    p1_lin  <- as.numeric(predict(model_objects$model_linear, newdata = nd1, type = "response"))

    p0_lny  <- exp(as.numeric(predict(model_objects$model_linear_lny, newdata = nd0, type = "link")))
    p1_lny  <- exp(as.numeric(predict(model_objects$model_linear_lny, newdata = nd1, type = "link")))

    p0_llog <- plogis(as.numeric(predict(model_objects$model_linear_logit, newdata = nd0, type = "response")))
    p1_llog <- plogis(as.numeric(predict(model_objects$model_linear_logit, newdata = nd1, type = "response")))

    p0_bin  <- if (!is.null(model_objects$model_binomial)) as.numeric(predict(model_objects$model_binomial, newdata = nd0, type = "response")) else NA_real_
    p1_bin  <- if (!is.null(model_objects$model_binomial)) as.numeric(predict(model_objects$model_binomial, newdata = nd1, type = "response")) else NA_real_

    p0_qb   <- if (!is.null(model_objects$model_quasibinomial)) as.numeric(predict(model_objects$model_quasibinomial, newdata = nd0, type = "response")) else NA_real_
    p1_qb   <- if (!is.null(model_objects$model_quasibinomial)) as.numeric(predict(model_objects$model_quasibinomial, newdata = nd1, type = "response")) else NA_real_

    # fit metrics table (quasibinomial appears only when !use_num_den)
    fit_metrics <- data.frame(
      model = c("linear", "linear ln(y)", if (use_num_den) "binomial" else NA, "linear logit(y)", if (!use_num_den) "quasibinomial" else NA),
      MAE = c(
        mae(df$health_indicator, df$pred_model_linear),
        mae(df$health_indicator, df$pred_model_linear_lny),
        if (use_num_den) if (!is.null(df$pred_model_binomial)) mae(df$health_indicator, df$pred_model_binomial) else NA_real_ else NA_real_,
        mae(df$health_indicator, df$pred_model_linear_logit),
        if (!is.null(model_objects$model_quasibinomial)) mae(df$health_indicator, df$pred_model_quasi) else NA_real_
      ),
      MSE = c(
        mse(df$health_indicator, df$pred_model_linear),
        mse(df$health_indicator, df$pred_model_linear_lny),
        if (use_num_den) if (!is.null(df$pred_model_binomial)) mse(df$health_indicator, df$pred_model_binomial) else NA_real_ else NA_real_,
        mse(df$health_indicator, df$pred_model_linear_logit),
        if (!is.null(model_objects$model_quasibinomial)) mse(df$health_indicator, df$pred_model_quasi) else NA_real_
      ),
      RMSE = c(
        rmse(df$health_indicator, df$pred_model_linear),
        rmse(df$health_indicator, df$pred_model_linear_lny),
        if (use_num_den) if (!is.null(df$pred_model_binomial)) rmse(df$health_indicator, df$pred_model_binomial) else NA_real_ else NA_real_,
        rmse(df$health_indicator, df$pred_model_linear_logit),
        if (!is.null(model_objects$model_quasibinomial)) rmse(df$health_indicator, df$pred_model_quasi) else NA_real_
      ),
      MAPE = c(
        mape(df$health_indicator, df$pred_model_linear),
        mape(df$health_indicator, df$pred_model_linear_lny),
        if (use_num_den) if (!is.null(df$pred_model_binomial)) mape(df$health_indicator, df$pred_model_binomial) else NA_real_ else NA_real_,
        mape(df$health_indicator, df$pred_model_linear_logit),
        if (!is.null(model_objects$model_quasibinomial)) mape(df$health_indicator, df$pred_model_quasi) else NA_real_
      ),
      AIC = c(
        model_objects$model_linear$aic,
        model_objects$model_linear_lny$aic,
        if (use_num_den) if (!is.null(model_objects$model_binomial)) model_objects$model_binomial$aic else NA_real_ else NA_real_,
        model_objects$model_linear_logit$aic,
        NA_real_  # quasibinomial does not have a meaningful AIC
      ),
      infeasible_estimates = c(
        any(df$pred_model_linear < 0 | df$pred_model_linear > 1, na.rm = TRUE),
        any(df$pred_model_linear_lny < 0 | df$pred_model_linear_lny > 1, na.rm = TRUE),
        if (use_num_den) if (!is.null(df$pred_model_binomial)) any(df$pred_model_binomial < 0 | df$pred_model_binomial > 1, na.rm = TRUE) else FALSE else FALSE,
        any(df$pred_model_linear_logit < 0 | df$pred_model_linear_logit > 1, na.rm = TRUE),
        if (!is.null(model_objects$model_quasibinomial)) any(df$pred_model_quasi < 0 | df$pred_model_quasi > 1, na.rm = TRUE) else FALSE
      ),
      pred_pos0 = c(p0_lin, p0_lny, p0_bin, p0_llog, p0_qb),
      pred_pos1 = c(p1_lin, p1_lny, p1_bin, p1_llog, p1_qb),
      indicator_has_zero = rep(any(df$health_indicator == 0, na.rm = TRUE), 5),
      stringsAsFactors = FALSE
    )

    fit_metrics1 <- fit_metrics %>% dplyr::filter(!is.na(model))

    fit_metrics2 <- fit_metrics1 %>%
      mutate(
        pos01_infeasible = (pred_pos0 < 0 | pred_pos0 > 1 | pred_pos1 < 0 | pred_pos1 > 1),
        model_viable = ifelse(pos01_infeasible | infeasible_estimates, FALSE, TRUE)
      )

    fit_metrics2 <- fit_metrics2 %>%
      rename(health_indicator_has_zero = indicator_has_zero) %>%
      select(-pos01_infeasible, -infeasible_estimates)

  }

  # If `models` (user selection) exists, keep only those rows in the fit table
  if (!is.null(models) && length(models) > 0) {
    models <- as.character(models)
    fit_metrics3 <- fit_metrics2 %>% dplyr::filter(model %in% models)
  }

  # ────────────────────────────────────────────────────────────────
  # Filter fit_metrics to user-requested models (support both 'models' and 'modelos')
  # ────────────────────────────────────────────────────────────────
  requested_models <- get0("models", ifnotfound = get0("modelos", ifnotfound = NULL))

  if (!is.null(requested_models)) {
    fit_metrics3 <- fit_metrics3 %>% dplyr::filter(model %in% requested_models)
  }

  # If none of the requested models are available -> informative error
  if (nrow(fit_metrics3) == 0) {
    stop("None of the requested models (in 'models' / 'modelos') could be evaluated. Check whether some models require numerator/denominator inputs.")
  }

  # ────────────────────────────────────────────────────
  # Select best model among viable candidates
  # ────────────────────────────────────────────────────
  candidates <- fit_metrics3 %>% dplyr::filter(model_viable == TRUE)

  if (nrow(candidates) == 0L) {
    stop(
      "No viable model found: predictions are inconsistent with the indicator type.\n",
      "- For 'rate' models: there are negative predictions for extreme social positions.\n",
      "- For 'proportion' models: there are predictions outside [0, 1] for extreme social positions.\n"
    )
  }

  # Determine which evaluation metric to use (support both English/Spanish param names)
  metric_name <- get0("model_eval_metric", ifnotfound = get0("metrica_evaluacion_modelo", ifnotfound = "MSE"))

  # Validate metric exists in table; fallback to "MSE" if missing
  if (!(metric_name %in% names(candidates))) {
    warning(sprintf("Requested metric '%s' not found in fit table. Falling back to 'MSE'.", metric_name))
    metric_name <- "MSE"
    if (!("MSE" %in% names(candidates))) {
      stop("Neither requested metric nor 'MSE' present in fit_metrics; cannot select best model.")
    }
  }

  # choose model with minimum metric (ignore NA values)
  best_model <- candidates %>%
    dplyr::filter(!!rlang::sym(metric_name) == min(!!rlang::sym(metric_name), na.rm = TRUE)) %>%
    dplyr::slice(1) %>%
    dplyr::pull(model)

  # ────────────────────────────────────────────────────────────────
  # Map chosen model name to fitted object and final predictions
  # ────────────────────────────────────────────────────────────────
  object_map <- list(
    "linear"            = model_objects$model_linear %||% NULL,
    "linear ln(y)"      = model_objects$model_linear_lny %||% NULL,
    "poisson"           = model_objects$model_poisson %||% NULL,
    "negative binomial" = model_objects$model_negbin %||% NULL,
    "binomial"          = model_objects$model_binomial %||% NULL,
    "linear logit(y)"   = model_objects$model_linear_logit %||% NULL,
    "quasibinomial"     = model_objects$model_quasibinomial %||% NULL
  )

  pred_col_map <- list(
    "linear"            = "pred_model_linear",
    "linear ln(y)"      = "pred_model_linear_lny",
    "poisson"           = "pred_model_poisson",
    "negative binomial" = "pred_model_negbin",
    "binomial"          = "pred_model_binomial",
    "linear logit(y)"   = "pred_model_linear_logit",
    "quasibinomial"     = "pred_model_quasi"
  )

  # assign selected model object (or NULL if not available)
  model <- if (!is.null(best_model) && best_model %in% names(object_map)) object_map[[best_model]] else NULL

  # assign final predictions column to df$final_model_pred (safe: check existence)
  final_pred_col <- if (!is.null(best_model) && best_model %in% names(pred_col_map)) pred_col_map[[best_model]] else NULL

  if (!is.null(final_pred_col) && final_pred_col %in% names(df)) {
    df$final_model_pred <- df[[final_pred_col]]
  } else {
    df$final_model_pred <- NA_real_
    if (!is.null(best_model)) {
      warning(sprintf("Selected model '%s' found but prediction column '%s' is missing. final_model_pred set to NA.", best_model, final_pred_col))
    } else {
      warning("No model was assigned; final_model_pred set to NA.")
    }
  }

  # ───────────────────────────────────────────────────────────────
  # Robust vcov (if a model was selected)
  # ───────────────────────────────────────────────────────────────
  robust_vcov <- tryCatch({
    if (!is.null(model)) sandwich::vcovHC(model, type = "HC1") else NULL
  }, error = function(e) NULL)

  # ───────────────────────────────────────────────────────────────
  # emmeans at extremes (social_position = 1 vs 0), contrast and RII + CI
  # ───────────────────────────────────────────────────────────────
  # Prepare df for emmeans: ensure denominator/factor exist (standardize exposure = 1)
  df_for_emm <- df %>%
    mutate(
      logit_y = (health_indicator) / (1 - health_indicator),
      health_denominator = 1,
      rate_scaling_factor = 1)

  emmeans_extremes <- tryCatch(
    emmeans::emmeans(
      model,
      specs = ~ social_position,
      at = list(social_position = c(1, 0)),
      type = "response",
      vcov. = robust_vcov,
      data = df_for_emm
    ),
    error = function(e) NULL
  )

  if (!is.null(best_model) && best_model %in% c("negative binomial", "poisson")) {
    rii_elements <- summary(
      pairs(emmeans_extremes, ratio = TRUE), type = "response",
      infer = c(TRUE, TRUE), level = conf_level
    )
  } else if ( !is.null(best_model) && best_model %in% c("linear" , "linear ln(y)")) {
    rii_elements <- summary(
      pairs(regrid(emmeans_extremes, transform = "log"), type = "response"),
      infer = c(TRUE, TRUE), level = conf_level
    )
  } else if (!is.null(best_model) && best_model %in% c("binomial", "quasibinomial", "linear logit(y)")) {
    rii_elements <- summary(
      pairs(regrid(emmeans_extremes, transform = "log"), type = "response"),
      infer = c(TRUE, TRUE), level = conf_level
    )
  }

  rii_est <- rii_elements$ratio
  se.formula <- rii_elements$SE
  rii_lwr <- if (is.null(rii_elements$lower.CL)) rii_elements$asymp.LCL else rii_elements$lower.CL
  rii_upr <- if (is.null(rii_elements$upper.CL)) rii_elements$asymp.UCL else rii_elements$upper.CL

  # Build result row
  summary_table <- data.frame(
    inequality_metric = "Relative index of inequality (RII)",
    selected_model = if (exists("best_model")) best_model else NA_character_,
    value = rii_est,
    ci_lower = rii_lwr,
    ci_upper = rii_upr,
    stringsAsFactors = FALSE
  )

  # Add social extremes & predictions
  has_0 <- any(df$social_position == 0)
  has_1 <- any(df$social_position == 1)
  df_extremes <- as.data.frame(emmeans_extremes)

  extreme_0 <- df_extremes %>%
    dplyr::filter(social_position == 0) %>%
    dplyr::pull(2)

  extreme_1 <- df_extremes %>%
    dplyr::filter(social_position == 1) %>%
    dplyr::pull(2)

  unit_analysis_ext <- c(if (!has_0) "-", df$unit_analysis, if (!has_1) "-")
  social_position_ext <- c(if (!has_0) 0, df$social_position, if (!has_1) 1)
  health_indicator_ext <- c(if (!has_0) NA_real_, df$health_indicator, if (!has_1) NA_real_)
  final_model_pred_ext <- c(if (!has_0) extreme_0, df$final_model_pred, if (!has_1)  extreme_1 )

  data_with_extremes <- data.frame(
    unit_analysis = unit_analysis_ext,
    social_position = social_position_ext,
    health_indicator = health_indicator_ext,
    final_model_pred = final_model_pred_ext
  )

  # ────────────────
  # Interpretations
  # ────────────────
  factor_legible <- format(rate_scaling_factor, scientific = FALSE, big.mark = ",")
  eps_legible <- format(eps, scientific = FALSE, big.mark = ",")

  if (health_indicator_type == "rate") {
    if (rii_est > 1) {
      part_main_english <- paste0(
        "The Relative Index of Inequality is ", round(rii_est, 2),
        ", ", conf_level*100, "%CI [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "This means that the Health Indicator is higher in the most advantaged group. ",
        "In relative terms, the most advantaged group has ", round(rii_est, 2),
        " times the number of events per ", factor_legible,
        " individuals in the reference population compared to the most disadvantaged group."
      )
    } else if (rii_est < 1) {
      part_main_english <- paste0(
        "The Relative Index of Inequality is ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%CI [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "This means that the Health Indicator is higher in the most disadvantaged group. ",
        "In relative terms, the most advantaged group has ", round(rii_est*100, 2),
        "% of the number of events per ", factor_legible,
        " individuals in the reference population compared to the most disadvantaged group."
      )
    } else if (rii_est == 1) {
      part_main_english <- paste0(
        "The Relative Index of Inequality is ", round(rii_est, 2),
        ", ", conf_level*100, "%CI [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "This indicates that, in relative terms, there is no difference in the number of events per ",
        factor_legible, " individuals in the reference population between the most advantaged and the most disadvantaged groups."
      )
    }
  } else if (health_indicator_type == "proportion") {
    if (rii_est > 1) {
      part_main_english <- paste0(
        "The Relative Index of Inequality is ", round(rii_est, 2),
        ", ", conf_level*100, "%CI [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "This means that the Health Indicator is higher in the most advantaged group. ",
        "In relative terms, the value in the most advantaged group is ", round(rii_est, 2),
        " times that of the most disadvantaged group."
      )
    } else if (rii_est < 1) {
      part_main_english <- paste0(
        "The Relative Index of Inequality is ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%CI [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "This means that the Health Indicator is higher in the most disadvantaged group. ",
        "In relative terms, the value in the most advantaged group represents ", round(rii_est*100, 2),
        "% of the value in the most disadvantaged group."
      )
    } else if (rii_est == 1) {
      part_main_english <- paste0(
        "The Relative Index of Inequality is ", round(rii_est, 2),
        ", ", conf_level*100, "%CI [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "This indicates that, in relative terms, there is no difference between the most advantaged and the most disadvantaged groups."
      )
    }
  }

  if (rii_lwr <= 1 & rii_upr >= 1) {
    part_sig_english <- " However, this relationship is not statistically significant."
  } else {
    part_sig_english <- " Moreover, this relationship is statistically significant."
  }

  interpretation_english <- paste0(part_main_english, part_sig_english)


  if (health_indicator_type == "rate") {
    if (rii_est > 1) {
      part_main_spanish <- paste0(
        "El Índice de Desigualdad Relativo es de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Esto significa que el Indicador de Salud es mayor en el extremo más favorecido. ",
        "En términos relativos, el extremo más favorecido tiene ", round(rii_est, 2),
        " veces el número de eventos por cada ", factor_legible,
        " personas de la población base en comparación con el extremo más desfavorecido."
      )
    } else if (rii_est < 1) {
      part_main_spanish <- paste0(
        "El Índice de Desigualdad Relativo es de ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Esto significa que el Indicador de Salud es mayor en el extremo más desfavorecido. ",
        "En términos relativos, el extremo más favorecido tiene ", round(rii_est*100, 2),
        "% del número de eventos por cada ", factor_legible,
        " personas de la población base en comparación con el extremo más desfavorecido."
      )
    } else if (rii_est == 1) {
      part_main_spanish <- paste0(
        "El Índice de Desigualdad Relativo es de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Esto indica que, en términos relativos, no hay diferencia en el número de eventos por cada ",
        factor_legible, " personas de la población base entre el extremo más favorecido y el extremo más desfavorecido."
      )
    }
  } else if (health_indicator_type == "proportion") {
    if (rii_est > 1) {
      part_main_spanish <- paste0(
        "El Índice de Desigualdad Relativo es de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Esto significa que el Indicador de Salud es mayor en el extremo más favorecido. ",
        "En términos relativos, el valor del extremo más favorecido es ", round(rii_est, 2),
        " veces el del extremo más desfavorecido."
      )
    } else if (rii_est < 1) {
      part_main_spanish <- paste0(
        "El Índice de Desigualdad Relativo es de ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Esto significa que el Indicador de Salud es mayor en el extremo más desfavorecido. ",
        "En términos relativos, el valor en el extremo más favorecido representa ", round(rii_est*100, 2),
        "% del valor del extremo más desfavorecido."
      )
    } else if (rii_est == 1) {
      part_main_spanish <- paste0(
        "El Índice de Desigualdad Relativo es de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Esto indica que, en términos relativos, no existe diferencia entre el extremo más favorecido y el extremo más desfavorecido."
      )
    }
  }

  if (rii_lwr <= 1 & rii_upr >= 1) {
    part_sig_spanish <- " Sin embargo, esta relación no es estadísticamente significativa."
  } else {
    part_sig_spanish <- " Además, esta relación es estadísticamente significativa."
  }

  interpretation_spanish <- paste0(part_main_spanish, part_sig_spanish)


  if (health_indicator_type == "rate") {
    if (rii_est > 1) {
      part_main_french <- paste0(
        "L'Indice d'Inégalité Relative est de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Cela signifie que l'Indicateur de Santé est plus élevé dans le groupe le plus avantagé. ",
        "En termes relatifs, le groupe le plus avantagé présente ", round(rii_est, 2),
        " fois le nombre d'événements pour ", factor_legible,
        " personnes de la population de base par rapport au groupe le plus défavorisé."
      )
    } else if (rii_est < 1) {
      part_main_french <- paste0(
        "L'Indice d'Inégalité Relative est de ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Cela signifie que l'Indicateur de Santé est plus élevé dans le groupe le plus défavorisé. ",
        "En termes relatifs, le groupe le plus avantagé présente ", round(rii_est*100, 2),
        "% du nombre d'événements pour ", factor_legible,
        " personnes de la population de base par rapport au groupe le plus défavorisé."
      )
    } else if (rii_est == 1) {
      part_main_french <- paste0(
        "L'Indice d'Inégalité Relative est de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Cela indique qu'il n'existe pas de différence relative dans le nombre d'événements pour ",
        factor_legible, " personnes de la population de base entre le groupe le plus avantagé et le groupe le plus défavorisé."
      )
    }
  } else if (health_indicator_type == "proportion") {
    if (rii_est > 1) {
      part_main_french <- paste0(
        "L'Indice d'Inégalité Relative est de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Cela signifie que l'Indicateur de Santé est plus élevé dans le groupe le plus avantagé. ",
        "En termes relatifs, la valeur dans le groupe le plus avantagé est ", round(rii_est, 2),
        " fois celle du groupe le plus défavorisé."
      )
    } else if (rii_est < 1) {
      part_main_french <- paste0(
        "L'Indice d'Inégalité Relative est de ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Cela signifie que l'Indicateur de Santé est plus élevé dans le groupe le plus défavorisé. ",
        "En termes relatifs, la valeur dans le groupe le plus avantagé représente ", round(rii_est*100, 2),
        "% de celle du groupe le plus défavorisé."
      )
    } else if (rii_est == 1) {
      part_main_french <- paste0(
        "L'Indice d'Inégalité Relative est de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Cela indique qu'il n'existe pas de différence relative entre le groupe le plus avantagé et le groupe le plus défavorisé."
      )
    }
  }

  if (rii_lwr <= 1 & rii_upr >= 1) {
    part_sig_french <- " Cependant, cette relation n’est pas statistiquement significative."
  } else {
    part_sig_french <- " De plus, cette relation est statistiquement significative."
  }

  interpretation_french <- paste0(part_main_french, part_sig_french)


  if (health_indicator_type == "rate") {
    if (rii_est > 1) {
      part_main_portuguese <- paste0(
        "O Índice de Desigualdade Relativa é de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Isso significa que o Indicador de Saúde é maior no extremo mais favorecido. ",
        "Em termos relativos, o extremo mais favorecido apresenta ", round(rii_est, 2),
        " vezes o número de eventos por cada ", factor_legible,
        " pessoas da população base em comparação com o extremo mais desfavorecido."
      )
    } else if (rii_est < 1) {
      part_main_portuguese <- paste0(
        "O Índice de Desigualdade Relativa é de ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Isso significa que o Indicador de Saúde é maior no extremo mais desfavorecido. ",
        "Em termos relativos, o extremo mais favorecido apresenta ", round(rii_est*100, 2),
        "% do número de eventos por cada ", factor_legible,
        " pessoas da população base em comparação com o extremo mais desfavorecido."
      )
    } else if (rii_est == 1) {
      part_main_portuguese <- paste0(
        "O Índice de Desigualdade Relativa é de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Isso indica que, em termos relativos, não há diferença no número de eventos por cada ",
        factor_legible, " pessoas da população base entre o extremo mais favorecido e o extremo mais desfavorecido."
      )
    }
  } else if (health_indicator_type == "proportion") {
    if (rii_est > 1) {
      part_main_portuguese <- paste0(
        "O Índice de Desigualdade Relativa é de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr, 2), "; ", round(rii_upr, 2), "]. ",
        "Isso significa que o Indicador de Saúde é maior no extremo mais favorecido. ",
        "Em termos relativos, o valor no extremo mais favorecido é ", round(rii_est, 2),
        " vezes o valor do extremo mais desfavorecido."
      )
    } else if (rii_est < 1) {
      part_main_portuguese <- paste0(
        "O Índice de Desigualdade Relativa é de ", round(rii_est*100, 2),
        "%, ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Isso significa que o Indicador de Saúde é maior no extremo mais desfavorecido. ",
        "Em termos relativos, o valor no extremo mais favorecido representa ", round(rii_est*100, 2),
        "% do valor do extremo mais desfavorecido."
      )
    } else if (rii_est == 1) {
      part_main_portuguese <- paste0(
        "O Índice de Desigualdade Relativa é de ", round(rii_est, 2),
        ", ", conf_level*100, "%IC [", round(rii_lwr*100, 2), "%; ", round(rii_upr*100, 2), "%]. ",
        "Isso indica que, em termos relativos, não há diferença entre o extremo mais favorecido e o extremo mais desfavorecido."
      )
    }
  }

  if (rii_lwr <= 1 & rii_upr >= 1) {
    part_sig_portuguese <- " No entanto, essa relação não é estatisticamente significativa."
  } else {
    part_sig_portuguese <- " Além disso, essa relação é estatisticamente significativa."
  }

  interpretation_portuguese <- paste0(part_main_portuguese, part_sig_portuguese)


  # Note based on indicator type
  note_english <- if (health_indicator_type == "rate") {
    paste0("Note: Events correspond to the numerator and the base population to the denominator. The disadvantaged extreme = posicion_social 0; advantaged = 1.")
  } else {
    paste0("Note: Indicator expressed as proportions (0–1). The disadvantaged extreme = posicion_social 0; advantaged = 1.")
  }
  alert_english <- NULL
  if (any(df$health_indicator == 0) && summary_table$selected_model == "linear ln(y)") {
    alert_english <- paste0("Health indicator contains at least one 0. To apply the 'linear ln(y)' model, zeros were replaced by ", eps_legible, ". Please assess whether this affects the validity of the model fit.")
    message(alert_english)
  }
  if (!is.null(alert_english)) {
    note_english <- paste(note_english, alert_english, sep = "\n\n")
  }


  note_spanish <- if (health_indicator_type == "rate") {
    paste0("Nota: Los eventos corresponden al numerador y la población base al denominador. El extremo desfavorecido = posicion_social 0; favorecido = 1.")
  } else {
    paste0("Nota: Indicador en proporciones (0–1). El extremo desfavorecido = posicion_social 0; favorecido = 1.")
  }
  alert_spanish <- NULL
  if (any(df$health_indicator == 0) && summary_table$selected_model == "linear ln(y)") {
    alert_spanish <- paste0("Health indicator contiene al menos un 0. Para aplicar el modelo 'linear ln(y)', los ceros fueron reemplazados por ", eps_legible, ". Evalúe si esto compromete la validez del ajuste.")
    message(alert_spanish)
  }
  if (!is.null(alert_spanish)) {note_spanish <- paste(note_spanish, alert_spanish, sep = "\n\n")}


  note_french <- if (health_indicator_type == "rate") {
    paste0("Remarque : Les événements correspondent au numérateur et la population de base au dénominateur. L'extrême défavorisé = position_sociale 0 ; favorisé = 1.")
  } else {
    paste0("Remarque : Indicateur exprimé en proportions (0–1). L'extrême défavorisé = position_sociale 0 ; favorisé = 1.")
  }
  alert_french <- NULL
  if (any(df$health_indicator == 0) && summary_table$selected_model == "linear ln(y)") {
    alert_french <- paste0("L'indicateur de santé contient au moins une valeur égale à 0. Pour appliquer le modèle 'linear ln(y)', les zéros ont été remplacés par ", eps_legible, ". Veuillez évaluer si cela compromet la validité de l'ajustement.")
    message(alert_french)
  }
  if (!is.null(alert_french)) {
    note_french <- paste(note_french, alert_french, sep = "\n\n")
  }


  note_portuguese <- if (health_indicator_type == "rate") {
    paste0("Nota: Os eventos correspondem ao numerador e a população base ao denominador. O extremo desfavorecido = posicion_social 0; favorecido = 1.")
  } else {
    paste0("Nota: Indicador em proporções (0–1). O extremo desfavorecido = posicion_social 0; favorecido = 1.")
  }
  alert_portuguese <- NULL
  if (any(df$health_indicator == 0) && summary_table$selected_model == "linear ln(y)") {
    alert_portuguese <- paste0("Health indicator contém pelo menos um valor igual a 0. Para aplicar o modelo 'linear ln(y)', os zeros foram substituídos por ", eps_legible, ". Avalie se isso compromete a validade do ajuste.")
    message(alert_portuguese)
  }
  if (!is.null(alert_portuguese)) {
    note_portuguese <- paste(note_portuguese, alert_portuguese, sep = "\n\n")
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
    global_health_mean = global_health_mean,
    dx_models_regression = fit_metrics2,
    data = df,
    data_with_extremes = data_with_extremes
  ))
}

