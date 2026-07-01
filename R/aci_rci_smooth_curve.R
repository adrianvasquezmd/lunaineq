#' Non-Parametric Smoothing Engine for Concentration Curves
#'
#' @description
#' `aci_rci_smooth_curve` computes generalized geometric projections of the discrete
#' empirical concentration curve for visualization purposes.
#'
#' The function fits optional continuous interpolations over the cumulative fractional
#' rank (x-axis) and cumulative health share (y-axis) space. These non-parametric
#' extensions (such as locally estimated scatterplot smoothing) are purely illustrative;
#' they strictly bypass the principal estimation pipeline and do not inform the core RCI
#' or ACI indices, bounds, or variance approximations.
#'
#' Implemented smoothing methods:
#' \describe{
#'   \item{loess}{
#'     Local polynomial regression with tricube weighting on the cumulative
#'     population share. This formulation provides high local flexibility, though
#'     asymptotic monotonicity is not formally guaranteed.
#'   }
#'   \item{ml}{
#'     Parametric maximum-likelihood/nonlinear least-squares style smoothing
#'     using a simple one-parameter concentration-curve family:
#'     L(p) = p^theta, with theta > 0. This preserves endpoints and monotonicity
#'     for theta > 0. The parameter is estimated by minimizing squared residuals
#'     on the empirical curve points.
#'   }
#' }
#'
#' Both smoothed curves are returned on the same empirical cumulative population
#' grid. Endpoint constraints are enforced so that the smoothed curves start at
#' (0, 0) and end at (1, 1). Values outside \code{[0, 1]} caused by smoothing are clipped
#' and reported through warnings.
#'
#' @param prepared A list returned by `.aci_rci_prepare_data()`. It must include
#'   `curve_data`.
#' @param config A validated configuration list returned by
#'   `.aci_rci_validate_inputs()`.
#' @param loess_span Numeric. Smoothing span used by LOESS. Defaults to 0.75.
#' @param tolerance Numeric. Numerical tolerance for endpoint and range checks.
#'
#' @return A list with:
#' \describe{
#'   \item{curve_data}{Curve data with empirical and optional smoothed columns.}
#'   \item{warnings}{Warnings generated while smoothing the curve.}
#'   \item{diagnostics}{Diagnostics describing requested and computed smoothing
#'     methods.}
#' }
#'
#' @keywords internal
#' @noRd
.aci_rci_smooth_curve <- function(
  prepared,
  config,
  loess_span = 0.75,
  tolerance = 1e-10
) {
  warnings_out <- character()
  diagnostics <- list()

  add_warning <- function(message) {
    if (isTRUE(config$verbose)) warning(message, call. = FALSE, immediate. = TRUE)
    warnings_out <<- c(warnings_out, message)
    invisible(NULL)
  }

  add_message <- function(message) {
    if (isTRUE(config$verbose)) message(message)
    invisible(NULL)
  }

  stop_smooth <- function(x) {
    stop(x, call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Basic checks
  # -------------------------------------------------------------------------

  if (!is.list(prepared)) {
    stop_smooth("`prepared` must be the list returned by `.aci_rci_prepare_data()`.")
  }

  if (!is.list(config)) {
    stop_smooth("`config` must be the list returned by `.aci_rci_validate_inputs()`.")
  }

  if (!is.numeric(loess_span) || length(loess_span) != 1L ||
      is.na(loess_span) || loess_span <= 0 || loess_span > 1) {
    stop_smooth("`loess_span` must be a single numeric value greater than 0 and less than or equal to 1.")
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) || tolerance <= 0) {
    stop_smooth("`tolerance` must be a single positive numeric value.")
  }

  if (is.null(prepared$curve_data) || !is.data.frame(prepared$curve_data)) {
    stop_smooth("`prepared$curve_data` must be a data frame.")
  }

  curve_data <- prepared$curve_data

  required_curve_cols <- c(
    ".cum_population",
    ".cum_health_empirical"
  )

  missing_curve_cols <- setdiff(required_curve_cols, names(curve_data))

  if (length(missing_curve_cols) > 0L) {
    stop_smooth(
      paste0(
        "`prepared$curve_data` is missing required column(s): ",
        paste(missing_curve_cols, collapse = ", "),
        "."
      )
    )
  }

  if (!is.numeric(curve_data$.cum_population) ||
      !is.numeric(curve_data$.cum_health_empirical)) {
    stop_smooth("Curve columns `.cum_population` and `.cum_health_empirical` must be numeric.")
  }

  if (any(!is.finite(curve_data$.cum_population)) ||
      any(!is.finite(curve_data$.cum_health_empirical))) {
    stop_smooth("Curve columns contain non-finite values.")
  }

  if (nrow(curve_data) < 2L) {
    stop_smooth("At least two concentration-curve points are required.")
  }

  p <- curve_data$.cum_population
  L <- curve_data$.cum_health_empirical

  if (any(diff(p) < -tolerance)) {
    stop_smooth("Cumulative population shares must be monotonically non-decreasing.")
  }

  if (abs(p[1L]) > tolerance || abs(L[1L]) > tolerance) {
    add_warning(
      "The empirical concentration curve does not start exactly at (0, 0). The curve will be kept as supplied, but smoothed curves will enforce the endpoint."
    )
  }

  if (abs(utils::tail(p, 1L) - 1) > tolerance ||
      abs(utils::tail(L, 1L) - 1) > tolerance) {
    add_warning(
      "The empirical concentration curve does not end exactly at (1, 1). The curve will be kept as supplied, but smoothed curves will enforce the endpoint."
    )
  }

  if (any(L < -tolerance | L > 1 + tolerance, na.rm = TRUE)) {
    add_warning(
      "The empirical cumulative health shares fall outside [0, 1]. Smoothed curves may be difficult to interpret."
    )
  }

  if (any(diff(L) < -tolerance, na.rm = TRUE)) {
    add_warning(
      "The empirical concentration curve is not monotone increasing. This usually reflects negative health values or unusual grouped data; smoothing is for visualization only."
    )
  }

  # Ensure standard public-facing empirical column exists. Keep legacy/internal
  # columns as well to preserve harmony with `.aci_rci_prepare_data()`.
  if (!"cum_population" %in% names(curve_data)) {
    curve_data$cum_population <- curve_data$.cum_population
  }

  if (!"cum_health_empirical" %in% names(curve_data)) {
    curve_data$cum_health_empirical <- curve_data$.cum_health_empirical
  }

  # -------------------------------------------------------------------------
  # Smoothing switches
  # -------------------------------------------------------------------------

  add_smoothed_curves <- isTRUE(config$add_smoothed_curves)
  smoothing_methods <- config$smoothing_methods

  if (is.null(smoothing_methods)) {
    smoothing_methods <- character()
  }

  smoothing_methods <- unique(tolower(smoothing_methods))

  allowed_methods <- c("ml", "loess")
  invalid_methods <- setdiff(smoothing_methods, allowed_methods)

  if (length(invalid_methods) > 0L) {
    stop_smooth(
      paste0(
        "Invalid smoothing method(s): ",
        paste(invalid_methods, collapse = ", "),
        ". Allowed methods are: ",
        paste(allowed_methods, collapse = ", "),
        "."
      )
    )
  }

  diagnostics$add_smoothed_curves <- add_smoothed_curves
  diagnostics$smoothing_methods_requested <- smoothing_methods
  diagnostics$smoothing_methods_computed <- character()
  diagnostics$smoothing_methods_failed <- character()
  diagnostics$smoothing_used_for_estimation <- FALSE
  diagnostics$loess_span <- loess_span

  if (!add_smoothed_curves || length(smoothing_methods) == 0L) {
    add_message(
      "`add_smoothed_curves = FALSE` or no smoothing methods were requested. Only the empirical concentration curve will be returned."
    )

    diagnostics$smoothing_methods_computed <- character()
    diagnostics$smoothing_methods_failed <- character()

    return(
      list(
        curve_data = curve_data,
        warnings = unique(warnings_out),
        diagnostics = diagnostics
      )
    )
  }

  add_message(
    "Smoothed concentration curves are computed for visualization only and are not used for ACI/RCI estimation, bounded corrections, standard errors, or confidence intervals."
  )

  # -------------------------------------------------------------------------
  # Utilities
  # -------------------------------------------------------------------------

  clip01 <- function(x, method_name) {
    if (any(x < -tolerance | x > 1 + tolerance, na.rm = TRUE)) {
      add_warning(
        paste0(
          "The ", method_name,
          "-smoothed curve produced values outside [0, 1]. Values were clipped to [0, 1]."
        )
      )
    }

    pmin(pmax(x, 0), 1)
  }

  enforce_endpoints <- function(x) {
    if (length(x) >= 1L) {
      x[1L] <- 0
      x[length(x)] <- 1
    }
    x
  }

  # -------------------------------------------------------------------------
  # LOESS smoothing
  # -------------------------------------------------------------------------

  if ("loess" %in% smoothing_methods) {
    loess_values <- tryCatch(
      {
        # LOESS needs enough distinct x values. This check is deliberately
        # conservative because ecological data often have few groups.
        if (length(unique(p)) < 4L) {
          stop(
            "LOESS smoothing requires at least four distinct cumulative population points.",
            call. = FALSE
          )
        }

        fit <- stats::loess(
          formula = L ~ p,
          span = loess_span,
          degree = 2,
          family = "gaussian",
          control = stats::loess.control(surface = "direct")
        )

        pred <- stats::predict(fit, newdata = data.frame(p = p))

        if (any(!is.finite(pred))) {
          stop("LOESS produced non-finite fitted values.", call. = FALSE)
        }

        pred <- clip01(pred, "LOESS")
        pred <- enforce_endpoints(pred)

        if (any(diff(pred) < -tolerance, na.rm = TRUE)) {
          add_warning(
            "The LOESS-smoothed concentration curve is not monotone increasing. This is expected in some data sets; the curve is for visualization only."
          )
        }

        pred
      },
      error = function(e) {
        add_warning(
          paste0(
            "LOESS smoothing failed: ",
            conditionMessage(e),
            " The `cum_health_loess` column will be returned as NA."
          )
        )
        rep(NA_real_, length(p))
      }
    )

    curve_data$cum_health_loess <- loess_values

    if (all(is.na(loess_values))) {
      diagnostics$smoothing_methods_failed <- c(
        diagnostics$smoothing_methods_failed,
        "loess"
      )
    } else {
      diagnostics$smoothing_methods_computed <- c(
        diagnostics$smoothing_methods_computed,
        "loess"
      )
    }
  }

  # -------------------------------------------------------------------------
  # ML / parametric smoothing
  # -------------------------------------------------------------------------

  if ("ml" %in% smoothing_methods) {
    ml_values <- tryCatch(
      {
        # Parametric curve family: Murray & Lopez (ML) smoothing method
        # L(p) = (exp(p / (k - p)) - 1) / (exp(1 / (k - 1)) - 1)
        # This preserves endpoints (0,0) and (1,1).
        
        func_lorenz <- function(x, k) {
          val <- (exp(x / (k - x)) - 1) / (exp(1 / (k - 1)) - 1)
          val[is.nan(val) | is.infinite(val)] <- x[is.nan(val) | is.infinite(val)]
          val
        }

        loss_lorenz <- function(k) {
          sum((L - func_lorenz(p, k))^2, na.rm = TRUE)
        }

        # Optimization splits search space to avoid k=1 discontinuity
        opt1 <- stats::optimize(loss_lorenz, interval = c(1.001, 100))
        opt2 <- stats::optimize(loss_lorenz, interval = c(-100, 0.999))

        if (opt1$objective < opt2$objective) {
          k_hat <- opt1$minimum
          obj_val <- opt1$objective
        } else {
          k_hat <- opt2$minimum
          obj_val <- opt2$objective
        }

        if (!is.finite(k_hat)) {
          stop("Estimated ML smoothing parameter is not finite.", call. = FALSE)
        }

        pred <- func_lorenz(p, k_hat)
        pred <- clip01(pred, "ML")
        pred <- enforce_endpoints(pred)

        if (any(diff(pred) < -tolerance, na.rm = TRUE)) {
          add_warning(
            "The ML-smoothed concentration curve is not monotone increasing. Check curve inputs."
          )
        }

        diagnostics$ml_k <- k_hat
        diagnostics$ml_objective_value <- obj_val

        pred
      },
      error = function(e) {
        add_warning(
          paste0(
            "ML smoothing failed: ",
            conditionMessage(e),
            " The `cum_health_ml` column will be returned as NA."
          )
        )
        diagnostics$ml_k <<- NA_real_
        diagnostics$ml_objective_value <<- NA_real_
        rep(NA_real_, length(p))
      }
    )

    curve_data$cum_health_ml <- ml_values

    if (all(is.na(ml_values))) {
      diagnostics$smoothing_methods_failed <- c(
        diagnostics$smoothing_methods_failed,
        "ml"
      )
    } else {
      diagnostics$smoothing_methods_computed <- c(
        diagnostics$smoothing_methods_computed,
        "ml"
      )
    }
  }

  # -------------------------------------------------------------------------
  # Final diagnostics
  # -------------------------------------------------------------------------

  diagnostics$smoothing_methods_computed <- unique(diagnostics$smoothing_methods_computed)
  diagnostics$smoothing_methods_failed <- unique(diagnostics$smoothing_methods_failed)
  diagnostics$n_curve_points <- nrow(curve_data)
  diagnostics$empirical_curve_monotone <- !any(diff(L) < -tolerance, na.rm = TRUE)
  diagnostics$smoothed_columns <- intersect(
    c("cum_health_ml", "cum_health_loess"),
    names(curve_data)
  )

  list(
    curve_data = curve_data,
    warnings = unique(warnings_out),
    diagnostics = diagnostics
  )
}
