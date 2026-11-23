# Phase of a fitted sinusoid  temp = YInt + s*sin(r) + c*cos(r).
#
#   "atan2" (default)  time of the fitted peak, anywhere in the cycle:
#                      period/(2*pi) * (atan2(s, c) mod 2*pi). Lags are wrapped
#                      into (-period/2, period/2].
#   "atan"             the PASTA / Diel PASTA R formulas, kept to reproduce them:
#                        annual  365/(2*pi) * (3*pi/2 - atan(c/s))
#                        diel     24/(2*pi) * (pi - atan(c/s))
#                      atan only sees c/s, so each phase covers half a cycle.

#' Phase methods
#'
#' `"atan2"` (default) gives the time of the fitted peak anywhere in the cycle
#' and wraps phase lags into plus or minus half a cycle. `"atan"` reproduces the
#' PASTA and Diel PASTA R formulas, which only see half a cycle (annual phase
#' 182.5-365 d, diel phase 6-18 h).
#' @export
PHASE_METHODS <- c("atan2", "atan")

#' @rdname PHASE_METHODS
#' @export
DEFAULT_PHASE_METHOD <- "atan2"

check_phase_method <- function(method) {
  if (length(method) != 1 || !(method %in% PHASE_METHODS)) {
    stop("phase_method must be one of ", paste(PHASE_METHODS, collapse = ", "), call. = FALSE)
  }
  method
}

#' Phase and phase lag of a fitted sinusoid
#'
#' @param s,c Sine and cosine coefficients of `temp ~ sin(r) + cos(r)`.
#' @param period 365 (days of year) or 24 (hours).
#' @param method `"atan2"` or `"atan"`, see [PHASE_METHODS].
#' @return `pasta_phase()`: phase in the units of `period`.
#' @export
pasta_phase <- function(s, c, period, method = DEFAULT_PHASE_METHOD) {
  k <- period / (2 * pi)
  if (check_phase_method(method) == "atan2") {
    return(k * (atan2(s, c) %% (2 * pi)))
  }
  offset <- if (period == 365) 1.5 * pi else pi
  k * (offset - atan(c / s))
}

#' @rdname pasta_phase
#' @param water_phase,air_phase Phases from `pasta_phase()`.
#' @return `pasta_lag()`: water minus air, wrapped into (-period/2, period/2] for atan2.
#' @export
pasta_lag <- function(water_phase, air_phase, period, method = DEFAULT_PHASE_METHOD) {
  d <- water_phase - air_phase
  if (check_phase_method(method) == "atan2") {
    half <- period / 2
    d <- half - ((half - d) %% period)
  }
  d
}
