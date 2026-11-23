# Diel signal: daily sinusoid fits and paired metrics (from
# DielPASTA_Share/R/PASTA_dielAnalysisFunctions.R and Diel_PASTA_Example.Rmd).

#' Paired stream and air records for a set of pairs
#' @param pairs List of c(site_id, stream, air), or a pair table.
#' @param loggers Named list of loggers.
#' @param water_min_C Stream records at or below this are dropped.
#' @return data.frame site_id, DateTime_r, str_temp, air_temp (the R `T.y`).
#' @export
build_paired_df <- function(pairs, loggers, water_min_C = 0.99) {
  pairs <- pair_tuples(pairs)
  rows <- lapply(pairs, function(p) {
    m <- merge(loggers[[p[2]]], loggers[[p[3]]], by = "timestamp", suffixes = c("_s", "_a"))
    m <- m[stats::complete.cases(m) & m$temp_C_s > water_min_C, , drop = FALSE]
    m <- m[order(m$timestamp), , drop = FALSE]
    if (!nrow(m)) return(NULL)
    df(site_id = p[1], DateTime_r = m$timestamp, str_temp = m$temp_C_s, air_temp = m$temp_C_a)
  })
  out <- rbind_all(rows)
  if (!nrow(out)) {
    out <- df(site_id = character(0),
              DateTime_r = as.POSIXct(numeric(0), origin = "1970-01-01", tz = .tz),
              str_temp = numeric(0), air_temp = numeric(0))
  }
  out
}

#' Daily sinusoid fits for every site, day and medium
#'
#' OLS fit `temp ~ sin(r) + cos(r)`, `r = 2*pi*h/24` with h = hour + minute/60
#' (as `rad_hour()` in R), solved for all days at once.
#' @param T_y Paired records from [build_paired_df()].
#' @param min_coverage Share of a full day of records needed to fit a day.
#' @param interval_min Logging interval in minutes.
#' @param phase_method `"atan2"` or `"atan"`, see [PHASE_METHODS].
#' @return data.frame amplitude_C, phase_h, AdjRsqr, RMSE, sinSlope, cosSlope, YInt,
#'   medium, site_id, date, year, n_obs.
#' @export
fit_daily_sinusoids <- function(T_y, min_coverage = 0.75, interval_min = 15,
                                phase_method = DEFAULT_PHASE_METHOD) {
  if (!nrow(T_y)) return(data.frame())
  lt <- as.POSIXlt(T_y$DateTime_r, tz = .tz)
  rd <- 2 * pi * (lt$hour + lt$min / 60) / 24
  s <- sin(rd)
  cc <- cos(rd)
  date <- day_of(T_y$DateTime_r)
  sites <- csort(unique(T_y$site_id))
  g <- interaction(factor(T_y$site_id, levels = sites), factor(as.numeric(date)),
                   drop = TRUE, lex.order = TRUE)
  per_day <- 24 * 60 / interval_min
  out <- list()
  for (medium in c("air", "water")) {
    y <- if (medium == "air") T_y$air_temp else T_y$str_temp
    S <- rowsum(cbind(n = 1, s = s, c = cc, y = y, ss = s * s, cc = cc * cc, sc = s * cc,
                      sy = s * y, cy = cc * y, yy = y * y), g, reorder = TRUE)
    n <- S[, "n"]
    # normal equations X'X b = X'y with X = [1, sin, cos], solved by Cramer's rule
    a11 <- n; a12 <- S[, "s"]; a13 <- S[, "c"]
    a22 <- S[, "ss"]; a23 <- S[, "sc"]; a33 <- S[, "cc"]
    b1 <- S[, "y"]; b2 <- S[, "sy"]; b3 <- S[, "cy"]
    det3 <- function(m11, m12, m13, m21, m22, m23, m31, m32, m33) {
      m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31)
    }
    D <- det3(a11, a12, a13, a12, a22, a23, a13, a23, a33)
    ok <- n >= max(4, min_coverage * per_day) & abs(D) > 1e-9
    B0 <- det3(b1, a12, a13, b2, a22, a23, b3, a23, a33) / D
    B1 <- det3(a11, b1, a13, a12, b2, a23, a13, b3, a33) / D
    B2 <- det3(a11, a12, b1, a12, a22, b2, a13, a23, b3) / D
    sse <- pmax(S[, "yy"] - (B0 * b1 + B1 * b2 + B2 * b3), 0)
    sst <- S[, "yy"] - S[, "y"]^2 / n
    r2 <- 1 - sse / sst
    keys <- do.call(rbind, strsplit(rownames(S), ".", fixed = TRUE))
    nk <- ncol(keys)
    res <- df(amplitude_C = sqrt(B1^2 + B2^2),
              phase_h = pasta_phase(B1, B2, 24, phase_method),
              AdjRsqr = 1 - (1 - r2) * (n - 1) / (n - 3),
              RMSE = sqrt(sse / n),
              sinSlope = B1, cosSlope = B2, YInt = B0, medium = medium,
              site_id = vapply(strsplit(rownames(S), ".", fixed = TRUE),
                               function(k) paste(k[-length(k)], collapse = "."), ""),
              date = as.Date(as.numeric(keys[, nk]), origin = "1970-01-01"),
              n_obs = as.integer(n))
    out[[medium]] <- res[ok, , drop = FALSE]
  }
  fits <- rbind_all(out)
  if (!nrow(fits)) return(data.frame())
  fits$year <- as.integer(format(fits$date, "%Y"))
  fits <- fits[corder(fits$site_id, fits$date, fits$medium), c(
    "amplitude_C", "phase_h", "AdjRsqr", "RMSE", "sinSlope", "cosSlope", "YInt",
    "medium", "site_id", "date", "year", "n_obs")]
  rownames(fits) <- NULL
  fits
}

#' Fitted daily sinusoid at given times
#' @param fit One row of [fit_daily_sinusoids()] output.
#' @param timestamps POSIXct times.
#' @export
diel_sinusoid_value <- function(fit, timestamps) {
  lt <- as.POSIXlt(timestamps, tz = .tz)
  rd <- 2 * pi * (lt$hour + lt$min / 60) / 24
  fit$sinSlope * sin(rd) + fit$cosSlope * cos(rd) + fit$YInt
}

#' Daily water/air metrics
#'
#' The metrics chunk of Diel_PASTA_Example.Rmd: keep days where both air and water
#' adjusted R\eqn{^2} exceed `r2_min`, then amplitude ratio, phase lag and mean ratios,
#' with centred rolling means.
#' @param fits Output of [fit_daily_sinusoids()].
#' @param r2_min Both fits must exceed this adjusted R\eqn{^2}.
#' @param roll_k Rolling-mean window (days, over kept rows).
#' @inheritParams fit_daily_sinusoids
#' @export
pasta_metrics <- function(fits, r2_min = 0.7, roll_k = 7, phase_method = DEFAULT_PHASE_METHOD) {
  if (!nrow(fits)) return(data.frame())
  a <- fits[fits$medium == "air", c("site_id", "date", "amplitude_C", "phase_h", "YInt", "AdjRsqr")]
  w <- fits[fits$medium == "water", c("site_id", "date", "amplitude_C", "phase_h", "YInt", "AdjRsqr")]
  m <- merge(a, w, by = c("site_id", "date"), suffixes = c("_air", "_water"))
  m <- m[!is.na(m$AdjRsqr_air) & !is.na(m$AdjRsqr_water) &
           m$AdjRsqr_air > r2_min & m$AdjRsqr_water > r2_min, , drop = FALSE]
  m <- m[corder(m$site_id, m$date), , drop = FALSE]
  out <- df(site_id = m$site_id, date = m$date,
            AmpRatio = np_round(m$amplitude_C_water / m$amplitude_C_air, 2),
            PhaseLag_h = np_round(pasta_lag(m$phase_h_water, m$phase_h_air, 24, phase_method), 2),
            Ratio_Mean_C = np_round(m$YInt_water / m$YInt_air, 5),
            Ratio_Mean_K = np_round((m$YInt_water + 273.15) / (m$YInt_air + 273.15), 5),
            Amp_air_C = np_round(m$amplitude_C_air, 3),
            Amp_water_C = np_round(m$amplitude_C_water, 3),
            Mean_air_C = np_round(m$YInt_air, 3),
            Mean_water_C = np_round(m$YInt_water, 3))
  out$year <- as.integer(format(out$date, "%Y"))
  out$DOY <- as.integer(format(out$date, "%j"))
  out$odd_year <- out$year %% 2 == 1
  out$AR_roll7 <- NA_real_
  out$PH_roll7 <- NA_real_
  for (sid in unique(out$site_id)) {
    i <- out$site_id == sid
    out$AR_roll7[i] <- roll_center(out$AmpRatio[i], roll_k)
    out$PH_roll7[i] <- roll_center(out$PhaseLag_h[i], roll_k)
  }
  rownames(out) <- NULL
  out
}
