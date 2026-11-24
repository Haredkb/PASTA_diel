# Annual signal (from the PASTA Shiny app: R/paired_annualT_signals.R,
# R/stream_thermal_sensitivity.R and therm_analysis() in R/rundata_v1.R).

#' Daily means of paired records
#' @param T_y Paired records (site_id, DateTime_r, str_temp, air_temp).
#' @param min_coverage Share of a full day of paired records needed for a daily mean.
#' @param interval_min Logging interval in minutes.
#' @return data.frame site_id, date, tavg_wat_C, tavg_air_C, n_obs.
#' @export
daily_means <- function(T_y, min_coverage = 0.75, interval_min = 15) {
  if (!nrow(T_y)) {
    return(df(site_id = character(0), date = as.Date(character(0)), tavg_wat_C = numeric(0),
              tavg_air_C = numeric(0), n_obs = integer(0)))
  }
  date <- day_of(T_y$DateTime_r)
  key <- paste(T_y$site_id, as.numeric(date), sep = "\r")
  first <- !duplicated(key)
  n <- tapply(T_y$str_temp, key, length)
  out <- df(site_id = T_y$site_id[first], date = date[first])
  k <- key[first]
  out$tavg_wat_C <- as.numeric(tapply(T_y$str_temp, key, mean)[k])
  out$tavg_air_C <- as.numeric(tapply(T_y$air_temp, key, mean)[k])
  out$n_obs <- as.integer(n[k])
  out <- out[corder(out$site_id, out$date), , drop = FALSE]
  out <- out[out$n_obs >= min_coverage * 24 * 60 / interval_min, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Water year (1 Oct - 30 Sep, named by the year it ends)
#' @param dates Dates.
#' @export
water_year <- function(dates) {
  y <- as.integer(format(as.Date(dates), "%Y"))
  m <- as.integer(format(as.Date(dates), "%m"))
  ifelse(m > 9, y + 1L, y)
}

#' Radian day 2*pi*d/365
#' @param dates Dates.
#' @param yr_type `"water"` (d = day of water year, 1 Oct = 1) or `"calendar"`.
#' @export
rad_day <- function(dates, yr_type = "water") {
  dates <- as.Date(dates)
  if (yr_type == "calendar") {
    d <- as.numeric(format(dates, "%j"))
  } else {
    start <- as.Date(paste0(water_year(dates) - 1, "-09-30"))
    d <- as.numeric(dates - start)
  }
  2 * pi * d / 365
}

ols <- function(X, y) {
  f <- stats::lm.fit(X, y)
  b <- f$coefficients
  n <- nrow(X)
  p <- ncol(X)
  sse <- sum(f$residuals^2)
  sst <- sum((y - mean(y))^2)
  r2 <- if (sst > 0) 1 - sse / sst else NaN
  adj <- if (n > p) 1 - (1 - r2) * (n - 1) / (n - p) else NaN
  list(b = unname(b), adj = adj, rmse = sqrt(sse / n))
}

#' Annual sinusoid fit
#'
#' `temp ~ sin(rday) + cos(rday)`; phase_d is the day of the fitted peak (`"atan2"`)
#' or the R formula (`"atan"`).
#' @param dates Dates.
#' @param temp Daily mean temperatures.
#' @param yr_type `"water"` or `"calendar"`.
#' @param phase_method `"atan2"` or `"atan"`, see [PHASE_METHODS].
#' @return List amplitude_C, phase_d, AdjRsqr, RMSE, sinSlope, cosSlope, YInt, or
#'   `NULL` with fewer than 4 days.
#' @export
fit_TAS <- function(dates, temp, yr_type = "water", phase_method = DEFAULT_PHASE_METHOD) {
  y <- as.numeric(temp)
  rd <- rad_day(dates, yr_type)
  ok <- is.finite(y)
  y <- y[ok]
  rd <- rd[ok]
  if (length(y) < 4) return(NULL)
  f <- ols(cbind(1, sin(rd), cos(rd)), y)
  b <- f$b
  list(amplitude_C = sqrt(b[2]^2 + b[3]^2), phase_d = pasta_phase(b[2], b[3], 365, phase_method),
       AdjRsqr = f$adj, RMSE = f$rmse, sinSlope = b[2], cosSlope = b[3], YInt = b[1])
}

#' Daily linear thermal sensitivity, water ~ air (water > 1 \enc{°}{deg}C)
#' @param a_temp,w_temp Daily mean air and water temperatures.
#' @export
fit_ThermalSens <- function(a_temp, w_temp) {
  a <- as.numeric(a_temp)
  w <- as.numeric(w_temp)
  ok <- is.finite(a) & is.finite(w) & w > 1
  if (sum(ok) < 3) return(list(TS__Slope = NaN, TS_AdjRsqr = NaN, TS_RMSE = NaN, TS_YInt = NaN))
  f <- ols(cbind(1, a[ok]), w[ok])
  list(TS__Slope = round(f$b[2], 2), TS_AdjRsqr = round(f$adj, 2), TS_RMSE = f$rmse,
       TS_YInt = round(f$b[1], 2))
}

#' Fitted annual sinusoid on given dates
#' @param fit One annual fit (list or data.frame row).
#' @param dates Dates.
#' @param yr_type `"water"` or `"calendar"`.
#' @export
annual_sinusoid_value <- function(fit, dates, yr_type = "water") {
  rd <- rad_day(dates, yr_type)
  fit$sinSlope * sin(rd) + fit$cosSlope * cos(rd) + fit$YInt
}

#' Data gap check per site
#' @param daily Daily means.
#' @return data.frame site_id, max_conseq_missing_days, count.
#' @export
data_gap_check <- function(daily) {
  d <- daily[!is.na(daily$tavg_wat_C) & !is.na(daily$tavg_air_C) & daily$tavg_wat_C < 70, ,
             drop = FALSE]
  rows <- lapply(unique(d$site_id), function(sid) {
    x <- sort(d$date[d$site_id == sid])
    dd <- as.numeric(diff(x))
    df(site_id = sid, max_conseq_missing_days = if (length(dd)) max(dd) else NaN,
       count = length(x))
  })
  out <- rbind_all(rows)
  if (!nrow(out)) out <- df(site_id = character(0), max_conseq_missing_days = numeric(0),
                            count = integer(0))
  out
}

#' Whole-record annual metrics per site
#' @param daily Daily means.
#' @inheritParams fit_TAS
#' @return list(metrics, fits); metrics: site_id, AmpRatio, PhaseLag_d, Ratio_Mean,
#'   TS__Slope, TS_AdjRsqr, TS_YInt.
#' @export
therm_analysis <- function(daily, yr_type = "water", phase_method = DEFAULT_PHASE_METHOD) {
  d <- daily[daily$tavg_wat_C > 1, , drop = FALSE]
  mets <- list()
  fits <- list()
  for (sid in unique(d$site_id)) {
    x <- d[d$site_id == sid, , drop = FALSE]
    fa <- fit_TAS(x$date, x$tavg_air_C, yr_type, phase_method)
    fw <- fit_TAS(x$date, x$tavg_wat_C, yr_type, phase_method)
    ts <- fit_ThermalSens(x$tavg_air_C, x$tavg_wat_C)
    for (med in c("air", "water")) {
      f <- if (med == "air") fa else fw
      if (!is.null(f)) {
        fits[[length(fits) + 1]] <- df(as.data.frame(f), medium = med, site_id = sid,
                                       n_days = nrow(x), yr_type = yr_type,
                                       phase_method = phase_method)
      }
    }
    if (!is.null(fa) && !is.null(fw)) {
      m <- df(site_id = sid, AmpRatio = round(fw$amplitude_C / fa$amplitude_C, 2),
              PhaseLag_d = round(pasta_lag(fw$phase_d, fa$phase_d, 365, phase_method), 2),
              Ratio_Mean = round(fw$YInt / fa$YInt, 2))
    } else {
      m <- df(site_id = sid, AmpRatio = NaN, PhaseLag_d = NaN, Ratio_Mean = NaN)
    }
    m$TS__Slope <- ts$TS__Slope
    m$TS_AdjRsqr <- ts$TS_AdjRsqr
    m$TS_YInt <- ts$TS_YInt
    mets[[length(mets) + 1]] <- m
  }
  metrics <- rbind_all(mets)
  if (!nrow(metrics)) {
    metrics <- df(site_id = character(0), AmpRatio = numeric(0), PhaseLag_d = numeric(0),
                  Ratio_Mean = numeric(0), TS__Slope = numeric(0), TS_AdjRsqr = numeric(0),
                  TS_YInt = numeric(0))
  }
  list(metrics = metrics, fits = rbind_all(fits))
}

fit_days <- function(daily) {
  d <- daily[!is.na(daily$tavg_wat_C) & !is.na(daily$tavg_air_C), , drop = FALSE]
  d[d$tavg_wat_C > 1 & d$tavg_wat_C < 70, , drop = FALSE]
}

.signal_cols <- c("AmpRatio", "PhaseLag_d", "Ratio_Mean")

#' Days in a water year (1 Oct of year-1 to 30 Sep) or calendar year: 365 or 366
#' @param year Year(s), named by the year they end.
#' @param yr_type `"water"` or `"calendar"`.
#' @export
days_in_year <- function(year, yr_type = "water") {
  y <- as.integer(year)
  if (yr_type == "water") {
    as.integer(as.Date(paste0(y, "-10-01")) - as.Date(paste0(y - 1, "-10-01")))
  } else {
    as.integer(as.Date(paste0(y + 1, "-01-01")) - as.Date(paste0(y, "-01-01")))
  }
}

#' Fit days as a share of a record span
#'
#' @param all_dates All paired days (cold days included) of the record or year.
#' @param fit_dates Days with a stream daily mean > 1 \enc{°}{deg}C (after the 1 \enc{°}{deg}C drop).
#' @return list(n_days_fit, span_days, pct_fit), span = first to last paired day.
#' @export
span_coverage <- function(all_dates, fit_dates) {
  if (!length(all_dates)) return(list(n_days_fit = 0L, span_days = 0L, pct_fit = 0))
  span <- as.integer(as.numeric(max(all_dates) - min(all_dates)) + 1)
  n <- length(fit_dates)
  list(n_days_fit = n, span_days = span, pct_fit = 100 * n / span)
}

#' Annual metrics for the whole record and each year, with the data rules
#'
#' Both rules are applied after the 1 \enc{°}{deg}C drop and count fit days: paired daily means
#' with stream > 1 \enc{°}{deg}C.
#' * yearly: `pct_year_fit = 100 * n_days_fit / n_days_in_year`, out of the full water
#'   (or calendar) year, 365 or 366 days. Below `min_year_pct` the year is dropped
#'   (`meets_year_pct = FALSE`): no AmpRatio, PhaseLag_d or Ratio_Mean; thermal
#'   sensitivity is still reported. A partial first or last year counts the days
#'   outside the record as missing.
#' * record: `pct_record_fit = 100 * n_days_fit / record_span_days`, the span from the
#'   first to the last paired day (cold days included). Below
#'   `min_record_pct` the signal metrics are NA (`meets_record_pct = FALSE`).
#' @param daily Daily means.
#' @param min_year_pct,min_record_pct The rules above (default 75 each).
#' @inheritParams fit_TAS
#' @return list(record, yearly, fits).
#' @export
annual_metrics <- function(daily, yr_type = "water", min_year_pct = 75, min_record_pct = 75,
                           phase_method = DEFAULT_PHASE_METHOD) {
  if (!nrow(daily)) return(list(record = data.frame(), yearly = data.frame(), fits = data.frame()))
  d <- daily
  d$date <- as.Date(d$date)
  d$year <- if (yr_type == "water") water_year(d$date) else as.integer(format(d$date, "%Y"))
  fd <- fit_days(d)
  all_by <- split(d$date, paste(d$site_id, d$year, sep = "\r"))
  fit_by <- split(fd$date, paste(fd$site_id, fd$year, sep = "\r"))
  year_check <- function(sid, yr) {
    k <- paste(sid, yr, sep = "\r")
    n <- if (k %in% names(fit_by)) length(fit_by[[k]]) else 0L
    tot <- days_in_year(yr, yr_type)
    pct <- 100 * n / tot
    list(n_days_fit = n, n_days_in_year = tot, pct_year_fit = round(pct, 1),
         meets_year_pct = pct >= min_year_pct)
  }

  # per year: fit every year, report signal metrics only for years that qualify
  yr_rows <- list()
  yr_fits <- list()
  for (yr in sort(unique(d$year))) {
    x <- d[d$year == yr, , drop = FALSE]
    ta <- therm_analysis(x, yr_type, phase_method)
    m <- merge(ta$metrics, data_gap_check(x), by = "site_id", all.x = TRUE, sort = FALSE)
    m <- m[match(ta$metrics$site_id, m$site_id), , drop = FALSE]
    if (nrow(m)) m$year <- as.integer(yr)
    f <- ta$fits
    if (nrow(f)) f$year <- as.character(yr)
    yr_rows[[length(yr_rows) + 1]] <- m
    yr_fits[[length(yr_fits) + 1]] <- f
  }
  yearly <- rbind_all(yr_rows)
  if (nrow(yearly)) {
    chk <- do.call(rbind, lapply(seq_len(nrow(yearly)), function(i) {
      as.data.frame(year_check(yearly$site_id[i], yearly$year[i]))
    }))
    yearly <- cbind(yearly, chk)
    yearly[!yearly$meets_year_pct, .signal_cols] <- NA_real_
    yearly$year_type <- yr_type
    yearly <- yearly[corder(yearly$site_id, yearly$year), , drop = FALSE]
    rownames(yearly) <- NULL
  }
  ok_years <- names(all_by)[vapply(names(all_by), function(k) {
    p <- strsplit(k, "\r", fixed = TRUE)[[1]]
    year_check(p[1], as.integer(p[2]))$meets_year_pct
  }, logical(1))]

  # whole record: all fit days; signal metrics only if they cover enough of the span
  ta <- therm_analysis(d, yr_type, phase_method)
  rec <- merge(ta$metrics, data_gap_check(d), by = "site_id", all.x = TRUE, sort = FALSE)
  rows <- lapply(unique(d$site_id), function(sid) {
    sc <- span_coverage(d$date[d$site_id == sid], fd$date[fd$site_id == sid])
    df(site_id = sid, n_days_fit = sc$n_days_fit, record_span_days = sc$span_days,
       pct_record_fit = round(sc$pct_fit, 1), meets_record_pct = sc$pct_fit >= min_record_pct)
  })
  base <- rbind_all(rows)
  rec <- merge(base, rec, by = "site_id", all.x = TRUE, sort = FALSE)
  rec <- rec[match(base$site_id, rec$site_id), , drop = FALSE]
  rec[!rec$meets_record_pct, .signal_cols] <- NA_real_
  cols <- c("site_id", "AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope", "TS_AdjRsqr",
            "TS_YInt", "max_conseq_missing_days", "count", "n_days_fit",
            "record_span_days", "pct_record_fit", "meets_record_pct")
  for (c in setdiff(cols, names(rec))) rec[[c]] <- NA
  rec <- rec[, cols]
  rownames(rec) <- NULL

  rf <- ta$fits
  if (nrow(rf)) {
    rf$year <- "all"
    rf$meets_criteria <- rec$meets_record_pct[match(rf$site_id, rec$site_id)] %in% TRUE
  }
  yf <- rbind_all(yr_fits)
  if (nrow(yf)) yf$meets_criteria <- paste(yf$site_id, yf$year, sep = "\r") %in% ok_years
  list(record = rec, yearly = yearly, fits = rbind_all(list(rf, yf)))
}
