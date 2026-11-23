# Main entry points: run_folder(), run_loggers(), run_pasta(), run_annual(), run_diel().
#
# The analysis functions take a data.frame whose first four columns are, in order
# (the PASTA R convention): site_id, date/time, stream temperature, air temperature.

#' Analyses
#' @export
ANALYSES <- c("annual", "diel")

standardize <- function(data) {
  if (!is.data.frame(data) || ncol(data) < 4) {
    stop("data must be a data.frame whose first four columns are site_id, date/time, ",
         "stream temperature, air temperature", call. = FALSE)
  }
  d <- df(site_id = as.character(data[[1]]), datetime = as_utc(data[[2]]),
          water = suppressWarnings(as.numeric(data[[3]])),
          air = suppressWarnings(as.numeric(data[[4]])))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d <- d[corder(d$site_id, d$datetime), , drop = FALSE]
  rownames(d) <- NULL
  d
}

is_daily <- function(d) {
  midnight <- all(format(d$datetime, "%H:%M:%S", tz = .tz) == "00:00:00")
  one <- !any(duplicated(d[, c("site_id", "datetime")]))
  midnight && one
}

infer_interval <- function(d) {
  st <- unlist(lapply(split(as.numeric(d$datetime), d$site_id), function(x) diff(x) / 60))
  st <- st[st > 0]
  if (length(st)) as.integer(round(stats::median(st))) else 15L
}

to_paired <- function(d, interval_min, water_min_C) {
  p <- df(site_id = d$site_id, DateTime_r = round_time(d$datetime, interval_min),
          str_temp = d$water, air_temp = d$air)
  p <- p[p$str_temp > water_min_C, , drop = FALSE]
  rownames(p) <- NULL
  p
}

#' Four-column input table from logger series
#' @param loggers Named list from [read_loggers()].
#' @param pairs Pair table from [pair_loggers()] (site_id, stream, air) or a list of
#'   c(site_id, stream, air).
#' @return data.frame site_id, DateTime_r, str_temp, air_temp (inner join on time).
#' @export
pair_data <- function(loggers, pairs) build_paired_df(pairs, loggers, water_min_C = -Inf)

new_result <- function(x, cls) structure(x, class = c(cls, "list"))

#' Diel PASTA
#'
#' Fits `temp ~ sin(2*pi*h/24) + cos(2*pi*h/24)` to air and stream for every site and
#' day, then water/air metrics for days where both fits pass `r2_min`.
#' @param data data.frame, first four columns site_id, date/time, stream, air
#'   (sub-daily records).
#' @param interval_min Logging interval for rounding times; inferred when `NULL`.
#' @param water_min_C Stream records at or below this are dropped (R: str_temp > 0.99).
#' @param min_coverage Share of a full day of records needed to fit a day.
#' @param r2_min Both fits must exceed this adjusted R\eqn{^2}.
#' @param roll_k Centred rolling-mean window (days).
#' @param phase_method `"atan2"` (default; phase_h is the hour of the fitted peak,
#'   lag wrapped into (-12, 12] h) or `"atan"` (the Diel PASTA R formula).
#' @return A `pasta_diel` list: paired, fits, metrics, summary, settings.
#' @export
run_diel <- function(data, interval_min = NULL, water_min_C = 0.99, min_coverage = 0.75,
                     r2_min = 0.7, roll_k = 7, phase_method = DEFAULT_PHASE_METHOD) {
  d <- standardize(data)
  if (is_daily(d)) stop("run_diel() needs sub-daily records; these look like daily means.",
                        call. = FALSE)
  if (is.null(interval_min)) interval_min <- infer_interval(d)
  paired <- to_paired(d, interval_min, water_min_C)
  check_phase_method(phase_method)
  fits <- fit_daily_sinusoids(paired, min_coverage, interval_min, phase_method)
  metrics <- pasta_metrics(fits, r2_min, as.integer(roll_k), phase_method)
  rows <- lapply(unique(paired$site_id), function(sid) {
    m <- if (nrow(metrics)) metrics[metrics$site_id == sid, , drop = FALSE] else metrics
    f <- if (nrow(fits)) fits[fits$site_id == sid, , drop = FALSE] else fits
    med <- function(col) if (nrow(m)) stats::median(m[[col]]) else NaN
    df(site_id = sid, diel_n_days_fitted = if (nrow(f)) length(unique(f$date)) else 0L,
       diel_n_days_kept = nrow(m), diel_median_AmpRatio = med("AmpRatio"),
       diel_median_PhaseLag_h = med("PhaseLag_h"), diel_median_Ratio_Mean_K = med("Ratio_Mean_K"))
  })
  settings <- list(interval_min = interval_min, water_min_C = water_min_C,
                   min_coverage = min_coverage, r2_min = r2_min, roll_k = roll_k,
                   phase_method = phase_method)
  new_result(list(paired = paired, fits = fits, metrics = metrics, summary = rbind_all(rows),
                  settings = settings), "pasta_diel")
}

#' Annual-signal PASTA
#'
#' Fits `temp ~ sin(2*pi*d/365) + cos(2*pi*d/365)` to daily mean air and stream
#' temperature for each site, for the whole record and each year, plus daily thermal
#' sensitivity (stream ~ air) and a data gap check.
#' @param data data.frame, first four columns site_id, date/time, stream, air.
#'   Sub-daily records are averaged to daily means; daily means are used as they are.
#' @param yr_type `"water"` (1 Oct - 30 Sep, named by the year it ends) or `"calendar"`.
#' @param max_gap_days A year with more than this many consecutive no-data days is
#'   dropped (no AmpRatio, PhaseLag_d, Ratio_Mean). Default 50. Counted after the
#'   1 \enc{°}{deg}C drop: a no-data day has no daily mean with stream > 1 \enc{°}{deg}C.
#' @param max_missing_pct A year whose no-data days are more than this % of the year
#'   (365/366) is also dropped. Default 75.
#' @param min_record_pct The whole record needs fit days on at least this % of the
#'   record span (first to last paired day). All fit days are used.
#' @param water_min_C Sub-daily stream records at or below this are dropped before
#'   averaging. Days with a stream mean <= 1 \enc{°}{deg}C are always dropped (as in PASTA).
#' @param min_coverage Share of a full day of paired records needed for a daily mean.
#' @param interval_min Logging interval for sub-daily data; inferred when `NULL`.
#' @param phase_method `"atan2"` (default; phase_d is the day of the fitted peak,
#'   lag wrapped into (-182.5, 182.5] d) or `"atan"` (the PASTA R formula).
#' @return A `pasta_annual` list: daily, record, yearly, fits, settings.
#' @export
run_annual <- function(data, yr_type = "water", max_gap_days = 50, max_missing_pct = 75,
                       min_record_pct = 75,
                       water_min_C = 0.99, min_coverage = 0.75, interval_min = NULL,
                       phase_method = DEFAULT_PHASE_METHOD) {
  if (!yr_type %in% c("water", "calendar")) stop("yr_type must be 'water' or 'calendar'", call. = FALSE)
  d <- standardize(data)
  if (is_daily(d)) {
    daily <- df(site_id = d$site_id, date = day_of(d$datetime), tavg_wat_C = d$water,
                tavg_air_C = d$air)
    interval_min <- NULL
  } else {
    if (is.null(interval_min)) interval_min <- infer_interval(d)
    daily <- daily_means(to_paired(d, interval_min, water_min_C), min_coverage, interval_min)
  }
  check_phase_method(phase_method)
  res <- annual_metrics(daily, yr_type, as.integer(max_gap_days), as.numeric(max_missing_pct),
                        as.numeric(min_record_pct), phase_method)
  settings <- list(yr_type = yr_type, max_gap_days = max_gap_days,
                   max_missing_pct = max_missing_pct,
                   min_record_pct = min_record_pct, water_min_C = water_min_C,
                   min_coverage = min_coverage, interval_min = interval_min,
                   phase_method = phase_method)
  new_result(list(daily = daily, record = res$record, yearly = res$yearly, fits = res$fits,
                  settings = settings), "pasta_annual")
}

.arg_names <- function(fn) setdiff(names(formals(fn)), "data")

#' Run the annual and/or diel analysis on the same data
#'
#' @param data data.frame, first four columns site_id, date/time, stream, air.
#' @param analyses `c("annual", "diel")` by default; one of them for one.
#' @param ... Settings for [run_annual()] / [run_diel()] (each takes the ones it knows).
#' @return A `pasta_result` list: annual, diel, summary, pairs, coverage, paired,
#'   date_range, run_ts, settings.
#' @export
run_pasta <- function(data, analyses = ANALYSES, ...) {
  kw <- list(...)
  analyses <- tolower(analyses)
  if (!length(analyses) || !all(analyses %in% ANALYSES)) {
    stop("analyses must be one or both of ", paste(ANALYSES, collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(names(kw), c(.arg_names(run_annual), .arg_names(run_diel)))
  if (length(unknown)) stop("Unknown settings: ", paste(unknown, collapse = ", "), call. = FALSE)
  pick <- function(fn) kw[names(kw) %in% .arg_names(fn)]
  out <- list(annual = NULL, diel = NULL, summary = NULL, pairs = NULL, coverage = NULL,
              paired = NULL, date_range = NULL, run_ts = Sys.time(), settings = list())
  if ("annual" %in% analyses) out$annual <- do.call(run_annual, c(list(data), pick(run_annual)))
  if ("diel" %in% analyses) out$diel <- do.call(run_diel, c(list(data), pick(run_diel)))
  d <- standardize(data)
  water_min <- if (!is.null(kw$water_min_C)) kw$water_min_C else 0.99
  interval <- kw$interval_min
  for (r in list(out$annual, out$diel)) {
    if (!is.null(r) && !is.null(r$settings$interval_min)) interval <- r$settings$interval_min
  }
  keep <- d[d$water > water_min, , drop = FALSE]
  out$paired <- df(site_id = keep$site_id, DateTime_r = keep$datetime, str_temp = keep$water,
                   air_temp = keep$air)
  s <- list(analyses = analyses, interval_min = interval, water_min_C = water_min)
  for (r in list(out$annual, out$diel)) {
    if (!is.null(r)) for (k in setdiff(names(r$settings), "interval_min")) s[[k]] <- r$settings[[k]]
  }
  out$settings <- s
  out <- new_result(out, "pasta_result")
  out$summary <- run_summary(out, d)
  out
}

#' Run summary, one row per site
#'
#' As in the app: run time, loggers, settings, then the annual whole-record metrics
#' and the diel medians.
#' @param res A `pasta_result`.
#' @param data The input data (or its standardized form).
#' @param pairs Pair table, to add stream_logger and air_logger.
#' @param date_range POSIXct c(start, end), end exclusive.
#' @export
run_summary <- function(res, data = NULL, pairs = NULL, date_range = NULL) {
  s <- res$settings
  d <- if (!is.null(data)) {
    if (all(c("datetime", "water") %in% names(data))) data else standardize(data)
  }
  sites <- if (!is.null(d)) unique(d$site_id) else character(0)
  pm <- NULL
  if (!is.null(pairs)) {
    pm <- pairs
    sites <- c(sites, setdiff(pm$site_id, sites))
  }
  if (!is.null(date_range)) {
    d0 <- date_range[1]
    d1 <- date_range[2] - 86400
  } else if (!is.null(d) && nrow(d)) {
    d0 <- min(d$datetime)
    d1 <- max(d$datetime)
  } else {
    d0 <- d1 <- NA
  }
  n_rec <- if (!is.null(res$paired)) table(res$paired$site_id) else integer(0)
  `%||%` <- function(a, b) if (is.null(a)) b else a
  rows <- lapply(sites, function(sid) {
    i <- if (!is.null(pm)) match(sid, pm$site_id) else NA
    row <- list(
      run_timestamp = if (!is.null(res$run_ts)) format(res$run_ts, "%Y-%m-%dT%H:%M:%S") else "",
      site_id = sid,
      stream_logger = if (!is.na(i)) pm$stream[i] else NA_character_,
      air_logger = if (!is.na(i)) pm$air[i] else NA_character_,
      analyses = paste(s$analyses, collapse = " + "),
      interval_min = s$interval_min %||% NA, water_min_C = s$water_min_C %||% NA,
      min_day_coverage_pct = if (!is.null(s$min_coverage)) 100 * s$min_coverage else NA,
      phase_method = s$phase_method %||% NA,
      date_from = if (!is.na(d0)) format(d0, "%Y-%m-%d") else "",
      date_to = if (!is.na(d1)) format(d1, "%Y-%m-%d") else "",
      n_paired_records = if (sid %in% names(n_rec)) as.integer(n_rec[[sid]]) else 0L)
    if (!is.null(res$annual)) {
      rec <- res$annual$record
      r <- if (nrow(rec)) rec[rec$site_id == sid, , drop = FALSE] else rec
      row$year_type <- s$yr_type
      row$max_gap_days_per_year <- s$max_gap_days
      row$max_missing_pct_per_year <- s$max_missing_pct
      row$min_record_pct <- s$min_record_pct
      row$annual_n_days <- if (nrow(r)) as.integer(r$n_days_fit[1]) else 0L
      for (c in c("AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope", "TS_AdjRsqr",
                  "max_conseq_missing_days", "record_span_days", "pct_record_fit",
                  "meets_record_pct")) {
        row[[paste0("annual_", c)]] <- if (nrow(r)) r[[c]][1] else NA
      }
    }
    if (!is.null(res$diel)) {
      ds <- res$diel$summary
      r <- if (nrow(ds)) ds[ds$site_id == sid, , drop = FALSE] else ds
      row$adjR2_min <- s$r2_min
      row$roll_k <- s$roll_k
      for (c in c("diel_n_days_fitted", "diel_n_days_kept", "diel_median_AmpRatio",
                  "diel_median_PhaseLag_h", "diel_median_Ratio_Mean_K")) {
        row[[c]] <- if (nrow(r)) r[[c]][1] else if (grepl("n_days", c)) 0L else NaN
      }
    }
    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  })
  rbind_all(rows)
}

pair_table <- function(pairs) {
  if (is.data.frame(pairs)) {
    t <- pairs
    if ("include" %in% names(t)) t <- t[as.logical(t$include), setdiff(names(t), "include"), drop = FALSE]
    t <- t[!is.na(t$air), , drop = FALSE]
    rownames(t) <- NULL
    return(t)
  }
  p <- pair_tuples(pairs)
  df(site_id = vapply(p, `[`, "", 1), stream = vapply(p, `[`, "", 2),
     air = vapply(p, `[`, "", 3), basis = "user")
}

#' The app's whole run in one call
#'
#' Logger series and pairs -> data coverage, date range, pairing and the annual
#' and/or diel analysis.
#' @param loggers Named list from [read_loggers()].
#' @param pairs Pair table (site_id, stream, air, optional include) or list of
#'   c(site_id, stream, air); suggested by [pair_loggers()] when `NULL`.
#' @param analyses `c("annual", "diel")` or one of them.
#' @param date_range `"all"` (each site's full record), `"common"` (only dates every
#'   site has paired data) or c(start, end) dates, inclusive.
#' @param interval_min Logging interval in minutes.
#' @param ... Settings for [run_annual()] / [run_diel()]: water_min_C, min_coverage,
#'   phase_method, yr_type, max_gap_days, max_missing_pct, min_record_pct, r2_min, roll_k.
#' @return A `pasta_result` with pairs, coverage (per site, in the range) and date_range.
#' @export
run_loggers <- function(loggers, pairs = NULL, analyses = ANALYSES, date_range = "all",
                        interval_min = 15, ...) {
  pt <- if (is.null(pairs)) pair_loggers(loggers, interval_min = interval_min) else pair_table(pairs)
  pt <- pt[!is.na(pt$air), , drop = FALSE]
  if (!nrow(pt)) stop("No air-stream pairs to run.", call. = FALSE)
  if (anyDuplicated(pt$site_id)) {
    stop("Site IDs must be unique: ", paste(unique(pt$site_id[duplicated(pt$site_id)]), collapse = ", "),
         call. = FALSE)
  }
  cov <- coverage(loggers, pt, interval_min)
  rng <- date_range(cov, date_range)   # R skips the non-function argument of the same name
  data <- pair_data(clip_loggers(loggers, rng[1], rng[2]), pt)
  res <- run_pasta(data, analyses = analyses, interval_min = interval_min, ...)
  res$pairs <- pt
  res$coverage <- coverage_in_range(cov, rng[1], rng[2])
  res$date_range <- rng
  res$summary <- run_summary(res, data, pairs = pt, date_range = rng)
  res
}

#' Read a folder of logger CSVs and run PASTA
#'
#' [read_loggers()] + [run_loggers()].
#' @param folder Folder with air and stream logger CSVs.
#' @param analyses `c("annual", "diel")` or one of them.
#' @param pairs Optional pairs to override the suggestions, e.g.
#'   `list(c("HF-003", "HF-003", "HF_103"))`.
#' @param interval_min Logging interval in minutes.
#' @param out_dir If given, every table is written there as `PASTA_<table>.csv`.
#' @param date_range `"all"`, `"common"` or c(start, end); see [run_loggers()].
#' @param ... Settings for [run_annual()] / [run_diel()].
#' @export
run_folder <- function(folder, analyses = ANALYSES, pairs = NULL, interval_min = 15,
                       out_dir = NULL, date_range = "all", ...) {
  loggers <- read_loggers(folder, interval_min = interval_min)
  res <- run_loggers(loggers, pairs = pairs, analyses = analyses, date_range = date_range,
                     interval_min = interval_min, ...)
  if (!is.null(out_dir)) write_pasta(res, out_dir)
  res
}

#' @export
print.pasta_result <- function(x, ...) {
  cat("PASTA result:", paste(x$settings$analyses, collapse = " + "), "\n")
  if (!is.null(x$date_range)) {
    cat("Date range:", format(x$date_range[1], "%Y-%m-%d"), "to",
        format(x$date_range[2] - 86400, "%Y-%m-%d"), "\n")
  }
  cat("Tables:", paste(names(pasta_tables(x)), collapse = ", "), "\n\n")
  print(x$summary[, intersect(c("site_id", "annual_AmpRatio", "annual_PhaseLag_d",
                                "annual_TS__Slope", "diel_median_AmpRatio",
                                "diel_median_PhaseLag_h"), names(x$summary)), drop = FALSE],
        row.names = FALSE)
  invisible(x)
}

#' @export
print.pasta_annual <- function(x, ...) {
  cat("Annual PASTA (", x$settings$yr_type, " year)\n", sep = "")
  print(x$record, row.names = FALSE)
  invisible(x)
}

#' @export
print.pasta_diel <- function(x, ...) {
  cat("Diel PASTA\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}
