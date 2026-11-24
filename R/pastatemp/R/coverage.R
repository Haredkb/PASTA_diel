# Data coverage and date range: when each logger and each air-stream pair has data,
# the period shared by all sites, and clipping loggers to a range.

pair_tuples <- function(pairs) {
  if (is.data.frame(pairs)) {
    if ("include" %in% names(pairs)) pairs <- pairs[as.logical(pairs$include), , drop = FALSE]
    pairs <- pairs[!is.na(pairs$air), , drop = FALSE]
    return(lapply(seq_len(nrow(pairs)), function(i) {
      c(as.character(pairs$site_id[i]), as.character(pairs$stream[i]), as.character(pairs$air[i]))
    }))
  }
  if (is.matrix(pairs)) return(lapply(seq_len(nrow(pairs)), function(i) as.character(pairs[i, ])))
  lapply(pairs, as.character)
}

#' Data coverage per site
#'
#' Start and end of each logger and of the paired record (both present).
#' @param loggers Named list from [read_loggers()].
#' @param pairs Pair table (site_id, stream, air) or list of c(site_id, stream, air).
#' @param interval_min Logging interval in minutes.
#' @return data.frame site_id, stream, air, stream_start, stream_end, air_start,
#'   air_end, paired_start, paired_end, paired_days, and a list column paired_index
#'   (the paired times; dropped by [coverage_in_range()]).
#' @export
coverage <- function(loggers, pairs, interval_min = 15) {
  per_day <- 24 * 60 / interval_min
  na_t <- as.POSIXct(NA_real_, origin = "1970-01-01", tz = .tz)
  rows <- lapply(pair_tuples(pairs), function(p) {
    gs <- loggers[[p[2]]]
    ga <- loggers[[p[3]]]
    gs <- gs[!is.na(gs$temp_C), , drop = FALSE]
    ga <- ga[!is.na(ga$temp_C), , drop = FALSE]
    both <- sort(gs$timestamp[gs$timestamp %in% ga$timestamp])
    r <- df(site_id = p[1], stream = p[2], air = p[3],
            stream_start = min(gs$timestamp), stream_end = max(gs$timestamp),
            air_start = min(ga$timestamp), air_end = max(ga$timestamp),
            paired_start = if (length(both)) min(both) else na_t,
            paired_end = if (length(both)) max(both) else na_t,
            paired_days = length(both) / per_day)
    r$paired_index <- I(list(both))
    r
  })
  rbind_all(rows)
}

#' Period shared by every site
#' @param cov Output of [coverage()].
#' @return POSIXct c(start, end), or `NULL` when the sites share no period.
#' @export
common_period <- function(cov) {
  c0 <- max(cov$paired_start)
  c1 <- min(cov$paired_end)
  if (is.na(c0) || is.na(c1) || c0 >= c1) return(NULL)
  c(c0, c1)
}

day_floor <- function(t) as.POSIXct(format(t, "%Y-%m-%d", tz = .tz), tz = .tz)

#' Analysis window as whole days
#'
#' @param cov Output of [coverage()].
#' @param mode `"all"` (first to last paired day of any site; each site uses its
#'   full record), `"common"` (only days every site has paired data), or two dates
#'   `c(start, end)`, inclusive.
#' @return POSIXct c(start, end), end exclusive.
#' @export
date_range <- function(cov, mode = "all") {
  if (length(mode) == 2) {
    m <- as_utc(mode)
    return(c(day_floor(m[1]), day_floor(m[2]) + 86400))
  }
  if (identical(mode, "all")) {
    c0 <- min(cov$paired_start)
    c1 <- max(cov$paired_end)
  } else if (identical(mode, "common")) {
    cp <- common_period(cov)
    if (is.null(cp)) stop("No period is shared by all sites; use date_range = 'all' or your own dates.",
                          call. = FALSE)
    c0 <- cp[1]
    c1 <- cp[2]
  } else {
    stop("date_range must be 'all', 'common' or c(start, end)", call. = FALSE)
  }
  c(day_floor(c0), day_floor(c1) + 86400)
}

#' Days from start to end of the longer of each site's stream and air records
#' @param cov Output of [coverage()].
#' @return Named numeric, by site_id.
#' @export
logger_span_days <- function(cov) {
  dd <- function(a, b) ifelse(is.na(a) | is.na(b), 0,
                              as.numeric(as.Date(b) - as.Date(a)) + 1)
  stats::setNames(pmax(dd(cov$stream_start, cov$stream_end), dd(cov$air_start, cov$air_end)),
                  cov$site_id)
}

#' Record span per site
#'
#' Used for '% of range' and for the annual whole-record rule (record_span_days).
#' @param cov Output of [coverage()].
#' @param start,end Range, end exclusive.
#' @param mode `"all"`: start to end of the longer of the site's stream and air
#'   records; otherwise the chosen range, end - start.
#' @return Named numeric, by site_id.
#' @export
record_spans <- function(cov, start, end, mode = "all") {
  if (is.character(mode) && length(mode) == 1 && mode == "all") return(logger_span_days(cov))
  n <- max(as.numeric(difftime(as_utc(end), as_utc(start), units = "days")), 1)
  stats::setNames(rep(n, nrow(cov)), cov$site_id)
}

#' Coverage of each site within a date range
#' @param cov Output of [coverage()].
#' @param start,end Range, end exclusive (e.g. from [date_range()]).
#' @param spans Optional named numeric from [record_spans()]: each site's range in
#'   days; default end - start.
#' @return The coverage table (without paired_days) plus days_in_range (paired days
#'   in the range), range_days, pct_of_range and status ('covers whole range',
#'   'starts ...', 'ends ...', 'gaps inside range', ...).
#' @export
coverage_in_range <- function(cov, start, end, spans = NULL) {
  out <- cov[, setdiff(names(cov), c("paired_index", "paired_days")), drop = FALSE]
  n_all <- max(as.numeric(difftime(end, start, units = "days")), 1)
  tol <- 86400
  days <- integer(nrow(cov))
  rng <- numeric(nrow(cov))
  status <- character(nrow(cov))
  for (i in seq_len(nrow(cov))) {
    sid <- cov$site_id[i]
    n_days <- if (!is.null(spans) && sid %in% names(spans)) max(spans[[sid]], 1) else n_all
    rng[i] <- n_days
    idx <- cov$paired_index[[i]]
    k <- idx[idx >= start & idx < end]
    days[i] <- length(unique(day_of(k)))
    ps <- cov$paired_start[i]
    pe <- cov$paired_end[i]
    status[i] <- if (days[i] == 0) {
      "no data in range"
    } else if (ps > start + tol && pe < end - 2 * tol) {
      "starts late and ends early"
    } else if (ps > start + tol) {
      paste("starts", format(ps, "%Y-%m-%d"))
    } else if (pe < end - 2 * tol) {
      paste("ends", format(pe, "%Y-%m-%d"))
    } else if (days[i] < 0.9 * n_days) {
      "gaps inside range"
    } else {
      "covers whole range"
    }
  }
  out$days_in_range <- days
  out$range_days <- rng
  out$pct_of_range <- 100 * days / rng
  out$status <- status
  out
}

#' TRUE when every site covers the whole range
#' @param cov_range Output of [coverage_in_range()].
#' @export
all_sites_cover <- function(cov_range) {
  nrow(cov_range) > 0 && all(cov_range$status == "covers whole range")
}

#' Keep records with start <= time < end in every logger
#' @param loggers Named list of loggers.
#' @param start,end POSIXct or date strings; `NULL` for no limit.
#' @export
clip_loggers <- function(loggers, start = NULL, end = NULL) {
  s <- if (!is.null(start)) as_utc(start)
  e <- if (!is.null(end)) as_utc(end)
  lapply(loggers, function(g) {
    keep <- rep(TRUE, nrow(g))
    if (!is.null(s)) keep <- keep & g$timestamp >= s
    if (!is.null(e)) keep <- keep & g$timestamp < e
    g[keep, , drop = FALSE]
  })
}
