# Sorting loggers into air vs stream and suggesting air-stream pairs.

# Name tokens (bounded by non-letters):
#   air:    air, atm, AT (air temperature)
#   stream: stream, str, ST (stream temperature), water, wat, WT, SW, river
.air_rx <- "(?i)(^|[^a-z])(air|atm|atmo|at)([^a-z]|$)"
.wat_rx <- "(?i)(^|[^a-z])(str|stream|st|wat|water|wt|sw|riv|river)([^a-z]|$)"

has_air_name <- function(x) grepl(.air_rx, x, perl = TRUE)
has_wat_name <- function(x) grepl(.wat_rx, x, perl = TRUE)

#' Features used to tell air from stream loggers
#' @param g Logger data.frame(timestamp, temp_C).
#' @return Named list: start, end, n_obs, mean_C, min_C, max_C,
#'   median_daily_range_C, pct_below_0, median_step_C.
#' @export
logger_features <- function(g) {
  g <- g[!is.na(g$temp_C), , drop = FALSE]
  day <- day_of(g$timestamp)
  cnt <- tapply(g$temp_C, day, length)
  full <- names(cnt)[cnt >= 0.75 * max(cnt)]
  rng <- (tapply(g$temp_C, day, max) - tapply(g$temp_C, day, min))[full]
  step <- abs(diff(g$temp_C))
  list(start = min(g$timestamp), end = max(g$timestamp), n_obs = nrow(g),
       mean_C = mean(g$temp_C), min_C = min(g$temp_C), max_C = max(g$temp_C),
       median_daily_range_C = if (length(rng)) stats::median(rng) else NaN,
       pct_below_0 = mean(g$temp_C < 0) * 100,
       median_step_C = if (length(step)) stats::median(step) else NaN)
}

#' Classify one logger as air or stream
#'
#' Stream water is buffered (small daily range, cannot sit far below 0 \enc{°}{deg}C); air has a
#' large daily range, noisy steps, and often freezes. A name containing an air or
#' stream code wins over the data.
#' @param name Logger name.
#' @param feats Output of [logger_features()].
#' @param range_threshold Median daily range (\enc{°}{deg}C) above which a logger looks like air.
#' @return Character vector: type, confidence, reason.
#' @export
classify_logger <- function(name, feats, range_threshold = 4) {
  if (has_air_name(name)) return(c("air", "high", "name contains an air keyword"))
  if (has_wat_name(name)) return(c("stream", "high", "name contains a water keyword"))
  rng <- feats$median_daily_range_C
  votes <- 0
  reasons <- character(0)
  if (isTRUE(rng > range_threshold)) {
    votes <- votes + 2
    reasons <- c(reasons, sprintf("daily range %.1f \u00b0C > %.1f", rng, range_threshold))
  } else {
    reasons <- c(reasons, sprintf("daily range %.1f \u00b0C \u2264 %.1f", rng, range_threshold))
  }
  if (isTRUE(feats$min_C < -1)) {
    votes <- votes + 1
    reasons <- c(reasons, sprintf("min %.1f \u00b0C (below freezing)", feats$min_C))
  }
  if (isTRUE(feats$median_step_C > 0.1)) {
    votes <- votes + 1
    reasons <- c(reasons, "noisy 15-min steps")
  }
  kind <- if (votes >= 2) "air" else "stream"
  margin <- abs(log(max(rng, 1e-3) / range_threshold))
  conf <- if (isTRUE(margin > 0.7)) "high" else if (isTRUE(margin > 0.3)) "medium" else "low"
  c(kind, conf, paste(reasons, collapse = "; "))
}

#' Sort every logger into air or stream
#' @param loggers Named list from [read_loggers()].
#' @inheritParams classify_logger
#' @return data.frame, one row per logger: logger, type, confidence, reason and the
#'   features used.
#' @export
classify_loggers <- function(loggers, range_threshold = 4) {
  rows <- lapply(names(loggers), function(n) {
    f <- logger_features(loggers[[n]])
    k <- classify_logger(n, f, range_threshold)
    df(logger = n, type = k[1], confidence = k[2], reason = k[3],
       start = f$start, end = f$end, n_obs = f$n_obs, mean_C = f$mean_C, min_C = f$min_C,
       max_C = f$max_C, median_daily_range_C = f$median_daily_range_C,
       pct_below_0 = f$pct_below_0, median_step_C = f$median_step_C)
  })
  rbind_all(rows)
}

# ---- pairing -------------------------------------------------------------------

digits_of <- function(x) regmatches(x, gregexpr("[0-9]+", x))[[1]]

type_stem <- function(name) {
  n <- gsub(.air_rx, "\\1\\3", name, perl = TRUE)
  n <- gsub(.wat_rx, "\\1\\3", n, perl = TRUE)
  n <- gsub("[_\\-\\s]+", "_", n, perl = TRUE)
  tolower(gsub("^_+|_+$", "", n))
}

letters_of <- function(stem) gsub("[^a-z]", "", stem)

#' Default site ID from a stream logger name
#'
#' `"LLC_LLC01_WT_corrected"` becomes `"LLC_LLC01"`, `"ST_002"` becomes `"site_002"`,
#' `"HF-003"` stays `"HF-003"`.
#' @param stream_name Stream logger name.
#' @export
site_label <- function(stream_name) {
  n <- gsub(.wat_rx, "\\1\\3", stream_name, perl = TRUE)
  n <- gsub("(?i)(^|[_\\-\\s])(corrected|raw|qc|final|clean|data)(?=$|[_\\-\\s])", "\\1", n,
            perl = TRUE)
  n <- gsub("[_\\-\\s]{2,}", "_", n, perl = TRUE)
  n <- gsub("^[_\\- ]+|[_\\- ]+$", "", n, perl = TRUE)
  if (grepl("^[0-9]+$", n)) n <- paste0("site_", n)
  if (nzchar(n)) n else stream_name
}

#' Name-match tiers, strongest first
#' @export
NAME_TIERS <- c(
  "1" = "same name except air/stream code",   # AT_002 <-> ST_002, LLC01_AT <-> LLC01_WT
  "2" = "same site name and number",          # Air_Brook_2 <-> Str_Brook_02
  "3" = "same number",                        # Air_2 <-> Brook_02
  "4" = "same last two digits",               # HF_003 <-> HF_103
  "5" = "same prefix"                         # AMTH_A01 <-> AMTH_air
)

#' Strongest name rule linking a stream and an air logger
#' @param stream_name,air_name Logger names.
#' @return Tier number 1-5 (see [NAME_TIERS]) or `NA`.
#' @export
name_match_level <- function(stream_name, air_name) {
  ss <- type_stem(stream_name)
  sa <- type_stem(air_name)
  typed <- has_air_name(air_name) || has_wat_name(stream_name)
  if (typed && nzchar(ss) && ss == sa) return(1L)
  sd <- digits_of(stream_name)
  ad <- digits_of(air_name)
  if (length(sd) && length(ad)) {
    ls <- sd[length(sd)]
    la <- ad[length(ad)]
    if (as.numeric(ls) == as.numeric(la)) {
      return(if (letters_of(ss) == letters_of(sa)) 2L else 3L)
    }
    if (nchar(ls) >= 2 && nchar(la) >= 2 &&
        substr(ls, nchar(ls) - 1, nchar(ls)) == substr(la, nchar(la) - 1, nchar(la))) {
      return(4L)
    }
  }
  s_pre <- tolower(strsplit(stream_name, "[_\\-\\s]", perl = TRUE)[[1]][1])
  a_pre <- tolower(strsplit(air_name, "[_\\-\\s]", perl = TRUE)[[1]][1])
  if (has_air_name(air_name) && identical(s_pre, a_pre) && nchar(s_pre) > 2) return(5L)
  NA_integer_
}

#' Diel cross-correlation of a stream and an air logger
#'
#' Maximum correlation of the 24 h detrended signals with the stream lagging air
#' by 0 to `max_lag_h` hours.
#' @param stream,air Logger data.frames.
#' @param interval_min Logging interval in minutes.
#' @param max_lag_h Largest lag tried, hours.
#' @return c(r, lag_h, overlap_days).
#' @export
diel_xcorr <- function(stream, air, interval_min = 15, max_lag_h = 6) {
  m <- merge(stream, air, by = "timestamp", suffixes = c("_w", "_a"))
  m <- m[stats::complete.cases(m), , drop = FALSE]
  m <- m[order(m$timestamp), , drop = FALSE]
  per_day <- round(24 * 60 / interval_min)
  if (nrow(m) < 3 * per_day) return(c(NaN, NaN, nrow(m) / per_day))
  w <- m$temp_C_w - roll_center(m$temp_C_w, per_day, per_day %/% 2)
  a <- m$temp_C_a - roll_center(m$temp_C_a, per_day, per_day %/% 2)
  best <- c(-1, 0)
  n <- length(w)
  for (k in 0:floor(max_lag_h * 60 / interval_min)) {
    ww <- if (k) w[(k + 1):n] else w
    aa <- if (k) a[1:(n - k)] else a
    ok <- is.finite(ww) & is.finite(aa)
    if (sum(ok) < per_day) next
    r <- stats::cor(ww[ok], aa[ok])
    if (is.finite(r) && r > best[1]) best <- c(r, k * interval_min / 60)
  }
  c(best, nrow(m) / per_day)
}

.score <- function(score, s, a) {
  v <- score[paste(s, a, sep = "\r")]
  if (length(v) == 0 || is.na(v)) NaN else unname(v)
}

#' Suggest one air logger for every stream logger (names first)
#'
#' 1. Names first: each stream takes the air logger(s) at its strongest name tier;
#'    several candidates are separated by the best diel correlation.
#' 2. Data last: streams with no name match get the best-scoring free air logger,
#'    then the best-scoring air logger overall.
#' @param streams,airs Logger names.
#' @param score Named numeric of diel r; each name is the stream and air logger
#'   names joined by a carriage return (`paste(stream, air, sep = "\\\\r")`).
#' @return Named list stream -> list(air, basis).
#' @export
suggest_pairs <- function(streams, airs, score) {
  out <- list()
  used <- character(0)
  key <- function(r) ifelse(is.finite(r), r, -9)
  for (s in streams) {
    lv <- vapply(airs, function(a) name_match_level(s, a), integer(1))
    lv <- lv[!is.na(lv)]
    if (!length(lv)) next
    best <- min(lv)
    hits <- names(lv)[lv == best]
    basis <- paste0("name: ", NAME_TIERS[[as.character(best)]])
    if (length(hits) > 1) {
      sc <- vapply(hits, function(a) .score(score, s, a), numeric(1))
      hits <- hits[order(-key(sc), method = "radix")]
      basis <- paste0(basis, " (", length(hits), " candidates; best diel r)")
    }
    out[[s]] <- list(air = hits[1], basis = basis)
    used <- c(used, hits[1])
  }
  remaining <- streams[!streams %in% names(out)]
  if (length(remaining) && length(airs)) {
    cs <- expand.grid(a = airs, s = remaining, stringsAsFactors = FALSE)[, c("s", "a")]
    cs$r <- mapply(function(s, a) .score(score, s, a), cs$s, cs$a)
    cs <- cs[order(-key(cs$r), method = "radix"), , drop = FALSE]
    for (i in seq_len(nrow(cs))) {
      s <- cs$s[i]
      a <- cs$a[i]
      if (s %in% names(out) || a %in% used || !is.finite(cs$r[i])) next
      out[[s]] <- list(air = a, basis = "data: best diel correlation")
      used <- c(used, a)
    }
  }
  for (s in remaining) {
    if (s %in% names(out)) next
    r <- vapply(airs, function(a) .score(score, s, a), numeric(1))
    ok <- is.finite(r)
    if (any(ok)) {
      top <- airs[ok][r[ok] == max(r[ok])]
      out[[s]] <- list(air = utils::tail(csort(top), 1),
                       basis = "data: best diel correlation (shared air)")
    } else {
      out[[s]] <- list(air = NA_character_, basis = "no name match and no overlapping data")
    }
  }
  out
}

#' Basis note for the pair table
#'
#' The basis, plus a 'confirm' flag when a data-based pick is uncertain.
#' @param stream,air,basis Stream, chosen air logger and basis from [suggest_pairs()].
#' @param airs Air logger names.
#' @param score Named numeric of diel r, as in [suggest_pairs()].
#' @export
basis_note <- function(stream, air, basis, airs, score) {
  if (startsWith(basis, "name") || is.na(air)) return(basis)
  rs <- vapply(airs, function(a) .score(score, stream, a), numeric(1))
  rs <- sort(rs[is.finite(rs)], decreasing = TRUE)
  if (length(rs) > 1 && rs[1] - rs[2] < 0.02) return(paste0(basis, " - confirm (air loggers score alike)"))
  if (length(rs) && rs[1] < 0.3) return(paste0(basis, " - confirm (weak diel link)"))
  basis
}

#' Suggest air-stream pairs
#'
#' @param loggers Named list from [read_loggers()].
#' @param types Named character `c(name = "air" | "stream" | "exclude")`; default
#'   from [classify_loggers()].
#' @param interval_min Logging interval in minutes.
#' @return data.frame site_id, stream, air, basis, diel_r, lag_h, overlap_days.
#'   'basis' says 'confirm ...' when the data cannot tell the air loggers apart.
#' @export
pair_loggers <- function(loggers, types = NULL, interval_min = 15) {
  if (is.null(types)) {
    ct <- classify_loggers(loggers)
    types <- stats::setNames(ct$type, ct$logger)
  }
  airs <- names(types)[types == "air"]
  streams <- names(types)[types == "stream"]
  sc <- list()
  for (s in streams) for (a in airs) {
    sc[[paste(s, a, sep = "\r")]] <- diel_xcorr(loggers[[s]], loggers[[a]], interval_min)
  }
  r_only <- vapply(sc, function(v) v[1], numeric(1))
  sugg <- suggest_pairs(streams, airs, r_only)
  rows <- lapply(streams, function(s) {
    a <- sugg[[s]]$air
    v <- if (!is.na(a)) sc[[paste(s, a, sep = "\r")]] else c(NaN, NaN, NaN)
    df(site_id = site_label(s), stream = s, air = a,
       basis = basis_note(s, a, sugg[[s]]$basis, airs, r_only),
       diel_r = v[1], lag_h = v[2], overlap_days = v[3])
  })
  out <- rbind_all(rows)
  if (!nrow(out)) {
    out <- df(site_id = character(0), stream = character(0), air = character(0),
              basis = character(0), diel_r = numeric(0), lag_h = numeric(0),
              overlap_days = numeric(0))
  }
  dup <- duplicated(out$site_id) | duplicated(out$site_id, fromLast = TRUE)
  out$site_id[dup] <- out$stream[dup]
  attr(out, "scores") <- sc
  out
}
