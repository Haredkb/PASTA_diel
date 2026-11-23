# Small helpers shared by the package. Times are POSIXct in "UTC" so that the
# logger clock time is kept exactly as written (no time-zone or DST shifts).

.tz <- "UTC"

# numpy / pandas .round(): scale, round half to even, unscale
np_round <- function(x, digits) round(x * 10^digits) / 10^digits

# pandas dt.round(): nearest multiple of the interval, ties to even
round_time <- function(t, interval_min) {
  step <- interval_min * 60
  as.POSIXct(round(as.numeric(t) / step) * step, origin = "1970-01-01", tz = .tz)
}

as_utc <- function(x) {
  if (inherits(x, "POSIXct")) return(as.POSIXct(format(x, "%Y-%m-%d %H:%M:%S"), tz = .tz))
  if (inherits(x, "Date")) return(as.POSIXct(format(x), tz = .tz))
  as.POSIXct(as.character(x), tz = .tz, tryFormats = c("%Y-%m-%d %H:%M:%OS", "%Y/%m/%d %H:%M:%OS",
                                                         "%Y-%m-%d %H:%M", "%Y/%m/%d %H:%M",
                                                         "%Y-%m-%d", "%Y/%m/%d"))
}

day_of <- function(t) as.Date(format(t, "%Y-%m-%d", tz = .tz))

# sort strings like Python (code-point order), not by locale
csort <- function(x) sort(x, method = "radix")
corder <- function(...) order(..., method = "radix")

# centred rolling mean over rows, as pandas rolling(k, center=TRUE, min_periods=mp)
roll_center <- function(x, k, min_periods = k) {
  n <- length(x)
  if (!n) return(numeric(0))
  ok <- is.finite(x)
  cs <- c(0, cumsum(ifelse(ok, x, 0)))
  cn <- c(0, cumsum(ok))
  lo_raw <- seq_len(n) - k %/% 2
  hi <- pmin(lo_raw + k - 1, n)
  lo <- pmax(lo_raw, 1)
  cnt <- cn[hi + 1] - cn[lo]
  out <- (cs[hi + 1] - cs[lo]) / cnt
  out[cnt < min_periods | cnt == 0] <- NA_real_
  out
}

# data.frame with no row names, strings kept
df <- function(...) data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)

rbind_all <- function(lst) {
  lst <- Filter(function(x) !is.null(x) && nrow(x) > 0, lst)
  if (!length(lst)) return(data.frame())
  out <- do.call(rbind, lst)
  rownames(out) <- NULL
  out
}
