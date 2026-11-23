# Reading logger files: HOBOware exports and simple date/time + temperature CSVs.
# A logger is a data.frame(timestamp = POSIXct "UTC", temp_C) on the logging grid.

.date_formats <- c(
  "%m/%d/%Y %H:%M:%S", "%m/%d/%y %H:%M:%S", "%m/%d/%Y %H:%M", "%m/%d/%y %H:%M",
  "%m/%d/%y %I:%M:%S %p", "%m/%d/%Y %I:%M:%S %p",
  "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%dT%H:%M:%S"
)

# full-match regex for a strptime format, so "%H:%M" does not accept "10:15:00 PM"
# (R's strptime ignores trailing text; Python's does not)
.fmt_regex <- function(fmt) {
  map <- c("%m" = "\\d{1,2}", "%d" = "\\d{1,2}", "%Y" = "\\d{4}", "%y" = "\\d{2}",
           "%H" = "\\d{1,2}", "%I" = "\\d{1,2}", "%M" = "\\d{1,2}", "%S" = "\\d{1,2}",
           "%p" = "(AM|PM|am|pm)")
  rx <- gsub("([.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1", fmt, perl = TRUE)
  for (k in names(map)) rx <- gsub(k, map[[k]], rx, fixed = TRUE)
  paste0("^", rx, "$")
}

parse_datetimes <- function(s) {
  s <- trimws(as.character(s))
  sample <- utils::head(s[!is.na(s) & nzchar(s)], 50)
  for (fmt in .date_formats) {
    if (!length(sample) || !all(grepl(.fmt_regex(fmt), sample, perl = TRUE))) next
    p <- as.POSIXct(strptime(sample, fmt, tz = .tz))
    if (anyNA(p)) next
    out <- as.POSIXct(strptime(s, fmt, tz = .tz))
    out[!grepl(.fmt_regex(fmt), s, perl = TRUE)] <- NA
    return(out)
  }
  out <- suppressWarnings(tryCatch(as_utc(s), error = function(e) NULL))
  if (is.null(out)) out <- as.POSIXct(rep(NA_real_, length(s)), origin = "1970-01-01", tz = .tz)
  out
}

#' Logger name from a file name
#'
#' Drops the HOBOware date stamp: `"HF_005 2026-07-22 12_44_20 EDT.csv"` becomes `"HF_005"`.
#' @param fname File name or path.
#' @export
logger_name_from_filename <- function(fname) {
  stem <- sub("\\.[^.]*$", "", basename(fname))
  out <- trimws(strsplit(stem, "\\s+\\d{4}-\\d{2}-\\d{2}", perl = TRUE)[[1]][1])
  if (is.na(out) || !nzchar(out)) stem else out
}

#' Parse one logger CSV
#'
#' Handles HOBOware exports (optional 'Plot Title' line, '#' index column,
#' 'Date-Time (EDT)' and 'Temperature (\enc{°}{deg}C/\enc{°}{deg}F)' columns, event-only rows)
#' and simple CSVs with a date/time column and a temperature column, or one numeric
#' column named after the logger (e.g. `Timestamp (UTC-05:00),LLC_LLC01_AT`).
#' Degrees F are converted to C.
#'
#' @param path Path to the CSV, or the file's lines when `text = TRUE`.
#' @param fname File name used in messages.
#' @param text `TRUE` when `path` holds the lines rather than a path.
#' @return data.frame(timestamp, temp_C), sorted, one row per time, with attribute
#'   `tz_label` (e.g. "EDT" or "UTC-05:00").
#' @export
parse_logger_csv <- function(path, fname = basename(path), text = FALSE) {
  lines <- if (text) path else readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines)) lines[1] <- sub("^\ufeff", "", lines[1])
  low <- tolower(utils::head(lines, 30))
  hit <- which((grepl("date", low, fixed = TRUE) | grepl("time", low, fixed = TRUE)) &
                 (grepl("temp", low, fixed = TRUE) | grepl("\u00b0c", low, fixed = TRUE) |
                    grepl("value", low, fixed = TRUE)))
  header_idx <- if (length(hit)) hit[1] else 1
  d <- utils::read.csv(text = lines[header_idx:length(lines)], check.names = FALSE,
                       colClasses = "character", strip.white = TRUE, fill = TRUE)
  cols <- names(d)
  lc <- tolower(cols)
  ts_col <- cols[grepl("date", lc, fixed = TRUE) | grepl("time", lc, fixed = TRUE)][1]
  t_col <- cols[grepl("temp", lc, fixed = TRUE)][1]
  if (is.na(t_col)) t_col <- cols[grepl("value", lc, fixed = TRUE)][1]
  if (is.na(t_col)) {
    num <- cols[vapply(cols, function(c) {
      !identical(c, ts_col) && trimws(c) != "#" &&
        mean(!is.na(suppressWarnings(as.numeric(d[[c]])))) > 0.5
    }, logical(1))]
    if (length(num) == 1) t_col <- num
  }
  if (is.na(ts_col) || is.na(t_col)) {
    stop("Could not find a date/time column and a temperature column in ", fname,
         ". Columns: ", paste(cols, collapse = ", "), call. = FALSE)
  }
  out <- df(timestamp = parse_datetimes(d[[ts_col]]),
            temp_C = suppressWarnings(as.numeric(d[[t_col]])))
  tl <- tolower(t_col)
  if (grepl("\u00b0f", tl, fixed = TRUE) || grepl("\\(\\s*f\\s*\\)|deg\\s*f", tl, perl = TRUE)) {
    out$temp_C <- (out$temp_C - 32) * 5 / 9
  }
  out <- out[stats::complete.cases(out), , drop = FALSE]
  out <- out[order(out$timestamp), , drop = FALSE]
  out <- out[!duplicated(out$timestamp), , drop = FALSE]
  rownames(out) <- NULL
  if (nrow(out) < 10) stop("Fewer than 10 valid temperature records in ", fname, call. = FALSE)
  m <- regmatches(ts_col, regexec("\\(([A-Z]{2,5}|(?:GMT|UTC)[^)]*)\\)", ts_col, perl = TRUE))[[1]]
  attr(out, "tz_label") <- if (length(m) > 1) m[2] else ""
  out
}

#' Round to the logging interval and average duplicates
#'
#' @param x data.frame(timestamp, temp_C) from [parse_logger_csv()].
#' @param interval_min Logging interval in minutes.
#' @return data.frame(timestamp, temp_C) on the grid (R: `round_date()`).
#' @export
to_grid <- function(x, interval_min = 15) {
  t <- round_time(x$timestamp, interval_min)
  key <- as.numeric(t)
  m <- tapply(x$temp_C, key, mean)
  out <- df(timestamp = as.POSIXct(as.numeric(names(m)), origin = "1970-01-01", tz = .tz),
            temp_C = as.numeric(m))
  out[order(out$timestamp), , drop = FALSE]
}

#' Read one logger CSV onto the logging grid
#' @inheritParams to_grid
#' @param path Path to the CSV.
#' @export
read_logger <- function(path, interval_min = 15) {
  to_grid(parse_logger_csv(path, basename(path)), interval_min)
}

#' Read several logger CSVs
#'
#' @param source A folder (every .csv in it is read) or a vector of file paths.
#' @param interval_min Logging interval in minutes.
#' @param on_error `"warn"` (default) skips files that are not logger files with a
#'   warning, as the app does; `"stop"` stops at the first bad file.
#' @return Named list of loggers, `name = data.frame(timestamp, temp_C)`. Names come
#'   from the file names. Skipped files are in `attr(x, "skipped")`.
#' @export
read_loggers <- function(source, interval_min = 15, on_error = c("warn", "stop")) {
  on_error <- match.arg(on_error)
  if (length(source) == 1 && dir.exists(source)) {
    paths <- list.files(source, pattern = "\\.csv$", ignore.case = TRUE, full.names = TRUE)
    paths <- paths[order(basename(paths), method = "radix")]
    if (!length(paths)) stop("No .csv files found in ", source, call. = FALSE)
  } else {
    paths <- source
  }
  out <- list()
  skipped <- character(0)
  for (p in paths) {
    name <- logger_name_from_filename(p)
    if (name %in% names(out)) name <- sub("\\.[^.]*$", "", basename(p))
    g <- tryCatch(read_logger(p, interval_min), error = function(e) e)
    if (inherits(g, "error")) {
      if (on_error == "stop") stop(g)
      skipped[basename(p)] <- conditionMessage(g)
      warning("Skipped ", basename(p), ": ", conditionMessage(g), call. = FALSE)
      next
    }
    out[[name]] <- g
  }
  if (!length(out)) {
    stop("None of the files could be read as logger files: ", paste(skipped, collapse = "; "),
         call. = FALSE)
  }
  attr(out, "skipped") <- skipped
  out
}
