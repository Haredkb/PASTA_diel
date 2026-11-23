# Saving results: all tables to a folder or into one zip, as the app's Download step.

#' The result tables, as the app downloads them
#'
#' run_summary, pairs, coverage, paired_df, annual_metrics, annual_metrics_yearly,
#' annual_fits, daily_means, diel_metrics, diel_fits (missing ones are skipped).
#' @param res A `pasta_result`, `pasta_annual` or `pasta_diel`.
#' @export
pasta_tables <- function(res) {
  t <- list()
  if (inherits(res, "pasta_annual")) {
    return(list(annual_metrics = res$record, annual_metrics_yearly = res$yearly,
                annual_fits = res$fits, daily_means = res$daily))
  }
  if (inherits(res, "pasta_diel")) {
    return(list(diel_metrics = res$metrics, diel_fits = res$fits, diel_summary = res$summary))
  }
  t$run_summary <- res$summary
  t$pairs <- res$pairs
  t$coverage <- res$coverage
  t$paired_df <- res$paired
  if (!is.null(res$annual)) t <- c(t, pasta_tables(res$annual))
  if (!is.null(res$diel)) t <- c(t, list(diel_metrics = res$diel$metrics, diel_fits = res$diel$fits))
  Filter(Negate(is.null), t)
}

#' File-name stem of the app: "PASTA_YYYYmmdd_HHMM"
#' @param res A `pasta_result`.
#' @export
default_prefix <- function(res) {
  ts <- if (!is.null(res$run_ts)) res$run_ts else Sys.time()
  paste0("PASTA_", format(ts, "%Y%m%d_%H%M"))
}

#' Clean a typed or pasted folder path
#'
#' Accepts quotes (Explorer's 'Copy as path'), ~ and environment variables.
#' @param raw Folder path as typed.
#' @export
clean_out_dir <- function(raw) {
  txt <- trimws(gsub("^[\"']|[\"']$", "", trimws(as.character(raw))))
  # %VAR% (Windows) and $VAR / ${VAR}
  for (v in regmatches(txt, gregexpr("%[A-Za-z_][A-Za-z0-9_]*%", txt))[[1]]) {
    val <- Sys.getenv(gsub("%", "", v), unset = NA)
    if (!is.na(val)) txt <- sub(v, val, txt, fixed = TRUE)
  }
  for (v in regmatches(txt, gregexpr("\\$\\{?[A-Za-z_][A-Za-z0-9_]*\\}?", txt))[[1]]) {
    val <- Sys.getenv(gsub("[${}]", "", v), unset = NA)
    if (!is.na(val)) txt <- sub(v, val, txt, fixed = TRUE)
  }
  if (!nzchar(txt)) stop("Type a folder path, e.g. C:\\Users\\you\\Documents\\PASTA_output", call. = FALSE)
  if (.Platform$OS.type != "windows" && grepl("^[A-Za-z]:", txt)) {
    stop("This is running on a server, not on your computer, so it cannot save to ", txt,
         ". Download the zip instead.", call. = FALSE)
  }
  path <- path.expand(txt)
  if (.Platform$OS.type == "windows" && grepl("^[A-Za-z]:", path) &&
      !dir.exists(paste0(substr(path, 1, 2), "/"))) {
    stop("Drive ", substr(path, 1, 2), " was not found.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_table <- function(x, path) {
  x <- as.data.frame(x)
  for (n in names(x)) if (is.list(x[[n]])) x[[n]] <- NULL
  for (n in names(x)) {
    if (inherits(x[[n]], "POSIXct")) {
      # as pandas: date only when every time in the column is midnight
      midnight <- all(is.na(x[[n]]) | format(x[[n]], "%H:%M:%S", tz = .tz) == "00:00:00")
      x[[n]] <- format(x[[n]], if (midnight) "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S", tz = .tz)
    }
  }
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

#' Write every table as `<prefix>_<table>.csv`
#' @param tables Named list of data.frames (e.g. [pasta_tables()]).
#' @param folder Output folder (created if needed; typed paths are cleaned).
#' @param prefix File-name prefix.
#' @return The written paths.
#' @export
save_tables <- function(tables, folder, prefix = "PASTA") {
  target <- clean_out_dir(folder)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  paths <- character(0)
  for (n in names(tables)) {
    if (is.null(tables[[n]])) next
    p <- file.path(target, paste0(prefix, "_", n, ".csv"))
    write_table(tables[[n]], p)
    paths <- c(paths, p)
  }
  paths
}

#' Write a result's tables to a folder
#' @param res A `pasta_result`, `pasta_annual` or `pasta_diel`.
#' @param folder Output folder.
#' @param prefix File-name prefix.
#' @export
write_pasta <- function(res, folder, prefix = "PASTA") save_tables(pasta_tables(res), folder, prefix)

#' Zip every table (the app's "Download all")
#' @param tables Named list of data.frames.
#' @param path Zip file to write.
#' @param prefix File-name prefix inside the zip.
#' @return `path`, invisibly. Uses the 'zip' package when installed, else `utils::zip()`.
#' @export
tables_to_zip <- function(tables, path, prefix = "PASTA") {
  tmp <- tempfile("pasta_zip")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  files <- basename(save_tables(tables, tmp, prefix))
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) file.remove(path)
  if (requireNamespace("zip", quietly = TRUE)) {
    zip::zip(path, files, root = tmp)
  } else {
    old <- setwd(tmp)
    on.exit(setwd(old), add = TRUE)
    utils::zip(path, files, flags = "-q")
  }
  invisible(path)
}

#' Zip a result's tables
#' @param res A `pasta_result`.
#' @param path Zip file to write; default `<prefix>.zip` in the working folder.
#' @param prefix File-name prefix; default the app's "PASTA_YYYYmmdd_HHMM".
#' @export
zip_pasta <- function(res, path = NULL, prefix = default_prefix(res)) {
  if (is.null(path)) path <- paste0(prefix, ".zip")
  tables_to_zip(pasta_tables(res), path, prefix)
}
