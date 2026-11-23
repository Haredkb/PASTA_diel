#' pastatemp: Paired Air and Stream Temperature Analysis (PASTA)
#'
#' The R twin of the 'pastatemp' Python package and its app: the same steps,
#' settings and output tables.
#'
#' | App step | Function |
#' |---|---|
#' | 1 Logger files | [read_loggers()] |
#' | 2 Air or stream? | [classify_loggers()] |
#' | 3 Pairs (names first) | [pair_loggers()] |
#' | 4 Date range and coverage | [coverage()], [common_period()], [date_range()], [coverage_in_range()], [clip_loggers()] |
#' | 4 Run | [run_loggers()] (all of step 4 in one call), [run_folder()] (steps 1-4) |
#' | 5 Results and figures | `res$summary`, `res$annual`, `res$diel`; [pasta_plots] |
#' | 6 Save | [write_pasta()], [zip_pasta()], [save_tables()] |
#'
#' The analysis functions [run_annual()], [run_diel()] and [run_pasta()] also take a
#' data.frame whose first four columns are site_id, date/time, stream and air
#' temperature (the PASTA R convention).
#'
#' Times are kept as the logger clock shows them (POSIXct in "UTC"), so there is no
#' time-zone or daylight-saving shift.
#'
#' @examples
#' \dontrun{
#' res <- run_folder("path/to/logger_csvs", date_range = "common")
#' res$summary
#' zip_pasta(res, "PASTA_results.zip")
#'
#' loggers <- read_loggers("path/to/logger_csvs")
#' pairs <- list(c("HF-003", "HF-003", "HF_103"), c("HF_005", "HF_005", "HF_101"))
#' res <- run_loggers(loggers, pairs, date_range = "common", max_gap_days = 50,
#'                    min_record_pct = 75)
#' plot_annual_fit(res$annual, "HF-003")
#' }
#' @keywords internal
"_PACKAGE"
