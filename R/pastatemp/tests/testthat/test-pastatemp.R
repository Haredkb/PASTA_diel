dat <- synthetic()

test_that("run_diel recovers ratio and lag", {
  res <- run_diel(dat, water_min_C = -99)
  m <- res$metrics
  expect_gt(nrow(m), 600)
  expect_equal(median(m$AmpRatio), DIEL_AMP_WAT / DIEL_AMP_AIR, tolerance = 0.02)
  expect_equal(median(m$PhaseLag_h), DIEL_LAG_H, tolerance = 0.1)
  expect_equal(res$settings$interval_min, 15)
})

test_that("run_annual recovers ratio and lag", {
  res <- run_annual(dat, yr_type = "water", water_min_C = -99)
  rec <- res$record[1, ]
  expect_equal(rec$AmpRatio, ANN_AMP_WAT / ANN_AMP_AIR, tolerance = 0.02)
  expect_equal(rec$PhaseLag_d, ANN_LAG_D, tolerance = 1)
  expect_equal(res$yearly$year, c(2025L, 2026L))
  expect_true(all(is.finite(res$yearly$AmpRatio)))
})

test_that("run_annual accepts daily means; run_diel refuses them", {
  dd <- aggregate(cbind(stream, air) ~ site + day, transform(dat, day = as.Date(time)), mean)
  res <- run_annual(dd)
  expect_null(res$settings$interval_min)
  expect_equal(res$record$AmpRatio[1], 0.5, tolerance = 0.02)
  expect_error(run_diel(dd), "sub-daily")
})

test_that("a partial year is judged against the full year", {
  short <- dat[dat$time < as.POSIXct("2024-12-31", tz = "UTC"), ]   # 91 days, all with data
  res <- run_annual(short)
  y <- res$yearly[1, ]
  expect_equal(c(y$n_days_in_year, y$pct_year_fit), c(365, 24.9))   # 91 / 365
  expect_false(y$meets_year_pct)
  expect_true(is.na(y$AmpRatio))
  expect_true(res$record$meets_record_pct[1])
  expect_equal(res$record$pct_record_fit[1], 100)
})

test_that("leap water years have 366 days", {
  expect_equal(days_in_year(c(2024, 2025), "water"), c(366L, 365L))
  expect_equal(days_in_year(2024, "calendar"), 366L)
})

wy2025 <- function(gap_days = 0, cold_days = 0) {
  days <- seq(as.Date("2024-10-01"), as.Date("2025-09-30"), by = 1)
  rd <- 2 * pi * (seq_along(days) - 1) / 365
  w <- 10 + 4 * sin(rd)
  if (cold_days) w[(201 - cold_days):200] <- 0.5
  d <- data.frame(s = "X", date = days, w = w, a = 10 + 8 * sin(rd))
  if (gap_days) d <- d[-(201:(200 + gap_days)), ]
  d
}

test_that("a year needs a minimum % of the year", {
  y <- run_annual(wy2025(91))$yearly[1, ]                    # 274 / 365 = 75.1 %
  expect_equal(c(y$n_days_in_year, y$n_days_fit, y$pct_year_fit), c(365, 274, 75.1))
  expect_true(y$meets_year_pct && is.finite(y$AmpRatio))
  y <- run_annual(wy2025(92))$yearly[1, ]                    # 273 / 365 = 74.8 %
  expect_false(y$meets_year_pct)
  expect_true(is.na(y$AmpRatio) && is.finite(y$TS__Slope))
  expect_true(run_annual(wy2025(92), min_year_pct = 70)$yearly$meets_year_pct[1])
})

test_that("year % counts days after the 1 C drop", {
  y <- run_annual(wy2025(50, cold_days = 50))$yearly[1, ]    # 265 / 365 = 72.6 %
  expect_equal(c(y$count, y$n_days_fit, y$pct_year_fit), c(315, 265, 72.6))
  expect_false(y$meets_year_pct)
})

test_that("record % counts days after the 1 C drop", {
  days <- seq(as.Date("2024-10-01"), by = 1, length.out = 200)
  rd <- 2 * pi * (0:119) / 365
  d <- data.frame(s = "X", date = days, w = c(rep(0.5, 80), 8 + 5 * sin(rd)),
                  a = c(rep(-5, 80), 10 + 10 * sin(rd)))
  res <- run_annual(d, min_year_pct = 30)                    # year: 120 / 365 = 32.9 %
  rec <- res$record[1, ]
  expect_equal(c(rec$count, rec$n_days_fit, rec$record_span_days), c(200, 120, 200))
  expect_equal(rec$pct_record_fit, 60)
  expect_false(rec$meets_record_pct)
  expect_true(is.na(rec$AmpRatio) && is.finite(rec$TS__Slope))
  expect_true(res$yearly$meets_year_pct[1])
  expect_true(run_annual(d, min_record_pct = 60)$record$meets_record_pct[1])
  expect_false(any(res$fits$meets_criteria[res$fits$year == "all"]))
})

test_that("run_pasta runs both or one, and rejects unknown settings", {
  both <- run_pasta(dat, water_min_C = -99)
  expect_false(is.null(both$annual) || is.null(both$diel))
  expect_true(all(c("annual_AmpRatio", "diel_median_AmpRatio") %in% names(both$summary)))
  one <- run_pasta(dat, analyses = "annual", water_min_C = -99)
  expect_null(one$diel)
  expect_error(run_pasta(dat, not_a_setting = 100), "Unknown settings")
})

test_that("name tiers", {
  cases <- list(c("ST_002", "AT_002", 1), c("LLC_LLC01_WT_corrected", "LLC_LLC01_AT_corrected", 1),
                c("Str_Brook_02", "Air_Brook_2", 2), c("SiteB_1_WT", "SiteA_1_AT", 3),
                c("HF-003", "HF_103", 4), c("ST_002", "AT_003", NA), c("HF_005", "HF_101", NA))
  for (k in cases) expect_identical(name_match_level(k[1], k[2]), as.integer(k[3]), info = k[1])
})

test_that("names beat data; data is the last resort", {
  sc <- c("ST_001\rAT_002" = 0.9, "ST_001\rAT_001" = 0.5, "ST_002\rAT_001" = 0.9, "ST_002\rAT_002" = 0.5)
  out <- suggest_pairs(c("ST_001", "ST_002"), c("AT_001", "AT_002"), sc)
  expect_equal(out$ST_001$air, "AT_001")
  expect_equal(out$ST_002$air, "AT_002")
  expect_true(startsWith(out$ST_001$basis, "name"))
  out <- suggest_pairs("Brook", c("AT_001", "AT_002"), c("Brook\rAT_002" = 0.8, "Brook\rAT_001" = 0.4))
  expect_equal(out$Brook, list(air = "AT_002", basis = "data: best diel correlation"))
})

test_that("site labels", {
  expect_equal(site_label("LLC_LLC01_WT_corrected"), "LLC_LLC01")
  expect_equal(site_label("ST_002"), "site_002")
  expect_equal(site_label("HF-003"), "HF-003")
})

test_that("phase helpers and atan2 default", {
  w <- 2 * pi * 150 / 365
  expect_equal(pasta_phase(sin(w), cos(w), 365, "atan2"), 150)
  expect_equal(pasta_phase(sin(w), cos(w), 365, "atan"), 332.5)
  w <- 2 * pi * 15 / 24
  expect_equal(pasta_phase(sin(w), cos(w), 24, "atan2"), 15)
  expect_equal(pasta_phase(sin(w), cos(w), 24, "atan"), 9)
  expect_equal(pasta_lag(1, 23, 24, "atan2"), 2)
  expect_equal(pasta_lag(1, 23, 24, "atan"), -22)
  expect_error(pasta_phase(1, 1, 24, "bad"))
  expect_equal(formals(run_annual)$phase_method, quote(DEFAULT_PHASE_METHOD))
  expect_equal(DEFAULT_PHASE_METHOD, "atan2")
})

phase_case <- function(peak_air_h, peak_wat_h, ann_air, ann_wat) {
  t <- seq(as.POSIXct("2023-10-01", tz = "UTC"), as.POSIXct("2025-09-30 23:45", tz = "UTC"), by = 900)
  doy <- as.POSIXlt(t)$yday + 1
  h <- as.POSIXlt(t)$hour + as.POSIXlt(t)$min / 60
  data.frame(s = "X", t = t,
             w = 12 + 6 * cos(2 * pi * (doy - ann_wat) / 365) + cos(2 * pi * (h - peak_wat_h) / 24),
             a = 10 + 12 * cos(2 * pi * (doy - ann_air) / 365) + 5 * cos(2 * pi * (h - peak_air_h) / 24))
}

test_that("atan2 fixes the late diel peak and annual peaks straddling July", {
  d <- phase_case(15, 1, 200, 210)
  expect_equal(median(run_diel(d, water_min_C = -99)$metrics$PhaseLag_h), 10, tolerance = 0.1)
  expect_equal(median(run_diel(d, water_min_C = -99, phase_method = "atan")$metrics$PhaseLag_h), -2,
               tolerance = 0.1)
  d <- phase_case(15, 17, 175, 190)
  a2 <- run_annual(d, yr_type = "calendar", water_min_C = -99)$record
  a1 <- run_annual(d, yr_type = "calendar", water_min_C = -99, phase_method = "atan")$record
  expect_equal(a2$PhaseLag_d, 15, tolerance = 1)
  expect_gt(abs(a1$PhaseLag_d - 15), 100)
  expect_equal(a1$AmpRatio, a2$AmpRatio)
})

test_that("coverage and date range", {
  lg <- two_site_loggers()
  pairs <- list(c("S1", "ST_001", "AT_001"), c("S2", "ST_002", "AT_002"))
  cov <- coverage(lg, pairs)
  expect_equal(common_period(cov)[1], as.POSIXct("2025-01-01", tz = "UTC"))
  r <- date_range(cov, "all")
  expect_equal(r, as.POSIXct(c("2024-10-01", "2025-10-01"), tz = "UTC"))
  inr <- coverage_in_range(cov, r[1], r[2])
  expect_false(all_sites_cover(inr))
  expect_true(startsWith(inr$status[inr$site_id == "S2"], "starts"))
  rc <- date_range(cov, "common")
  expect_true(all_sites_cover(coverage_in_range(cov, rc[1], rc[2])))
  expect_equal(date_range(cov, c("2025-02-01", "2025-02-28")),
               as.POSIXct(c("2025-02-01", "2025-03-01"), tz = "UTC"))
})

test_that("run_loggers: common period, summary and outputs", {
  lg <- two_site_loggers()
  res <- run_loggers(lg, date_range = "common", water_min_C = -99)
  expect_equal(res$pairs$site_id, c("site_001", "site_002"))
  expect_equal(res$date_range[1], as.POSIXct("2025-01-01", tz = "UTC"))
  expect_true(min(res$paired$DateTime_r) >= as.POSIXct("2025-01-01", tz = "UTC"))
  s <- res$summary
  expect_equal(s$date_from[s$site_id == "site_001"], "2025-01-01")
  expect_equal(s$stream_logger[s$site_id == "site_001"], "ST_001")
  expect_true(all(c("annual_AmpRatio", "annual_pct_record_fit", "diel_median_AmpRatio",
                    "min_year_pct") %in% names(s)))
  expect_equal(names(pasta_tables(res)),
               c("run_summary", "pairs", "coverage", "paired_df", "annual_metrics",
                 "annual_metrics_yearly", "annual_fits", "daily_means", "diel_metrics", "diel_fits"))
  tmp <- withr::local_tempdir()
  z <- zip_pasta(res, file.path(tmp, "r.zip"), prefix = "X")
  expect_true(file.exists(z))
  expect_length(save_tables(pasta_tables(res), paste0('"', file.path(tmp, "out"), '"'), "X"), 10)
})

test_that("reading HOBO-style and one-column logger files, folder run", {
  tmp <- withr::local_tempdir()
  d <- synthetic(start = "2025-04-01 00:00", end = "2025-06-30 23:45")
  stamp <- format(d$time, "%m/%d/%Y %H:%M")
  write.csv(data.frame(`Timestamp (UTC-05:00)` = stamp, SITE_A01_AT = d$air, check.names = FALSE),
            file.path(tmp, "SITE_A01_AT_corrected.csv"), row.names = FALSE)
  write.csv(data.frame(`Timestamp (UTC-05:00)` = stamp, SITE_A01_WT = d$stream, check.names = FALSE),
            file.path(tmp, "SITE_A01_WT_corrected.csv"), row.names = FALSE)
  writeLines(c("a,b", "1,2"), file.path(tmp, "notes.csv"))
  expect_warning(lg <- read_loggers(tmp), "Skipped notes.csv")
  expect_equal(names(attr(lg, "skipped")), "notes.csv")
  ct <- classify_loggers(lg)
  expect_equal(setNames(ct$type, ct$logger),
               c(SITE_A01_AT_corrected = "air", SITE_A01_WT_corrected = "stream"))
  out <- file.path(tmp, "out")
  res <- suppressWarnings(run_folder(tmp, out_dir = out))
  expect_equal(res$pairs$site_id, "SITE_A01")
  expect_equal(res$pairs$basis, "name: same name except air/stream code")
  expect_equal(median(res$diel$metrics$AmpRatio), 0.2, tolerance = 0.02)
  expect_true(all(c("PASTA_run_summary.csv", "PASTA_diel_metrics.csv") %in% list.files(out)))
})

test_that("HOBOware export with title line and degrees F", {
  tmp <- withr::local_tempfile(fileext = ".csv")
  t <- format(seq(as.POSIXct("2025-07-01", tz = "UTC"), by = 900, length.out = 50), "%m/%d/%y %I:%M:%S %p")
  writeLines(c('"Plot Title: HF_005"',
               paste0('"#","Date Time, GMT-04:00","Temp, °F (LGR S/N: 1)"'),
               paste0(seq_along(t), ",", t, ",", 68)), tmp, useBytes = FALSE)
  x <- parse_logger_csv(tmp)
  expect_equal(nrow(x), 50)
  expect_equal(x$temp_C[1], 20)
  expect_equal(format(x$timestamp[2], "%H:%M"), "00:15")
})

test_that("figures draw without error", {
  lg <- two_site_loggers()
  res <- run_loggers(lg, water_min_C = -99)
  f <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  cov <- coverage(lg, res$pairs)
  expect_silent({
    plot_loggers(lg, c(ST_001 = "stream", AT_001 = "air", ST_002 = "stream", AT_002 = "air"))
    plot_pair(lg, "ST_001", "AT_001")
    plot_coverage(cov, res$date_range[1], res$date_range[2])
    plot_annual_fit(res$annual, "site_001")
    plot_thermal_sensitivity(res$annual, "site_001")
    plot_annual_by_year(res$annual)
    plot_diel_metrics(res$diel)
    plot_diel_day(res$diel, "site_001", "2025-06-01")
  })
})
