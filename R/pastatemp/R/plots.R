# The app's figures in base R graphics (no extra packages).

.air_col <- "#d62728"
.wat_col <- "#1f77b4"
.site_cols <- c("#3F4570", "#BD9C8F", "#7C83BC", "#8ABD6A", "#557044",
                "#C06C84", "#355C7D", "#F8B195", "#6C5B7B", "#99B898")

thin <- function(x, n = 2000) if (nrow(x) <= n) x else x[seq(1, nrow(x), by = max(1, nrow(x) %/% n)), , drop = FALSE]

#' Figures
#'
#' Base-graphics versions of the app's plots.
#' * `plot_loggers()`: every logger (air thin, stream thick).
#' * `plot_pair()`: the air and stream logger of one pair.
#' * `plot_coverage()`: air, stream and paired record per site, with the date range.
#' * `plot_annual_fit()`: daily means with the whole-record annual sinusoids.
#' * `plot_thermal_sensitivity()`: stream vs air daily means, 1:1 line and slope.
#' * `plot_annual_by_year()`: amplitude ratio vs phase lag for each site and year.
#' * `plot_diel_metrics()`: daily amplitude ratio, phase lag and mean ratio.
#' * `plot_diel_day()`: one day's records and fitted sinusoids.
#' @param loggers Named list of loggers.
#' @param types Named character c(name = "air" | "stream").
#' @name pasta_plots
NULL

#' @rdname pasta_plots
#' @export
plot_loggers <- function(loggers, types = NULL) {
  nm <- names(loggers)
  if (!is.null(types)) nm <- nm[!(types[nm] %in% "exclude")]
  t <- do.call(c, lapply(loggers[nm], function(g) g$timestamp))
  y <- unlist(lapply(loggers[nm], function(g) g$temp_C))
  graphics::plot(range(t), range(y, na.rm = TRUE), type = "n", xlab = "", ylab = "Temperature (\u00b0C)")
  cols <- grDevices::hcl.colors(max(length(nm), 2), "Dark 3")
  for (i in seq_along(nm)) {
    g <- thin(loggers[[nm[i]]])
    air <- !is.null(types) && identical(unname(types[nm[i]]), "air")
    graphics::lines(g$timestamp, g$temp_C, col = grDevices::adjustcolor(cols[i], if (air) 0.55 else 1),
                    lwd = if (air) 1 else 1.6)
  }
  graphics::legend("topleft", legend = nm, col = cols[seq_along(nm)], lwd = 1.5, bty = "n", cex = 0.8)
  invisible(NULL)
}

#' @rdname pasta_plots
#' @param stream,air Logger names.
#' @export
plot_pair <- function(loggers, stream, air) {
  a <- thin(loggers[[air]])
  s <- thin(loggers[[stream]])
  graphics::plot(a$timestamp, a$temp_C, type = "l", col = .air_col, xlab = "",
                 ylab = "Temperature (\u00b0C)", ylim = range(c(a$temp_C, s$temp_C), na.rm = TRUE))
  graphics::lines(s$timestamp, s$temp_C, col = .wat_col, lwd = 1.2)
  graphics::legend("topleft", c(paste("air:", air), paste("stream:", stream)),
                   col = c(.air_col, .wat_col), lwd = 1.5, bty = "n")
  invisible(NULL)
}

#' @rdname pasta_plots
#' @param cov Output of [coverage()].
#' @param start,end Selected date range (shaded).
#' @export
plot_coverage <- function(cov, start = NULL, end = NULL) {
  n <- nrow(cov)
  xr <- range(c(cov$air_start, cov$air_end, cov$stream_start, cov$stream_end), na.rm = TRUE)
  graphics::plot(xr, c(0.5, n + 0.5), type = "n", yaxt = "n", xlab = "", ylab = "")
  graphics::axis(2, at = seq_len(n), labels = cov$site_id, las = 1, cex.axis = 0.8)
  if (!is.null(start) && !is.null(end)) {
    graphics::rect(start, 0, end, n + 1, col = grDevices::adjustcolor(.wat_col, 0.08), border = NA)
  }
  for (i in seq_len(n)) {
    graphics::segments(cov$air_start[i], i + 0.18, cov$air_end[i], i + 0.18, col = .air_col, lwd = 4)
    graphics::segments(cov$paired_start[i], i, cov$paired_end[i], i, col = "#333333", lwd = 9)
    graphics::segments(cov$stream_start[i], i - 0.18, cov$stream_end[i], i - 0.18, col = .wat_col, lwd = 4)
  }
  graphics::legend("top", c("air", "paired", "stream"), col = c(.air_col, "#333333", .wat_col),
                   lwd = c(4, 9, 4), horiz = TRUE, bty = "n", cex = 0.8, xpd = TRUE, inset = -0.12)
  invisible(NULL)
}

#' @rdname pasta_plots
#' @param annual A `pasta_annual` (e.g. `res$annual`).
#' @param site Site ID.
#' @export
plot_annual_fit <- function(annual, site) {
  yt <- annual$settings$yr_type
  d <- annual$daily[annual$daily$site_id == site, , drop = FALSE]
  f <- annual$fits
  f <- f[f$site_id == site & f$year == "all", , drop = FALSE]
  graphics::plot(d$date, d$tavg_air_C, pch = 16, cex = 0.5, col = grDevices::adjustcolor(.air_col, 0.5),
                 ylim = range(c(d$tavg_air_C, d$tavg_wat_C), na.rm = TRUE), xlab = "",
                 ylab = "Daily mean temperature (\u00b0C)")
  graphics::points(d$date, d$tavg_wat_C, pch = 16, cex = 0.5, col = grDevices::adjustcolor(.wat_col, 0.5))
  for (med in c("air", "water")) {
    fr <- f[f$medium == med, , drop = FALSE]
    if (nrow(fr)) graphics::lines(d$date, annual_sinusoid_value(fr[1, ], d$date, yt),
                                  col = if (med == "air") .air_col else .wat_col, lwd = 2)
  }
  graphics::legend("topleft", c("air", "water"), col = c(.air_col, .wat_col), lwd = 2, bty = "n")
  invisible(NULL)
}

#' @rdname pasta_plots
#' @export
plot_thermal_sensitivity <- function(annual, site) {
  d <- annual$daily[annual$daily$site_id == site, , drop = FALSE]
  rec <- annual$record[annual$record$site_id == site, , drop = FALSE]
  graphics::plot(d$tavg_air_C, d$tavg_wat_C, pch = 16, cex = 0.5,
                 col = grDevices::adjustcolor(.wat_col, 0.5), xlab = "Air (\u00b0C)", ylab = "Stream (\u00b0C)")
  graphics::abline(0, 1, lty = 2, col = "gray")
  if (nrow(rec) && is.finite(rec$TS__Slope[1])) {
    graphics::abline(rec$TS_YInt[1], rec$TS__Slope[1], lwd = 2)
    graphics::legend("topleft", sprintf("slope %.2f", rec$TS__Slope[1]), lwd = 2, bty = "n")
  }
  invisible(NULL)
}

#' @rdname pasta_plots
#' @export
plot_annual_by_year <- function(annual) {
  y <- annual$yearly
  y <- y[is.finite(y$AmpRatio) & is.finite(y$PhaseLag_d), , drop = FALSE]
  if (!nrow(y)) return(invisible(NULL))
  sites <- unique(y$site_id)
  col <- .site_cols[(match(y$site_id, sites) - 1) %% length(.site_cols) + 1]
  graphics::plot(y$PhaseLag_d, y$AmpRatio, pch = 16, col = col, xlab = "Phase lag (d)",
                 ylab = "Amplitude ratio", main = "Annual signal by year")
  graphics::text(y$PhaseLag_d, y$AmpRatio, y$year, pos = 3, cex = 0.7)
  graphics::legend("topright", sites, col = .site_cols[seq_along(sites)], pch = 16, bty = "n", cex = 0.8)
  invisible(NULL)
}

#' @rdname pasta_plots
#' @param diel A `pasta_diel` (e.g. `res$diel`).
#' @param x `"date"` or `"DOY"`.
#' @param max_amp_ratio Days with AmpRatio at or above this are hidden (`NULL` shows all).
#' @export
plot_diel_metrics <- function(diel, x = "date", max_amp_ratio = 1.2) {
  m <- diel$metrics
  if (!is.null(max_amp_ratio) && nrow(m)) m <- m[m$AmpRatio < max_amp_ratio, , drop = FALSE]
  op <- graphics::par(mfrow = c(3, 1), mar = c(2, 4, 2, 1))
  on.exit(graphics::par(op))
  sites <- unique(m$site_id)
  for (p in list(c("AmpRatio", "AR_roll7", "Amplitude ratio (water / air)"),
                 c("PhaseLag_h", "PH_roll7", "Phase lag (h, water - air)"),
                 c("Ratio_Mean_K", NA, "Mean ratio (K, water / air)"))) {
    graphics::plot(m[[x]], m[[p[1]]], type = "n", xlab = "", ylab = "", main = p[3])
    for (i in seq_along(sites)) {
      d <- m[m$site_id == sites[i], , drop = FALSE]
      col <- .site_cols[(i - 1) %% length(.site_cols) + 1]
      graphics::points(d[[x]], d[[p[1]]], pch = 16, cex = 0.6, col = grDevices::adjustcolor(col, 0.7))
      if (!is.na(p[2])) graphics::lines(d[[x]], d[[p[2]]], col = col, lwd = 2)
    }
  }
  invisible(NULL)
}

#' @rdname pasta_plots
#' @param day A date.
#' @export
plot_diel_day <- function(diel, site, day) {
  day <- as.Date(day)
  p <- diel$paired
  obs <- p[p$site_id == site & day_of(p$DateTime_r) == day, , drop = FALSE]
  fs <- diel$fits[diel$fits$site_id == site & diel$fits$date == day, , drop = FALSE]
  graphics::plot(obs$DateTime_r, obs$air_temp, pch = 16, cex = 0.5, col = .air_col,
                 ylim = range(c(obs$air_temp, obs$str_temp), na.rm = TRUE), xlab = "",
                 ylab = "Temperature (\u00b0C)")
  graphics::points(obs$DateTime_r, obs$str_temp, pch = 16, cex = 0.5, col = .wat_col)
  lab <- character(0)
  for (med in c("air", "water")) {
    fr <- fs[fs$medium == med, , drop = FALSE]
    if (nrow(fr)) {
      graphics::lines(obs$DateTime_r, diel_sinusoid_value(fr[1, ], obs$DateTime_r),
                      col = if (med == "air") .air_col else .wat_col, lwd = 2)
      lab <- c(lab, sprintf("%s: amp %.2f \u00b0C, phase %.2f h, adj R\u00b2 %.2f", med,
                            fr$amplitude_C[1], fr$phase_h[1], fr$AdjRsqr[1]))
    }
  }
  if (length(lab)) graphics::legend("topleft", lab, col = c(.air_col, .wat_col)[seq_along(lab)],
                                    lwd = 2, bty = "n", cex = 0.8)
  invisible(NULL)
}
