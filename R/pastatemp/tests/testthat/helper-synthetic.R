# Air and stream records with known annual and diel amplitudes and lags
# (the same signals as the Python tests).
ANN_AMP_AIR <- 12; ANN_AMP_WAT <- 6; ANN_LAG_D <- 10      # annual AmpRatio 0.5, lag 10 d
DIEL_AMP_AIR <- 5; DIEL_AMP_WAT <- 1; DIEL_LAG_H <- 2     # diel AmpRatio 0.2, lag 2 h
AIR_PEAK_H <- 14

synthetic <- function(start = "2024-10-01 00:00", end = "2026-09-30 23:45", step_min = 15, seed = 0) {
  t <- seq(as.POSIXct(start, tz = "UTC"), as.POSIXct(end, tz = "UTC"), by = step_min * 60)
  doy <- as.numeric(as.Date(t) - as.Date("2024-09-30")) %% 365
  lt <- as.POSIXlt(t)
  hour <- lt$hour + lt$min / 60
  ann <- function(amp, lag) amp * cos(2 * pi * (doy - 290 - lag) / 365)
  dl <- function(amp, peak) amp * cos(2 * pi * (hour - peak) / 24)
  set.seed(seed)
  air <- 10 + ann(ANN_AMP_AIR, 0) + dl(DIEL_AMP_AIR, AIR_PEAK_H) + stats::rnorm(length(t), 0, 0.2)
  wat <- 12 + ann(ANN_AMP_WAT, ANN_LAG_D) + dl(DIEL_AMP_WAT, AIR_PEAK_H + DIEL_LAG_H) +
    stats::rnorm(length(t), 0, 0.05)
  data.frame(site = "S1", time = t, stream = wat, air = air, stringsAsFactors = FALSE)
}

two_site_loggers <- function() {
  d <- synthetic(end = "2025-09-30 23:45")
  late <- d$time >= as.POSIXct("2025-01-01", tz = "UTC")
  lg <- function(v, i) data.frame(timestamp = d$time[i], temp_C = v[i])
  all <- rep(TRUE, nrow(d))
  list(ST_001 = lg(d$stream, all), AT_001 = lg(d$air, all),
       ST_002 = lg(d$stream, late), AT_002 = lg(d$air, late))
}
