"""
Annual signal (Python port of the PASTA Shiny app's
R/paired_annualT_signals.R, R/stream_thermal_sensitivity.R and
therm_analysis() in R/rundata_v1.R).

Works on daily mean air and stream temperature per site:
    - annual sinusoid fit to air and water (fit_TAS)
    - AmpRatio, PhaseLag_d, Ratio_Mean (therm_analysis)
    - linear thermal sensitivity, water ~ air (fit_ThermalSens, Kelleher et al. 2012)
    - data gap check (data_gap_check)
for the whole record and for each water or calendar year (TMy_output).
"""

import numpy as np
import pandas as pd

from . import phase as _ph


# ==============================================================================
#  Daily means
# ==============================================================================

def daily_means(T_y, min_coverage=0.75, interval_min=15):
    """
    Sub-daily paired df (site_id, DateTime_r, str_temp, air_temp) ->
    daily df (site_id, date, tavg_wat_C, tavg_air_C) in PASTA's column order.
    Days with less than `min_coverage` of a full day of paired records are dropped.
    """
    if T_y.empty:
        return pd.DataFrame(columns=["site_id", "date", "tavg_wat_C", "tavg_air_C", "n_obs"])
    d = T_y.assign(date=pd.to_datetime(T_y["DateTime_r"]).dt.normalize())
    g = d.groupby(["site_id", "date"], sort=True)
    out = g.agg(tavg_wat_C=("str_temp", "mean"), tavg_air_C=("air_temp", "mean"),
                n_obs=("str_temp", "size")).reset_index()
    per_day = 24 * 60 / interval_min
    return out[out["n_obs"] >= min_coverage * per_day].reset_index(drop=True)


# ==============================================================================
#  Radian day and water year
# ==============================================================================

def water_year(dates):
    dates = pd.to_datetime(pd.Series(dates))
    return np.where(dates.dt.month > 9, dates.dt.year + 1, dates.dt.year)


def rad_day(dates, yr_type="water"):
    """2*pi*d/365, d = day of water year (Oct 1 = 1) or day of calendar year."""
    dates = pd.to_datetime(pd.Series(dates)).reset_index(drop=True)
    if yr_type == "calendar":
        d = dates.dt.dayofyear.to_numpy(dtype=float)
    else:
        wy = water_year(dates)
        start = pd.to_datetime(pd.Series(wy - 1).astype(str) + "-09-30")
        d = (dates - start).dt.days.to_numpy(dtype=float)
    return 2 * np.pi * d / 365.0


# ==============================================================================
#  Fits
# ==============================================================================

def _ols(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    n, p = X.shape
    sse = float(resid @ resid)
    sst = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - sse / sst if sst > 0 else np.nan
    adj = 1 - (1 - r2) * (n - 1) / (n - p) if n > p else np.nan
    return b, adj, np.sqrt(sse / n)


def fit_TAS(dates, temp, yr_type="water", phase_method=_ph.DEFAULT_PHASE_METHOD):
    """Annual sinusoid fit: temp ~ sin(rday) + cos(rday). phase_d is the day of the
    fitted peak ('atan2') or the R formula ('atan'); see pastatemp.phase."""
    y = np.asarray(temp, dtype=float)
    rd = rad_day(dates, yr_type)
    ok = np.isfinite(y)
    y, rd = y[ok], rd[ok]
    if len(y) < 4:
        return None
    X = np.column_stack([np.ones_like(rd), np.sin(rd), np.cos(rd)])
    b, adj, rmse = _ols(X, y)
    phase = _ph.phase(b[1], b[2], 365, phase_method)
    return {"amplitude_C": float(np.hypot(b[1], b[2])), "phase_d": float(phase),
            "AdjRsqr": adj, "RMSE": rmse,
            "sinSlope": float(b[1]), "cosSlope": float(b[2]), "YInt": float(b[0])}


def fit_ThermalSens(a_temp, w_temp):
    """Daily linear thermal sensitivity, water ~ air (water > 1 °C)."""
    a = np.asarray(a_temp, dtype=float)
    w = np.asarray(w_temp, dtype=float)
    ok = np.isfinite(a) & np.isfinite(w) & (w > 1)
    if ok.sum() < 3:
        return {"TS__Slope": np.nan, "TS_AdjRsqr": np.nan, "TS_RMSE": np.nan, "TS_YInt": np.nan}
    X = np.column_stack([np.ones(ok.sum()), a[ok]])
    b, adj, rmse = _ols(X, w[ok])
    return {"TS__Slope": round(float(b[1]), 2), "TS_AdjRsqr": round(float(adj), 2),
            "TS_RMSE": rmse, "TS_YInt": round(float(b[0]), 2)}


def sinusoid_value(fit, dates, yr_type="water"):
    rd = rad_day(dates, yr_type)
    return fit["sinSlope"] * np.sin(rd) + fit["cosSlope"] * np.cos(rd) + fit["YInt"]


# ==============================================================================
#  Metrics
# ==============================================================================

def data_gap_check(daily):
    """Per site: longest gap between consecutive days and number of days."""
    rows = []
    d = daily.dropna(subset=["tavg_wat_C", "tavg_air_C"])
    d = d[d["tavg_wat_C"] < 70]
    for sid, x in d.groupby("site_id", sort=False):
        dd = pd.to_datetime(x["date"]).sort_values().diff().dt.days
        rows.append({"site_id": sid,
                     "max_conseq_missing_days": float(dd.max()) if dd.notna().any() else np.nan,
                     "count": int(len(x))})
    return pd.DataFrame(rows, columns=["site_id", "max_conseq_missing_days", "count"])


def therm_analysis(daily, yr_type="water", phase_method=_ph.DEFAULT_PHASE_METHOD):
    """
    Whole-record annual metrics per site. Returns (metrics, sine_fits).
    metrics: site_id, AmpRatio, PhaseLag_d, Ratio_Mean, TS__Slope, TS_AdjRsqr, TS_YInt
    """
    d = daily[daily["tavg_wat_C"] > 1]
    mets, fits = [], []
    for sid, x in d.groupby("site_id", sort=False):
        fa = fit_TAS(x["date"], x["tavg_air_C"], yr_type, phase_method)
        fw = fit_TAS(x["date"], x["tavg_wat_C"], yr_type, phase_method)
        ts = fit_ThermalSens(x["tavg_air_C"], x["tavg_wat_C"])
        for f, med in ((fa, "air"), (fw, "water")):
            if f:
                fits.append({**f, "medium": med, "site_id": sid, "n_days": int(len(x)),
                             "yr_type": yr_type, "phase_method": phase_method})
        if fa and fw:
            m = {"site_id": sid,
                 "AmpRatio": round(fw["amplitude_C"] / fa["amplitude_C"], 2),
                 "PhaseLag_d": round(float(_ph.lag(fw["phase_d"], fa["phase_d"], 365,
                                                  phase_method)), 2),
                 "Ratio_Mean": round(fw["YInt"] / fa["YInt"], 2)}
        else:
            m = {"site_id": sid, "AmpRatio": np.nan, "PhaseLag_d": np.nan, "Ratio_Mean": np.nan}
        m.update({"TS__Slope": ts["TS__Slope"], "TS_AdjRsqr": ts["TS_AdjRsqr"],
                  "TS_YInt": ts["TS_YInt"]})
        mets.append(m)
    cols = ["site_id", "AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope", "TS_AdjRsqr", "TS_YInt"]
    return pd.DataFrame(mets, columns=cols), pd.DataFrame(fits)


def _fit_days(daily):
    """Days that enter the sinusoid fit: paired daily means with stream > 1 °C."""
    d = daily.dropna(subset=["tavg_wat_C", "tavg_air_C"])
    return d[(d["tavg_wat_C"] > 1) & (d["tavg_wat_C"] < 70)]


SIGNAL_COLS = ["AmpRatio", "PhaseLag_d", "Ratio_Mean"]


def days_in_year(year, yr_type="water"):
    """Days in a water year (1 Oct of year-1 to 30 Sep) or calendar year: 365 or 366."""
    y = int(year)
    if yr_type == "water":
        return (pd.Timestamp(y, 10, 1) - pd.Timestamp(y - 1, 10, 1)).days
    return (pd.Timestamp(y + 1, 1, 1) - pd.Timestamp(y, 1, 1)).days


def span_coverage(all_dates, fit_dates):
    """(n_days_fit, span_days, pct_fit): fit days as a % of the span from the first
    to the last paired day (cold days included in the span)."""
    if len(all_dates) == 0:
        return 0, 0, 0.0
    span = int((all_dates.max() - all_dates.min()).days + 1)
    n = int(len(fit_dates))
    return n, span, 100.0 * n / span


def annual_metrics(daily, yr_type="water", min_year_pct=75, min_record_pct=75,
                   phase_method=_ph.DEFAULT_PHASE_METHOD):
    """
    Returns dict of DataFrames:
      'record' : whole-record metrics + gap check per site (PASTA metric_table)
      'yearly' : per water/calendar year (TMy_output)
      'fits'   : annual sinusoid coefficients (year 'all' = whole record, and per year)

    Both rules are applied after the 1 °C drop and count fit days: paired daily means
    with stream > 1 °C.
      yearly : pct_year_fit = 100 * n_days_fit / n_days_in_year, out of the full water
               (or calendar) year, 365 or 366 days. Below `min_year_pct` the year is
               dropped (meets_year_pct = False): no AmpRatio / PhaseLag_d / Ratio_Mean;
               its thermal sensitivity is still reported, as in PASTA. A partial
               first or last year counts the days outside the record as missing.
      record : pct_record_fit = 100 * n_days_fit / record_span_days, the span from the
               first to the last paired day (cold days included). Below
               `min_record_pct` the signal metrics are NA (meets_record_pct = False).
               All fit days are used; no year is dropped.
    'count' (all paired days, as in the R data_gap_check) is reported alongside.
    """
    if daily.empty:
        e = pd.DataFrame()
        return {"record": e, "yearly": e, "fits": e}

    d = daily.copy()
    d["date"] = pd.to_datetime(d["date"])
    d["year"] = water_year(d["date"]) if yr_type == "water" else d["date"].dt.year
    fd = _fit_days(d)
    all_by = {k: g["date"] for k, g in d.groupby(["site_id", "year"])}
    n_fit = fd.groupby(["site_id", "year"]).size()

    def _year_check(sid, yr):
        """(n_days_fit, n_days_in_year, pct_year_fit, meets_year_pct)"""
        n = int(n_fit.get((sid, yr), 0))
        tot = days_in_year(yr, yr_type)
        pct = 100.0 * n / tot
        return n, tot, round(pct, 1), bool(pct >= min_year_pct)

    # ---- per year: fit every year, report signal metrics only for years that qualify
    yr_rows, yr_fits = [], []
    for yr, x in d.groupby("year", sort=True):
        m, f = therm_analysis(x, yr_type, phase_method)
        m = m.merge(data_gap_check(x), on="site_id", how="left")
        m["year"] = int(yr)
        f["year"] = str(int(yr))
        yr_rows.append(m)
        yr_fits.append(f)
    yearly = pd.concat(yr_rows, ignore_index=True) if yr_rows else pd.DataFrame()
    if not yearly.empty:
        chk = pd.DataFrame([_year_check(s, y) for s, y in zip(yearly["site_id"], yearly["year"])],
                           columns=["n_days_fit", "n_days_in_year", "pct_year_fit",
                                    "meets_year_pct"])
        yearly = pd.concat([yearly, chk], axis=1)
        yearly.loc[~yearly["meets_year_pct"], SIGNAL_COLS] = np.nan
        yearly["year_type"] = yr_type
        yearly = yearly.sort_values(["site_id", "year"]).reset_index(drop=True)
    ok_years = {k for k in all_by if _year_check(*k)[-1]}

    # ---- whole record: all fit days; signal metrics only if they cover enough of the span
    rec, rec_fits = therm_analysis(d, yr_type, phase_method)
    rec = rec.merge(data_gap_check(d), on="site_id", how="left")
    rows = []
    for sid, x in d.groupby("site_id", sort=False):
        n, span, pct = span_coverage(x["date"], fd.loc[fd["site_id"] == sid, "date"])
        rows.append({"site_id": sid, "n_days_fit": n,
                     "record_span_days": span, "pct_record_fit": round(pct, 1),
                     "meets_record_pct": bool(pct >= min_record_pct)})
    rec = pd.DataFrame(rows).merge(rec, on="site_id", how="left")
    rec.loc[~rec["meets_record_pct"], SIGNAL_COLS] = np.nan
    cols = ["site_id", "AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope", "TS_AdjRsqr",
            "TS_YInt", "max_conseq_missing_days", "count", "n_days_fit",
            "record_span_days", "pct_record_fit", "meets_record_pct"]
    rec = rec.reindex(columns=cols)

    rec_fits["year"] = "all"
    if len(rec_fits):
        ok_rec = dict(zip(rec["site_id"], rec["meets_record_pct"]))
        rec_fits["meets_criteria"] = [bool(ok_rec.get(s, False)) for s in rec_fits["site_id"]]
    yf = pd.concat(yr_fits, ignore_index=True) if yr_fits else pd.DataFrame()
    if len(yf):
        yf["meets_criteria"] = [(s, int(y)) in ok_years for s, y in zip(yf["site_id"], yf["year"])]
    fits = pd.concat([rec_fits, yf], ignore_index=True)
    return {"record": rec, "yearly": yearly, "fits": fits}
