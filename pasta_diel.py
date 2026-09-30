"""
pasta_diel.py -- core Diel PASTA functions (Python port of
DielPASTA_Share/R/PASTA_dielAnalysisFunctions.R + the metrics chunk of
Diel_PASTA_Example.Rmd).

Everything here is pure pandas/numpy so it can be cached by Streamlit and
tested without the UI.

Sections
    1. Logger file parsing (HOBOware exports + simple 2-column CSVs)
    2. Air / stream classification
    3. Air-stream pairing suggestions
    4. Daily sinusoid fits (fit_TDAS_hour / TMd_output, vectorised)
    5. Paired metrics (AmpRatio, PhaseLag_h, mean ratios, 7-day rolls)
"""

import io
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ==============================================================================
#  1. Parsing
# ==============================================================================

_DATE_FORMATS = (
    "%m/%d/%Y %H:%M:%S", "%m/%d/%y %H:%M:%S", "%m/%d/%Y %H:%M", "%m/%d/%y %H:%M",
    "%m/%d/%y %I:%M:%S %p", "%m/%d/%Y %I:%M:%S %p",
    "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%dT%H:%M:%S",
)


def _parse_datetimes(s):
    """Try common logger formats first (fast), then fall back to inference."""
    s = s.astype(str).str.strip()
    sample = s.dropna().head(50)
    for fmt in _DATE_FORMATS:
        try:
            pd.to_datetime(sample, format=fmt)
        except (ValueError, TypeError):
            continue
        return pd.to_datetime(s, format=fmt, errors="coerce")
    return pd.to_datetime(s, errors="coerce")


def logger_name_from_filename(fname):
    """'HF_005 2026-07-22 12_44_20 EDT.csv' -> 'HF_005'."""
    stem = Path(str(fname)).stem
    stem = re.split(r"\s+\d{4}-\d{2}-\d{2}", stem)[0]   # drop HOBOware date stamp
    return stem.strip() or Path(str(fname)).stem


def parse_logger_csv(file_bytes, fname):
    """
    Parse one logger CSV into DataFrame['timestamp', 'temp_C'].

    Handles HOBOware exports (optional 'Plot Title' line, '#' index column,
    'Date-Time (EDT)' + 'Temperature (°C/°F)' columns, event-only rows at the
    end) and simple CSVs with a date/time column and a temperature column.
    """
    text = file_bytes.decode("utf-8-sig", errors="ignore")
    lines = text.splitlines()

    header_idx = None
    for i, line in enumerate(lines[:30]):
        low = line.lower()
        if ("date" in low or "time" in low) and ("temp" in low or "°c" in low or "value" in low):
            header_idx = i
            break
    if header_idx is None:
        header_idx = 0

    df = pd.read_csv(io.StringIO("\n".join(lines[header_idx:])), low_memory=False)
    cols = [str(c) for c in df.columns]

    ts_col = next((c for c in cols if "date" in c.lower() or "time" in c.lower()), None)
    t_col = next((c for c in cols if "temp" in c.lower()), None)
    if t_col is None:
        t_col = next((c for c in cols if "value" in c.lower()), None)
    if ts_col is None or t_col is None:
        raise ValueError(
            "Could not find a date/time column and a temperature column in "
            + str(fname) + ". Columns: " + str(cols)
        )

    out = pd.DataFrame({
        "timestamp": _parse_datetimes(df[ts_col]),
        "temp_C": pd.to_numeric(df[t_col], errors="coerce"),
    })
    if "°f" in t_col.lower() or re.search(r"\(\s*f\s*\)|deg\s*f", t_col.lower()):
        out["temp_C"] = (out["temp_C"] - 32.0) * 5.0 / 9.0

    out = (out.dropna()
              .sort_values("timestamp")
              .drop_duplicates(subset="timestamp")
              .reset_index(drop=True))
    if len(out) < 10:
        raise ValueError("Fewer than 10 valid temperature records in " + str(fname))

    tz = re.search(r"\(([A-Z]{2,5}|GMT[^)]*)\)", ts_col)
    out.attrs["tz_label"] = tz.group(1) if tz else ""
    return out


def to_grid(df, interval_min=15):
    """Round timestamps to the logging interval and average duplicates (R: round_date)."""
    g = df.copy()
    g["timestamp"] = g["timestamp"].dt.round(str(int(interval_min)) + "min")
    return g.groupby("timestamp", as_index=True)["temp_C"].mean()


# ==============================================================================
#  2. Classification
# ==============================================================================

_AIR_NAME = re.compile(r"(^|[^a-z])(air|atm|atmo)([^a-z]|$)", re.IGNORECASE)
_WAT_NAME = re.compile(r"(^|[^a-z])(str|stream|wat|water|sw|riv|river)([^a-z]|$)",
                       re.IGNORECASE)


def logger_features(series):
    """Summary features used to tell air from stream loggers."""
    s = series.dropna()
    daily = s.groupby(s.index.normalize())
    counts = daily.count()
    full = counts[counts >= 0.75 * counts.max()].index
    rng = (daily.max() - daily.min()).loc[full]
    step = s.diff().abs()
    return {
        "start": s.index.min(),
        "end": s.index.max(),
        "n_obs": int(len(s)),
        "mean_C": float(s.mean()),
        "min_C": float(s.min()),
        "max_C": float(s.max()),
        "median_daily_range_C": float(rng.median()) if len(rng) else float("nan"),
        "pct_below_0": float((s < 0).mean() * 100),
        "median_step_C": float(step.median()),
    }


def classify_logger(name, feats, range_threshold=4.0):
    """
    Return (type, confidence, reason).

    Stream water is buffered: small diel range and it cannot sit far below 0 °C.
    Air has a large diel range and noisy 15-min steps, and often freezes.
    A name containing 'air' / 'stream' etc. wins over the data.
    """
    if _AIR_NAME.search(name):
        return "air", "high", "name contains an air keyword"
    if _WAT_NAME.search(name):
        return "stream", "high", "name contains a water keyword"

    rng = feats["median_daily_range_C"]
    votes_air = 0
    reasons = []
    if rng > range_threshold:
        votes_air += 2
        reasons.append("daily range {:.1f} °C > {:.1f}".format(rng, range_threshold))
    else:
        reasons.append("daily range {:.1f} °C ≤ {:.1f}".format(rng, range_threshold))
    if feats["min_C"] < -1.0:
        votes_air += 1
        reasons.append("min {:.1f} °C (below freezing)".format(feats["min_C"]))
    if feats["median_step_C"] > 0.1:
        votes_air += 1
        reasons.append("noisy 15-min steps")

    kind = "air" if votes_air >= 2 else "stream"
    margin = abs(np.log(max(rng, 1e-3) / range_threshold))
    conf = "high" if margin > 0.7 else ("medium" if margin > 0.3 else "low")
    return kind, conf, "; ".join(reasons)


# ==============================================================================
#  3. Pairing
# ==============================================================================

def _digits(name):
    return re.findall(r"\d+", name)


def name_match(stream_name, air_name):
    """True when the loggers look like a named pair (HF_003 <-> HF_103, AMTH_A01 <-> AMTH_air)."""
    sd, ad = _digits(stream_name), _digits(air_name)
    if sd and ad and len(sd[-1]) >= 2 and len(ad[-1]) >= 2 and sd[-1][-2:] == ad[-1][-2:]:
        return True
    s_pre = re.split(r"[_\-\s]", stream_name)[0].lower()
    a_pre = re.split(r"[_\-\s]", air_name)[0].lower()
    return bool(_AIR_NAME.search(air_name)) and s_pre == a_pre and len(s_pre) > 2


def diel_xcorr(stream_s, air_s, interval_min=15, max_lag_h=6):
    """
    Max correlation of the high-passed (24 h detrended) signals with the
    stream lagging air by 0..max_lag_h hours. Returns (r, lag_h, overlap_days).
    """
    df = pd.concat([stream_s.rename("w"), air_s.rename("a")], axis=1, join="inner").dropna()
    per_day = int(round(24 * 60 / interval_min))
    if len(df) < 3 * per_day:
        return float("nan"), float("nan"), len(df) / per_day
    hp = df - df.rolling(per_day, center=True, min_periods=per_day // 2).mean()
    w = hp["w"].to_numpy()
    a = hp["a"].to_numpy()
    best = (-1.0, 0.0)
    for k in range(0, int(max_lag_h * 60 / interval_min) + 1):
        ww, aa = (w[k:], a[:len(a) - k]) if k else (w, a)
        ok = np.isfinite(ww) & np.isfinite(aa)
        if ok.sum() < per_day:
            continue
        r = np.corrcoef(ww[ok], aa[ok])[0, 1]
        if r > best[0]:
            best = (r, k * interval_min / 60)
    return float(best[0]), float(best[1]), len(df) / per_day


def suggest_pairs(streams, airs, score):
    """
    streams, airs: lists of names; score: dict[(stream, air)] -> r.
    Name matches first, then the best-scoring unused air, then the best air
    overall (one air logger may serve several streams).
    Returns dict stream -> (air or None, basis).
    """
    out, used = {}, set()
    for s in streams:
        hits = [a for a in airs if name_match(s, a)]
        if len(hits) == 1:
            out[s] = (hits[0], "name match")
            used.add(hits[0])
    remaining = [s for s in streams if s not in out]
    cands = sorted(
        ((score.get((s, a), float("nan")), s, a) for s in remaining for a in airs),
        key=lambda t: -np.nan_to_num(t[0], nan=-9),
    )
    for r, s, a in cands:
        if s in out or a in used or not np.isfinite(r):
            continue
        out[s] = (a, "best diel correlation")
        used.add(a)
    for s in remaining:
        if s in out:
            continue
        opts = [(score.get((s, a), float("nan")), a) for a in airs]
        opts = [o for o in opts if np.isfinite(o[0])]
        out[s] = (max(opts)[1], "best diel correlation (shared air)") if opts else (None, "no overlap")
    return out


# ==============================================================================
#  4. Daily sinusoid fits -- vectorised fit_TDAS_hour over every day
# ==============================================================================

def build_paired_df(pairs, grids, water_min_C=0.99):
    """
    pairs: list of (site_id, stream_name, air_name); grids: name -> gridded Series.
    Returns long df: site_id, DateTime_r, str_temp, air_temp  (R: T.y).
    """
    frames = []
    for site_id, s_name, a_name in pairs:
        d = pd.concat([grids[s_name].rename("str_temp"), grids[a_name].rename("air_temp")],
                      axis=1, join="inner").dropna()
        d = d[d["str_temp"] > water_min_C]
        d = d.rename_axis("DateTime_r").reset_index()
        d.insert(0, "site_id", site_id)
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=["site_id", "DateTime_r", "str_temp", "air_temp"])
    return pd.concat(frames, ignore_index=True)


def fit_daily_sinusoids(T_y, min_coverage=0.75, interval_min=15):
    """
    Per site/day/medium OLS fit  temp ~ sin(rday) + cos(rday), rday = 2*pi*h/24
    with h = hour + minute/60 (as rad_hour() in R). Returns the same columns as
    TMd_output(): amplitude_C, phase_h, AdjRsqr, RMSE, sinSlope, cosSlope,
    YInt, medium, site_id, date, year, plus n_obs.
    """
    if T_y.empty:
        return pd.DataFrame()
    dt = pd.to_datetime(T_y["DateTime_r"])
    h = dt.dt.hour + dt.dt.minute / 60.0
    rd = 2 * np.pi * h.to_numpy() / 24.0
    base = pd.DataFrame({
        "site_id": T_y["site_id"].to_numpy(),
        "date": dt.dt.date.to_numpy(),
        "s": np.sin(rd),
        "c": np.cos(rd),
    })
    per_day = 24 * 60 / interval_min
    out = []
    for medium, col in (("air", "air_temp"), ("water", "str_temp")):
        d = base.copy()
        d["y"] = T_y[col].to_numpy(dtype=float)
        d["ss"], d["cc"], d["sc"] = d.s * d.s, d.c * d.c, d.s * d.c
        d["sy"], d["cy"], d["yy"] = d.s * d.y, d.c * d.y, d.y * d.y
        g = d.groupby(["site_id", "date"], sort=True)
        S = g[["s", "c", "y", "ss", "cc", "sc", "sy", "cy", "yy"]].sum()
        n = g.size().to_numpy(dtype=float)

        # Normal equations X'X b = X'y with X = [1, sin, cos]  (batched 3x3)
        XtX = np.empty((len(S), 3, 3))
        XtX[:, 0, 0] = n
        XtX[:, 0, 1] = XtX[:, 1, 0] = S["s"]
        XtX[:, 0, 2] = XtX[:, 2, 0] = S["c"]
        XtX[:, 1, 1] = S["ss"]
        XtX[:, 1, 2] = XtX[:, 2, 1] = S["sc"]
        XtX[:, 2, 2] = S["cc"]
        Xty = np.stack([S["y"], S["sy"], S["cy"]], axis=1)
        ok = (n >= max(4, min_coverage * per_day)) & (np.abs(np.linalg.det(XtX)) > 1e-9)
        b = np.full((len(S), 3), np.nan)
        if ok.any():
            b[ok] = np.linalg.solve(XtX[ok], Xty[ok][..., None])[..., 0]

        sse = S["yy"].to_numpy() - np.einsum("ij,ij->i", b, Xty)
        sst = S["yy"].to_numpy() - S["y"].to_numpy() ** 2 / n
        sse = np.clip(sse, 0, None)
        with np.errstate(divide="ignore", invalid="ignore"):
            r2 = 1 - sse / sst
            adj = 1 - (1 - r2) * (n - 1) / (n - 3)
            rmse = np.sqrt(sse / n)
            phase = (24 / (2 * np.pi)) * (np.pi - np.arctan(b[:, 2] / b[:, 1]))
        amp = np.sqrt(b[:, 1] ** 2 + b[:, 2] ** 2)

        res = S.index.to_frame(index=False)
        res = res.assign(
            amplitude_C=amp, phase_h=phase, AdjRsqr=adj, RMSE=rmse,
            sinSlope=b[:, 1], cosSlope=b[:, 2], YInt=b[:, 0],
            medium=medium, n_obs=n.astype(int),
        )
        out.append(res[ok])
    fits = pd.concat(out, ignore_index=True)
    fits["date"] = pd.to_datetime(fits["date"])
    fits["year"] = fits["date"].dt.year
    cols = ["amplitude_C", "phase_h", "AdjRsqr", "RMSE", "sinSlope", "cosSlope", "YInt",
            "medium", "site_id", "date", "year", "n_obs"]
    return fits[cols].sort_values(["site_id", "date", "medium"]).reset_index(drop=True)


def sinusoid_value(fit_row, timestamps):
    """fitted_sinValue_hour() for one fit row."""
    ts = pd.to_datetime(pd.Series(timestamps))
    rd = 2 * np.pi * (ts.dt.hour + ts.dt.minute / 60.0) / 24.0
    return fit_row["sinSlope"] * np.sin(rd) + fit_row["cosSlope"] * np.cos(rd) + fit_row["YInt"]


# ==============================================================================
#  5. Paired metrics
# ==============================================================================

def pasta_metrics(fits, r2_min=0.7, roll_k=7):
    """
    Metrics chunk of Diel_PASTA_Example.Rmd: keep days where both air and
    water AdjRsqr > r2_min, then water/air ratios and lags.
    """
    if fits.empty:
        return pd.DataFrame()
    w = fits.pivot_table(index=["site_id", "date"], columns="medium",
                         values=["amplitude_C", "phase_h", "YInt", "AdjRsqr"])
    w.columns = [a + "_" + b for a, b in w.columns]
    w = w.dropna(subset=["AdjRsqr_air", "AdjRsqr_water"])
    w = w[(w["AdjRsqr_air"] > r2_min) & (w["AdjRsqr_water"] > r2_min)].reset_index()

    m = pd.DataFrame({
        "site_id": w["site_id"],
        "date": w["date"],
        "AmpRatio": (w["amplitude_C_water"] / w["amplitude_C_air"]).round(2),
        "PhaseLag_h": (w["phase_h_water"] - w["phase_h_air"]).round(2),
        "Ratio_Mean_C": (w["YInt_water"] / w["YInt_air"]).round(5),
        "Ratio_Mean_K": ((w["YInt_water"] + 273.15) / (w["YInt_air"] + 273.15)).round(5),
        "Amp_air_C": w["amplitude_C_air"].round(3),
        "Amp_water_C": w["amplitude_C_water"].round(3),
        "Mean_air_C": w["YInt_air"].round(3),
        "Mean_water_C": w["YInt_water"].round(3),
    }).sort_values(["site_id", "date"]).reset_index(drop=True)

    m["year"] = m["date"].dt.year
    m["DOY"] = m["date"].dt.dayofyear
    m["odd_year"] = m["year"] % 2 == 1
    # zoo::rollmean(k = 7, fill = NA): centred, over rows within each site
    g = m.groupby("site_id")
    m["AR_roll7"] = g["AmpRatio"].transform(lambda x: x.rolling(roll_k, center=True).mean())
    m["PH_roll7"] = g["PhaseLag_h"].transform(lambda x: x.rolling(roll_k, center=True).mean())
    return m
