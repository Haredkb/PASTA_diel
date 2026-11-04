"""
Diel signal: daily sinusoid fits and paired metrics (Python port of
DielPASTA_Share/R/PASTA_dielAnalysisFunctions.R and the metrics chunk of
Diel_PASTA_Example.Rmd).
"""

import numpy as np
import pandas as pd


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
