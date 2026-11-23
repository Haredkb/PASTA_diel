"""
The app's figures as Plotly functions. Each returns a plotly Figure
(fig.show() in a notebook, fig.write_html("x.html") to save).

Needs plotly:  pip install "pastatemp[plots]"
"""

import numpy as np
import pandas as pd

from . import annual as _annual
from . import diel as _diel

AIR_COLOR = "#d62728"
WATER_COLOR = "#1f77b4"
SITE_COLORS = ["#3F4570", "#BD9C8F", "#7C83BC", "#8ABD6A", "#557044",
               "#C06C84", "#355C7D", "#F8B195", "#6C5B7B", "#99B898"]
MAX_PLOT_POINTS = 2_000
_MARGIN = {"l": 20, "r": 20, "t": 20, "b": 20}


def _go():
    try:
        import plotly.graph_objects as go
    except ImportError as e:
        raise ImportError('Plots need plotly: pip install "pastatemp[plots]"') from e
    return go


def downsample(obj, n=MAX_PLOT_POINTS):
    """Every k-th row so that len(result) <= n (keeps plots responsive)."""
    if len(obj) <= n:
        return obj
    return obj.iloc[:: max(1, len(obj) // n)]


def plot_loggers(loggers, types=None):
    """Overview of every logger. types: {name: 'air'|'stream'} (e.g. from
    classify_loggers); air lines are thinner and faded."""
    go = _go()
    fig = go.Figure()
    for n, g in loggers.items():
        kind = (types or {}).get(n, "")
        if kind == "exclude":
            continue
        s = downsample(g)
        is_air = kind == "air"
        fig.add_trace(go.Scattergl(x=s.index, y=s.to_numpy(), mode="lines",
                                   name=n + (" ({})".format(kind) if kind else ""),
                                   line={"width": 1.0 if is_air else 1.6},
                                   opacity=0.55 if is_air else 1.0,
                                   legendgroup=kind or None))
    fig.update_layout(height=320, yaxis_title="Temperature (°C)", margin=_MARGIN)
    return fig


def plot_pair(loggers, stream, air):
    """Air and stream logger of one pair on one axis."""
    go = _go()
    fig = go.Figure()
    for n, col, nm in ((air, AIR_COLOR, "air: "), (stream, WATER_COLOR, "stream: ")):
        s = downsample(loggers[n])
        fig.add_trace(go.Scattergl(x=s.index, y=s.to_numpy(), mode="lines", name=nm + n,
                                   line={"width": 1.2, "color": col}))
    fig.update_layout(height=300, yaxis_title="Temperature (°C)", margin=_MARGIN)
    return fig


def plot_coverage(cov, start=None, end=None):
    """Timeline of air, stream and paired records per site (coverage() output),
    with the selected date range shaded."""
    go = _go()
    fig = go.Figure()
    if start is not None and end is not None:
        fig.add_vrect(x0=start, x1=end, fillcolor="rgba(31,119,180,0.08)", line_width=0,
                      annotation_text="selected range", annotation_position="top left")
    for i, r in enumerate(cov.itertuples()):
        for st, en, col, nm in ((r.air_start, r.air_end, AIR_COLOR, "air"),
                                (r.paired_start, r.paired_end, "#333333", "paired"),
                                (r.stream_start, r.stream_end, WATER_COLOR, "stream")):
            if pd.isna(st):
                continue
            fig.add_trace(go.Scatter(x=[st, en], y=[r.site_id, r.site_id], mode="lines",
                                     line={"color": col, "width": 9 if nm == "paired" else 4},
                                     name=nm, legendgroup=nm, showlegend=(i == 0),
                                     hovertemplate=nm + ": %{x|%Y-%m-%d}<extra>" + str(r.site_id) + "</extra>"))
    fig.update_layout(height=max(160, 46 * len(cov) + 80),
                      margin={"l": 20, "r": 20, "t": 30, "b": 20},
                      legend={"orientation": "h", "y": 1.12})
    return fig


def plot_annual_fit(annual, site):
    """Daily means of air and stream with the whole-record annual sinusoids.
    annual: AnnualResult (run_annual) or PastaResult.annual."""
    go = _go()
    yr_type = annual.settings.get("yr_type", "water")
    d = annual.daily[annual.daily["site_id"] == site]
    f = annual.fits
    f = f[(f["site_id"] == site) & (f["year"] == "all")] if len(f) else f
    fig = go.Figure()
    for med, col, ycol in (("air", AIR_COLOR, "tavg_air_C"), ("water", WATER_COLOR, "tavg_wat_C")):
        fig.add_trace(go.Scattergl(x=d["date"], y=d[ycol], mode="markers", name=med + " daily mean",
                                   marker={"size": 4, "color": col, "opacity": 0.5}))
        fr = f[f["medium"] == med] if len(f) else f
        if len(fr) and len(d):
            fig.add_trace(go.Scatter(x=d["date"], y=_annual.sinusoid_value(fr.iloc[0], d["date"], yr_type),
                                     mode="lines", name=med + " annual fit",
                                     line={"color": col, "width": 2}))
    fig.update_layout(height=340, yaxis_title="Daily mean temperature (°C)", margin=_MARGIN)
    return fig


def plot_thermal_sensitivity(annual, site):
    """Daily stream vs air means with the 1:1 line and the thermal-sensitivity fit."""
    go = _go()
    d = annual.daily[annual.daily["site_id"] == site]
    rec = annual.record[annual.record["site_id"] == site] if len(annual.record) else annual.record
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=d["tavg_air_C"], y=d["tavg_wat_C"], mode="markers", name="days",
                               marker={"size": 4, "color": WATER_COLOR, "opacity": 0.5}))
    if len(d):
        lo, hi = float(d["tavg_air_C"].min()), float(d["tavg_air_C"].max())
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", name="1:1",
                                 line={"dash": "dash", "color": "gray", "width": 1}))
        if len(rec) and np.isfinite(rec["TS__Slope"].iloc[0]):
            sl, yi = rec["TS__Slope"].iloc[0], rec["TS_YInt"].iloc[0]
            fig.add_trace(go.Scatter(x=[lo, hi], y=[yi + sl * lo, yi + sl * hi], mode="lines",
                                     name="slope {:.2f}".format(sl),
                                     line={"color": "black", "width": 2}))
    fig.update_layout(height=340, xaxis_title="Air (°C)", yaxis_title="Stream (°C)",
                      margin=_MARGIN, legend={"orientation": "h", "y": -0.25})
    return fig


def plot_annual_by_year(annual):
    """AmpRatio vs phase lag for every site and qualifying year."""
    go = _go()
    yl = annual.yearly.dropna(subset=["AmpRatio", "PhaseLag_d"]) if len(annual.yearly) else annual.yearly
    fig = go.Figure()
    for i, (sid, d) in enumerate(yl.groupby("site_id", sort=False)):
        fig.add_trace(go.Scatter(x=d["PhaseLag_d"], y=d["AmpRatio"], mode="markers+text",
                                 text=d["year"].astype(str), textposition="top center", name=sid,
                                 marker={"size": 9, "color": SITE_COLORS[i % len(SITE_COLORS)]}))
    fig.update_layout(height=340, xaxis_title="Phase lag (d)", yaxis_title="Amplitude ratio",
                      margin={"l": 20, "r": 20, "t": 30, "b": 20},
                      title={"text": "Annual signal by year", "font": {"size": 14}})
    return fig


def plot_diel_metrics(diel, x="date", max_amp_ratio=1.2):
    """Daily AmpRatio, phase lag and mean ratio per site, with rolling means.
    x: 'date' or 'DOY'. Days with AmpRatio >= max_amp_ratio are hidden (None = show all)."""
    go = _go()
    from plotly.subplots import make_subplots
    m = diel.metrics
    if max_amp_ratio is not None and len(m):
        m = m[m["AmpRatio"] < max_amp_ratio]
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.05,
                        subplot_titles=["Amplitude ratio (water / air)",
                                        "Phase lag (h, water − air)",
                                        "Mean ratio (K, water / air)"])
    for i, (sid, d) in enumerate(m.groupby("site_id", sort=False)):
        col = SITE_COLORS[i % len(SITE_COLORS)]
        for row, y, roll in ((1, "AmpRatio", "AR_roll7"), (2, "PhaseLag_h", "PH_roll7"),
                             (3, "Ratio_Mean_K", None)):
            fig.add_trace(go.Scattergl(x=d[x], y=d[y], mode="markers", name=sid, legendgroup=sid,
                                       showlegend=(row == 1),
                                       marker={"size": 5, "color": col, "opacity": 0.7}),
                          row=row, col=1)
            if roll:
                fig.add_trace(go.Scatter(x=d[x], y=d[roll], mode="lines", name=sid + " roll",
                                         legendgroup=sid, showlegend=False,
                                         line={"color": col, "width": 2}), row=row, col=1)
    fig.update_layout(height=720, margin={"l": 20, "r": 20, "t": 40, "b": 20})
    return fig


def plot_diel_day(diel, site, day):
    """Records and fitted daily sinusoids of air and stream for one site and day."""
    go = _go()
    day = pd.Timestamp(day).date()
    p = diel.paired
    obs = p[(p["site_id"] == site) & (pd.to_datetime(p["DateTime_r"]).dt.date == day)]
    fs = diel.fits[diel.fits["site_id"] == site] if len(diel.fits) else diel.fits
    fig = go.Figure()
    for medium, col, ycol in (("air", AIR_COLOR, "air_temp"), ("water", WATER_COLOR, "str_temp")):
        fig.add_trace(go.Scatter(x=obs["DateTime_r"], y=obs[ycol], mode="markers", name=medium,
                                 marker={"size": 4, "color": col, "opacity": 0.6}))
        fr = fs[(fs["date"].dt.date == day) & (fs["medium"] == medium)] if len(fs) else fs
        if len(fr):
            r = fr.iloc[0]
            fig.add_trace(go.Scatter(x=obs["DateTime_r"], y=_diel.sinusoid_value(r, obs["DateTime_r"]),
                                     mode="lines",
                                     name="{} fit (amp {:.2f} °C, phase {:.2f} h, adj R² {:.2f})".format(
                                         medium, r["amplitude_C"], r["phase_h"], r["AdjRsqr"]),
                                     line={"color": col, "width": 2}))
    fig.update_layout(height=320, yaxis_title="Temperature (°C)", margin=_MARGIN)
    return fig
