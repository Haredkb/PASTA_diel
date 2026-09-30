"""
Diel PASTA -- Streamlit Application
Paired Air and Stream Temperature Analysis at the diel (daily) scale.

Upload air and stream temperature logger files (HOBOware exports or simple
date/time + temperature CSVs); the app sorts loggers into air vs stream,
suggests air-stream pairs, fits a daily sinusoid to each record and returns
the paired metrics (amplitude ratio, phase lag, mean ratio).

Usage:
    streamlit run diel_pasta_app.py        (or double-click run_diel_pasta.bat)

Core functions live in pasta_diel.py (Python port of
DielPASTA_Share/R/PASTA_dielAnalysisFunctions.R).
"""

# -- Auto-install missing packages from requirements_diel_pasta_app.txt --------
import importlib
import subprocess
import sys
import os
from pathlib import Path

_REQUIREMENTS_FILE = Path(__file__).parent / "requirements_diel_pasta_app.txt"


def _check_and_install():
    """Install any missing packages; return True if the app should re-exec."""
    if not _REQUIREMENTS_FILE.exists():
        return False
    missing = []
    for line in _REQUIREMENTS_FILE.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        pkg = line.split(">=")[0].split("==")[0].split("<=")[0].strip()
        try:
            importlib.import_module(pkg.lower().replace("-", "_"))
        except ImportError:
            missing.append(line)
    if not missing:
        return False
    print("[diel_pasta_app] Installing missing packages: " + str(missing), flush=True)
    any_failed = False
    for spec in missing:
        r = subprocess.run([sys.executable, "-m", "pip", "install", "--prefer-binary", spec])
        any_failed = any_failed or r.returncode != 0
    return not any_failed


if _check_and_install():
    os.execv(sys.executable, [sys.executable] + sys.argv)
# -- End auto-install ----------------------------------------------------------

import io
import hashlib
import datetime as _dt

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pasta_diel as pdl

pio.templates.default = "plotly_white"

# Maximum points in any single Plotly trace (browser stays responsive).
_MAX_PLOT_POINTS = 2_000

AIR_COLOR = "#d62728"
WATER_COLOR = "#1f77b4"
SITE_COLORS = ["#3F4570", "#BD9C8F", "#7C83BC", "#8ABD6A", "#557044",
               "#C06C84", "#355C7D", "#F8B195", "#6C5B7B", "#99B898"]

# Pre-fill the folder box by setting DIEL_PASTA_DATA to your data folder.
_DEFAULT_FOLDER = os.environ.get("DIEL_PASTA_DATA", "")

# st.fragment (Streamlit >= 1.37) reruns only the decorated block when its own
# widgets change -- used for the plot selectors so the rest of the page is not
# recomputed. Falls back to a no-op decorator on older Streamlit.
_fragment = getattr(st, "fragment", None) or (lambda f: f)


def _downsample(obj, n=_MAX_PLOT_POINTS):
    """Return every k-th row so that len(result) <= n."""
    if len(obj) <= n:
        return obj
    return obj.iloc[:: max(1, len(obj) // n)]


st.set_page_config(
    page_title="Diel PASTA",
    page_icon="🌡️",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ==============================================================================
#  Module-level cache helpers
# ==============================================================================

@st.cache_data(show_spinner=False, max_entries=200)
def _load_logger(file_bytes, fname, interval_min):
    """Parse + grid one logger file and compute its features (cached on bytes)."""
    raw = pdl.parse_logger_csv(file_bytes, fname)
    grid = pdl.to_grid(raw, interval_min)
    feats = pdl.logger_features(grid)
    feats["tz_label"] = raw.attrs.get("tz_label", "")
    return grid, feats


@st.cache_data(show_spinner=False)
def _read_folder(folder_str):
    """Read every CSV in a folder once; returns [(name, bytes)]."""
    return [(p.name, p.read_bytes()) for p in sorted(Path(folder_str).glob("*.csv"))]


@st.cache_data(show_spinner=False)
def _pair_scores(sig, _grids, streams, airs, interval_min):
    """Diel cross-correlation for every stream x air combination.
    `sig` identifies the inputs; `_grids` is not hashed (leading underscore)."""
    out = {}
    for s in streams:
        for a in airs:
            out[(s, a)] = pdl.diel_xcorr(_grids[s], _grids[a], interval_min)
    return out


@st.cache_data(show_spinner=False)
def _overview_frame(sig, _grids, names):
    """Downsampled long table of every logger for the overview plot."""
    parts = []
    for n in names:
        s = _downsample(_grids[n])
        parts.append(pd.DataFrame({"timestamp": s.index, "temp_C": s.to_numpy(), "logger": n}))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


@st.cache_data(show_spinner=False)
def _df_to_csv(df):
    """Convert DataFrame to UTF-8 CSV bytes (cached)."""
    return df.to_csv(index=False).encode("utf-8")


def _sig(*parts):
    return hashlib.md5(repr(parts).encode()).hexdigest()


# ==============================================================================
#  UI
# ==============================================================================

st.title("Diel PASTA -- Paired Air & Stream Temperature Analysis")
st.caption(
    "Daily sinusoid fits to paired air and stream temperature  |  "
    "amplitude ratio · phase lag · mean ratio  |  port of PASTA_dielAnalysisFunctions.R"
)
st.markdown("---")

# ---- 1: Load logger files ----------------------------------------------------
st.subheader("1 - Logger files")

c_src, c_int = st.columns([3, 1])
with c_src:
    source_mode = st.radio(
        "Data source",
        ["📤 Upload CSV files", "📂 Local folder path"],
        horizontal=True,
        key="src_mode",
        help="**Upload** works everywhere. **Local folder** reads every CSV in a folder "
             "and only works when the app runs on your own computer.",
    )
with c_int:
    interval_min = st.number_input(
        "Logging interval (min)", value=15, min_value=1, max_value=120, step=1,
        help="Timestamps are rounded to this interval before air and stream are joined "
             "(R: round_date(..., '15 mins')).",
    )

files = []  # [(fname, bytes)]
if source_mode == "📤 Upload CSV files":
    st.caption(
        "Select all air and stream files at once (Ctrl/Cmd-click to multi-select). "
        "HOBOware exports work as-is; other CSVs need a date/time column and a "
        "temperature column."
    )
    ups = st.file_uploader("Logger CSV files", type=["csv"], accept_multiple_files=True,
                           key="uploader", label_visibility="collapsed")
    files = [(u.name, u.getvalue()) for u in (ups or [])]
else:
    folder = st.text_input("Folder with logger CSVs", value=_DEFAULT_FOLDER, key="folder",
                           placeholder=r"C:\path\to\logger_csvs")
    folder = folder.strip().strip('"')
    if not folder:
        pass
    elif Path(folder).is_dir():
        files = _read_folder(folder)
    else:
        st.warning("Folder not found: `" + folder + "`")

if not files:
    st.info("⬆️ Add air and stream temperature files to get started.")
    st.stop()

grids, feats_rows, errors = {}, [], []
for fname, b in files:
    name = pdl.logger_name_from_filename(fname)
    if name in grids:  # two files with the same logger name -> keep both, disambiguate
        name = Path(fname).stem
    try:
        g, f = _load_logger(b, fname, int(interval_min))
    except Exception as exc:
        errors.append("**" + fname + "**: " + str(exc))
        continue
    grids[name] = g
    kind, conf, why = pdl.classify_logger(name, f)
    feats_rows.append({"logger": name, "file": fname, **f,
                       "auto_type": kind, "confidence": conf, "reason": why})
for e in errors:
    st.error(e)
if not grids:
    st.stop()

feats = pd.DataFrame(feats_rows)
tzs = sorted(set(feats["tz_label"]) - {""})
if len(tzs) > 1:
    st.warning("Files use different time zones (" + ", ".join(tzs) + "). "
               "Re-export them in one zone before pairing, or the phase lag will be off.")

# ---- 2: Air or stream --------------------------------------------------------
st.markdown("---")
st.subheader("2 - Air or stream?")
st.caption(
    "Each logger is sorted from its data: stream water has a small daily range "
    "(typically < 4 °C) and stays near or above 0 °C; air swings widely and freezes. "
    "A name containing *air* or *stream/water* overrides the data. "
    "Change **Type** to correct a logger, or set it to *exclude*."
)

cls_view = feats[["logger", "auto_type", "confidence", "median_daily_range_C", "min_C",
                  "mean_C", "max_C", "start", "end", "n_obs", "reason"]].copy()
cls_view.insert(1, "type", cls_view["auto_type"])
cls_edit = st.data_editor(
    cls_view,
    column_config={
        "logger": st.column_config.TextColumn("Logger", disabled=True),
        "type": st.column_config.SelectboxColumn("Type", options=["air", "stream", "exclude"],
                                                 required=True, width="small"),
        "auto_type": st.column_config.TextColumn("Suggested", disabled=True, width="small"),
        "confidence": st.column_config.TextColumn("Confidence", disabled=True, width="small"),
        "median_daily_range_C": st.column_config.NumberColumn("Median daily range (°C)",
                                                              format="%.2f", disabled=True),
        "min_C": st.column_config.NumberColumn("Min (°C)", format="%.2f", disabled=True),
        "mean_C": st.column_config.NumberColumn("Mean (°C)", format="%.2f", disabled=True),
        "max_C": st.column_config.NumberColumn("Max (°C)", format="%.2f", disabled=True),
        "start": st.column_config.DatetimeColumn("Start", format="YYYY-MM-DD HH:mm", disabled=True),
        "end": st.column_config.DatetimeColumn("End", format="YYYY-MM-DD HH:mm", disabled=True),
        "n_obs": st.column_config.NumberColumn("Records", disabled=True),
        "reason": st.column_config.TextColumn("Why", disabled=True),
    },
    hide_index=True,
    width="stretch",
    num_rows="fixed",
    key="cls_editor_" + _sig(tuple(feats["logger"])),
)

types = dict(zip(cls_edit["logger"], cls_edit["type"]))
airs = [n for n, t in types.items() if t == "air"]
streams = [n for n, t in types.items() if t == "stream"]

with st.expander("Overview plot: all loggers", expanded=True):
    ov_names = tuple(airs + streams)
    ov = _overview_frame(_sig(tuple(sorted((f, len(b)) for f, b in files)), interval_min),
                         grids, ov_names)
    fig_ov = go.Figure()
    for n in ov_names:
        d = ov[ov["logger"] == n]
        is_air = types[n] == "air"
        fig_ov.add_trace(go.Scattergl(
            x=d["timestamp"], y=d["temp_C"], mode="lines",
            name=n + (" (air)" if is_air else " (stream)"),
            line={"width": 1.0 if is_air else 1.6},
            opacity=0.55 if is_air else 1.0,
            legendgroup="air" if is_air else "stream",
        ))
    fig_ov.update_layout(height=320, yaxis_title="Temperature (°C)",
                         margin={"l": 20, "r": 20, "t": 20, "b": 20})
    st.plotly_chart(fig_ov, width="stretch")

if not airs or not streams:
    st.info("Need at least one **air** and one **stream** logger to pair.")
    st.stop()

# ---- 3: Pairs ----------------------------------------------------------------
st.markdown("---")
st.subheader("3 - Air-stream pairs")
st.caption(
    "Each stream logger is paired with one air logger. Suggestions use matching "
    "numbers in the names first (HF_003 ↔ HF_103), then the diel correlation: how "
    "closely the stream's daily cycle follows the air's, allowing a 0–6 h lag. "
    "One air logger can serve several streams. Change **Air logger** to override."
)

score_sig = _sig(tuple(sorted((f, len(b)) for f, b in files)), interval_min)
scores = _pair_scores(score_sig, grids, tuple(streams), tuple(airs), int(interval_min))
suggest = pdl.suggest_pairs(streams, airs, {k: v[0] for k, v in scores.items()})

pair_rows = []
for s in streams:
    a, basis = suggest[s]
    rs = sorted((scores[(s, x)][0] for x in airs if np.isfinite(scores[(s, x)][0])), reverse=True)
    if basis == "name match":
        check = "name match"
    elif len(rs) > 1 and rs[0] - rs[1] < 0.02:
        check = "confirm (air loggers score alike)"
    elif rs and rs[0] < 0.3:
        check = "confirm (weak diel link)"
    else:
        check = "data"
    r, lag, days = scores.get((s, a), (np.nan, np.nan, np.nan)) if a else (np.nan, np.nan, np.nan)
    pair_rows.append({"include": a is not None, "site_id": s, "stream": s, "air": a,
                      "basis": check, "diel_r": r, "lag_h": lag, "overlap_days": days})
pairs_df = pd.DataFrame(pair_rows)

pairs_edit = st.data_editor(
    pairs_df,
    column_config={
        "include": st.column_config.CheckboxColumn("Run", width="small"),
        "site_id": st.column_config.TextColumn("Site ID", help="Label used in the outputs"),
        "stream": st.column_config.TextColumn("Stream logger", disabled=True),
        "air": st.column_config.SelectboxColumn("Air logger", options=airs, required=True),
        "basis": st.column_config.TextColumn("Suggested by", disabled=True),
        "diel_r": st.column_config.NumberColumn("Diel r (suggested)", format="%.3f", disabled=True),
        "lag_h": st.column_config.NumberColumn("Lag (h)", format="%.2f", disabled=True),
        "overlap_days": st.column_config.NumberColumn("Overlap (d)", format="%.1f", disabled=True),
    },
    hide_index=True,
    width="stretch",
    num_rows="fixed",
    key="pair_editor_" + _sig(tuple(streams), tuple(airs)),
)

with st.expander("Diel correlation of every stream × air combination"):
    mat = pd.DataFrame({a: [scores[(s, a)][0] for s in streams] for a in airs}, index=streams)
    st.dataframe(mat.round(3), width="stretch")
    st.caption("Rows are stream loggers, columns air loggers. Air loggers at one site "
               "often score within a few thousandths of each other; then the pairing "
               "should come from field notes rather than the data.")

sel_pairs = pairs_edit[pairs_edit["include"] & pairs_edit["air"].notna()]
dup_ids = sel_pairs["site_id"][sel_pairs["site_id"].duplicated()].tolist()
if dup_ids:
    st.error("Site IDs must be unique: " + ", ".join(dup_ids))
    st.stop()
if sel_pairs.empty:
    st.info("Tick at least one pair to run.")
    st.stop()


@_fragment
def _pair_check_plot():
    opts = sel_pairs["site_id"].tolist()
    pick = st.selectbox("Check a pair", opts, key="pair_check")
    row = sel_pairs[sel_pairs["site_id"] == pick].iloc[0]
    fig = go.Figure()
    for n, col, nm in ((row["air"], AIR_COLOR, "air: "), (row["stream"], WATER_COLOR, "stream: ")):
        s = _downsample(grids[n])
        fig.add_trace(go.Scattergl(x=s.index, y=s.to_numpy(), mode="lines", name=nm + n,
                                   line={"width": 1.2, "color": col}))
    fig.update_layout(height=300, yaxis_title="Temperature (°C)",
                      margin={"l": 20, "r": 20, "t": 20, "b": 20})
    st.plotly_chart(fig, width="stretch")


with st.expander("Pair check plot", expanded=False):
    _pair_check_plot()

# ---- 4: Settings & run -------------------------------------------------------
st.markdown("---")
st.subheader("4 - Run Diel PASTA")

with st.expander("Settings", expanded=True):
    s1, s2, s3 = st.columns(3)
    with s1:
        water_min = st.number_input(
            "Drop stream records ≤ (°C)", value=0.99, step=0.01, format="%.2f",
            help="R: filter(str_temp > 0.99). Removes ice-affected / out-of-water records.")
        all_times = pd.concat([feats["start"], feats["end"]])
        d0, d1 = all_times.min().date(), all_times.max().date()
        date_range = st.date_input("Date range", value=(d0, d1), min_value=d0, max_value=d1)
    with s2:
        min_cov = st.slider(
            "Minimum day coverage (%)", 0, 100, 75, 5,
            help="Days with fewer records than this share of a full day are not fitted. "
                 "Set 0 to fit every day, as the R script does (partial deploy/retrieval "
                 "days then usually fail the R² filter).")
        r2_min = st.slider("Minimum adjusted R² (air and water)", 0.0, 1.0, 0.70, 0.05,
                           help="Days kept only when both fits exceed this (R: AdjRsqr > 0.7).")
    with s3:
        roll_k = st.number_input("Rolling mean window (days)", value=7, min_value=1, step=2,
                                 help="Centred rolling mean over kept days (R: zoo::rollmean k = 7).")
        st.caption("Phase comes from the fitted sinusoid using the logger clock time "
                   "as written. The R `rad_hour()` converts to the computer's time zone, which shifts phase by 1 h across DST.")

if isinstance(date_range, (list, tuple)) and len(date_range) == 2:
    dr0, dr1 = pd.Timestamp(date_range[0]), pd.Timestamp(date_range[1]) + pd.Timedelta(days=1)
else:
    dr0, dr1 = pd.Timestamp(d0), pd.Timestamp(d1) + pd.Timedelta(days=1)

pair_list = [(r.site_id, r.stream, r.air) for r in sel_pairs.itertuples()]
run_sig = _sig(score_sig, tuple(pair_list), water_min, min_cov, r2_min, roll_k,
               str(dr0), str(dr1))

if st.button("Run Diel PASTA", type="primary"):
    with st.spinner("Fitting daily sinusoids..."):
        try:
            clipped = {n: g[(g.index >= dr0) & (g.index < dr1)] for n, g in grids.items()}
            T_y = pdl.build_paired_df(pair_list, clipped, water_min_C=water_min)
            fits = pdl.fit_daily_sinusoids(T_y, min_coverage=min_cov / 100,
                                           interval_min=int(interval_min))
            metrics = pdl.pasta_metrics(fits, r2_min=r2_min, roll_k=int(roll_k))
            run_ts = pd.Timestamp.now()
            summary = pd.DataFrame([{
                "run_timestamp": run_ts.isoformat(timespec="seconds"),
                "site_id": sid, "stream_logger": s, "air_logger": a,
                "interval_min": int(interval_min), "water_min_C": water_min,
                "min_day_coverage_pct": min_cov, "adjR2_min": r2_min, "roll_k": int(roll_k),
                "date_from": dr0.date().isoformat(),
                "date_to": (dr1 - pd.Timedelta(days=1)).date().isoformat(),
                "n_paired_records": int((T_y["site_id"] == sid).sum()),
                "n_days_fitted": int(fits.loc[fits["site_id"] == sid, "date"].nunique()) if len(fits) else 0,
                "n_days_kept": int((metrics["site_id"] == sid).sum()) if len(metrics) else 0,
                "median_AmpRatio": metrics.loc[metrics["site_id"] == sid, "AmpRatio"].median() if len(metrics) else np.nan,
                "median_PhaseLag_h": metrics.loc[metrics["site_id"] == sid, "PhaseLag_h"].median() if len(metrics) else np.nan,
                "median_Ratio_Mean_K": metrics.loc[metrics["site_id"] == sid, "Ratio_Mean_K"].median() if len(metrics) else np.nan,
            } for sid, s, a in pair_list])
            st.session_state["res"] = {"T_y": T_y, "fits": fits, "metrics": metrics,
                                       "summary": summary, "sig": run_sig, "run_ts": run_ts}
        except Exception as exc:
            st.error("Run failed: " + str(exc))
            import traceback
            with st.expander("Traceback"):
                st.code(traceback.format_exc())

if "res" not in st.session_state:
    st.stop()

res = st.session_state["res"]
T_y, fits, metrics, summary = res["T_y"], res["fits"], res["metrics"], res["summary"]
if res["sig"] != run_sig:
    st.warning("Pairs or settings changed since the last run. Press **Run Diel PASTA** "
               "to update the results below.")

# ---- 5: Results --------------------------------------------------------------
st.markdown("---")
st.subheader("5 - Results")
st.caption("Run time: " + res["run_ts"].strftime("%Y-%m-%d %H:%M:%S"))

st.dataframe(
    summary[["site_id", "stream_logger", "air_logger", "n_days_fitted", "n_days_kept",
             "median_AmpRatio", "median_PhaseLag_h", "median_Ratio_Mean_K"]]
    .rename(columns={"n_days_fitted": "days fitted", "n_days_kept": "days kept (both R² ok)",
                     "median_AmpRatio": "median AmpRatio",
                     "median_PhaseLag_h": "median phase lag (h)",
                     "median_Ratio_Mean_K": "median mean ratio (K)"}),
    hide_index=True, width="stretch",
)

if metrics.empty:
    st.warning("No days passed the R² filter. Lower the minimum adjusted R² or day coverage.")
else:
    @_fragment
    def _metric_plots():
        c1, c2 = st.columns([3, 1])
        with c2:
            x_mode = st.radio("x axis", ["Date", "Day of year"], key="xmode")
            ar_cap = st.checkbox("Hide AmpRatio ≥ 1.2", value=True, key="arcap",
                                 help="Values this high are usually poor data (Rmd plot filter).")
        xcol = "date" if x_mode == "Date" else "DOY"
        m = metrics[metrics["AmpRatio"] < 1.2] if ar_cap else metrics
        fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.05,
                            subplot_titles=["Amplitude ratio (water / air)",
                                            "Phase lag (h, water − air)",
                                            "Mean ratio (K, water / air)"])
        for i, (sid, d) in enumerate(m.groupby("site_id", sort=False)):
            col = SITE_COLORS[i % len(SITE_COLORS)]
            for row, y, roll in ((1, "AmpRatio", "AR_roll7"), (2, "PhaseLag_h", "PH_roll7"),
                                 (3, "Ratio_Mean_K", None)):
                fig.add_trace(go.Scattergl(x=d[xcol], y=d[y], mode="markers", name=sid,
                                           legendgroup=sid, showlegend=(row == 1),
                                           marker={"size": 5, "color": col, "opacity": 0.7}),
                              row=row, col=1)
                if roll:
                    fig.add_trace(go.Scatter(x=d[xcol], y=d[roll], mode="lines", name=sid + " roll",
                                             legendgroup=sid, showlegend=False,
                                             line={"color": col, "width": 2}), row=row, col=1)
        fig.update_layout(height=720, margin={"l": 20, "r": 20, "t": 40, "b": 20})
        with c1:
            st.plotly_chart(fig, width="stretch")

    _metric_plots()


@_fragment
def _day_fit_plot():
    st.markdown("#### Daily fit check")
    c1, c2 = st.columns(2)
    with c1:
        sid = st.selectbox("Site", summary["site_id"].tolist(), key="fit_site")
    fs = fits[fits["site_id"] == sid]
    days = sorted(fs["date"].dt.date.unique())
    if not days:
        st.info("No fitted days for this site.")
        return
    with c2:
        day = st.select_slider("Day", options=days, value=days[len(days) // 2], key="fit_day")
    obs = T_y[(T_y["site_id"] == sid) & (pd.to_datetime(T_y["DateTime_r"]).dt.date == day)]
    fig = go.Figure()
    txt = []
    for medium, col, ycol in (("air", AIR_COLOR, "air_temp"), ("water", WATER_COLOR, "str_temp")):
        fr = fs[(fs["date"].dt.date == day) & (fs["medium"] == medium)]
        fig.add_trace(go.Scatter(x=obs["DateTime_r"], y=obs[ycol], mode="markers", name=medium,
                                 marker={"size": 4, "color": col, "opacity": 0.6}))
        if len(fr):
            r = fr.iloc[0]
            fig.add_trace(go.Scatter(x=obs["DateTime_r"], y=pdl.sinusoid_value(r, obs["DateTime_r"]),
                                     mode="lines", name=medium + " fit",
                                     line={"color": col, "width": 2}))
            txt.append("{}: amp {:.2f} °C, phase {:.2f} h, adj R² {:.2f}".format(
                medium, r["amplitude_C"], r["phase_h"], r["AdjRsqr"]))
    fig.update_layout(height=320, yaxis_title="Temperature (°C)",
                      margin={"l": 20, "r": 20, "t": 20, "b": 20})
    st.plotly_chart(fig, width="stretch")
    st.caption("  |  ".join(txt))


_day_fit_plot()

with st.expander("Daily metrics table"):
    st.dataframe(metrics.round(4), width="stretch")
with st.expander("Daily sinusoid fits (air and water)"):
    st.dataframe(fits.round(4), width="stretch")

# ---- 6: Download -------------------------------------------------------------
st.markdown("---")
st.subheader("6 - Download")
stem = "dielPASTA_" + res["run_ts"].strftime("%Y%m%d_%H%M")
outputs = {
    stem + "_metrics.csv": metrics,
    stem + "_fits.csv": fits,
    stem + "_paired_df.csv": T_y,
    stem + "_run_summary.csv": summary,
}
dcols = st.columns(4)
labels = ["Metrics CSV", "Fits CSV", "Paired data CSV", "Run summary CSV"]
for col, (fn, df), lab in zip(dcols, outputs.items(), labels):
    with col:
        st.download_button("⬇️ " + lab, data=_df_to_csv(df), file_name=fn, mime="text/csv")

o1, o2 = st.columns([3, 1])
with o1:
    default_out = str(Path(folder) / "dielPASTA_output") if source_mode == "📂 Local folder path" \
        else str(Path(__file__).parent / "output")
    out_dir = st.text_input("Or save all four to a folder", value=default_out, key="out_dir")
with o2:
    st.write("")
    st.write("")
    if st.button("💾 Save to folder"):
        try:
            Path(out_dir).mkdir(parents=True, exist_ok=True)
            for fn, df in outputs.items():
                df.to_csv(Path(out_dir) / fn, index=False)
            st.success("Saved " + str(len(outputs)) + " files to `" + out_dir + "`")
        except Exception as exc:
            st.error("Could not save: " + str(exc))
