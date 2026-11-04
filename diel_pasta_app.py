"""
PASTA -- Streamlit Application
Paired Air and Stream Temperature Analysis at the annual and diel scales.

Upload air and stream temperature logger files (HOBOware exports or simple
date/time + temperature CSVs); the app sorts loggers into air vs stream,
suggests air-stream pairs, then runs either or both analyses:
  - Annual signal: one sinusoid per year fitted to daily means (AmpRatio,
    PhaseLag_d, Ratio_Mean, thermal sensitivity), whole record and per year.
  - Diel: one sinusoid per day fitted to the logger records (AmpRatio,
    PhaseLag_h, mean ratios) for every day.

Usage:
    streamlit run diel_pasta_app.py        (or double-click run_diel_pasta.bat)

Analysis code lives in the pasta_signals package (run_annual, run_diel); see README.
"""

# -- Auto-install missing packages from requirements.txt --------
import importlib
import subprocess
import sys
import os
from pathlib import Path

_REQUIREMENTS_FILE = Path(__file__).parent / "requirements.txt"


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
import zipfile
import datetime as _dt

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pasta_signals as ps
from pasta_signals import annual as an, diel as dl, pairing as pr, readers as rd

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
    page_title="PASTA",
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
    raw = rd.parse_logger_csv(file_bytes, fname)
    grid = rd.to_grid(raw, interval_min)
    feats = pr.logger_features(grid)
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
            out[(s, a)] = pr.diel_xcorr(_grids[s], _grids[a], interval_min)
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
def _coverage(sig, _grids, pairs, interval_min):
    """Start / end of each logger and of the paired (both present) record, per site."""
    per_day = 24 * 60 / interval_min
    rows = []
    for sid, s, a in pairs:
        gs, ga = _grids[s].dropna(), _grids[a].dropna()
        both = gs.index.intersection(ga.index)
        rows.append({"site_id": sid, "stream": s, "air": a,
                     "stream_start": gs.index.min(), "stream_end": gs.index.max(),
                     "air_start": ga.index.min(), "air_end": ga.index.max(),
                     "paired_start": both.min() if len(both) else pd.NaT,
                     "paired_end": both.max() if len(both) else pd.NaT,
                     "paired_days": len(both) / per_day,
                     "paired_index": both})
    return pd.DataFrame(rows)


def _coverage_in_range(cov, dr0, dr1):
    """Per site: paired days inside [dr0, dr1) and whether the site spans the range."""
    out = cov.drop(columns=["paired_index"]).copy()
    n_days = max((dr1 - dr0).days, 1)
    tol = pd.Timedelta(days=1)
    in_range, status = [], []
    for r in cov.itertuples():
        idx = r.paired_index
        k = idx[(idx >= dr0) & (idx < dr1)]
        days = len(pd.unique(k.normalize())) if len(k) else 0
        in_range.append(days)
        if days == 0:
            status.append("no data in range")
        elif r.paired_start > dr0 + tol and r.paired_end < dr1 - 2 * tol:
            status.append("starts late and ends early")
        elif r.paired_start > dr0 + tol:
            status.append("starts {}".format(r.paired_start.date()))
        elif r.paired_end < dr1 - 2 * tol:
            status.append("ends {}".format(r.paired_end.date()))
        elif days < 0.9 * n_days:
            status.append("gaps inside range")
        else:
            status.append("covers whole range")
    out["days_in_range"] = in_range
    out["pct_of_range"] = [100 * d / n_days for d in in_range]
    out["status"] = status
    return out


def _coverage_plot(cov, dr0, dr1):
    fig = go.Figure()
    fig.add_vrect(x0=dr0, x1=dr1, fillcolor="rgba(31,119,180,0.08)", line_width=0,
                  annotation_text="selected range", annotation_position="top left")
    for r in cov.itertuples():
        for y_off, start, end, col, nm in ((0.18, r.air_start, r.air_end, AIR_COLOR, "air"),
                                           (0.0, r.paired_start, r.paired_end, "#333333", "paired"),
                                           (-0.18, r.stream_start, r.stream_end, WATER_COLOR, "stream")):
            if pd.isna(start):
                continue
            fig.add_trace(go.Scatter(x=[start, end], y=[r.site_id, r.site_id], mode="lines",
                                     line={"color": col, "width": 9 if nm == "paired" else 4},
                                     name=nm, legendgroup=nm, showlegend=(r.Index == 0),
                                     hovertemplate=nm + ": %{x|%Y-%m-%d}<extra>" + r.site_id + "</extra>"))
    fig.update_layout(height=max(160, 46 * len(cov) + 80),
                      margin={"l": 20, "r": 20, "t": 30, "b": 20},
                      legend={"orientation": "h", "y": 1.12})
    st.plotly_chart(fig, width="stretch")


def _clean_out_dir(raw):
    """User-typed folder -> absolute Path. Accepts quotes (Copy as path), ~, %VARS%."""
    txt = os.path.expandvars(str(raw).strip().strip('"').strip("'").strip())
    if not txt:
        raise ValueError("Type a folder path, e.g. C:\\Users\\you\\Documents\\PASTA_output")
    if os.name != "nt" and len(txt) > 1 and txt[1] == ":":
        raise ValueError("This app is running on a server, not on your computer, so it cannot "
                         "save to " + txt + ". Use **Download all (.zip)** instead.")
    path = Path(txt).expanduser()
    if os.name == "nt" and path.drive and not Path(path.drive + "\\").exists():
        raise ValueError("Drive " + path.drive + " was not found.")
    return path.resolve()


@st.cache_data(show_spinner=False)
def _df_to_csv(df):
    """Convert DataFrame to UTF-8 CSV bytes (cached)."""
    return df.to_csv(index=False).encode("utf-8")


def _r4(df):
    """Round numeric columns for display."""
    return df.round({c: 4 for c in df.select_dtypes("number").columns})


def _sig(*parts):
    return hashlib.md5(repr(parts).encode()).hexdigest()


# ==============================================================================
#  UI
# ==============================================================================

st.title("PASTA -- Paired Air & Stream Temperature Analysis")
st.caption(
    "Annual and diel sinusoid fits to paired air and stream temperature  |  "
    "amplitude ratio · phase lag · mean ratio · thermal sensitivity  |  "
    "port of the PASTA and Diel PASTA R code"
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
    name = rd.logger_name_from_filename(fname)
    if name in grids:  # two files with the same logger name -> keep both, disambiguate
        name = Path(fname).stem
    try:
        g, f = _load_logger(b, fname, int(interval_min))
    except Exception as exc:
        errors.append("**" + fname + "**: " + str(exc))
        continue
    grids[name] = g
    kind, conf, why = pr.classify_logger(name, f)
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
    "Each stream logger is paired with one air logger. **Names are checked first**, "
    "strongest rule first: the same name apart from the air/stream code "
    "(AT_002 ↔ ST_002, LLC01_AT ↔ LLC01_WT), then the same site name and number, "
    "the same number, the same last two digits (HF_003 ↔ HF_103), and the same prefix. "
    "The data is used only as a last resort: for streams with no name match, or to choose "
    "between several air loggers that match equally well. It uses the diel correlation, "
    "i.e. how closely the stream's daily cycle follows the air's (0–6 h lag). "
    "One air logger can serve several streams. Change **Air logger** to override."
)

score_sig = _sig(tuple(sorted((f, len(b)) for f, b in files)), interval_min)
scores = _pair_scores(score_sig, grids, tuple(streams), tuple(airs), int(interval_min))
suggest = pr.suggest_pairs(streams, airs, {k: v[0] for k, v in scores.items()})

pair_rows = []
for s in streams:
    a, basis = suggest[s]
    check = pr.basis_note(s, a, basis, airs, {k: v[0] for k, v in scores.items()})
    r, lag, days = scores.get((s, a), (np.nan, np.nan, np.nan)) if a else (np.nan, np.nan, np.nan)
    pair_rows.append({"include": a is not None, "site_id": pr.site_label(s), "stream": s, "air": a,
                      "basis": check, "diel_r": r, "lag_h": lag, "overlap_days": days})
pairs_df = pd.DataFrame(pair_rows)
_dup = pairs_df["site_id"].duplicated(keep=False)          # shortened labels collided
pairs_df.loc[_dup, "site_id"] = pairs_df.loc[_dup, "stream"]

pairs_edit = st.data_editor(
    pairs_df,
    column_config={
        "include": st.column_config.CheckboxColumn("Run", width="small"),
        "site_id": st.column_config.TextColumn("Site ID", help="Label used in the outputs"),
        "stream": st.column_config.TextColumn("Stream logger", disabled=True),
        "air": st.column_config.SelectboxColumn("Air logger", options=airs, required=True),
        "basis": st.column_config.TextColumn("Suggested by", disabled=True),
        "diel_r": st.column_config.NumberColumn("Diel r", format="%.3f", disabled=True,
                                                help="Diel correlation of the suggested pair"),
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
st.subheader("4 - Run PASTA")

analyses = st.multiselect(
    "Analyses", ["Annual signal", "Diel"], default=["Annual signal", "Diel"],
    help="**Annual signal**: one sinusoid per year fitted to daily means "
         "(PASTA Shiny app). **Diel**: one sinusoid per day fitted to the logger "
         "records (Diel PASTA).",
)
do_annual = "Annual signal" in analyses
do_diel = "Diel" in analyses
if not analyses:
    st.info("Select at least one analysis.")
    st.stop()

pair_list = [(r.site_id, r.stream, r.air) for r in sel_pairs.itertuples()]
cov = _coverage(_sig(score_sig, tuple(pair_list)), grids, tuple(pair_list), int(interval_min))

st.markdown("#### Date range")
range_mode = st.radio(
    "Date range", ["All data", "Common period (all sites)", "Custom"], horizontal=True,
    key="range_mode", label_visibility="collapsed",
    help="**All data**: every site uses its full paired record. **Common period**: only "
         "dates where every selected site has paired air and stream data. **Custom**: pick dates.")
c_first, c_last = cov["paired_start"].min(), cov["paired_end"].max()
common0, common1 = cov["paired_start"].max(), cov["paired_end"].min()
if range_mode == "All data":
    dr0, dr1 = c_first.normalize(), c_last.normalize() + pd.Timedelta(days=1)
elif range_mode == "Common period (all sites)":
    if common0 >= common1:
        st.error("No period is shared by all selected sites. Choose **All data**, "
                 "**Custom**, or untick sites in the pair table.")
        st.stop()
    dr0, dr1 = common0.normalize(), common1.normalize() + pd.Timedelta(days=1)
else:
    d0, d1 = c_first.date(), c_last.date()
    date_range = st.date_input("Custom dates", value=(d0, d1), min_value=d0, max_value=d1,
                               key="custom_dates")
    if isinstance(date_range, (list, tuple)) and len(date_range) == 2:
        dr0, dr1 = pd.Timestamp(date_range[0]), pd.Timestamp(date_range[1]) + pd.Timedelta(days=1)
    else:
        dr0, dr1 = pd.Timestamp(d0), pd.Timestamp(d1) + pd.Timedelta(days=1)

cov_view = _coverage_in_range(cov, dr0, dr1)
n_full = int((cov_view["status"] == "covers whole range").sum())
range_txt = "{} to {} ({} days)".format(dr0.date(), (dr1 - pd.Timedelta(days=1)).date(),
                                        (dr1 - dr0).days)
if n_full == len(cov_view):
    st.success("All {} sites have paired data across the whole range: {}.".format(
        len(cov_view), range_txt))
else:
    short = cov_view[cov_view["status"] != "covers whole range"]
    st.warning(
        "The range {} is **not covered by all sites**: {} of {} sites have data for only "
        "part of it ({}). Sites will be analysed over different periods, so their annual "
        "metrics may not be comparable. Choose **Common period (all sites)** to use only the "
        "dates every site shares.".format(
            range_txt, len(short), len(cov_view),
            ", ".join("{} ({})".format(r.site_id, r.status) for r in short.itertuples())))
with st.expander("Data coverage by site", expanded=n_full < len(cov_view)):
    _coverage_plot(cov, dr0, dr1)
    st.dataframe(cov_view, hide_index=True, width="stretch",
                 column_config={
                     "stream_start": st.column_config.DatetimeColumn("Stream start", format="YYYY-MM-DD"),
                     "stream_end": st.column_config.DatetimeColumn("Stream end", format="YYYY-MM-DD"),
                     "air_start": st.column_config.DatetimeColumn("Air start", format="YYYY-MM-DD"),
                     "air_end": st.column_config.DatetimeColumn("Air end", format="YYYY-MM-DD"),
                     "paired_start": st.column_config.DatetimeColumn("Paired start", format="YYYY-MM-DD"),
                     "paired_end": st.column_config.DatetimeColumn("Paired end", format="YYYY-MM-DD"),
                     "paired_days": st.column_config.NumberColumn("Paired days (all)", format="%.0f"),
                     "days_in_range": st.column_config.NumberColumn("Paired days in range", format="%.0f"),
                     "pct_of_range": st.column_config.NumberColumn("% of range", format="%.0f"),
                 })
    st.caption("Paired = times with both the stream and the air record. Common period "
               "for all sites: {}.".format(
                   "{} to {}".format(common0.date(), common1.date()) if common0 < common1
                   else "none"))

with st.expander("Settings", expanded=True):
    s1, s2, s3 = st.columns(3)
    with s1:
        st.markdown("**Both**")
        water_min = st.number_input(
            "Drop stream records ≤ (°C)", value=0.99, step=0.01, format="%.2f",
            help="R: filter(str_temp > 0.99). Removes ice-affected / out-of-water records. "
                 "The annual analysis also drops days with a stream mean ≤ 1 °C, as PASTA does.")
        min_cov = st.slider(
            "Minimum day coverage (%)", 0, 100, 75, 5,
            help="Days with fewer paired records than this share of a full day are not "
                 "fitted (diel) and get no daily mean (annual). Set 0 to use every day.")
    with s2:
        st.markdown("**Annual signal**")
        yr_type = st.radio("Year type", ["water", "calendar"], horizontal=True,
                           disabled=not do_annual,
                           help="Water year runs 1 Oct – 30 Sep and is named by the year it ends.")
        min_days = st.number_input("Minimum days per year", value=100, min_value=1, step=10,
                                   disabled=not do_annual,
                                   help="Years with fewer daily means get no AmpRatio, phase lag "
                                        "or mean ratio (R: count < 100).")
    with s3:
        st.markdown("**Diel**")
        r2_min = st.slider("Minimum adjusted R² (air and water)", 0.0, 1.0, 0.70, 0.05,
                           disabled=not do_diel,
                           help="Days kept only when both fits exceed this (R: AdjRsqr > 0.7).")
        roll_k = st.number_input("Rolling mean window (days)", value=7, min_value=1, step=2,
                                 disabled=not do_diel,
                                 help="Centred rolling mean over kept days (R: zoo::rollmean k = 7).")
        st.caption("Diel phase uses the logger clock time as written. The R `rad_hour()` "
                   "converts to the computer's time zone, which shifts phase by 1 h across DST.")

run_sig = _sig(score_sig, tuple(pair_list), tuple(analyses), water_min, min_cov, yr_type,
               min_days, r2_min, roll_k, str(dr0), str(dr1))


def _site_med(df, sid, col):
    return df.loc[df["site_id"] == sid, col].median() if len(df) else np.nan


if st.button("Run PASTA", type="primary"):
    with st.spinner("Fitting sinusoids..."):
        try:
            clipped = {n: g[(g.index >= dr0) & (g.index < dr1)] for n, g in grids.items()}
            data = ps.pair_data(clipped, pair_list)
            T_y = data[data["str_temp"] > water_min].reset_index(drop=True)
            res = {"T_y": T_y, "analyses": list(analyses), "yr_type": yr_type,
                   "min_days": int(min_days),
                   "pairs": sel_pairs.drop(columns=["include"]).reset_index(drop=True),
                   "coverage": cov_view}

            if do_diel:
                dres = ps.run_diel(data, interval_min=int(interval_min), water_min_C=water_min,
                                   min_coverage=min_cov / 100, r2_min=r2_min, roll_k=int(roll_k))
                fits, metrics = dres.fits, dres.metrics
                res.update(fits=fits, metrics=metrics)
            if do_annual:
                ares = ps.run_annual(data, yr_type=yr_type, min_days=int(min_days),
                                     water_min_C=water_min, min_coverage=min_cov / 100,
                                     interval_min=int(interval_min))
                ann = {"record": ares.record, "yearly": ares.yearly, "fits": ares.fits}
                res.update(daily=ares.daily, ann=ann)

            run_ts = pd.Timestamp.now()
            rows = []
            for sid, s, a in pair_list:
                row = {"run_timestamp": run_ts.isoformat(timespec="seconds"),
                       "site_id": sid, "stream_logger": s, "air_logger": a,
                       "analyses": " + ".join(analyses),
                       "interval_min": int(interval_min), "water_min_C": water_min,
                       "min_day_coverage_pct": min_cov,
                       "date_from": dr0.date().isoformat(),
                       "date_to": (dr1 - pd.Timedelta(days=1)).date().isoformat(),
                       "n_paired_records": int((T_y["site_id"] == sid).sum())}
                if do_annual:
                    rec = ann["record"]
                    r = rec[rec["site_id"] == sid] if len(rec) else rec
                    row.update({"year_type": yr_type, "min_days_per_year": int(min_days),
                                "annual_n_days": int(r["count"].iloc[0]) if len(r) else 0})
                    for c in ["AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope",
                              "TS_AdjRsqr", "max_conseq_missing_days"]:
                        row["annual_" + c] = r[c].iloc[0] if len(r) else np.nan
                if do_diel:
                    row.update({
                        "adjR2_min": r2_min, "roll_k": int(roll_k),
                        "diel_n_days_fitted": int(fits.loc[fits["site_id"] == sid, "date"].nunique()) if len(fits) else 0,
                        "diel_n_days_kept": int((metrics["site_id"] == sid).sum()) if len(metrics) else 0,
                        "diel_median_AmpRatio": _site_med(metrics, sid, "AmpRatio"),
                        "diel_median_PhaseLag_h": _site_med(metrics, sid, "PhaseLag_h"),
                        "diel_median_Ratio_Mean_K": _site_med(metrics, sid, "Ratio_Mean_K"),
                    })
                rows.append(row)
            res.update(summary=pd.DataFrame(rows), sig=run_sig, run_ts=run_ts)
            st.session_state["res"] = res
        except Exception as exc:
            st.error("Run failed: " + str(exc))
            import traceback
            with st.expander("Traceback"):
                st.code(traceback.format_exc())

if "res" not in st.session_state:
    st.stop()

res = st.session_state["res"]
if "analyses" not in res:          # results from an older app version in this session
    st.session_state.pop("res")
    st.stop()
T_y, summary = res["T_y"], res["summary"]
ran_annual = "Annual signal" in res["analyses"]
ran_diel = "Diel" in res["analyses"]
if res["sig"] != run_sig:
    st.warning("Analyses, pairs or settings changed since the last run. Press **Run PASTA** "
               "to update the results below.")

# ---- 5: Results --------------------------------------------------------------
st.markdown("---")
st.subheader("5 - Results")
st.caption("Run time: " + res["run_ts"].strftime("%Y-%m-%d %H:%M:%S"))

sum_cols = {"site_id": "site", "stream_logger": "stream", "air_logger": "air"}
if ran_annual:
    sum_cols.update({"annual_AmpRatio": "annual AmpRatio", "annual_PhaseLag_d": "annual phase lag (d)",
                     "annual_Ratio_Mean": "annual mean ratio", "annual_TS__Slope": "thermal sensitivity"})
if ran_diel:
    sum_cols.update({"diel_n_days_kept": "diel days kept",
                     "diel_median_AmpRatio": "diel median AmpRatio",
                     "diel_median_PhaseLag_h": "diel median phase lag (h)"})
st.dataframe(summary[list(sum_cols)].rename(columns=sum_cols), hide_index=True, width="stretch")

tab_names = (["Annual signal"] if ran_annual else []) + (["Diel"] if ran_diel else [])
tabs = dict(zip(tab_names, st.tabs(tab_names)))

# -- Annual results --
if ran_annual:
    daily, ann = res["daily"], res["ann"]
    ytype = res["yr_type"]

    @_fragment
    def _annual_plots():
        sid = st.selectbox("Site", summary["site_id"].tolist(), key="ann_site")
        d = daily[daily["site_id"] == sid]
        f = ann["fits"]
        f = f[(f["site_id"] == sid) & (f["year"] == "all")] if len(f) else f
        if d.empty:
            st.info("No daily means for this site.")
            return
        c1, c2 = st.columns([2, 1])
        fig = go.Figure()
        for med, col, ycol in (("air", AIR_COLOR, "tavg_air_C"), ("water", WATER_COLOR, "tavg_wat_C")):
            fig.add_trace(go.Scattergl(x=d["date"], y=d[ycol], mode="markers", name=med + " daily mean",
                                       marker={"size": 4, "color": col, "opacity": 0.5}))
            fr = f[f["medium"] == med] if len(f) else f
            if len(fr):
                fig.add_trace(go.Scatter(x=d["date"], y=an.sinusoid_value(fr.iloc[0], d["date"], ytype),
                                         mode="lines", name=med + " annual fit",
                                         line={"color": col, "width": 2}))
        fig.update_layout(height=340, yaxis_title="Daily mean temperature (°C)",
                          margin={"l": 20, "r": 20, "t": 20, "b": 20})
        with c1:
            st.plotly_chart(fig, width="stretch")
        rec = ann["record"][ann["record"]["site_id"] == sid] if len(ann["record"]) else ann["record"]
        fig2 = go.Figure()
        fig2.add_trace(go.Scattergl(x=d["tavg_air_C"], y=d["tavg_wat_C"], mode="markers",
                                    name="days", marker={"size": 4, "color": WATER_COLOR, "opacity": 0.5}))
        lo, hi = float(d["tavg_air_C"].min()), float(d["tavg_air_C"].max())
        fig2.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", name="1:1",
                                  line={"dash": "dash", "color": "gray", "width": 1}))
        if len(rec) and np.isfinite(rec["TS__Slope"].iloc[0]):
            sl, yi = rec["TS__Slope"].iloc[0], rec["TS_YInt"].iloc[0]
            fig2.add_trace(go.Scatter(x=[lo, hi], y=[yi + sl * lo, yi + sl * hi], mode="lines",
                                      name="slope {:.2f}".format(sl),
                                      line={"color": "black", "width": 2}))
        fig2.update_layout(height=340, xaxis_title="Air (°C)", yaxis_title="Stream (°C)",
                           margin={"l": 20, "r": 20, "t": 20, "b": 20},
                           legend={"orientation": "h", "y": -0.25})
        with c2:
            st.plotly_chart(fig2, width="stretch")

    with tabs["Annual signal"]:
        if len(daily):
            span = daily.groupby("site_id")["date"].agg(lambda x: (x.max() - x.min()).days + 1)
            short_sites = span[span < 300]
            if len(short_sites):
                st.warning(
                    "Less than ~300 days of record for: "
                    + ", ".join("{} ({} d)".format(k, v) for k, v in short_sites.items())
                    + ". An annual sinusoid fitted to part of a year is poorly constrained; "
                      "treat these annual metrics with caution.")
        st.markdown("#### Whole record")
        st.dataframe(ann["record"], hide_index=True, width="stretch")
        st.caption("AmpRatio and mean ratio are water / air; phase lag is water − air in days. "
                   "TS__Slope is the daily thermal sensitivity (stream ~ air, Kelleher et al. 2012).")
        st.markdown("#### By " + ytype + " year")
        st.dataframe(ann["yearly"], hide_index=True, width="stretch")
        st.caption("Years with fewer than {} daily means show no signal metrics.".format(res["min_days"]))
        _annual_plots()
        yl = ann["yearly"].dropna(subset=["AmpRatio", "PhaseLag_d"]) if len(ann["yearly"]) else ann["yearly"]
        if len(yl):
            fig3 = go.Figure()
            for i, (sid, d) in enumerate(yl.groupby("site_id", sort=False)):
                fig3.add_trace(go.Scatter(x=d["PhaseLag_d"], y=d["AmpRatio"], mode="markers+text",
                                          text=d["year"].astype(str), textposition="top center",
                                          name=sid, marker={"size": 9, "color": SITE_COLORS[i % len(SITE_COLORS)]}))
            fig3.update_layout(height=340, xaxis_title="Phase lag (d)", yaxis_title="Amplitude ratio",
                               margin={"l": 20, "r": 20, "t": 30, "b": 20},
                               title={"text": "Annual signal by year", "font": {"size": 14}})
            st.plotly_chart(fig3, width="stretch")
        with st.expander("Annual sinusoid fits"):
            st.dataframe(_r4(ann["fits"]), hide_index=True, width="stretch")

# -- Diel results --
if ran_diel:
    fits, metrics = res["fits"], res["metrics"]

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

    @_fragment
    def _day_fit_plot():
        st.markdown("#### Daily fit check")
        c1, c2 = st.columns(2)
        with c1:
            sid = st.selectbox("Site", summary["site_id"].tolist(), key="fit_site")
        fs = fits[fits["site_id"] == sid] if len(fits) else fits
        days = sorted(fs["date"].dt.date.unique()) if len(fs) else []
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
                fig.add_trace(go.Scatter(x=obs["DateTime_r"], y=dl.sinusoid_value(r, obs["DateTime_r"]),
                                         mode="lines", name=medium + " fit",
                                         line={"color": col, "width": 2}))
                txt.append("{}: amp {:.2f} °C, phase {:.2f} h, adj R² {:.2f}".format(
                    medium, r["amplitude_C"], r["phase_h"], r["AdjRsqr"]))
        fig.update_layout(height=320, yaxis_title="Temperature (°C)",
                          margin={"l": 20, "r": 20, "t": 20, "b": 20})
        st.plotly_chart(fig, width="stretch")
        st.caption("  |  ".join(txt))

    with tabs["Diel"]:
        if metrics.empty:
            st.warning("No days passed the R² filter. Lower the minimum adjusted R² or day coverage.")
        else:
            _metric_plots()
        _day_fit_plot()
        with st.expander("Daily metrics table"):
            st.dataframe(_r4(metrics), width="stretch")
        with st.expander("Daily sinusoid fits (air and water)"):
            st.dataframe(_r4(fits), width="stretch")

# ---- 6: Download -------------------------------------------------------------
st.markdown("---")
st.subheader("6 - Download")
stem = "PASTA_" + res["run_ts"].strftime("%Y%m%d_%H%M")
outputs = {stem + "_run_summary.csv": ("Run summary", summary),
           stem + "_pairs.csv": ("Pairs used", res.get("pairs", pd.DataFrame())),
           stem + "_coverage.csv": ("Data coverage", res.get("coverage", pd.DataFrame())),
           stem + "_paired_df.csv": ("Paired data", T_y)}
if ran_annual:
    outputs.update({
        stem + "_annual_metrics.csv": ("Annual metrics (record)", res["ann"]["record"]),
        stem + "_annual_metrics_yearly.csv": ("Annual metrics (by year)", res["ann"]["yearly"]),
        stem + "_annual_fits.csv": ("Annual fits", res["ann"]["fits"]),
        stem + "_daily_means.csv": ("Daily means", res["daily"]),
    })
if ran_diel:
    outputs.update({
        stem + "_diel_metrics.csv": ("Diel metrics", res["metrics"]),
        stem + "_diel_fits.csv": ("Diel fits", res["fits"]),
    })


@st.cache_data(show_spinner=False)
def _zip_outputs(sig, _outputs):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for fn, (_, df) in _outputs.items():
            z.writestr(fn, df.to_csv(index=False))
    return buf.getvalue()


st.download_button(
    "⬇️ Download all ({} files, .zip)".format(len(outputs)),
    data=_zip_outputs(res["sig"] + stem + str(len(outputs)), outputs),
    file_name=stem + ".zip", mime="application/zip", type="primary", key="dl_zip")

with st.expander("Download single files"):
    items = list(outputs.items())
    for i in range(0, len(items), 4):
        cols = st.columns(4)
        for col, (fn, (lab, df)) in zip(cols, items[i:i + 4]):
            with col:
                st.download_button("⬇️ " + lab, data=_df_to_csv(df), file_name=fn,
                                   mime="text/csv", key="dl_" + fn)

st.markdown("**Or save all files to a folder on this computer**")
default_out = str(Path(folder) / "PASTA_output") if source_mode == "📂 Local folder path" else ""
o1, o2 = st.columns([4, 1])
with o1:
    out_dir = st.text_input(
        "Output folder", value=default_out, key="out_dir", label_visibility="collapsed",
        placeholder=r"e.g. C:\Users\you\Documents\PASTA_output",
        help="Type or paste any folder, e.g. C:\\Users\\you\\Documents\\PASTA_output. "
             "Quotes from Explorer's *Copy as path* are fine. The folder is created if needed.")
with o2:
    save_clicked = st.button("💾 Save to folder", width="stretch")
if save_clicked:
    try:
        target = _clean_out_dir(out_dir)
        target.mkdir(parents=True, exist_ok=True)
        written = []
        for fn, (_, df) in outputs.items():
            df.to_csv(target / fn, index=False)
            written.append(fn)
        st.success("Saved {} files to `{}`".format(len(written), target))
        st.caption("  \n".join(written))
    except PermissionError:
        st.error("No permission to write to `{}`. Choose a folder you can write to, "
                 "e.g. inside C:\\Users\\<you>\\Documents.".format(out_dir))
    except Exception as exc:
        st.error("Could not save: " + str(exc))
