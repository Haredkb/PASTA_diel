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

Every step is a function in the pastatemp package (run_loggers, run_annual,
run_diel, coverage, plots, export); this file is only the user interface. See README.
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

import hashlib

import streamlit as st
import pandas as pd
import numpy as np
import plotly.io as pio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pastatemp as ps


def _reload_package_if_changed():
    """Streamlit (and Streamlit Cloud after a git push) reruns this script with new
    code but keeps imported modules in memory, so an old pastatemp could be used
    with a new app ("unexpected keyword argument"). Reload the package whenever its
    files have changed since it was imported."""
    pkg_dir = Path(ps.__file__).parent
    stamp = max(p.stat().st_mtime for p in pkg_dir.glob("*.py"))
    if getattr(ps, "_app_loaded_stamp", None) == stamp:
        return ps
    order = ["phase", "readers", "pairing", "annual", "diel", "coverage", "export", "core", "plots"]
    for name in order:
        mod = sys.modules.get("pastatemp." + name)
        if mod is not None:
            importlib.reload(mod)
    pkg = importlib.reload(sys.modules["pastatemp"])
    pkg._app_loaded_stamp = stamp
    return pkg


ps = _reload_package_if_changed()
from pastatemp import pairing as pr, plots as pp, readers as rd  # noqa: E402

pio.templates.default = "plotly_white"

# Pre-fill the folder box by setting DIEL_PASTA_DATA to your data folder.
_DEFAULT_FOLDER = os.environ.get("DIEL_PASTA_DATA", "")

# st.fragment (Streamlit >= 1.37) reruns only the decorated block when its own
# widgets change -- used for the plot selectors so the rest of the page is not
# recomputed. Falls back to a no-op decorator on older Streamlit.
_fragment = getattr(st, "fragment", None) or (lambda f: f)

ANALYSIS_LABELS = {"Annual signal": "annual", "Diel": "diel"}
FILE_LABELS = {"run_summary": "Run summary", "pairs": "Pairs used", "coverage": "Data coverage",
               "paired_df": "Paired data", "annual_metrics": "Annual metrics (record)",
               "annual_metrics_yearly": "Annual metrics (by year)", "annual_fits": "Annual fits",
               "daily_means": "Daily means", "diel_metrics": "Diel metrics", "diel_fits": "Diel fits"}


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
def _coverage(sig, _grids, pairs, interval_min):
    """ps.coverage(), cached; `_grids` is not hashed (leading underscore)."""
    return ps.coverage(_grids, list(pairs), interval_min)


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
    st.plotly_chart(pp.plot_loggers({n: grids[n] for n in airs + streams}, types), width="stretch")

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
    st.plotly_chart(pp.plot_pair(grids, row["stream"], row["air"]), width="stretch")


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
pkg_analyses = tuple(ANALYSIS_LABELS[a] for a in analyses)

pair_list = [(r.site_id, r.stream, r.air) for r in sel_pairs.itertuples()]
cov = _coverage(_sig(score_sig, tuple(pair_list)), grids, tuple(pair_list), int(interval_min))

st.markdown("#### Date range")
range_mode = st.radio(
    "Date range", ["All data", "Common period (all sites)", "Custom"], horizontal=True,
    key="range_mode", label_visibility="collapsed",
    help="**All data**: every site uses its full paired record. **Common period**: only "
         "dates where every selected site has paired air and stream data. **Custom**: pick dates.")
c_first, c_last = cov["paired_start"].min(), cov["paired_end"].max()
common0, common1 = ps.common_period(cov)
if range_mode == "All data":
    dr0, dr1 = ps.date_range(cov, "all")
elif range_mode == "Common period (all sites)":
    if common0 is None:
        st.error("No period is shared by all selected sites. Choose **All data**, "
                 "**Custom**, or untick sites in the pair table.")
        st.stop()
    dr0, dr1 = ps.date_range(cov, "common")
else:
    d0, d1 = c_first.date(), c_last.date()
    date_range = st.date_input("Custom dates", value=(d0, d1), min_value=d0, max_value=d1,
                               key="custom_dates")
    if isinstance(date_range, (list, tuple)) and len(date_range) == 2:
        dr0, dr1 = pd.Timestamp(date_range[0]), pd.Timestamp(date_range[1]) + pd.Timedelta(days=1)
    else:
        dr0, dr1 = pd.Timestamp(d0), pd.Timestamp(d1) + pd.Timedelta(days=1)

cov_view = ps.coverage_in_range(cov, dr0, dr1)
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
    st.plotly_chart(pp.plot_coverage(cov, dr0, dr1), width="stretch")
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
                   "{} to {}".format(common0.date(), common1.date()) if common0 is not None
                   else "none"))

with st.expander("Settings", expanded=True):
    s1, s2, s3 = st.columns(3)
    with s1:
        st.markdown("**Both**")
        water_min = st.number_input(
            "Drop stream records ≤ (°C)", value=0.99, step=0.01, format="%.2f",
            help="R: filter(str_temp > 0.99). Removes ice-affected / out-of-water records. "
                 "The annual analysis also drops days with a stream mean ≤ 1 °C, as PASTA does.")
        phase_method = st.radio(
            "Phase calculation", ["atan2", "atan"], horizontal=True, key="phase_method",
            format_func=lambda m: {"atan2": "atan2: fitted peak time (default)",
                                   "atan": "atan: match the R scripts"}[m],
            help="**atan2** returns the time of the fitted peak anywhere in the cycle and "
                 "wraps phase lags into ±half a cycle. **atan** reproduces the PASTA and "
                 "Diel PASTA R formulas, which only see half a cycle: annual phase 182.5–365 d "
                 "(a peak earlier in the year comes out half a year late) and diel phase 6–18 h "
                 "(the peak hour − 6 h; a stream peaking after midnight jumps by 12 h).")
        min_cov = st.slider(
            "Minimum day coverage (%)", 0, 100, 75, 5,
            help="Days with fewer paired records than this share of a full day are not "
                 "fitted (diel) and get no daily mean (annual). Set 0 to use every day.")
    with s2:
        st.markdown("**Annual signal**")
        yr_type = st.radio("Year type", ["water", "calendar"], horizontal=True,
                           disabled=not do_annual,
                           help="Water year runs 1 Oct – 30 Sep and is named by the year it ends.")
        min_year_pct = st.slider(
            "Single year: minimum % of the year", 0, 100, 75, 5,
            disabled=not do_annual,
            help="A year is fitted only if daily means with stream > 1 °C (after the 1 °C "
                 "drop) cover at least this share of the full water or calendar year "
                 "(365/366 days). Days before the record starts or after it ends count as "
                 "missing. Other years are dropped: no AmpRatio, phase lag or mean ratio.")
        min_rec_pct = st.slider(
            "Whole record: minimum % of record", 0, 100, 75, 5,
            disabled=not do_annual,
            help="The whole-record fit is reported only if daily means with stream > 1 °C "
                 "(after the 1 °C drop) cover at least this share of the record span, "
                 "first to last paired day.")
    with s3:
        st.markdown("**Diel**")
        r2_min = st.slider("Minimum adjusted R² (air and water)", 0.0, 1.0, 0.70, 0.05,
                           disabled=not do_diel,
                           help="Days kept only when both fits exceed this (R: AdjRsqr > 0.7).")
        roll_k = st.number_input("Rolling mean window (days)", value=7, min_value=1, step=2,
                                 disabled=not do_diel,
                                 help="Centred rolling mean over kept days (R: zoo::rollmean k = 7).")
        st.caption("Diel phase uses the logger clock time as written. The R `rad_hour()` "
                   "converts to the computer's time zone, which shifts phase by 1 h across DST. "
                   "With atan2, phase_h is the hour of the fitted daily peak.")

run_sig = _sig(score_sig, tuple(pair_list), tuple(analyses), water_min, min_cov, yr_type,
               phase_method,
               min_year_pct, min_rec_pct, r2_min, roll_k, str(dr0), str(dr1))


if st.button("Run PASTA", type="primary"):
    with st.spinner("Fitting sinusoids..."):
        try:
            # the same call a script would make: pastatemp.run_loggers()
            result = ps.run_loggers(
                grids, sel_pairs.reset_index(drop=True), analyses=pkg_analyses,
                date_range=(dr0, dr1 - pd.Timedelta(days=1)), interval_min=int(interval_min),
                water_min_C=water_min, min_coverage=min_cov / 100, phase_method=phase_method,
                yr_type=yr_type, min_year_pct=float(min_year_pct),
                min_record_pct=float(min_rec_pct),
                r2_min=r2_min, roll_k=int(roll_k))
            res = {"result": result, "T_y": result.paired, "summary": result.summary,
                   "analyses": list(analyses), "yr_type": yr_type,
                   "min_year_pct": int(min_year_pct), "min_rec_pct": int(min_rec_pct),
                   "sig": run_sig, "run_ts": result.run_ts}
            if result.diel is not None:
                res.update(fits=result.diel.fits, metrics=result.diel.metrics)
            if result.annual is not None:
                a = result.annual
                res.update(daily=a.daily, ann={"record": a.record, "yearly": a.yearly, "fits": a.fits})
            st.session_state["res"] = res
        except Exception as exc:
            st.error("Run failed: " + str(exc))
            import traceback
            with st.expander("Traceback"):
                st.code(traceback.format_exc())

if "res" not in st.session_state:
    st.stop()

res = st.session_state["res"]
if "result" not in res:            # results from an older app version in this session
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
    ares = res["result"].annual

    @_fragment
    def _annual_plots():
        sid = st.selectbox("Site", summary["site_id"].tolist(), key="ann_site")
        if daily[daily["site_id"] == sid].empty:
            st.info("No daily means for this site.")
            return
        c1, c2 = st.columns([2, 1])
        with c1:
            st.plotly_chart(pp.plot_annual_fit(ares, sid), width="stretch")
        with c2:
            st.plotly_chart(pp.plot_thermal_sensitivity(ares, sid), width="stretch")

    with tabs["Annual signal"]:
        rec_t = ann["record"]
        if len(rec_t):
            low = rec_t[~rec_t["meets_record_pct"].fillna(False).astype(bool)]
            if len(low):
                st.error(
                    "Days in the fit cover less than {}% of the record for: {}. These sites have "
                    "no whole-record AmpRatio, phase lag or mean ratio; lower **Whole record: "
                    "minimum % of record** to include them.".format(
                        res["min_rec_pct"], ", ".join("{} ({:.0f}%)".format(r.site_id, r.pct_record_fit)
                                                     for r in low.itertuples())))
            yl0 = ann["yearly"]
            if len(yl0):
                fail = yl0[~yl0["meets_year_pct"]]
                if len(fail):
                    st.info(
                        "Years with fit days on less than {}% of the year were dropped (no "
                        "signal metrics): {}.".format(
                            res["min_year_pct"], "; ".join(
                                "{} {} ({} of {} d, {:.0f}%)".format(
                                    r.site_id, r.year, r.n_days_fit, r.n_days_in_year,
                                    r.pct_year_fit)
                                for r in fail.itertuples())))
            short_sites = rec_t[(rec_t["record_span_days"] > 0) & (rec_t["record_span_days"] < 300)]
            if len(short_sites):
                st.warning(
                    "Less than ~300 days of record in the fit for: "
                    + ", ".join("{} ({} d)".format(r.site_id, r.record_span_days)
                                for r in short_sites.itertuples())
                    + ". An annual sinusoid fitted to part of a year is poorly constrained; "
                      "treat these annual metrics with caution.")
        st.markdown("#### Whole record")
        st.dataframe(ann["record"], hide_index=True, width="stretch")
        st.caption("Fitted to every day with stream > 1 °C (n_days_fit). Signal metrics are "
                   "reported only when those days cover at least {}% of the record span "
                   "(pct_record_fit = n_days_fit / record_span_days). AmpRatio and mean ratio are water / air; "
                   .format(res["min_rec_pct"]) +
                   "phase lag is water − air in days. TS__Slope is the daily thermal sensitivity "
                   "(stream ~ air, Kelleher et al. 2012).")
        st.markdown("#### By " + res["yr_type"] + " year")
        st.dataframe(ann["yearly"], hide_index=True, width="stretch")
        st.caption("pct_year_fit = n_days_fit / n_days_in_year: days with stream > 1 °C as a "
                   "share of the full year (365/366 days). A year below {}% is dropped "
                   "(meets_year_pct = False): no AmpRatio, phase lag or mean ratio; its thermal "
                   "sensitivity is still reported, as in PASTA. count is all paired days, "
                   "including cold days.".format(res["min_year_pct"]))
        _annual_plots()
        yl = ann["yearly"].dropna(subset=["AmpRatio", "PhaseLag_d"]) if len(ann["yearly"]) else ann["yearly"]
        if len(yl):
            st.plotly_chart(pp.plot_annual_by_year(ares), width="stretch")
        with st.expander("Annual sinusoid fits"):
            st.dataframe(_r4(ann["fits"]), hide_index=True, width="stretch")

# -- Diel results --
if ran_diel:
    fits, metrics = res["fits"], res["metrics"]
    dres = res["result"].diel

    @_fragment
    def _metric_plots():
        c1, c2 = st.columns([3, 1])
        with c2:
            x_mode = st.radio("x axis", ["Date", "Day of year"], key="xmode")
            ar_cap = st.checkbox("Hide AmpRatio ≥ 1.2", value=True, key="arcap",
                                 help="Values this high are usually poor data (Rmd plot filter).")
        fig = pp.plot_diel_metrics(dres, x="date" if x_mode == "Date" else "DOY",
                                   max_amp_ratio=1.2 if ar_cap else None)
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
        st.plotly_chart(pp.plot_diel_day(dres, sid, day), width="stretch")

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
result = res["result"]
stem = result.default_prefix()
tables = result.tables()            # same tables as result.to_csv() / to_zip() in a script
outputs = {"{}_{}.csv".format(stem, k): (FILE_LABELS.get(k, k), df) for k, df in tables.items()}


@st.cache_data(show_spinner=False)
def _zip_outputs(sig, _result):
    return _result.to_zip()


st.download_button(
    "⬇️ Download all ({} files, .zip)".format(len(outputs)),
    data=_zip_outputs(res["sig"] + stem + str(len(outputs)), result),
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
        written = ps.save_tables(tables, out_dir, stem)
        st.success("Saved {} files to `{}`".format(len(written), written[0].parent))
        st.caption("  \n".join(p.name for p in written))
    except PermissionError:
        st.error("No permission to write to `{}`. Choose a folder you can write to, "
                 "e.g. inside C:\\Users\\<you>\\Documents.".format(out_dir))
    except Exception as exc:
        st.error("Could not save: " + str(exc))
