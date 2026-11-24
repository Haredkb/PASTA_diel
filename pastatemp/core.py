"""
Main entry points: run_annual(), run_diel() and run_pasta().

All three take a DataFrame whose first four columns are, in this order
(the same convention as the PASTA R code):

    site_id, date/time, stream (water) temperature °C, air temperature °C

Column names do not matter. Sub-daily records (e.g. 15-min logger data) work
for both analyses; daily means work for run_annual() only.

Logger files can be turned into that table with read_loggers() ->
pair_loggers() -> pair_data(), or all at once with run_folder().
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import annual as _annual
from . import diel as _diel
# import the functions, not the module: the package exposes a function named
# `coverage`, which hides the submodule of the same name on reload
from .coverage import clip_loggers as _clip_loggers
from .coverage import coverage as _coverage
from .coverage import coverage_in_range as _coverage_in_range
from .coverage import date_range as _date_range
from .coverage import record_spans as _record_spans
from .phase import DEFAULT_PHASE_METHOD, check_method
from .pairing import classify_loggers, pair_loggers
from .readers import read_loggers

ANALYSES = ("annual", "diel")


# ==============================================================================
#  Results
# ==============================================================================

def _write(tables, folder, prefix):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, df in tables.items():
        if df is None:
            continue
        p = folder / "{}_{}.csv".format(prefix, name)
        df.to_csv(p, index=False)
        paths.append(p)
    return paths


@dataclass
class DielResult:
    """Output of run_diel()."""
    paired: pd.DataFrame        # site_id, DateTime_r, str_temp, air_temp (after filters)
    fits: pd.DataFrame          # one row per site / day / medium
    metrics: pd.DataFrame       # one row per site / day where both fits pass r2_min
    summary: pd.DataFrame       # one row per site
    settings: dict = field(default_factory=dict)

    def tables(self):
        return {"diel_metrics": self.metrics, "diel_fits": self.fits,
                "diel_summary": self.summary}

    def to_csv(self, folder, prefix="PASTA"):
        """Write diel_metrics, diel_fits and diel_summary CSVs; returns the paths."""
        return _write(self.tables(), folder, prefix)


@dataclass
class AnnualResult:
    """Output of run_annual()."""
    daily: pd.DataFrame         # site_id, date, tavg_wat_C, tavg_air_C
    record: pd.DataFrame        # whole-record metrics, one row per site
    yearly: pd.DataFrame        # metrics per site and water/calendar year
    fits: pd.DataFrame          # annual sinusoid coefficients (year 'all' = whole record)
    settings: dict = field(default_factory=dict)

    def tables(self):
        return {"annual_metrics": self.record, "annual_metrics_yearly": self.yearly,
                "annual_fits": self.fits, "daily_means": self.daily}

    def to_csv(self, folder, prefix="PASTA"):
        """Write annual_metrics, annual_metrics_yearly, annual_fits and daily_means CSVs."""
        return _write(self.tables(), folder, prefix)


@dataclass
class PastaResult:
    """Output of run_pasta() / run_loggers() / run_folder()."""
    annual: AnnualResult = None
    diel: DielResult = None
    summary: pd.DataFrame = None    # one row per site: settings, annual and diel side by side
    pairs: pd.DataFrame = None      # pairing table used (run_loggers / run_folder)
    coverage: pd.DataFrame = None   # data coverage per site in the date range (run_loggers)
    paired: pd.DataFrame = None     # site_id, DateTime_r, str_temp, air_temp (stream > water_min_C)
    date_range: tuple = None        # (start, end) analysed, end exclusive (run_loggers)
    run_ts: pd.Timestamp = None
    settings: dict = field(default_factory=dict)

    def tables(self):
        """The same tables the app downloads, in the same order: run_summary, pairs,
        coverage, paired_df, annual_metrics, annual_metrics_yearly, annual_fits,
        daily_means, diel_metrics, diel_fits (missing ones are skipped)."""
        t = {"run_summary": self.summary, "pairs": self.pairs, "coverage": self.coverage,
             "paired_df": self.paired}
        if self.annual is not None:
            t.update(self.annual.tables())
        if self.diel is not None:
            t.update({"diel_metrics": self.diel.metrics, "diel_fits": self.diel.fits})
        return {k: v for k, v in t.items() if v is not None}

    def default_prefix(self):
        """'PASTA_YYYYmmdd_HHMM', the app's file-name stem."""
        ts = self.run_ts if self.run_ts is not None else pd.Timestamp.now()
        return "PASTA_" + ts.strftime("%Y%m%d_%H%M")

    def to_csv(self, folder, prefix="PASTA"):
        """Write every table as <prefix>_<table>.csv in `folder`; returns the paths."""
        return _write(self.tables(), folder, prefix)

    def to_zip(self, path=None, prefix=None):
        """All tables in one zip (the app's *Download all*). Returns the bytes, and
        also writes them to `path` when given."""
        from .export import tables_to_zip
        return tables_to_zip(self.tables(), prefix or self.default_prefix(), path)


# ==============================================================================
#  Input handling
# ==============================================================================

def _standardize(data):
    """First four columns -> site_id, datetime, water, air (R convention)."""
    if not isinstance(data, pd.DataFrame) or data.shape[1] < 4:
        raise ValueError("data must be a DataFrame whose first four columns are "
                         "site_id, date/time, stream temperature, air temperature")
    d = data.iloc[:, :4].copy()
    d.columns = ["site_id", "datetime", "water", "air"]
    d["site_id"] = d["site_id"].astype(str)
    d["datetime"] = pd.to_datetime(d["datetime"])
    d["water"] = pd.to_numeric(d["water"], errors="coerce")
    d["air"] = pd.to_numeric(d["air"], errors="coerce")
    return d.dropna().sort_values(["site_id", "datetime"]).reset_index(drop=True)


def _is_daily(d):
    t = d["datetime"]
    midnight = (t == t.dt.normalize()).all()
    one_per_day = not d.duplicated(["site_id", "datetime"]).any()
    return bool(midnight and one_per_day)


def _infer_interval(d):
    """Median logging interval in minutes (across sites)."""
    steps = d.groupby("site_id")["datetime"].diff().dt.total_seconds().div(60).dropna()
    steps = steps[steps > 0]
    return int(round(steps.median())) if len(steps) else 15


def _to_paired(d, interval_min, water_min_C):
    """Round to the logging interval and drop cold stream records. Every record is
    kept (as in R's round_date + lm), so finer data rounded to a coarser interval
    still contributes all its points to the fits."""
    p = pd.DataFrame({"site_id": d["site_id"],
                      "DateTime_r": d["datetime"].dt.round(str(int(interval_min)) + "min"),
                      "str_temp": d["water"], "air_temp": d["air"]})
    return p[p["str_temp"] > water_min_C].reset_index(drop=True)


def pair_data(loggers, pairs):
    """
    Build the four-column input table from logger series.

    loggers: {name: gridded Series} from read_loggers().
    pairs:   DataFrame from pair_loggers() (columns site_id, stream, air), or a
             list of (site_id, stream_logger, air_logger) tuples.
    Returns DataFrame site_id, DateTime_r, str_temp, air_temp (inner join on time).
    """
    if isinstance(pairs, pd.DataFrame):
        pairs = [(r.site_id, r.stream, r.air) for r in pairs.itertuples() if pd.notna(r.air)]
    return _diel.build_paired_df(pairs, loggers, water_min_C=-np.inf)


# ==============================================================================
#  Main functions
# ==============================================================================

def run_diel(data, interval_min=None, water_min_C=0.99, min_coverage=0.75,
             r2_min=0.7, roll_k=7, phase_method=DEFAULT_PHASE_METHOD):
    """
    Diel PASTA: fit temp ~ sin(2πh/24) + cos(2πh/24) to air and stream for every
    site and day, then compute water/air metrics for days where both fits pass r2_min.

    data:          DataFrame, first four columns site_id, date/time, stream °C, air °C
                   (sub-daily records).
    interval_min:  logging interval for rounding timestamps; inferred when None.
    water_min_C:   stream records at or below this are dropped (R: str_temp > 0.99).
    min_coverage:  share of a full day of records needed to fit a day (0 = every day, as in R).
    r2_min:        both fits must exceed this adjusted R² for the day to be kept.
    roll_k:        centred rolling-mean window (days) for AR_roll7 / PH_roll7.
    phase_method:  'atan2' (default): phase_h = hour of the fitted daily peak, lag
                   wrapped into (-12, 12] h. 'atan': the Diel PASTA R formula, to
                   reproduce R results (phase confined to 6-18 h; see pastatemp.phase).

    Returns DielResult(paired, fits, metrics, summary, settings).
    """
    d = _standardize(data)
    if _is_daily(d):
        raise ValueError("run_diel() needs sub-daily records; these look like daily means.")
    interval_min = interval_min or _infer_interval(d)
    paired = _to_paired(d, interval_min, water_min_C)
    check_method(phase_method)
    fits = _diel.fit_daily_sinusoids(paired, min_coverage=min_coverage, interval_min=interval_min,
                                     phase_method=phase_method)
    metrics = _diel.pasta_metrics(fits, r2_min=r2_min, roll_k=int(roll_k),
                                  phase_method=phase_method)

    rows = []
    for sid in paired["site_id"].unique():
        m = metrics[metrics["site_id"] == sid] if len(metrics) else metrics
        f = fits[fits["site_id"] == sid] if len(fits) else fits
        rows.append({"site_id": sid,
                     "diel_n_days_fitted": int(f["date"].nunique()) if len(f) else 0,
                     "diel_n_days_kept": int(len(m)),
                     "diel_median_AmpRatio": m["AmpRatio"].median() if len(m) else np.nan,
                     "diel_median_PhaseLag_h": m["PhaseLag_h"].median() if len(m) else np.nan,
                     "diel_median_Ratio_Mean_K": m["Ratio_Mean_K"].median() if len(m) else np.nan})
    settings = {"interval_min": interval_min, "water_min_C": water_min_C,
                "min_coverage": min_coverage, "r2_min": r2_min, "roll_k": roll_k,
                "phase_method": phase_method}
    return DielResult(paired, fits, metrics, pd.DataFrame(rows), settings)


def run_annual(data, yr_type="water", min_year_pct=25, max_gap_days=50, min_record_pct=75,
               water_min_C=0.99, interval_min=None, phase_method=DEFAULT_PHASE_METHOD,
               record_span_days=None):
    """
    Annual-signal PASTA: fit temp ~ sin(2πd/365) + cos(2πd/365) to daily mean air
    and stream temperature for each site, for the whole record and for each year,
    plus daily thermal sensitivity (stream ~ air) and a data gap check.

    data:          DataFrame, first four columns site_id, date/time, stream °C, air °C.
                   Sub-daily records are averaged to daily means first; daily means are
                   used as they are.
    yr_type:       'water' (1 Oct – 30 Sep, named by the year it ends) or 'calendar'.
    min_year_pct:  a year needs fit days (daily means with stream > 1 °C, i.e. after
                   the 1 °C drop) on at least this % of the full water (or calendar)
                   year, 365/366 days (default 25). Other years get no AmpRatio /
                   PhaseLag_d / Ratio_Mean in `yearly`.
    max_gap_days:  a year is also dropped when its longest run of days without a fit
                   day (between fit days) is longer than this (default 50 d, following
                   Johnson et al. 2020; see pastatemp.annual.JOHNSON_2020).
    min_record_pct: the whole record needs fit days on at least this % of the record
                   span (default 75). Otherwise `record` gets no AmpRatio / PhaseLag_d /
                   Ratio_Mean.
    record_span_days: {site_id: days} record span for that rule; run_loggers() passes
                   the longer of the stream and air records ('all') or the chosen date
                   range. When None: first to last paired day.
    water_min_C:   sub-daily stream records at or below this are dropped before averaging.
                   Days with a stream mean ≤ 1 °C are always dropped (as in PASTA).
                   Every day with records gets a daily mean (no day-coverage rule;
                   min_coverage applies to the diel analysis only).
    interval_min:  logging interval for sub-daily data; inferred when None.
    phase_method:  'atan2' (default): phase_d = day of the fitted peak, lag wrapped
                   into (-182.5, 182.5] d. 'atan': the PASTA R formula, to reproduce
                   R results (phase confined to 182.5-365 d; see pastatemp.phase).

    Returns AnnualResult(daily, record, yearly, fits, settings).
    """
    if yr_type not in ("water", "calendar"):
        raise ValueError("yr_type must be 'water' or 'calendar'")
    d = _standardize(data)
    if _is_daily(d):
        daily = d.rename(columns={"datetime": "date", "water": "tavg_wat_C",
                                  "air": "tavg_air_C"})
        interval_min = None
    else:
        interval_min = interval_min or _infer_interval(d)
        paired = _to_paired(d, interval_min, water_min_C)
        daily = _annual.daily_means(paired, min_coverage=0, interval_min=interval_min)
    check_method(phase_method)
    res = _annual.annual_metrics(daily, yr_type=yr_type, min_year_pct=float(min_year_pct),
                                 max_gap_days=float(max_gap_days),
                                 min_record_pct=float(min_record_pct),
                                 phase_method=phase_method, record_span_days=record_span_days)
    rec = res["record"]
    settings = {"yr_type": yr_type, "min_year_pct": min_year_pct, "max_gap_days": max_gap_days,
                "min_record_pct": min_record_pct, "water_min_C": water_min_C,
                "interval_min": interval_min, "phase_method": phase_method}
    return AnnualResult(daily, rec, res["yearly"], res["fits"], settings)


def run_pasta(data, analyses=ANALYSES, **kwargs):
    """
    Run run_annual() and/or run_diel() on the same data.

    analyses: ('annual', 'diel') by default; pass ('annual',) or ('diel',) for one.
    kwargs:   passed to whichever function accepts them (e.g. yr_type, r2_min).

    Returns PastaResult(annual, diel, summary) with one summary row per site.
    """
    analyses = tuple(a.lower() for a in analyses)
    bad = set(analyses) - set(ANALYSES)
    if bad or not analyses:
        raise ValueError("analyses must be one or both of " + str(ANALYSES))

    def _kw(fn):
        return {k: v for k, v in kwargs.items() if k in _argnames(fn)}

    unknown = set(kwargs) - set(_argnames(run_annual)) - set(_argnames(run_diel))
    if unknown:
        raise TypeError("Unknown settings: " + ", ".join(sorted(unknown)))

    out = PastaResult(run_ts=pd.Timestamp.now())
    if "annual" in analyses:
        out.annual = run_annual(data, **_kw(run_annual))
    if "diel" in analyses:
        out.diel = run_diel(data, **_kw(run_diel))

    d = _standardize(data)
    water_min = kwargs.get("water_min_C", 0.99)
    interval = kwargs.get("interval_min")
    for r in (out.annual, out.diel):
        if r is not None and r.settings.get("interval_min"):
            interval = r.settings["interval_min"]
    keep = d[d["water"] > water_min]
    out.paired = pd.DataFrame({"site_id": keep["site_id"], "DateTime_r": keep["datetime"],
                               "str_temp": keep["water"], "air_temp": keep["air"]}).reset_index(drop=True)
    out.settings = {"analyses": analyses, "interval_min": interval, "water_min_C": water_min,
                    **{k: v for r in (out.annual, out.diel) if r is not None
                       for k, v in r.settings.items() if k != "interval_min"}}
    out.summary = run_summary(out, d)
    return out


def _argnames(fn):
    return fn.__code__.co_varnames[:fn.__code__.co_argcount]


def run_summary(res, data=None, pairs=None, date_range=None):
    """
    One row per site, as in the app's run summary: run time, loggers, settings,
    then the annual whole-record metrics and the diel medians.
    """
    s = res.settings
    d = _standardize(data) if data is not None and "water" not in getattr(data, "columns", []) else data
    sites = list(pd.unique(d["site_id"])) if d is not None else []
    pmap = {}
    if pairs is not None:
        pmap = {r.site_id: (r.stream, r.air) for r in pairs.itertuples()}
        sites += [p for p in pmap if p not in sites]
    if date_range is not None:
        d0, d1 = date_range[0], date_range[1] - pd.Timedelta(days=1)
    elif d is not None and len(d):
        d0, d1 = d["datetime"].min(), d["datetime"].max()
    else:
        d0 = d1 = pd.NaT
    n_rec = res.paired["site_id"].value_counts() if res.paired is not None else pd.Series(dtype=int)
    rows = []
    for sid in sites:
        row = {"run_timestamp": res.run_ts.isoformat(timespec="seconds") if res.run_ts is not None else "",
               "site_id": sid,
               "stream_logger": pmap.get(sid, (None, None))[0],
               "air_logger": pmap.get(sid, (None, None))[1],
               "analyses": " + ".join(s.get("analyses", ())),
               "interval_min": s.get("interval_min"), "water_min_C": s.get("water_min_C"),
               "min_day_coverage_pct": 100 * s["min_coverage"] if "min_coverage" in s else np.nan,
               "phase_method": s.get("phase_method"),
               "date_from": d0.date().isoformat() if pd.notna(d0) else "",
               "date_to": d1.date().isoformat() if pd.notna(d1) else "",
               "n_paired_records": int(n_rec.get(sid, 0))}
        if res.annual is not None:
            rec = res.annual.record
            r = rec[rec["site_id"] == sid] if len(rec) else rec
            row.update({"year_type": s.get("yr_type"),
                        "min_year_pct": s.get("min_year_pct"),
                        "max_gap_days": s.get("max_gap_days"),
                        "min_record_pct": s.get("min_record_pct"),
                        "annual_n_days": int(r["n_days_fit"].iloc[0]) if len(r) else 0})
            for c in ["AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope", "TS_AdjRsqr",
                      "max_conseq_missing_days", "record_span_days", "pct_record_fit",
                      "meets_record_pct"]:
                row["annual_" + c] = r[c].iloc[0] if len(r) else np.nan
        if res.diel is not None:
            ds = res.diel.summary
            r = ds[ds["site_id"] == sid] if len(ds) else ds
            row.update({"adjR2_min": s.get("r2_min"), "roll_k": s.get("roll_k")})
            for c in ["diel_n_days_fitted", "diel_n_days_kept", "diel_median_AmpRatio",
                      "diel_median_PhaseLag_h", "diel_median_Ratio_Mean_K"]:
                row[c] = r[c].iloc[0] if len(r) else (0 if "n_days" in c else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def _pair_table(pairs):
    """List of (site_id, stream, air) or a pairing DataFrame -> DataFrame."""
    if isinstance(pairs, pd.DataFrame):
        t = pairs.copy()
        if "include" in t.columns:
            t = t[t["include"].astype(bool)].drop(columns=["include"])
        return t[t["air"].notna()].reset_index(drop=True)
    t = pd.DataFrame([tuple(p) for p in pairs], columns=["site_id", "stream", "air"])
    t["basis"] = "user"
    return t


def run_loggers(loggers, pairs=None, analyses=ANALYSES, date_range="all", interval_min=15,
                **kwargs):
    """
    The app's whole run as one call: logger series + pairs -> data coverage, date
    range, pairing, and the annual and/or diel analysis.

    loggers:    {name: gridded Series} from read_loggers().
    pairs:      pairing table (site_id, stream, air[, include]) or list of
                (site_id, stream, air); suggested by pair_loggers() when None.
    date_range: 'all' (each site's full record), 'common' (only dates every site
                has paired data) or (start, end) dates, inclusive.
    kwargs:     settings for run_annual() / run_diel() (water_min_C, min_coverage (diel),
                phase_method, yr_type, min_year_pct, max_gap_days, min_record_pct, r2_min,
                roll_k).

    The annual whole-record span (record_span_days) is the longer of each site's
    stream and air records for 'all', else the chosen date range.

    Returns PastaResult with pairs, coverage (per site in the range) and date_range set.
    """
    pair_table = pair_loggers(loggers, interval_min=interval_min) if pairs is None else _pair_table(pairs)
    pair_table = pair_table[pair_table["air"].notna()].reset_index(drop=True)
    if pair_table.empty:
        raise ValueError("No air-stream pairs to run.")
    if pair_table["site_id"].duplicated().any():
        raise ValueError("Site IDs must be unique: "
                         + ", ".join(pair_table["site_id"][pair_table["site_id"].duplicated()]))
    cov = _coverage(loggers, pair_table, interval_min)
    start, end = _date_range(cov, date_range)
    spans = _record_spans(cov, start, end, date_range)
    data = pair_data(_clip_loggers(loggers, start, end), pair_table)
    if "annual" in tuple(a.lower() for a in analyses):
        kwargs.setdefault("record_span_days", spans)
    res = run_pasta(data, analyses=analyses, interval_min=interval_min, **kwargs)
    res.pairs = pair_table
    res.coverage = _coverage_in_range(cov, start, end, spans)
    res.date_range = (start, end)
    res.summary = run_summary(res, data, pairs=pair_table, date_range=(start, end))
    return res


def run_folder(folder, analyses=ANALYSES, pairs=None, interval_min=15, out_dir=None,
               date_range="all", **kwargs):
    """
    Read every logger CSV in a folder, sort air vs stream, pair them, and run PASTA
    (read_loggers() + run_loggers()).

    pairs:      optional list of (site_id, stream_logger, air_logger) to override the
                suggested pairs (e.g. from field notes).
    date_range: 'all', 'common' or (start, end); see run_loggers().
    out_dir:    if given, all tables are written there as CSVs (PASTA_<table>.csv).

    Returns PastaResult with .pairs set to the pairing table that was used.
    """
    loggers = read_loggers(folder, interval_min=interval_min)
    res = run_loggers(loggers, pairs=pairs, analyses=analyses, date_range=date_range,
                      interval_min=interval_min, **kwargs)
    if out_dir:
        res.to_csv(out_dir)
    return res
