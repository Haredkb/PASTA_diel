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
    """Output of run_pasta() / run_folder()."""
    annual: AnnualResult = None
    diel: DielResult = None
    summary: pd.DataFrame = None    # one row per site, both analyses side by side
    pairs: pd.DataFrame = None      # set by run_folder()

    def to_csv(self, folder, prefix="PASTA"):
        tables = {"run_summary": self.summary}
        if self.pairs is not None:
            tables["pairs"] = self.pairs
        for r in (self.annual, self.diel):
            if r is not None:
                tables.update(r.tables())
        return _write(tables, folder, prefix)


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
             r2_min=0.7, roll_k=7):
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

    Returns DielResult(paired, fits, metrics, summary, settings).
    """
    d = _standardize(data)
    if _is_daily(d):
        raise ValueError("run_diel() needs sub-daily records; these look like daily means.")
    interval_min = interval_min or _infer_interval(d)
    paired = _to_paired(d, interval_min, water_min_C)
    fits = _diel.fit_daily_sinusoids(paired, min_coverage=min_coverage, interval_min=interval_min)
    metrics = _diel.pasta_metrics(fits, r2_min=r2_min, roll_k=int(roll_k))

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
                "min_coverage": min_coverage, "r2_min": r2_min, "roll_k": roll_k}
    return DielResult(paired, fits, metrics, pd.DataFrame(rows), settings)


def run_annual(data, yr_type="water", min_days=100, water_min_C=0.99,
               min_coverage=0.75, interval_min=None):
    """
    Annual-signal PASTA: fit temp ~ sin(2πd/365) + cos(2πd/365) to daily mean air
    and stream temperature for each site, for the whole record and for each year,
    plus daily thermal sensitivity (stream ~ air) and a data gap check.

    data:          DataFrame, first four columns site_id, date/time, stream °C, air °C.
                   Sub-daily records are averaged to daily means first; daily means are
                   used as they are.
    yr_type:       'water' (1 Oct – 30 Sep, named by the year it ends) or 'calendar'.
    min_days:      years with fewer daily means get no AmpRatio / PhaseLag_d / Ratio_Mean.
    water_min_C:   sub-daily stream records at or below this are dropped before averaging.
                   Days with a stream mean ≤ 1 °C are always dropped (as in PASTA).
    min_coverage:  share of a full day of paired records needed for a daily mean.
    interval_min:  logging interval for sub-daily data; inferred when None.

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
        daily = _annual.daily_means(paired, min_coverage=min_coverage, interval_min=interval_min)
    res = _annual.annual_metrics(daily, yr_type=yr_type, min_days=int(min_days))
    span = daily.groupby("site_id")["date"].agg(lambda x: (x.max() - x.min()).days + 1)
    rec = res["record"]
    if len(rec):
        rec = rec.merge(span.rename("record_span_days").reset_index(), on="site_id", how="left")
    settings = {"yr_type": yr_type, "min_days": min_days, "water_min_C": water_min_C,
                "min_coverage": min_coverage, "interval_min": interval_min}
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
        names = fn.__code__.co_varnames[:fn.__code__.co_argcount]
        return {k: v for k, v in kwargs.items() if k in names}

    out = PastaResult()
    parts = []
    if "annual" in analyses:
        out.annual = run_annual(data, **_kw(run_annual))
        rec = out.annual.record
        if len(rec):
            cols = ["AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope", "TS_AdjRsqr",
                    "max_conseq_missing_days", "count", "record_span_days"]
            parts.append(rec[["site_id"] + cols].rename(
                columns={c: "annual_" + c for c in cols}))
    if "diel" in analyses:
        out.diel = run_diel(data, **_kw(run_diel))
        parts.append(out.diel.summary)
    summary = None
    for p in parts:
        summary = p if summary is None else summary.merge(p, on="site_id", how="outer")
    out.summary = summary if summary is not None else pd.DataFrame()
    return out


def run_folder(folder, analyses=ANALYSES, pairs=None, interval_min=15, out_dir=None,
               **kwargs):
    """
    Read every logger CSV in a folder, sort air vs stream, pair them, and run PASTA.

    pairs:   optional list of (site_id, stream_logger, air_logger) to override the
             suggested pairs (e.g. from field notes).
    out_dir: if given, all tables are written there as CSVs.

    Returns PastaResult with .pairs set to the pairing table that was used.
    """
    loggers = read_loggers(folder, interval_min=interval_min)
    if pairs is None:
        pair_table = pair_loggers(loggers, interval_min=interval_min)
    else:
        pair_table = pd.DataFrame(pairs, columns=["site_id", "stream", "air"])
        pair_table["basis"] = "user"
    data = pair_data(loggers, pair_table)
    res = run_pasta(data, analyses=analyses, interval_min=interval_min, **kwargs)
    res.pairs = pair_table
    if out_dir:
        res.to_csv(out_dir)
    return res
