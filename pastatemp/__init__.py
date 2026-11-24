"""
pastatemp -- Paired Air and Stream Temperature Analysis (PASTA)
at the annual and diel scales. Every step of the PASTA app is a function here.

Main functions
    run_folder(folder, ...)                   read logger CSVs, pair them, run both
    run_loggers(loggers, pairs, ...)          the app's run: coverage, date range, analyses
    run_pasta(data, analyses=("annual", "diel"), ...)   both, with a combined summary
    run_annual(data, yr_type="water", ...)    annual sinusoid of daily means
    run_diel(data, ...)                       daily sinusoid of sub-daily records

`data` is a DataFrame whose first four columns are
site_id, date/time, stream temperature (°C), air temperature (°C).

Logger files                          (app steps 1-3)
    read_loggers(folder_or_paths)     -> {name: Series}
    classify_loggers(loggers)         -> air / stream table
    pair_loggers(loggers, types)      -> suggested air-stream pairs (by name)
    pair_data(loggers, pairs)         -> four-column table for run_*()

Date range and coverage               (app step 4)
    coverage(loggers, pairs), coverage_in_range(cov, start, end),
    common_period(cov), date_range(cov, "all" | "common" | (start, end)), clip_loggers()

Saving                                (app step 6)
    result.to_csv(folder), result.to_zip(path), save_tables(), tables_to_zip()

Figures (needs plotly)                (app step 5)
    pastatemp.plots: plot_loggers, plot_pair, plot_coverage, plot_annual_fit,
    plot_thermal_sensitivity, plot_annual_by_year, plot_diel_metrics, plot_diel_day

Example
    >>> import pastatemp as pt
    >>> res = pt.run_folder("path/to/logger_csvs", date_range="common")
    >>> res.summary
    >>> res.to_zip("PASTA_results.zip")
"""

from .core import (ANALYSES, AnnualResult, DielResult, PastaResult, pair_data,
                   run_annual, run_diel, run_folder, run_loggers, run_pasta, run_summary)
from .coverage import (all_sites_cover, clip_loggers, common_period, coverage,
                       coverage_in_range, date_range, record_spans)
from .export import clean_out_dir, save_tables, tables_to_zip
from .pairing import classify_loggers, pair_loggers
from .phase import DEFAULT_PHASE_METHOD, PHASE_METHODS
from .readers import parse_logger_csv, read_logger, read_loggers

__version__ = "0.2.0"

__all__ = [
    "run_folder", "run_loggers", "run_pasta", "run_annual", "run_diel", "run_summary",
    "read_loggers", "read_logger", "parse_logger_csv",
    "classify_loggers", "pair_loggers", "pair_data",
    "coverage", "coverage_in_range", "common_period", "date_range", "record_spans", "clip_loggers",
    "all_sites_cover",
    "save_tables", "tables_to_zip", "clean_out_dir",
    "AnnualResult", "DielResult", "PastaResult", "ANALYSES",
    "PHASE_METHODS", "DEFAULT_PHASE_METHOD",
]
