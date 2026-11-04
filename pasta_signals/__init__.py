"""
pasta_signals -- Paired Air and Stream Temperature Analysis (PASTA)
at the annual and diel scales.

Main functions
    run_annual(data, yr_type="water", ...)   annual sinusoid of daily means
    run_diel(data, ...)                      daily sinusoid of sub-daily records
    run_pasta(data, analyses=("annual", "diel"), ...)   both, with a combined summary
    run_folder(folder, ...)                  read logger CSVs, pair them, run both

`data` is a DataFrame whose first four columns are
site_id, date/time, stream temperature (°C), air temperature (°C).

Helpers for logger files
    read_loggers(folder_or_paths)     -> {name: Series}
    classify_loggers(loggers)         -> air / stream table
    pair_loggers(loggers, types)      -> suggested air-stream pairs
    pair_data(loggers, pairs)         -> four-column table for run_*()

Example
    >>> import pasta_signals as ps
    >>> res = ps.run_folder("path/to/logger_csvs", out_dir="path/to/output")
    >>> res.summary
"""

from .core import (ANALYSES, AnnualResult, DielResult, PastaResult, pair_data,
                   run_annual, run_diel, run_folder, run_pasta)
from .pairing import classify_loggers, pair_loggers
from .readers import parse_logger_csv, read_logger, read_loggers

__version__ = "0.1.0"

__all__ = [
    "run_annual", "run_diel", "run_pasta", "run_folder",
    "read_loggers", "read_logger", "parse_logger_csv",
    "classify_loggers", "pair_loggers", "pair_data",
    "AnnualResult", "DielResult", "PastaResult", "ANALYSES",
]
