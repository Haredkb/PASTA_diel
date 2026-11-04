"""
Reading logger files: HOBOware exports and simple date/time + temperature CSVs.
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
    if t_col is None:
        # e.g. 'Timestamp (UTC-05:00),LLC_LLC01_AT': the one numeric column
        # that is not the timestamp or a HOBO '#' row counter
        num = [c for c in cols if c != ts_col and c.strip() != "#"
               and pd.to_numeric(df[c], errors="coerce").notna().mean() > 0.5]
        if len(num) == 1:
            t_col = num[0]
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

    tz = re.search(r"\(([A-Z]{2,5}|(?:GMT|UTC)[^)]*)\)", ts_col)
    out.attrs["tz_label"] = tz.group(1) if tz else ""
    return out


def to_grid(df, interval_min=15):
    """Round timestamps to the logging interval and average duplicates (R: round_date)."""
    g = df.copy()
    g["timestamp"] = g["timestamp"].dt.round(str(int(interval_min)) + "min")
    return g.groupby("timestamp", as_index=True)["temp_C"].mean()


def read_logger(path, interval_min=15):
    """Read one logger CSV from disk -> gridded temperature Series (index = timestamp)."""
    path = Path(path)
    return to_grid(parse_logger_csv(path.read_bytes(), path.name), interval_min)


def read_loggers(source, interval_min=15):
    """
    Read several logger CSVs.

    source: a folder (every *.csv in it is read) or a list of file paths.
    Returns {logger_name: gridded Series}. Logger names come from the file
    names ('HF_005 2026-07-22 12_44_20 EDT.csv' -> 'HF_005').
    """
    if isinstance(source, (str, Path)) and Path(source).is_dir():
        paths = sorted(Path(source).glob("*.csv"))
    else:
        paths = [Path(p) for p in ([source] if isinstance(source, (str, Path)) else source)]
    out = {}
    for p in paths:
        name = logger_name_from_filename(p.name)
        if name in out:
            name = p.stem
        out[name] = read_logger(p, interval_min)
    return out
