"""
Data coverage and date range: when each logger and each air-stream pair has data,
the period shared by all sites, and clipping loggers to a range.
"""

import numpy as np
import pandas as pd

DATE_RANGE_MODES = ("all", "common")


def _pair_tuples(pairs):
    """DataFrame (site_id, stream, air) or list of tuples -> list of tuples."""
    if isinstance(pairs, pd.DataFrame):
        return [(r.site_id, r.stream, r.air) for r in pairs.itertuples()]
    return [tuple(p) for p in pairs]


def coverage(loggers, pairs, interval_min=15):
    """
    Start / end of each logger and of the paired record (both present), per site.

    loggers: {name: gridded Series} from read_loggers().
    pairs:   pairing table (site_id, stream, air) or list of (site_id, stream, air).
    Returns one row per site with stream_/air_/paired_ start and end, paired_days
    and paired_index (the paired timestamps; dropped by coverage_in_range()).
    """
    per_day = 24 * 60 / interval_min
    rows = []
    for sid, s, a in _pair_tuples(pairs):
        gs, ga = loggers[s].dropna(), loggers[a].dropna()
        both = gs.index.intersection(ga.index)
        rows.append({"site_id": sid, "stream": s, "air": a,
                     "stream_start": gs.index.min(), "stream_end": gs.index.max(),
                     "air_start": ga.index.min(), "air_end": ga.index.max(),
                     "paired_start": both.min() if len(both) else pd.NaT,
                     "paired_end": both.max() if len(both) else pd.NaT,
                     "paired_days": len(both) / per_day,
                     "paired_index": both})
    return pd.DataFrame(rows)


def common_period(cov):
    """(start, end) of the paired period shared by every site, or (None, None)."""
    c0, c1 = cov["paired_start"].max(), cov["paired_end"].min()
    if pd.isna(c0) or pd.isna(c1) or c0 >= c1:
        return None, None
    return c0, c1


def date_range(cov, mode="all"):
    """
    Analysis window [start, end) as whole days.

    mode: 'all'    - first to last paired day of any site (each site uses its full record)
          'common' - only the days every site has paired data (raises if there are none)
          (start, end) - your own dates, inclusive; strings or Timestamps
    """
    if isinstance(mode, (tuple, list)) and len(mode) == 2:
        return pd.Timestamp(mode[0]).normalize(), pd.Timestamp(mode[1]).normalize() + pd.Timedelta(days=1)
    if mode == "all":
        c0, c1 = cov["paired_start"].min(), cov["paired_end"].max()
    elif mode == "common":
        c0, c1 = common_period(cov)
        if c0 is None:
            raise ValueError("No period is shared by all sites; use date_range='all' or your own dates.")
    else:
        raise ValueError("date_range must be 'all', 'common' or (start, end)")
    return c0.normalize(), c1.normalize() + pd.Timedelta(days=1)


def coverage_in_range(cov, start, end):
    """
    Per site: paired days inside [start, end), % of the range, and a status
    ('covers whole range', 'starts …', 'ends …', 'gaps inside range', ...).
    """
    out = cov.drop(columns=["paired_index"], errors="ignore").copy()
    n_days = max((end - start).days, 1)
    tol = pd.Timedelta(days=1)
    in_range, status = [], []
    for r in cov.itertuples():
        idx = r.paired_index
        k = idx[(idx >= start) & (idx < end)]
        days = len(pd.unique(k.normalize())) if len(k) else 0
        in_range.append(days)
        if days == 0:
            status.append("no data in range")
        elif r.paired_start > start + tol and r.paired_end < end - 2 * tol:
            status.append("starts late and ends early")
        elif r.paired_start > start + tol:
            status.append("starts {}".format(r.paired_start.date()))
        elif r.paired_end < end - 2 * tol:
            status.append("ends {}".format(r.paired_end.date()))
        elif days < 0.9 * n_days:
            status.append("gaps inside range")
        else:
            status.append("covers whole range")
    out["days_in_range"] = in_range
    out["pct_of_range"] = [100 * d / n_days for d in in_range]
    out["status"] = status
    return out


def all_sites_cover(cov_range):
    """True when every site covers the whole range (coverage_in_range() output)."""
    return bool(len(cov_range)) and bool((cov_range["status"] == "covers whole range").all())


def clip_loggers(loggers, start=None, end=None):
    """Keep records with start <= time < end in every logger."""
    out = {}
    for n, g in loggers.items():
        m = np.ones(len(g), dtype=bool)
        if start is not None:
            m &= g.index >= pd.Timestamp(start)
        if end is not None:
            m &= g.index < pd.Timestamp(end)
        out[n] = g[m]
    return out
