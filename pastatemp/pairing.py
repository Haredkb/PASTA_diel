"""
Sorting loggers into air vs stream and suggesting air-stream pairs.
"""

import re

import numpy as np
import pandas as pd


# ==============================================================================
#  2. Classification
# ==============================================================================

# Name tokens (bounded by non-letters):
#   air:    air, atm, AT (air temperature)
#   stream: stream, str, ST (stream temperature), water, wat, WT (water temperature), SW, river
_AIR_NAME = re.compile(r"(^|[^a-z])(air|atm|atmo|at)([^a-z]|$)", re.IGNORECASE)
_WAT_NAME = re.compile(r"(^|[^a-z])(str|stream|st|wat|water|wt|sw|riv|river)([^a-z]|$)",
                       re.IGNORECASE)


def logger_features(series):
    """Summary features used to tell air from stream loggers."""
    s = series.dropna()
    daily = s.groupby(s.index.normalize())
    counts = daily.count()
    full = counts[counts >= 0.75 * counts.max()].index
    rng = (daily.max() - daily.min()).loc[full]
    step = s.diff().abs()
    return {
        "start": s.index.min(),
        "end": s.index.max(),
        "n_obs": int(len(s)),
        "mean_C": float(s.mean()),
        "min_C": float(s.min()),
        "max_C": float(s.max()),
        "median_daily_range_C": float(rng.median()) if len(rng) else float("nan"),
        "pct_below_0": float((s < 0).mean() * 100),
        "median_step_C": float(step.median()),
    }


def classify_logger(name, feats, range_threshold=4.0):
    """
    Return (type, confidence, reason).

    Stream water is buffered: small diel range and it cannot sit far below 0 °C.
    Air has a large diel range and noisy 15-min steps, and often freezes.
    A name containing 'air' / 'stream' etc. wins over the data.
    """
    if _AIR_NAME.search(name):
        return "air", "high", "name contains an air keyword"
    if _WAT_NAME.search(name):
        return "stream", "high", "name contains a water keyword"

    rng = feats["median_daily_range_C"]
    votes_air = 0
    reasons = []
    if rng > range_threshold:
        votes_air += 2
        reasons.append("daily range {:.1f} °C > {:.1f}".format(rng, range_threshold))
    else:
        reasons.append("daily range {:.1f} °C ≤ {:.1f}".format(rng, range_threshold))
    if feats["min_C"] < -1.0:
        votes_air += 1
        reasons.append("min {:.1f} °C (below freezing)".format(feats["min_C"]))
    if feats["median_step_C"] > 0.1:
        votes_air += 1
        reasons.append("noisy 15-min steps")

    kind = "air" if votes_air >= 2 else "stream"
    margin = abs(np.log(max(rng, 1e-3) / range_threshold))
    conf = "high" if margin > 0.7 else ("medium" if margin > 0.3 else "low")
    return kind, conf, "; ".join(reasons)


# ==============================================================================
#  3. Pairing
# ==============================================================================

def _digits(name):
    return re.findall(r"\d+", name)


def _type_stem(name):
    """Name with the air/stream code removed and separators normalised:
    'AT_002' -> '002', 'LLC_LLC01_WT_corrected' -> 'llc_llc01_corrected'."""
    n = _AIR_NAME.sub(lambda m: m.group(1) + m.group(3), name)
    n = _WAT_NAME.sub(lambda m: m.group(1) + m.group(3), n)
    return re.sub(r"[_\-\s]+", "_", n).strip("_").lower()


def _letters(stem):
    return re.sub(r"[^a-z]", "", stem)


def site_label(stream_name):
    """Default site ID: 'LLC_LLC01_WT_corrected' -> 'LLC_LLC01'; 'ST_002' -> 'site_002';
    'HF-003' stays 'HF-003'."""
    n = _WAT_NAME.sub(lambda m: m.group(1) + m.group(3), stream_name)
    n = re.sub(r"(^|[_\-\s])(corrected|raw|qc|final|clean|data)(?=$|[_\-\s])", r"\1", n,
               flags=re.IGNORECASE)
    n = re.sub(r"[_\-\s]{2,}", "_", n).strip("_- ")
    if n.isdigit():
        n = "site_" + n
    return n or stream_name


# Name-match tiers, strongest first. A stream is paired at its strongest tier.
# Pairing uses names only, never the logger data.
NAME_TIERS = {
    1: "same name except air/stream code",   # AT_002 <-> ST_002, LLC01_AT <-> LLC01_WT
    2: "same site name and number",          # Air_Brook_2 <-> Str_Brook_02
    3: "same number",                        # Air_2 <-> Brook_02
    4: "same last two digits",               # HF_003 <-> HF_103
    5: "same prefix",                        # AMTH_A01 <-> AMTH_air
}


def name_match_level(stream_name, air_name):
    """Strongest NAME_TIERS level that links the two names, or None."""
    ss, sa = _type_stem(stream_name), _type_stem(air_name)
    typed = _AIR_NAME.search(air_name) or _WAT_NAME.search(stream_name)
    if typed and ss and ss == sa:
        return 1
    sd, ad = _digits(stream_name), _digits(air_name)
    if sd and ad and int(sd[-1]) == int(ad[-1]):
        return 2 if _letters(ss) == _letters(sa) else 3
    if sd and ad and len(sd[-1]) >= 2 and len(ad[-1]) >= 2 and sd[-1][-2:] == ad[-1][-2:]:
        return 4
    s_pre = re.split(r"[_\-\s]", stream_name)[0].lower()
    a_pre = re.split(r"[_\-\s]", air_name)[0].lower()
    if _AIR_NAME.search(air_name) and s_pre == a_pre and len(s_pre) > 2:
        return 5
    return None


def name_match(stream_name, air_name):
    """True when any name rule links the two loggers (see NAME_TIERS)."""
    return name_match_level(stream_name, air_name) is not None


def overlap_days(stream_s, air_s, interval_min=15):
    """Days of records where both loggers have data (paired records / records per day)."""
    both = stream_s.dropna().index.intersection(air_s.dropna().index)
    return len(both) / (24 * 60 / interval_min)


def suggest_pairs(streams, airs):
    """
    Pair each stream logger with an air logger by name only (NAME_TIERS); the
    logger data is not used.

    Each stream takes the air logger(s) at its strongest name level. One
    candidate -> that one; several -> the first by name, flagged 'confirm'.
    A stream with no name match gets the only air logger when there is just
    one, else None (choose it yourself). One air logger may serve several streams.

    Returns dict stream -> (air or None, basis).
    """
    out = {}
    for s in streams:
        levels = {a: name_match_level(s, a) for a in airs}
        levels = {a: lv for a, lv in levels.items() if lv is not None}
        if levels:
            best = min(levels.values())
            hits = sorted(a for a, lv in levels.items() if lv == best)
            basis = "name: " + NAME_TIERS[best]
            if len(hits) > 1:
                basis += " - confirm ({} air loggers match)".format(len(hits))
            out[s] = (hits[0], basis)
        elif len(airs) == 1:
            out[s] = (airs[0], "only air logger - confirm")
        else:
            out[s] = (None, "no name match - choose the air logger")
    return out


def classify_loggers(loggers, range_threshold=4.0):
    """
    loggers: {name: gridded Series}. Returns a DataFrame, one row per logger, with
    the suggested type ('air' / 'stream'), confidence, reason and the features used.
    """
    rows = []
    for name, s in loggers.items():
        f = logger_features(s)
        kind, conf, why = classify_logger(name, f, range_threshold)
        rows.append({"logger": name, "type": kind, "confidence": conf, "reason": why, **f})
    return pd.DataFrame(rows)


def pair_loggers(loggers, types=None, interval_min=15):
    """
    Suggest one air logger for every stream logger, by name only (suggest_pairs()).

    loggers: {name: gridded Series}; types: {name: 'air' | 'stream' | 'exclude'}
    (default: classify_loggers()). Returns a DataFrame with columns
    site_id, stream, air, basis, overlap_days. 'basis' says 'confirm' when the
    pick should be checked, and air is None when no name rule applies.
    """
    if types is None:
        types = dict(zip(*classify_loggers(loggers)[["logger", "type"]].T.values))
    airs = [n for n, t in types.items() if t == "air"]
    streams = [n for n, t in types.items() if t == "stream"]
    sugg = suggest_pairs(streams, airs)
    rows = []
    for s in streams:
        a, basis = sugg[s]
        days = overlap_days(loggers[s], loggers[a], interval_min) if a else np.nan
        rows.append({"site_id": site_label(s), "stream": s, "air": a, "basis": basis,
                     "overlap_days": days})
    out = pd.DataFrame(rows, columns=["site_id", "stream", "air", "basis", "overlap_days"])
    dup = out["site_id"].duplicated(keep=False)        # shortened labels collided
    out.loc[dup, "site_id"] = out.loc[dup, "stream"]
    return out
