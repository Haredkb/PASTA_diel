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


# Name-match tiers, strongest first. A stream is paired at its strongest tier;
# the logger data only decides between several air loggers in that tier, or
# when no name rule applies at all.
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


def diel_xcorr(stream_s, air_s, interval_min=15, max_lag_h=6):
    """
    Max correlation of the high-passed (24 h detrended) signals with the
    stream lagging air by 0..max_lag_h hours. Returns (r, lag_h, overlap_days).
    """
    df = pd.concat([stream_s.rename("w"), air_s.rename("a")], axis=1, join="inner").dropna()
    per_day = int(round(24 * 60 / interval_min))
    if len(df) < 3 * per_day:
        return float("nan"), float("nan"), len(df) / per_day
    hp = df - df.rolling(per_day, center=True, min_periods=per_day // 2).mean()
    w = hp["w"].to_numpy()
    a = hp["a"].to_numpy()
    best = (-1.0, 0.0)
    for k in range(0, int(max_lag_h * 60 / interval_min) + 1):
        ww, aa = (w[k:], a[:len(a) - k]) if k else (w, a)
        ok = np.isfinite(ww) & np.isfinite(aa)
        if ok.sum() < per_day:
            continue
        r = np.corrcoef(ww[ok], aa[ok])[0, 1]
        if r > best[0]:
            best = (r, k * interval_min / 60)
    return float(best[0]), float(best[1]), len(df) / per_day


def suggest_pairs(streams, airs, score):
    """
    streams, airs: lists of names; score: dict[(stream, air)] -> diel r.

    1. Names first: each stream takes the air logger(s) at its strongest
       NAME_TIERS level. One candidate -> that one; several -> the best diel r
       among them. Name-matched air loggers may serve several streams.
    2. Data last: streams with no name match get the best-scoring air logger
       not already used, then the best-scoring air logger overall.

    Returns dict stream -> (air or None, basis).
    """
    out, used = {}, set()
    for s in streams:
        levels = {a: name_match_level(s, a) for a in airs}
        levels = {a: lv for a, lv in levels.items() if lv is not None}
        if not levels:
            continue
        best = min(levels.values())
        hits = [a for a, lv in levels.items() if lv == best]
        basis = "name: " + NAME_TIERS[best]
        if len(hits) > 1:
            hits.sort(key=lambda a: -np.nan_to_num(score.get((s, a), np.nan), nan=-9))
            basis += " ({} candidates; best diel r)".format(len(hits))
        out[s] = (hits[0], basis)
        used.add(hits[0])

    remaining = [s for s in streams if s not in out]
    cands = sorted(
        ((score.get((s, a), float("nan")), s, a) for s in remaining for a in airs),
        key=lambda t: -np.nan_to_num(t[0], nan=-9),
    )
    for r, s, a in cands:
        if s in out or a in used or not np.isfinite(r):
            continue
        out[s] = (a, "data: best diel correlation")
        used.add(a)
    for s in remaining:
        if s in out:
            continue
        opts = [(score.get((s, a), float("nan")), a) for a in airs]
        opts = [o for o in opts if np.isfinite(o[0])]
        out[s] = ((max(opts)[1], "data: best diel correlation (shared air)") if opts
                  else (None, "no name match and no overlapping data"))
    return out


def basis_note(stream, air, basis, airs, score):
    """Note for the pair table: the basis, plus a 'confirm' flag when a
    data-based pick is uncertain."""
    if basis.startswith("name") or air is None:
        return basis
    rs = sorted((score[(stream, x)] for x in airs
                 if np.isfinite(score.get((stream, x), np.nan))), reverse=True)
    if len(rs) > 1 and rs[0] - rs[1] < 0.02:
        return basis + " - confirm (air loggers score alike)"
    if rs and rs[0] < 0.3:
        return basis + " - confirm (weak diel link)"
    return basis


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
    Suggest one air logger for every stream logger.

    loggers: {name: gridded Series}; types: {name: 'air' | 'stream' | 'exclude'}
    (default: classify_loggers()). Returns a DataFrame with columns
    site_id, stream, air, basis, diel_r, lag_h, overlap_days. 'basis' says
    'confirm ...' when the data cannot tell the air loggers apart.
    """
    if types is None:
        types = dict(zip(*classify_loggers(loggers)[["logger", "type"]].T.values))
    airs = [n for n, t in types.items() if t == "air"]
    streams = [n for n, t in types.items() if t == "stream"]
    scores = {(s, a): diel_xcorr(loggers[s], loggers[a], interval_min)
              for s in streams for a in airs}
    sugg = suggest_pairs(streams, airs, {k: v[0] for k, v in scores.items()})
    rows = []
    for s in streams:
        a, basis = sugg[s]
        check = basis_note(s, a, basis, airs, {k: v[0] for k, v in scores.items()})
        r, lag, days = scores[(s, a)] if a else (np.nan, np.nan, np.nan)
        rows.append({"site_id": site_label(s), "stream": s, "air": a, "basis": check,
                     "diel_r": r, "lag_h": lag, "overlap_days": days})
    out = pd.DataFrame(rows, columns=["site_id", "stream", "air", "basis", "diel_r",
                                      "lag_h", "overlap_days"])
    dup = out["site_id"].duplicated(keep=False)        # shortened labels collided
    out.loc[dup, "site_id"] = out.loc[dup, "stream"]
    out.attrs["scores"] = scores
    return out
