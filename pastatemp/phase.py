"""
Phase of a fitted sinusoid  temp = YInt + s*sin(r) + c*cos(r).

Two methods:

  "atan2" (default)  peak time of the fitted sinusoid, anywhere in the cycle:
                     period/(2π) · (atan2(s, c) mod 2π). Phase lags are wrapped
                     into (−period/2, +period/2].

  "atan"             the formulas in the PASTA / Diel PASTA R code, kept for
                     reproducing those results:
                       annual  365/(2π) · (3π/2 − atan(c/s))
                       diel     24/(2π) · (π − atan(c/s))
                     atan only sees the ratio c/s, so each phase is confined to
                     half a cycle: annual 182.5–365 days (the peak day when the
                     peak falls in that half, otherwise peak + 182.5), diel 6–18 h
                     (the peak hour − 6 h for peaks between noon and midnight).
                     Lags are differences of these values, unwrapped.
"""

import numpy as np

PHASE_METHODS = ("atan2", "atan")
DEFAULT_PHASE_METHOD = "atan2"


def check_method(method):
    if method not in PHASE_METHODS:
        raise ValueError("phase_method must be one of " + str(PHASE_METHODS))
    return method


def phase(s, c, period, method=DEFAULT_PHASE_METHOD):
    """Phase in the units of `period` (365 for days of year, 24 for hours)."""
    s, c = np.asarray(s, dtype=float), np.asarray(c, dtype=float)
    k = period / (2 * np.pi)
    with np.errstate(divide="ignore", invalid="ignore"):
        if check_method(method) == "atan2":
            return k * np.mod(np.arctan2(s, c), 2 * np.pi)
        offset = 1.5 * np.pi if period == 365 else np.pi      # R annual vs diel formula
        return k * (offset - np.arctan(c / s))


def lag(water_phase, air_phase, period, method=DEFAULT_PHASE_METHOD):
    """water − air; wrapped into (−period/2, period/2] for atan2."""
    d = np.asarray(water_phase, dtype=float) - np.asarray(air_phase, dtype=float)
    if check_method(method) == "atan2":
        half = period / 2
        d = half - np.mod(half - d, period)
    return d
