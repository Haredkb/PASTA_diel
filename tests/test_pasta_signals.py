"""
Synthetic checks: build air and stream records with known annual and diel
amplitudes and lags, and make sure the package recovers them.

Run:  python -m pytest tests
"""

import numpy as np
import pandas as pd
import pytest

import pasta_signals as ps

# Known signals
ANN_AMP_AIR, ANN_AMP_WAT, ANN_LAG_D = 12.0, 6.0, 10.0      # annual AmpRatio 0.5, lag 10 d
DIEL_AMP_AIR, DIEL_AMP_WAT, DIEL_LAG_H = 5.0, 1.0, 2.0     # diel AmpRatio 0.2, lag 2 h
AIR_PEAK_H = 14.0


def _synthetic(start="2024-10-01", end="2026-09-30 23:45", freq="15min", seed=0):
    t = pd.date_range(start, end, freq=freq)
    doy = (t - pd.Timestamp("2024-09-30")).days.to_numpy() % 365
    hour = t.hour.to_numpy() + t.minute.to_numpy() / 60
    ann = lambda amp, lag: amp * np.cos(2 * np.pi * (doy - 290 - lag) / 365)
    dl = lambda amp, peak: amp * np.cos(2 * np.pi * (hour - peak) / 24)
    rng = np.random.default_rng(seed)
    air = 10 + ann(ANN_AMP_AIR, 0) + dl(DIEL_AMP_AIR, AIR_PEAK_H) + rng.normal(0, 0.2, len(t))
    wat = 12 + ann(ANN_AMP_WAT, ANN_LAG_D) + dl(DIEL_AMP_WAT, AIR_PEAK_H + DIEL_LAG_H) \
        + rng.normal(0, 0.05, len(t))
    return t, wat, air


@pytest.fixture(scope="module")
def data():
    t, wat, air = _synthetic()
    return pd.DataFrame({"site": "S1", "time": t, "stream": wat, "air": air})


def test_run_diel_recovers_ratio_and_lag(data):
    res = ps.run_diel(data, water_min_C=-99)
    m = res.metrics
    assert len(m) > 600
    assert m["AmpRatio"].median() == pytest.approx(DIEL_AMP_WAT / DIEL_AMP_AIR, abs=0.02)
    assert m["PhaseLag_h"].median() == pytest.approx(DIEL_LAG_H, abs=0.1)
    assert res.settings["interval_min"] == 15


def test_run_annual_recovers_ratio_and_lag(data):
    res = ps.run_annual(data, yr_type="water", water_min_C=-99)
    rec = res.record.iloc[0]
    assert rec["AmpRatio"] == pytest.approx(ANN_AMP_WAT / ANN_AMP_AIR, abs=0.02)
    assert rec["PhaseLag_d"] == pytest.approx(ANN_LAG_D, abs=1.0)
    assert list(res.yearly["year"]) == [2025, 2026]
    assert res.yearly["AmpRatio"].notna().all()


def test_annual_accepts_daily_means(data):
    d = (data.assign(time=data["time"].dt.normalize())
             .groupby(["site", "time"], as_index=False)[["stream", "air"]].mean())
    res = ps.run_annual(d)
    assert res.settings["interval_min"] is None
    assert res.record.iloc[0]["AmpRatio"] == pytest.approx(0.5, abs=0.02)


def test_diel_rejects_daily_means(data):
    d = (data.assign(time=data["time"].dt.normalize())
             .groupby(["site", "time"], as_index=False)[["stream", "air"]].mean())
    with pytest.raises(ValueError):
        ps.run_diel(d)


def test_short_years_get_no_signal_metrics(data):
    short = data[data["time"] < "2024-12-31"]      # 91 days
    res = ps.run_annual(short, min_days=100)
    assert res.yearly["AmpRatio"].isna().all()
    assert not res.yearly["meets_min_days"].any()
    # no qualifying year -> no whole-record metrics either
    rec = res.record.iloc[0]
    assert np.isnan(rec["AmpRatio"]) and rec["years_used"] == "" and rec["years_dropped"] == "2025"


def test_whole_record_drops_short_years(data):
    # WY2025 complete; WY2026 cut to 60 days -> WY2026 must be left out of the record fit
    d = data[data["time"] < "2025-11-30"]
    res = ps.run_annual(d, min_days=100, water_min_C=-99)
    rec = res.record.iloc[0]
    assert rec["years_used"] == "2025" and rec["years_dropped"] == "2026"
    full_wy2025 = ps.run_annual(data[data["time"] < "2025-10-01"], water_min_C=-99).record.iloc[0]
    for c in ("AmpRatio", "PhaseLag_d", "Ratio_Mean", "TS__Slope", "n_days_fit"):
        assert rec[c] == full_wy2025[c]                 # identical to fitting WY2025 alone
    yr = res.yearly.set_index("year")
    assert yr.loc[2025, "meets_min_days"] and not yr.loc[2026, "meets_min_days"]
    assert np.isnan(yr.loc[2026, "AmpRatio"]) and np.isfinite(yr.loc[2026, "TS__Slope"])
    f = res.fits
    assert set(f.loc[f["year"] == "2026", "meets_min_days"]) == {False}


def test_min_days_counts_days_in_fit_not_cold_days():
    # 150 paired days in WY2025, but 80 have stream <= 1 °C -> only 70 days enter the fit
    days = pd.date_range("2024-10-01", periods=150, freq="D")
    wat = np.r_[np.full(80, 0.5), np.linspace(3, 15, 70)]
    air = np.r_[np.full(80, -5.0), np.linspace(5, 25, 70)]
    df = pd.DataFrame({"s": "X", "date": days, "w": wat, "a": air})
    res = ps.run_annual(df, min_days=100)
    y = res.yearly.iloc[0]
    assert y["count"] == 150 and y["n_days_fit"] == 70
    assert not y["meets_min_days"] and np.isnan(y["AmpRatio"])
    assert res.record.iloc[0]["years_dropped"] == "2025"


def test_run_pasta_both_and_one(data):
    both = ps.run_pasta(data, water_min_C=-99)
    assert both.annual is not None and both.diel is not None
    assert {"annual_AmpRatio", "diel_median_AmpRatio"} <= set(both.summary.columns)
    only = ps.run_pasta(data, analyses=("diel",), water_min_C=-99)
    assert only.annual is None and "annual_AmpRatio" not in only.summary.columns


def test_run_folder_reads_pairs_and_writes(tmp_path):
    t, wat, air = _synthetic(start="2025-04-01", end="2025-06-30 23:45")
    stamp = t.strftime("%m/%d/%Y %H:%M")
    pd.DataFrame({"Timestamp (UTC-05:00)": stamp, "SITE_A01_AT": air}).to_csv(
        tmp_path / "SITE_A01_AT_corrected.csv", index=False)
    pd.DataFrame({"Timestamp (UTC-05:00)": stamp, "SITE_A01_WT": wat}).to_csv(
        tmp_path / "SITE_A01_WT_corrected.csv", index=False)

    types = ps.classify_loggers(ps.read_loggers(tmp_path))
    assert dict(zip(types["logger"], types["type"])) == {
        "SITE_A01_AT_corrected": "air", "SITE_A01_WT_corrected": "stream"}

    out = tmp_path / "out"
    res = ps.run_folder(tmp_path, out_dir=out)
    p = res.pairs.iloc[0]
    assert (p["site_id"], p["air"]) == ("SITE_A01", "SITE_A01_AT_corrected")
    assert p["basis"] == "name: same name except air/stream code"
    assert res.diel.metrics["AmpRatio"].median() == pytest.approx(0.2, abs=0.02)
    written = sorted(f.name for f in out.glob("*.csv"))
    assert "PASTA_run_summary.csv" in written and "PASTA_diel_metrics.csv" in written


@pytest.mark.parametrize("stream, air, level", [
    ("ST_002", "AT_002", 1),                                       # same name except code
    ("LLC_LLC01_WT_corrected", "LLC_LLC01_AT_corrected", 1),
    ("Str_Brook_02", "Air_Brook_2", 2),                            # same site name and number
    ("SiteB_1_WT", "SiteA_1_AT", 3),                               # same number only
    ("HF-003", "HF_103", 4),                                       # same last two digits
    ("ST_002", "AT_003", None),
    ("HF_005", "HF_101", None),
])
def test_name_match_levels(stream, air, level):
    from pasta_signals.pairing import name_match_level
    assert name_match_level(stream, air) == level


def test_names_beat_data():
    from pasta_signals.pairing import suggest_pairs
    # data prefers the wrong air logger for every stream; names must still win
    score = {("ST_001", "AT_002"): 0.9, ("ST_001", "AT_001"): 0.5,
             ("ST_002", "AT_001"): 0.9, ("ST_002", "AT_002"): 0.5}
    out = suggest_pairs(["ST_001", "ST_002"], ["AT_001", "AT_002"], score)
    assert out["ST_001"][0] == "AT_001" and out["ST_002"][0] == "AT_002"
    assert out["ST_001"][1].startswith("name")
    # no name match -> data as last resort
    out = suggest_pairs(["Brook"], ["AT_001", "AT_002"], {("Brook", "AT_002"): 0.8,
                                                          ("Brook", "AT_001"): 0.4})
    assert out["Brook"] == ("AT_002", "data: best diel correlation")


# ---- phase method -------------------------------------------------------------

def _phase_case(peak_air_h, peak_wat_h, ann_peak_air, ann_peak_wat):
    t = pd.date_range("2023-10-01", "2025-09-30 23:45", freq="15min")
    doy = t.dayofyear.to_numpy()
    h = t.hour.to_numpy() + t.minute.to_numpy() / 60
    air = 10 + 12 * np.cos(2 * np.pi * (doy - ann_peak_air) / 365) \
        + 5 * np.cos(2 * np.pi * (h - peak_air_h) / 24)
    wat = 12 + 6 * np.cos(2 * np.pi * (doy - ann_peak_wat) / 365) \
        + 1 * np.cos(2 * np.pi * (h - peak_wat_h) / 24)
    return pd.DataFrame({"s": "X", "t": t, "w": wat, "a": air})


def test_phase_helpers():
    from pasta_signals.phase import lag, phase
    w = 2 * np.pi * 150 / 365                                     # peak on day 150
    assert phase(np.sin(w), np.cos(w), 365, "atan2") == pytest.approx(150)
    assert phase(np.sin(w), np.cos(w), 365, "atan") == pytest.approx(332.5)   # R: half a year late
    w = 2 * np.pi * 15 / 24                                       # peak at 15:00
    assert phase(np.sin(w), np.cos(w), 24, "atan2") == pytest.approx(15)
    assert phase(np.sin(w), np.cos(w), 24, "atan") == pytest.approx(9)        # R: peak - 6 h
    assert lag(1, 23, 24, "atan2") == pytest.approx(2)            # wrapped across midnight
    assert lag(1, 23, 24, "atan") == pytest.approx(-22)
    with pytest.raises(ValueError):
        phase(1, 1, 24, "bad")


def test_default_is_atan2():
    import inspect
    for fn in (ps.run_annual, ps.run_diel):
        assert inspect.signature(fn).parameters["phase_method"].default == "atan2"


def test_atan2_fixes_late_diel_peak():
    df = _phase_case(15, 1, 200, 210)                  # stream peaks 01:00, 10 h after air
    lag2 = ps.run_diel(df, water_min_C=-99).metrics["PhaseLag_h"].median()
    lag1 = ps.run_diel(df, water_min_C=-99, phase_method="atan").metrics["PhaseLag_h"].median()
    assert lag2 == pytest.approx(10, abs=0.1)
    assert lag1 == pytest.approx(-2, abs=0.1)          # the R formula's 12 h jump


def test_atan2_fixes_annual_peaks_straddling_july():
    df = _phase_case(15, 17, 175, 190)                 # calendar peaks either side of 2 July
    a2 = ps.run_annual(df, yr_type="calendar", water_min_C=-99).record.iloc[0]
    a1 = ps.run_annual(df, yr_type="calendar", water_min_C=-99, phase_method="atan").record.iloc[0]
    assert a2["PhaseLag_d"] == pytest.approx(15, abs=1)
    assert abs(a1["PhaseLag_d"] - 15) > 100            # the R formula's half-year jump
    assert a1["AmpRatio"] == a2["AmpRatio"]            # amplitude is unaffected
