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
