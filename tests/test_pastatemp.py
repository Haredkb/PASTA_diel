"""
Synthetic checks: build air and stream records with known annual and diel
amplitudes and lags, and make sure the package recovers them.

Run:  python -m pytest tests
"""

import numpy as np
import pandas as pd
import pytest

import pastatemp as ps

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


def test_partial_year_judged_against_full_year(data):
    short = data[data["time"] < "2024-12-31"]      # 91 days, all with data
    res = ps.run_annual(short)
    y = res.yearly.iloc[0]
    assert y["n_days_in_year"] == 365 and y["pct_year_fit"] == 24.9   # 91 / 365
    assert not y["meets_year_pct"] and np.isnan(y["AmpRatio"])
    # the whole record is judged on its own span
    rec = res.record.iloc[0]
    assert rec["meets_record_pct"] and rec["pct_record_fit"] == 100.0


def test_leap_water_year_has_366_days():
    from pastatemp.annual import days_in_year
    assert days_in_year(2024, "water") == 366 and days_in_year(2025, "water") == 365
    assert days_in_year(2024, "calendar") == 366


def _wy2025(gap_days=0, cold_days=0):
    """WY2025, 1 Oct - 30 Sep (365 days): `gap_days` days removed from day 200,
    and `cold_days` days with stream <= 1 °C just before them."""
    days = pd.date_range("2024-10-01", "2025-09-30", freq="D")
    rd = 2 * np.pi * np.arange(len(days)) / 365
    w = 10 + 4 * np.sin(rd)
    w[200 - cold_days:200] = 0.5
    df = pd.DataFrame({"s": "X", "date": days, "w": w, "a": 10 + 8 * np.sin(rd)})
    return df.drop(index=range(200, 200 + gap_days))


def test_year_needs_min_pct_of_the_year():
    y = ps.run_annual(_wy2025(91), min_year_pct=75).yearly.iloc[0]  # 274 / 365 = 75.1 %
    assert y["n_days_in_year"] == 365 and y["n_days_fit"] == 274 and y["pct_year_fit"] == 75.1
    assert y["meets_year_pct"] and np.isfinite(y["AmpRatio"])
    y = ps.run_annual(_wy2025(92), min_year_pct=75).yearly.iloc[0]          # 273 / 365 = 74.8 %
    assert not y["meets_year_pct"] and np.isnan(y["AmpRatio"])
    assert np.isfinite(y["TS__Slope"])                       # thermal sensitivity still reported
    assert ps.run_annual(_wy2025(92), min_year_pct=70).yearly.iloc[0]["meets_year_pct"]


def test_year_pct_default_is_25():
    assert ps.run_annual(_wy2025().iloc[:92]).yearly.iloc[0]["meets_year_pct"]      # 25.2 %
    assert not ps.run_annual(_wy2025().iloc[:91]).yearly.iloc[0]["meets_year_pct"]  # 24.9 %


def test_year_pct_counts_after_1C_drop():
    # 50 days with no records + 50 days with stream <= 1 °C: 265 / 365 = 72.6 %
    y = ps.run_annual(_wy2025(50, cold_days=50), min_year_pct=75).yearly.iloc[0]
    assert y["count"] == 315 and y["n_days_fit"] == 265 and y["pct_year_fit"] == 72.6
    assert not y["meets_year_pct"]


def test_record_pct_counts_days_after_1C_drop():
    # 200 paired days; 80 have stream <= 1 °C -> 120 fit days = 60 % of the record
    days = pd.date_range("2024-10-01", periods=200, freq="D")
    rd = 2 * np.pi * np.arange(120) / 365
    wat = np.r_[np.full(80, 0.5), 8 + 5 * np.sin(rd)]
    air = np.r_[np.full(80, -5.0), 10 + 10 * np.sin(rd)]
    df = pd.DataFrame({"s": "X", "date": days, "w": wat, "a": air})
    res = ps.run_annual(df, min_year_pct=30)                # record: default 75 %
    rec = res.record.iloc[0]
    assert rec["count"] == 200 and rec["n_days_fit"] == 120 and rec["record_span_days"] == 200
    assert rec["pct_record_fit"] == 60.0 and not rec["meets_record_pct"]
    assert np.isnan(rec["AmpRatio"]) and np.isfinite(rec["TS__Slope"])
    y = res.yearly.iloc[0]
    assert y["pct_year_fit"] == 32.9 and y["meets_year_pct"]   # 120 / 365 of the year
    rec = ps.run_annual(df, min_record_pct=60).record.iloc[0]
    assert rec["meets_record_pct"] and np.isfinite(rec["AmpRatio"])
    f = res.fits
    assert set(f.loc[f["year"] == "all", "meets_criteria"]) == {False}


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
    from pastatemp.pairing import name_match_level
    assert name_match_level(stream, air) == level


def test_names_beat_data():
    from pastatemp.pairing import suggest_pairs
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
    from pastatemp.phase import lag, phase
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


# ---- app steps as package functions ---------------------------------------------

def _two_site_loggers():
    """Site A: Oct 2024 - Sep 2025; site B starts 1 Jan 2025 -> common period from 2025-01-01."""
    t, wat, air = _synthetic(end="2025-09-30 23:45")
    s = lambda v, idx: pd.Series(v, index=idx)
    late = t >= "2025-01-01"
    return {"ST_001": s(wat, t), "AT_001": s(air, t),
            "ST_002": s(wat[late], t[late]), "AT_002": s(air[late], t[late])}


def test_coverage_and_date_range():
    lg = _two_site_loggers()
    pairs = [("S1", "ST_001", "AT_001"), ("S2", "ST_002", "AT_002")]
    cov = ps.coverage(lg, pairs)
    c0, c1 = ps.common_period(cov)
    assert c0 == pd.Timestamp("2025-01-01")
    start, end = ps.date_range(cov, "all")
    assert (start, end) == (pd.Timestamp("2024-10-01"), pd.Timestamp("2025-10-01"))
    inr = ps.coverage_in_range(cov, start, end)
    assert not ps.all_sites_cover(inr) and inr.set_index("site_id").loc["S2", "status"].startswith("starts")
    assert ps.all_sites_cover(ps.coverage_in_range(cov, *ps.date_range(cov, "common")))
    assert ps.date_range(cov, ("2025-02-01", "2025-02-28")) == (pd.Timestamp("2025-02-01"),
                                                                 pd.Timestamp("2025-03-01"))


def test_run_loggers_common_period_and_outputs(tmp_path):
    lg = _two_site_loggers()
    res = ps.run_loggers(lg, date_range="common", water_min_C=-99)
    assert list(res.pairs["site_id"]) == ["site_001", "site_002"]
    assert res.date_range[0] == pd.Timestamp("2025-01-01")
    assert res.paired["DateTime_r"].min() >= pd.Timestamp("2025-01-01")
    s = res.summary.set_index("site_id")
    assert s.loc["site_001", "date_from"] == "2025-01-01"
    assert s.loc["site_001", "stream_logger"] == "ST_001" and s.loc["site_001", "air_logger"] == "AT_001"
    assert {"annual_AmpRatio", "annual_pct_record_fit", "diel_median_AmpRatio",
            "min_year_pct"} <= set(res.summary.columns)
    assert list(res.tables()) == ["run_summary", "pairs", "coverage", "paired_df", "annual_metrics",
                                  "annual_metrics_yearly", "annual_fits", "daily_means",
                                  "diel_metrics", "diel_fits"]
    import zipfile, io
    z = zipfile.ZipFile(io.BytesIO(res.to_zip(tmp_path / "r.zip", prefix="X")))
    assert len(z.namelist()) == 10 and "X_run_summary.csv" in z.namelist()
    assert (tmp_path / "r.zip").exists()
    written = ps.save_tables(res.tables(), '"' + str(tmp_path / "out") + '"', "X")
    assert len(written) == 10


def test_run_pasta_rejects_unknown_setting(data):
    with pytest.raises(TypeError):
        ps.run_pasta(data, not_a_setting=100)


def test_plots_build():
    pytest.importorskip("plotly")
    from pastatemp import plots as pp
    lg = _two_site_loggers()
    res = ps.run_loggers(lg, water_min_C=-99)
    cov = ps.coverage(lg, res.pairs)
    figs = [pp.plot_loggers(lg, {"ST_001": "stream", "AT_001": "air"}),
            pp.plot_pair(lg, "ST_001", "AT_001"), pp.plot_coverage(cov, *res.date_range),
            pp.plot_annual_fit(res.annual, "site_001"),
            pp.plot_thermal_sensitivity(res.annual, "site_001"),
            pp.plot_annual_by_year(res.annual), pp.plot_diel_metrics(res.diel),
            pp.plot_diel_day(res.diel, "site_001", "2025-06-01")]
    assert all(len(f.data) for f in figs)


def test_old_name_still_imports():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import pasta_signals
        from pasta_signals import annual
    assert pasta_signals.run_annual is ps.run_annual and annual is ps.annual
