# Diel PASTA app

Streamlit app for **Paired Air and Stream Temperature Analysis** at the diel scale.
A Python port of `DielPASTA_Share/R/PASTA_dielAnalysisFunctions.R` and the metrics chunk of
`Diel_PASTA_Example.Rmd`, laid out like the Tempest 1D app.

## Run

Double-click `run_diel_pasta.bat`. The first run creates `%USERPROFILE%\diel_pasta_venv` and installs
packages (no admin). After that it launches straight away at http://localhost:8501.

Set the environment variable `DIEL_PASTA_DATA` to a folder path to pre-fill the folder box.

Or, from any Python with the packages in `requirements_diel_pasta_app.txt`:

```
streamlit run diel_pasta_app.py
```

## Workflow

1. **Logger files**: upload CSVs, or point at a folder. HOBOware exports work as-is (the
   `#, Date-Time (EDT), Temperature (°C), ...` layout; °F is converted). Any other CSV works if it has a
   date/time column and a temperature column. The logger name is the filename up to the HOBOware
   date stamp (`HF_005 2026-07-22 12_44_20 EDT.csv` becomes `HF_005`).
2. **Air or stream?**: each logger is sorted automatically. A median daily range above 4 °C, a minimum
   below −1 °C and noisy 15-min steps point to air. A name containing `air` or `stream/water/str` overrides
   the data. You can change or exclude any logger in the table.
3. **Air-stream pairs**: each stream gets one air logger. Matching trailing numbers come first
   (HF_003 ↔ HF_103). After that, the app picks the highest diel cross-correlation: the 24-h-detrended
   signals, with the stream lagging 0–6 h. The **Suggested by** column says
   *confirm (air loggers score alike)* when two air loggers score within 0.02 of each other. At one site
   they usually do, and then the pairing should come from field notes. Override in the table.
4. **Run**: daily OLS fit `temp ~ sin(rday) + cos(rday)` for air and water on every day.
5. **Results**: amplitude ratio, phase lag and mean ratio (°C and K) for days where both fits pass the
   adjusted R² threshold, with 7-day centred rolling means, plus a single-day fit check plot.
6. **Download**: metrics, fits, paired data and a run summary (settings + pairs), or save all four to a folder.

## Differences from the R script

- **Time zone**: `rad_hour()` in R calls `as.POSIXlt(dateTime, "")`, which converts timestamps to the
  computer's local zone. For EST data read on a machine in EDT, every phase comes out 1 h late (phase
  *lag* cancels, except on the DST-change day). The app uses logger clock time as written.
- **Minimum day coverage** (default 75 %): partial days are not fitted. Set it to 0 to fit every day, as R does.
- The phase formula is kept exactly: `24/(2π) · (π − atan(cos/sin))`. That formula can only return 6–18 h.

Checked against `DielPASTA_Share/output`: with coverage 0, amplitude, R², RMSE and mean match
`EG_wma_dielPASTA_fits.csv`, and phase differs by exactly the 1 h offset above. AmpRatio and
Ratio_Mean_K match `EG_dielPASTA_metrics.csv` on all 649 rows.

## CPU-saving techniques (same as Tempest)

- `st.cache_data` on file parsing and gridding (keyed on file bytes), the pair-score matrix, the overview
  plot frame and CSV encoding, so widget changes do not re-parse files.
- Every plot trace is downsampled to ≤ 2,000 points and drawn with `Scattergl`.
- The analysis runs only when **Run Diel PASTA** is pressed. Results live in `st.session_state`, and
  the app warns when settings have changed since the last run.
- Daily fits are solved for all days at once (batched 3×3 normal equations) instead of one `lm()` per day.
- Plot selectors run inside `st.fragment`, so changing a site or day redraws only that plot.

## Files

| File | Purpose |
|---|---|
| `diel_pasta_app.py` | Streamlit UI |
| `pasta_diel.py` | Parsing, classification, pairing, sinusoid fits, metrics (no Streamlit) |
| `requirements_diel_pasta_app.txt` | Packages (auto-installed on launch if missing) |
| `run_diel_pasta.bat` | Windows launcher with its own venv |
