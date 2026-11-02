# PASTA app (annual + diel)

Streamlit app for **Paired Air and Stream Temperature Analysis** at the annual and diel scales.
A Python port of two R code bases, laid out like the Tempest 1D app:

- **Annual signal**: the PASTA Shiny app (`900_SupportScripts/PASTA`: `R/paired_annualT_signals.R`,
  `R/stream_thermal_sensitivity.R`, `therm_analysis()` in `R/rundata_v1.R`).
- **Diel**: `DielPASTA_Share/R/PASTA_dielAnalysisFunctions.R` and the metrics chunk of `Diel_PASTA_Example.Rmd`.

Both run by default; either can be switched off.

## Run

Double-click `run_diel_pasta.bat`. The first run creates `%USERPROFILE%\diel_pasta_venv` and installs
packages (no admin). After that it launches straight away at http://localhost:8501.

Set the environment variable `DIEL_PASTA_DATA` to a folder path to pre-fill the folder box.

Or, from any Python with the packages in `requirements.txt`:

```
streamlit run diel_pasta_app.py
```

## Share online (Streamlit Community Cloud)

1. Go to https://share.streamlit.io and sign in with GitHub.
2. **Create app** → *Deploy a public app from GitHub* → repository `Haredkb/PASTA_diel`, branch `main`,
   main file `diel_pasta_app.py` → **Deploy**.
3. Share the `*.streamlit.app` link. Because the GitHub repo is private, the app is private too: invite
   viewers by email under **Share** in the app's top-right menu, or make the repo public to allow anyone.

On the cloud only **Upload CSV files** works (the server can't see your local folders), and results come
back through the download buttons. **Save to folder** writes to the server, not your computer.

## Workflow

1. **Logger files**: upload CSVs, or point at a folder. HOBOware exports work as-is (the
   `#, Date-Time (EDT), Temperature (°C), ...` layout; °F is converted). Any other CSV works if it has a
   date/time column and a temperature column. The logger name is the filename up to the HOBOware
   date stamp (`HF_005 2026-07-22 12_44_20 EDT.csv` becomes `HF_005`). Files with a timestamp column and one
   temperature column named after the logger also work (e.g. `Timestamp (UTC-05:00),LLC_LLC01_AT`).
2. **Air or stream?**: each logger is sorted automatically. A median daily range above 4 °C, a minimum
   below −1 °C and noisy 15-min steps point to air. A name containing `air`/`AT` or `stream/water/str/WT`
   overrides the data. You can change or exclude any logger in the table.
3. **Air-stream pairs**: each stream gets one air logger. Matching names come first: the same name apart
   from the air/water token (`LLC_LLC01_WT_corrected` ↔ `LLC_LLC01_AT_corrected`), or matching trailing
   numbers (HF_003 ↔ HF_103). After that, the app picks the highest diel cross-correlation: the 24-h-detrended
   signals, with the stream lagging 0–6 h. The **Suggested by** column says
   *confirm (air loggers score alike)* when two air loggers score within 0.02 of each other. At one site
   they usually do, and then the pairing should come from field notes. Override in the table.
4. **Run PASTA**: choose *Annual signal*, *Diel* or both (default).
   - **Annual signal**: paired records are averaged to daily means (days below the coverage threshold are
     dropped, as are days with a stream mean ≤ 1 °C). Each site gets a sinusoid `temp ~ sin(rday) + cos(rday)`
     with `rday = 2π·d/365` (d = day of water year or calendar year), for the whole record and for each year.
     Outputs: AmpRatio, PhaseLag_d, Ratio_Mean, daily thermal sensitivity (stream ~ air slope;
     Kelleher et al. 2012), longest data gap and day count. Years with fewer than 100 days (adjustable) get
     no signal metrics. A warning appears when a site has less than ~300 days of record.
   - **Diel**: the same sinusoid with `rday = 2π·h/24`, fitted to every day for air and stream. Days where
     both fits pass the adjusted R² threshold give AmpRatio, PhaseLag_h and mean ratios (°C and K), with
     7-day centred rolling means.
5. **Results**: a combined summary per site, then an *Annual signal* tab (record and yearly tables, daily
   means with the annual fits, stream vs air with the thermal-sensitivity line, AmpRatio vs phase lag by year)
   and a *Diel* tab (daily metrics through time and a single-day fit check).
6. **Download**: run summary and paired data, plus annual metrics (record and by year), annual fits and daily
   means, and diel metrics and fits. Or save them all to a folder.

## Differences from the R code

Annual signal:

- **Year type**: in R, `therm_analysis()` never passes `yr_type` on to `fit_TAS()`, so radian days are always
  water-year days even when *calendar* is chosen (calendar only changes the yearly grouping). The app uses
  the chosen year type for both. AmpRatio is unaffected and phase lag only changes when a phase wraps.
- **Thermal-sensitivity columns** are named `TS__Slope`, `TS_AdjRsqr` and `TS_YInt` (R: `TS__Slope`,
  `AdjRsqr`, `YInt`) so they are not confused with the sinusoid fit's R² and intercept.
- Daily means are computed from the paired sub-daily records, so air and stream means cover the same times.

Diel:

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
- The analysis runs only when **Run PASTA** is pressed. Results live in `st.session_state`, and
  the app warns when settings have changed since the last run.
- Diel fits are solved for all days at once (batched 3×3 normal equations) instead of one `lm()` per day.
- Plot selectors run inside `st.fragment`, so changing a site or day redraws only that plot.

## Files

| File | Purpose |
|---|---|
| `diel_pasta_app.py` | Streamlit UI |
| `pasta_diel.py` | Parsing, classification, pairing, diel sinusoid fits and metrics (no Streamlit) |
| `pasta_annual.py` | Daily means, annual sinusoid fits, thermal sensitivity, gap check, annual metrics (no Streamlit) |
| `requirements.txt` | Packages (auto-installed on launch if missing) |
| `run_diel_pasta.bat` | Windows launcher with its own venv |
