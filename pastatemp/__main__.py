"""
Command line:  pastatemp FOLDER [--out OUT] [--zip] [--only annual|diel]
               [--date-range all|common|START:END] [--year-type water|calendar]
               [--pair SITE STREAM AIR ...]
(or python -m pastatemp ...)
"""

import argparse
import sys
from pathlib import Path

from .core import ANALYSES, run_folder
from .export import save_tables


def _date_range(txt):
    if txt in ("all", "common"):
        return txt
    if ":" in txt:
        a, b = txt.split(":", 1)
        return (a, b)
    raise argparse.ArgumentTypeError("use all, common or START:END, e.g. 2024-10-01:2025-09-30")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pastatemp",
                                 description="Annual and diel paired air-stream temperature analysis "
                                             "(PASTA) for a folder of logger CSVs.")
    ap.add_argument("folder", help="folder with air and stream logger CSVs")
    ap.add_argument("--out", help="output folder (default: FOLDER/PASTA_output)")
    ap.add_argument("--zip", action="store_true", help="also write every table into one .zip")
    ap.add_argument("--only", choices=ANALYSES, help="run only one analysis (default: both)")
    ap.add_argument("--date-range", type=_date_range, default="all",
                    help="all (default; each site's full record), common (dates every site "
                         "shares) or START:END")
    ap.add_argument("--year-type", choices=["water", "calendar"], default="water")
    ap.add_argument("--phase-method", choices=["atan2", "atan"], default="atan2",
                    help="atan2 (default): phase = fitted peak time; atan: the PASTA R formula")
    ap.add_argument("--interval", type=int, default=15, help="logging interval, minutes")
    ap.add_argument("--water-min", type=float, default=0.99,
                    help="drop stream records at or below this, °C (default 0.99)")
    ap.add_argument("--min-coverage", type=float, default=75,
                    help="minimum %% of a full day of records (default 75)")
    ap.add_argument("--min-year-pct", type=float, default=25,
                    help="annual, single year: %% of the full water/calendar year (365/366 "
                         "days) with stream > 1 °C days needed (default 25)")
    ap.add_argument("--min-record-pct", type=float, default=75,
                    help="annual, whole record: %% of the record span with stream > 1 °C days "
                         "needed (default 75)")
    ap.add_argument("--r2-min", type=float, default=0.7, help="diel: minimum adjusted R² (default 0.7)")
    ap.add_argument("--roll-k", type=int, default=7, help="diel: rolling mean window, days (default 7)")
    ap.add_argument("--pair", nargs=3, action="append", metavar=("SITE", "STREAM", "AIR"),
                    help="set a pair yourself (repeat for each site); overrides suggestions")
    a = ap.parse_args(argv)

    res = run_folder(a.folder, analyses=(a.only,) if a.only else ANALYSES, pairs=a.pair,
                     interval_min=a.interval, date_range=a.date_range, yr_type=a.year_type,
                     phase_method=a.phase_method, water_min_C=a.water_min,
                     min_coverage=a.min_coverage / 100, min_year_pct=a.min_year_pct,
                     min_record_pct=a.min_record_pct, r2_min=a.r2_min, roll_k=a.roll_k)
    out = a.out or str(Path(a.folder) / "PASTA_output")
    prefix = res.default_prefix()
    written = save_tables(res.tables(), out, prefix)
    if a.zip:
        res.to_zip(Path(out) / (prefix + ".zip"))
    print("Pairs:")
    print(res.pairs[["site_id", "stream", "air", "basis"]].to_string(index=False))
    print("\nDate range: {} to {}".format(res.date_range[0].date(),
                                         (res.date_range[1] - __import__("pandas").Timedelta(days=1)).date()))
    print("\nSummary:")
    print(res.summary.to_string(index=False))
    print("\nWrote {} files to {}".format(len(written) + int(a.zip), out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
