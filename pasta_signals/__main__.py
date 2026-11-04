"""
Command line:  python -m pasta_signals FOLDER [--out OUT] [--only annual|diel]
               [--year-type water|calendar] [--pair SITE STREAM AIR ...]
"""

import argparse
import sys

from .core import ANALYSES, run_folder


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m pasta_signals",
                                 description="Annual and diel paired air-stream temperature analysis "
                                             "for a folder of logger CSVs.")
    ap.add_argument("folder", help="folder with air and stream logger CSVs")
    ap.add_argument("--out", help="output folder (default: FOLDER/PASTA_output)")
    ap.add_argument("--only", choices=ANALYSES, help="run only one analysis (default: both)")
    ap.add_argument("--year-type", choices=["water", "calendar"], default="water")
    ap.add_argument("--interval", type=int, default=15, help="logging interval, minutes")
    ap.add_argument("--pair", nargs=3, action="append", metavar=("SITE", "STREAM", "AIR"),
                    help="set a pair yourself (repeat for each site); overrides suggestions")
    a = ap.parse_args(argv)

    out = a.out or str(__import__("pathlib").Path(a.folder) / "PASTA_output")
    res = run_folder(a.folder, analyses=(a.only,) if a.only else ANALYSES, pairs=a.pair,
                     interval_min=a.interval, out_dir=out, yr_type=a.year_type)
    print("Pairs:")
    print(res.pairs[["site_id", "stream", "air", "basis"]].to_string(index=False))
    print("\nSummary:")
    print(res.summary.to_string(index=False))
    print("\nWrote CSVs to " + out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
