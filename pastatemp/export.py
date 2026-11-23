"""
Saving results: all tables to a folder or into one zip, as the app's Download step does.
"""

import io
import os
import zipfile
from pathlib import Path


def clean_out_dir(raw):
    """
    User-typed folder -> absolute Path. Accepts quotes (Explorer's 'Copy as path'),
    ~ and %VARS% / $VARS. Raises ValueError for an empty path, a Windows path on a
    non-Windows machine (e.g. a cloud server), or a missing drive.
    """
    txt = os.path.expandvars(str(raw).strip().strip('"').strip("'").strip())
    if not txt:
        raise ValueError("Type a folder path, e.g. C:\\Users\\you\\Documents\\PASTA_output")
    if os.name != "nt" and len(txt) > 1 and txt[1] == ":":
        raise ValueError("This is running on a server, not on your computer, so it cannot "
                         "save to " + txt + ". Download the zip instead.")
    path = Path(txt).expanduser()
    if os.name == "nt" and path.drive and not Path(path.drive + "\\").exists():
        raise ValueError("Drive " + path.drive + " was not found.")
    return path.resolve()


def output_files(tables, prefix="PASTA"):
    """{table: df} -> {'<prefix>_<table>.csv': df}."""
    return {"{}_{}.csv".format(prefix, k): v for k, v in tables.items() if v is not None}


def tables_to_zip(tables, prefix="PASTA", path=None):
    """Zip every table as <prefix>_<table>.csv. Returns the zip bytes; also writes
    them to `path` if given."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for fn, df in output_files(tables, prefix).items():
            z.writestr(fn, df.to_csv(index=False))
    data = buf.getvalue()
    if path is not None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_bytes(data)
    return data


def save_tables(tables, folder, prefix="PASTA"):
    """Write every table to `folder` (created if needed; typed paths are cleaned with
    clean_out_dir). Returns the written paths."""
    target = clean_out_dir(folder)
    target.mkdir(parents=True, exist_ok=True)
    out = []
    for fn, df in output_files(tables, prefix).items():
        df.to_csv(target / fn, index=False)
        out.append(target / fn)
    return out
