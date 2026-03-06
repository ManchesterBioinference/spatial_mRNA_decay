"""
Download and process Berrocal_2020 supplementary data.

Downloads the eLife 61635 supplementary zip, which contains a nested
``Berrocal_2020_File_S1.tgz``.  Extracts only the Data/ CSVs and the
Berrocal_2020.ipynb notebook from that tarball (Figures/ and Movies/ are
discarded to save disk space), then executes the first 17 cells of the
notebook plus a new save cell to produce:

    data/Berrocal_2020/Data/eve_data_longform_w_nuclei_060520_FILTERED.csv

The notebook is modified in-place (cells 18+ removed, save cell appended) so
the file on disk records exactly what was executed.

Usage
-----
    uv run scripts/00_download_berrocal_data.py
    uv run scripts/00_download_berrocal_data.py --output-dir path/to/Berrocal_2020
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "pandas==2.3.3",
#   "numpy==2.3.5",
#   "matplotlib==3.10.8",
#   "scipy==1.17.0",
#   "seaborn==0.13.2",
#   "statsmodels==0.14.6",
#   "nbformat==5.10.4",
#   "nbclient==0.10.4",
#   "ipykernel==7.2.0",
# ]
# ///

import argparse
import io
import tarfile
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import nbclient
import nbformat

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DOWNLOAD_URL = "https://cdn.elifesciences.org/articles/61635/elife-61635-supp1-v2.zip"
DEFAULT_OUTPUT_DIR = Path("data/Berrocal_2020")

# Number of notebook cells to keep (cells 1-17, i.e. indices 0-16).
# Cell 18 (index 17) onward is removed; a new save cell is appended as cell 18.
N_CELLS_TO_KEEP = 17

# Output path relative to the notebook working directory (data/Berrocal_2020/)
FILTERED_CSV_REL = "Data/eve_data_longform_w_nuclei_060520_FILTERED.csv"

# The save expression appended as the new cell 18
SAVE_CELL_SOURCE = f"eve.to_csv('{FILTERED_CSV_REL}', index=False)"

# Replacement source for cell 1 (the original imports block).
# The original cell pulls in h5py, sklearn, matplotlib.animation, etc. which
# are not needed to produce the filtered CSV.
MINIMAL_IMPORTS_CELL = """\
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance
""".strip()

# Name of the tarball inside the outer zip
INNER_TGZ_NAME = "Berrocal_2020_File_S1.tgz"

# Prefix inside the tarball that all Berrocal files live under
TGZ_PREFIX = "Berrocal_2020/"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def download_zip(url: str, dest: Path) -> None:
    """Download *url* to *dest*, printing progress."""
    print(f"Downloading:\n  {url}", flush=True)
    urllib.request.urlretrieve(url, dest)
    size_mb = dest.stat().st_size / (1024 * 1024)
    print(f"  -> saved to {dest} ({size_mb:.1f} MB)", flush=True)


def extract_selective(zip_path: Path, output_dir: Path) -> None:
    """
    Extract only Data/ CSVs and the notebook from the nested archive.

    The outer zip contains a single ``Berrocal_2020_File_S1.tgz``.  That
    tarball holds a ``Berrocal_2020/`` directory.  This function:

    1. Reads the ``.tgz`` directly from the zip (no temp file needed).
    2. Selectively extracts ``Data/`` CSVs and ``Berrocal_2020.ipynb``.
    3. Strips the ``Berrocal_2020/`` prefix so files land in *output_dir*.

    Figures/ and Movies/ are skipped.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as zf:
        # Find the inner tgz (tolerate minor name variations)
        tgz_names = [n for n in zf.namelist() if n.endswith(".tgz") and "MACOSX" not in n]
        if not tgz_names:
            raise FileNotFoundError(
                f"No .tgz found inside {zip_path}. "
                f"Contents: {zf.namelist()}"
            )
        tgz_name = tgz_names[0]
        print(f"  Found inner archive: {tgz_name}", flush=True)

        # Read tgz bytes into memory so we don't need a second temp file
        tgz_bytes = io.BytesIO(zf.read(tgz_name))

    with tarfile.open(fileobj=tgz_bytes, mode="r:gz") as tf:
        for member in tf.getmembers():
            rel = member.name  # e.g. "Berrocal_2020/Data/eve_data_key.csv"

            # Only keep items inside Data/ or the top-level notebook
            is_data = rel.startswith(f"{TGZ_PREFIX}Data/")
            is_notebook = rel == f"{TGZ_PREFIX}Berrocal_2020.ipynb"
            if not (is_data or is_notebook):
                continue

            # Strip the "Berrocal_2020/" prefix → relative path inside output_dir
            stripped = rel[len(TGZ_PREFIX):]
            if not stripped:
                continue

            dest = output_dir / stripped
            if member.isdir():
                dest.mkdir(parents=True, exist_ok=True)
            else:
                dest.parent.mkdir(parents=True, exist_ok=True)
                fobj = tf.extractfile(member)
                if fobj is not None:
                    with open(dest, "wb") as tgt:
                        tgt.write(fobj.read())
                    print(f"  extracted: {dest}", flush=True)


def build_modified_notebook(
    nb_path: Path,
    n_keep: int,
    save_source: str,
) -> nbformat.NotebookNode:
    """
    Read *nb_path*, retain the first *n_keep* cells, replace the imports in
    cell 1 with the minimal set, and append a new save cell.

    The modified notebook is written back to *nb_path* so the file on disk
    documents exactly what was executed.
    """
    nb = nbformat.read(nb_path, as_version=4)

    original_n = len(nb.cells)
    nb.cells = nb.cells[:n_keep]
    print(
        f"  Notebook: kept {n_keep}/{original_n} cells, removed {original_n - n_keep}",
        flush=True,
    )

    # Replace the first cell with minimal imports (drops h5py, sklearn, etc.)
    nb.cells[0].source = MINIMAL_IMPORTS_CELL
    print("  Replaced cell 1 with minimal imports", flush=True)

    save_cell = nbformat.v4.new_code_cell(source=save_source)
    nb.cells.append(save_cell)
    print(f"  Appended save cell: {save_source}", flush=True)

    nbformat.write(nb, nb_path)
    print(f"  Modified notebook written to {nb_path}", flush=True)
    return nb


def execute_notebook(nb: nbformat.NotebookNode, working_dir: Path) -> None:
    """
    Execute *nb* with the kernel working directory set to *working_dir*.

    Uses nbclient so no display / X11 is required.  The ``MPLBACKEND=Agg``
    environment variable should be set by the caller (or the Snakemake shell
    command) to suppress matplotlib GUI backends.
    """
    print(f"Executing notebook (working dir: {working_dir}) ...", flush=True)
    client = nbclient.NotebookClient(
        nb,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(working_dir)}},
    )
    client.execute()
    print("  -> notebook execution complete", flush=True)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination directory (default: data/Berrocal_2020)",
    )
    parser.add_argument(
        "--url",
        default=DOWNLOAD_URL,
        help="Download URL for the supplementary zip",
    )
    args = parser.parse_args()

    output_dir: Path = args.output_dir.resolve()

    # --- Download ---
    with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        download_zip(args.url, tmp_path)
        extract_selective(tmp_path, output_dir)
    finally:
        tmp_path.unlink(missing_ok=True)

    # --- Modify notebook ---
    nb_path = output_dir / "Berrocal_2020.ipynb"
    if not nb_path.exists():
        raise FileNotFoundError(
            f"Notebook not found after extraction: {nb_path}\n"
            "Check that the zip structure matches the expected layout."
        )
    nb = build_modified_notebook(nb_path, N_CELLS_TO_KEEP, SAVE_CELL_SOURCE)

    # --- Execute notebook ---
    execute_notebook(nb, output_dir)

    # --- Verify output ---
    filtered_csv = output_dir / FILTERED_CSV_REL
    if filtered_csv.exists():
        import pandas as pd

        df = pd.read_csv(filtered_csv)
        print(
            f"\nSUCCESS: {filtered_csv}\n"
            f"  rows={len(df):,}  cols={df.shape[1]}",
            flush=True,
        )
    else:
        raise FileNotFoundError(
            f"Expected output not found after notebook execution:\n  {filtered_csv}\n"
            "Check the log for notebook errors."
        )


if __name__ == "__main__":
    main()
