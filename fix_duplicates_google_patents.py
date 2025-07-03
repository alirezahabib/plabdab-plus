#!/usr/bin/env python3
"""fix_duplicates_google_patents.py

Utility script to organise Google Patent HTML pages downloaded for antibody sequence
rows listed in `data.csv`.

Behaviour (for rows 0-999):
1. Load `./data.csv` (expects the header described by the user).
2. For rows whose `url` starts with "https://patents.google.com" **and** whose `crawl` column contains
   the string `google_patents`, check that a file named `<num>.html` exists inside
   `./data_google_patents/`.
3. If the file exists and its destination `<patent_id>.html` is **not already** present inside
   `./data_google_patents/patents/` **and** hasn't been processed during the current run,
   move (rename) the file there.
4. If the destination already exists (either prior or during the run) simply delete the
   redundant `<num>.html` file.
5. Print a descriptive error if the expected source file does not exist.

Run:  python fix_duplicates_google_patents.py
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
from pathlib import Path
from typing import Set

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATA_CSV_PATH = Path("data.csv")
SRC_DIR = Path("data_google_patents")
DEST_DIR = SRC_DIR / "patents"

# Regular expression to extract the patent ID between "/patent/" and the next slash
PATENT_ID_RE = re.compile(r"/patent/([^/]+)/")


def extract_patent_id(url: str) -> str | None:
    """Return the patent identifier from a Google Patents URL or ``None`` if not found."""
    match = PATENT_ID_RE.search(url)
    return match.group(1) if match else None


def load_rows(csv_path: Path, limit: int | None = None) -> list[dict[str, str]]:
    """Read CSV into a list of dictionaries. If *limit* is provided, stop after that many rows."""
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = []
        for i, row in enumerate(reader):
            if limit is not None and i >= limit:
                break
            rows.append(row)
        return rows


def gather_existing_patent_ids(dest_dir: Path) -> Set[str]:
    """Return a set of patent IDs already present in *dest_dir*."""
    if not dest_dir.exists():
        return set()
    return {
        p.stem  # filename without .html
        for p in dest_dir.glob("*.html")
        if p.is_file()
    }


def main() -> None:
    # Sanity checks
    if not DATA_CSV_PATH.exists():
        raise FileNotFoundError(f"Expected CSV file not found: {DATA_CSV_PATH}")

    # Ensure destination directory exists
    DEST_DIR.mkdir(parents=True, exist_ok=True)

    parser = argparse.ArgumentParser(description="Organise Google Patent HTML files")
    parser.add_argument("start", type=int, help="Starting num (inclusive)")
    parser.add_argument("end", type=int, help="Ending num (inclusive)")
    args = parser.parse_args()

    start_num, end_num = args.start, args.end
    if start_num > end_num:
        parser.error("start must be <= end")

    existing_patent_ids = gather_existing_patent_ids(DEST_DIR)
    processed_patent_ids: Set[str] = set()

    # Read full CSV; we will filter by num values.
    rows = load_rows(DATA_CSV_PATH)

    for row in rows:
        url = (row.get("url") or "").strip()
        num = (row.get("num") or "").strip()

        # Validate numeric range
        try:
            num_int = int(num)
        except ValueError:
            print(f"Warning: Non-integer num '{num}' encountered: skipping row")
            continue

        if num_int < start_num or num_int > end_num:
            continue  # Outside desired range

        # Filter rows of interest
        if not url.startswith("https://patents.google.com"):
            continue

        patent_id = extract_patent_id(url)
        if not patent_id:
            print(f"Warning: Could not parse patent ID from URL '{url}' (num={num})")
            continue

        src_file = SRC_DIR / f"{num}.html"
        dest_file = DEST_DIR / f"{patent_id}.html"

        if not src_file.exists():
            print(f"Error: Missing source HTML for num={num} (expected {src_file})")
            continue

        # Skip if already present (either before the run or processed during this run)
        if patent_id in existing_patent_ids or patent_id in processed_patent_ids:
            try:
                src_file.unlink()
                print(f"Info: Duplicate patent '{patent_id}' for num={num}; source file removed.")
            except OSError as e:
                print(f"Warning: Could not delete duplicate source file {src_file}: {e}")
            continue

        # Move (rename) the file to destination
        try:
            shutil.move(str(src_file), str(dest_file))
            processed_patent_ids.add(patent_id)
            print(f"Moved: {src_file} -> {dest_file}")
        except (OSError, shutil.Error) as e:
            print(f"Error: Failed to move {src_file} to {dest_file}: {e}")


if __name__ == "__main__":
    main()
