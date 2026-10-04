#!/usr/bin/env python3
"""Combine SRA RunInfo / SraRunTable CSVs into a single tab-separated run list.

Usage:
    python3 scripts/00_build_run_list.py <runinfo_dir> <output.tsv> [selected_runs.txt]

The optional third argument is a file with one SRR accession per line;
only those runs are kept.
"""
import csv
import glob
import os
import sys

import socket
if socket.gethostname().startswith("iv"):
    sys.exit("ERROR: you are on the login node. Request a compute node first:\n"
             "  salloc --time=24:00:00 --cpus-per-task=4 --mem=8G --partition bio")

# output column -> possible column names in RunInfo / SraRunTable files
COLUMNS = {
    "run": ["Run"],
    "bioproject": ["BioProject"],
    "library_name": ["LibraryName", "Library Name"],
    "library_strategy": ["LibraryStrategy", "Assay Type"],
    "library_layout": ["LibraryLayout"],
    "platform": ["Platform"],
    "model": ["Model", "Instrument"],
    "spots": ["spots"],
    "avg_length": ["avgLength", "AvgSpotLen"],
    "sample_name": ["SampleName", "Sample Name"],
}


def get(row, names):
    for n in names:
        value = (row.get(n) or "").strip()
        if value:
            return value
    return "NA"


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    indir, out = sys.argv[1], sys.argv[2]
    keep = None
    if len(sys.argv) > 3:
        with open(sys.argv[3]) as fh:
            keep = {line.strip() for line in fh if line.strip()}

    files = sorted(glob.glob(os.path.join(indir, "*.csv")))
    if not files:
        sys.exit(f"No CSV files found in {indir}")

    rows, seen = [], set()
    for path in files:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                run = get(row, COLUMNS["run"])
                # skip empty lines, repeated headers and duplicates
                if run in ("NA", "Run") or run in seen:
                    continue
                if keep is not None and run not in keep:
                    continue
                seen.add(run)
                rows.append([get(row, names) for names in COLUMNS.values()]
                            + [os.path.basename(path)])

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(list(COLUMNS) + ["source_file"])
        writer.writerows(rows)

    print(f"{len(rows)} runs from {len(files)} files written to {out}")
    if keep is not None:
        missing = keep - seen
        if missing:
            print(f"WARNING: {len(missing)} selected runs not found, e.g. {sorted(missing)[:5]}")


if __name__ == "__main__":
    main()
