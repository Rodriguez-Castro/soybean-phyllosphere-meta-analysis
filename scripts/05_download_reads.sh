#!/usr/bin/env bash
# Step 05: download the complete reads of every run in the sample sheet.
#
# Usage:
#   bash scripts/05_download_reads.sh <sample_sheet.tsv> <raw_dir> [threads]
#
# Writes <raw_dir>/<run>_1.fastq.gz and <run>_2.fastq.gz (or <run>.fastq.gz if single-end)
# and a <run>.done marker. Resumable: runs with a .done marker are skipped.
# Failed runs are listed in logs/05_download_failed.tsv.
# Requires: sra-toolkit (fasterq-dump), run from the repository root.

set -uo pipefail

SHEET=${1:?sample sheet required}
RAW=${2:?output directory required}
THREADS=${3:-4}

mkdir -p "$RAW/tmp" logs
GZ=$(command -v pigz > /dev/null && echo "pigz -p $THREADS" || echo gzip)

total=$(($(wc -l < "$SHEET") - 1)); i=0

tail -n +2 "$SHEET" | cut -f1 | while read -r run; do
  i=$((i + 1))
  [[ -f "$RAW/${run}.done" ]] && continue
  echo "[$i/$total] $run"
  rm -f "$RAW/${run}"*.fastq "$RAW/${run}"*.fastq.gz

  if fasterq-dump --split-files --skip-technical --threads "$THREADS" --temp "$RAW/tmp" \
       -O "$RAW" "$run" < /dev/null > /dev/null 2>> logs/05_download_errors.log \
     && ls "$RAW/${run}"*.fastq > /dev/null 2>&1; then
    $GZ "$RAW/${run}"*.fastq && touch "$RAW/${run}.done"
  else
    echo -e "$run\tdownload_failed" >> logs/05_download_failed.tsv
    echo "  download failed"
  fi
done

rm -rf "$RAW/tmp"
echo "Done: $(ls "$RAW"/*.done 2> /dev/null | wc -l) of $total runs downloaded in $RAW"
