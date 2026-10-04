#!/usr/bin/env bash
# Step 05: download the complete reads of every run in the sample sheet.
#
# Usage:
#   bash scripts/05_download_reads.sh <sample_sheet.tsv> <raw_dir> [threads]
#
# Writes <raw_dir>/<run>_1.fastq.gz and <run>_2.fastq.gz and a <run>.done marker.
# --split-3 keeps R1 and R2 properly paired: reads without mate go to <run>.fastq.gz (not used).
# (--split-files must not be used: unmated reads would shift the pairing between _1 and _2.) Resumable: runs with a .done marker are skipped.
# Failed runs are listed in logs/05_download_failed.tsv.
# Requires: sra-toolkit (fasterq-dump), run from the repository root.

set -uo pipefail

# refuse to run on the login node (lab rule: no computation on iv12; use salloc --partition bio)
if [[ "$(hostname -s)" == iv* ]]; then
  echo "ERROR: you are on the login node ($(hostname -s)). Request a compute node first:" >&2
  echo "  salloc --time=24:00:00 --cpus-per-task=4 --mem=8G --partition bio" >&2
  exit 1
fi

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

  if fasterq-dump --split-3 --skip-technical --threads "$THREADS" --temp "$RAW/tmp" \
       -O "$RAW" "$run" < /dev/null > /dev/null 2>> logs/05_download_errors.log \
     && ls "$RAW/${run}"*.fastq > /dev/null 2>&1; then
    $GZ "$RAW/${run}"*.fastq
    if [[ -f "$RAW/${run}_2.fastq.gz" ]]; then
      n1=$(( $(zcat "$RAW/${run}_1.fastq.gz" | wc -l) / 4 ))
      n2=$(( $(zcat "$RAW/${run}_2.fastq.gz" | wc -l) / 4 ))
      if [[ "$n1" != "$n2" ]]; then
        echo -e "$run\tunpaired_${n1}_${n2}" >> logs/05_download_failed.tsv
        echo "  R1/R2 read numbers differ ($n1 vs $n2)"; continue
      fi
    fi
    touch "$RAW/${run}.done"
  else
    echo -e "$run\tdownload_failed" >> logs/05_download_failed.tsv
    echo "  download failed"
  fi
done

rm -rf "$RAW/tmp"
echo "Done: $(ls "$RAW"/*.done 2> /dev/null | wc -l) of $total runs downloaded in $RAW"
