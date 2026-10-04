#!/usr/bin/env bash
# Step 05b: check that R1 and R2 of every downloaded run contain the same number of reads.
#
# Usage:
#   bash scripts/05b_check_pairs.sh <sample_sheet.tsv> <raw_dir> [--reset]
#
# Writes results/05_download/pair_check.tsv. With --reset, the files and .done marker of
# mismatched runs are deleted, so that 05_download_reads.sh downloads them again.

set -uo pipefail

SHEET=${1:?sample sheet required}
RAW=${2:?raw reads directory required}
RESET=${3:-}
OUT=results/05_download/pair_check.tsv
mkdir -p results/05_download

echo -e "run\tr1_reads\tr2_reads\tstatus" > "$OUT"
tail -n +2 "$SHEET" | cut -f1 | while read -r run; do
  f1="$RAW/${run}_1.fastq.gz"; f2="$RAW/${run}_2.fastq.gz"
  n1=NA; n2=NA
  [[ -f "$f1" ]] && n1=$(( $(zcat "$f1" | wc -l) / 4 ))
  [[ -f "$f2" ]] && n2=$(( $(zcat "$f2" | wc -l) / 4 ))
  if [[ "$n1" != NA && "$n1" == "$n2" ]]; then st=ok; else st=mismatch; fi
  echo -e "$run\t$n1\t$n2\t$st" >> "$OUT"
  if [[ "$st" == mismatch && "$RESET" == "--reset" ]]; then
    rm -f "$RAW/${run}.done" "$RAW/${run}"*.fastq.gz
  fi
done

echo "$(awk -F'\t' '$4 == "mismatch"' "$OUT" | wc -l) of $(($(wc -l < "$OUT") - 1)) runs with R1/R2 mismatch (see $OUT)"
[[ "$RESET" == "--reset" ]] && echo "Mismatched runs deleted: run 05_download_reads.sh again to re-download them"
