#!/usr/bin/env bash
# Verify primers and sequencing run (instrument / run / flow cell) of each SRA run
# directly from its reads.
#
# Usage:
#   bash scripts/01_verify_runs.sh <run_list.tsv> <output.tsv> [n_reads]
#
# <run_list.tsv>: output of 00_build_run_list.py (columns: run, bioproject, ...)
# The script is resumable: runs already present in <output.tsv> are skipped.
# Requires: sra-toolkit (fastq-dump), run from the repository root.

set -uo pipefail

RUNLIST=${1:?run list required}
OUT=${2:?output file required}
NREADS=${3:-1000}
PRIMERS=config/primers.tsv

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# IUPAC degenerate bases -> regular expression
iupac2re() {
  echo "$1" | sed -e 's/R/[AG]/g; s/Y/[CT]/g; s/S/[CG]/g; s/W/[AT]/g; s/K/[GT]/g; s/M/[AC]/g' \
                  -e 's/B/[CGT]/g; s/D/[AGT]/g; s/H/[ACT]/g; s/V/[ACG]/g; s/N/[ACGT]/g'
}

# load primer catalogue
PNAME=(); PRE=()
while IFS=$'\t' read -r name seq _; do
  [[ -z "$name" || "$name" == "name" || "$name" == \#* ]] && continue
  PNAME+=("$name")
  PRE+=("$(iupac2re "$seq")")
done < "$PRIMERS"
echo "Loaded ${#PNAME[@]} primers from $PRIMERS"

mkdir -p "$(dirname "$OUT")"
[[ -s "$OUT" ]] || echo -e "run\tbioproject\tinstrument\tseq_run\tflowcell\tread\tn_reads\tprimer\tcount" > "$OUT"

total=$(($(wc -l < "$RUNLIST") - 1)); i=0

tail -n +2 "$RUNLIST" | while IFS=$'\t' read -r run bioproject _; do
  i=$((i + 1))
  if grep -q -P "^${run}\t" "$OUT"; then continue; fi
  echo "[$i/$total] $run ($bioproject)"

  if ! fastq-dump -X "$NREADS" --split-files -O "$TMP" "$run" < /dev/null > /dev/null 2>&1; then
    echo -e "$run\t$bioproject\tNA\tNA\tNA\tNA\t0\tDOWNLOAD_FAILED\t0" >> "$OUT"
    echo "  download failed"
    continue
  fi

  # sequencing run from the original read name, e.g. M03152:573:000000000-J8KJ7:1:1101:15366:1768
  first=$(ls "$TMP/${run}"*.fastq | head -1)
  hdr=$(head -1 "$first" | awk '{print $2}')
  if [[ "$hdr" == *:*:*:* ]]; then
    IFS=: read -r inst srun fc _ <<< "$hdr"
  else
    inst=NA; srun=NA; fc=NA   # original read names not kept by SRA
  fi

  for f in "$TMP/${run}"*.fastq; do
    case "$f" in
      *_1.fastq) rd=R1 ;;
      *_2.fastq) rd=R2 ;;
      *)         rd=SE ;;
    esac
    awk 'NR % 4 == 2 { print substr($0, 1, 30) }' "$f" > "$TMP/starts.txt"
    n=$(wc -l < "$TMP/starts.txt")
    for k in "${!PNAME[@]}"; do
      c=$(grep -c -E "${PRE[$k]}" "$TMP/starts.txt")
      echo -e "$run\t$bioproject\t$inst\t$srun\t$fc\t$rd\t$n\t${PNAME[$k]}\t$c"
    done >> "$OUT"
  done

  rm -f "$TMP/${run}"*.fastq
done

echo "Done: $OUT"
