#!/usr/bin/env bash
# Step 06: remove primers (where present) and orient all reads so that file _1 starts at the
# forward primer site and file _2 at the reverse primer site. Driven by metadata/sample_sheet.tsv.
#
# Usage:
#   bash scripts/06_trim_orient.sh <sample_sheet.tsv> <raw_dir> <out_dir> [threads]
#
# Per run, from the sample sheet:
#   primers present  -> cutadapt: 5' primers required (--discard-untrimmed, as in the predoc;
#                       non-anchored, so spacers/barcodes before the primer are removed too),
#                       then reverse-complemented primers removed from 3' ends (read-through)
#   primers removed  -> reads copied unchanged
#   R1_reverse       -> R1/R2 swapped
#   mixed            -> two cutadapt passes; the pairs found in each orientation are written as
#                       separate sets (<run>_A and <run>_B) so that DADA2 can learn their error
#                       models separately (set B forward reads come from the original R2)
# Extra cutadapt options per batch come from config/dada2_params.tsv (column cutadapt_extra).
#
# Output: <out_dir>/<batch>/<run>[_A|_B]_1.fastq.gz / _2.fastq.gz and
#         results/06_trim/trim_manifest.tsv (reads in/out per run). Resumable.
# Requires: cutadapt, run from the repository root.

set -uo pipefail

SHEET=${1:?sample sheet required}
RAW=${2:?raw reads directory required}
OUT=${3:?output directory required}
THREADS=${4:-4}
MANIFEST=results/06_trim/trim_manifest.tsv
mkdir -p "$OUT" results/06_trim logs

rc() { echo "$1" | rev | tr 'ACGTRYKMBDHVN' 'TGCAYRMKVHDBN'; }
nreads() { if [[ -f "$1" ]]; then echo $(( $(zcat "$1" | wc -l) / 4 )); else echo 0; fi; }

# extra cutadapt options per batch
# (tabs are whitespace for "read", so empty fields would collapse: use "|" as separator)
declare -A EXTRA
while IFS='|' read -r batch extra; do
  EXTRA[$batch]="$extra"
done < <(awk -F'\t' 'NR > 1 { print $1 "|" $3 }' config/dada2_params.tsv)

[[ -s "$MANIFEST" ]] || echo -e "run\tbatch\tset\taction\treads_in\treads_out\tr1\tr2" > "$MANIFEST"

# fixed column order regardless of the sample sheet layout
awk -F'\t' -v OFS='|' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
  { print $h["run"], $h["batch"], $h["orientation"], $h["fwd_seq"], $h["fwd_status"], $h["rev_seq"], $h["rev_status"] }' "$SHEET" \
| while IFS='|' read -r run batch orient F fstat R rstat; do

  grep -q -P "^${run}\t" "$MANIFEST" && continue
  dir="$OUT/$batch"; mkdir -p "$dir"
  in1="$RAW/${run}_1.fastq.gz"; in2="$RAW/${run}_2.fastq.gz"
  extra=${EXTRA[$batch]:-}
  log="logs/06_cutadapt_${run}.log"

  if [[ ! -f "$in1" || ! -f "$in2" ]]; then
    echo -e "$run\t$batch\tNA\tmissing_paired_input\t0\t0\tNA\tNA" >> "$MANIFEST"
    echo "$run: paired input not found, skipped"; continue
  fi
  n_in=$(nreads "$in1")

  # R1_reverse: the reverse primer side is in R1, so swap
  if [[ "$orient" == "R1_reverse" ]]; then a1="$in2"; a2="$in1"; else a1="$in1"; a2="$in2"; fi

  write_row() {  # set action r1 r2
    echo -e "$run\t$batch\t$1\t$2\t$n_in\t$(nreads "$3")\t$3\t$4" >> "$MANIFEST"
  }
  trim3() {      # remove read-through into the opposite primer: in1 in2 out1 out2
    cutadapt -j "$THREADS" -a "$(rc "$R")" -A "$(rc "$F")" --minimum-length 50 \
      -o "$3" -p "$4" "$1" "$2" >> "$log" 2>&1
  }

  if [[ "$fstat" == "removed" && "$rstat" == "removed" ]]; then
    o1="$dir/${run}_1.fastq.gz"; o2="$dir/${run}_2.fastq.gz"
    cp "$a1" "$o1"; cp "$a2" "$o2"
    write_row main "copy$([[ "$orient" == "R1_reverse" ]] && echo "_swapped")" "$o1" "$o2"

  elif [[ "$orient" == "mixed" ]]; then
    t="$dir/tmp_${run}"; mkdir -p "$t"
    # set A: forward primer in R1
    cutadapt -j "$THREADS" -g "$F" -G "$R" $extra \
      --untrimmed-output "$t/u1.fq.gz" --untrimmed-paired-output "$t/u2.fq.gz" \
      -o "$t/a1.fq.gz" -p "$t/a2.fq.gz" "$in1" "$in2" > "$log" 2>&1
    # set B: forward primer in R2 (swap the untrimmed pairs)
    cutadapt -j "$THREADS" -g "$F" -G "$R" --discard-untrimmed $extra \
      -o "$t/b1.fq.gz" -p "$t/b2.fq.gz" "$t/u2.fq.gz" "$t/u1.fq.gz" >> "$log" 2>&1
    trim3 "$t/a1.fq.gz" "$t/a2.fq.gz" "$dir/${run}_A_1.fastq.gz" "$dir/${run}_A_2.fastq.gz"
    trim3 "$t/b1.fq.gz" "$t/b2.fq.gz" "$dir/${run}_B_1.fastq.gz" "$dir/${run}_B_2.fastq.gz"
    rm -rf "$t"
    write_row A cutadapt_mixed "$dir/${run}_A_1.fastq.gz" "$dir/${run}_A_2.fastq.gz"
    write_row B cutadapt_mixed "$dir/${run}_B_1.fastq.gz" "$dir/${run}_B_2.fastq.gz"

  else
    t="$dir/tmp_${run}"; mkdir -p "$t"
    cutadapt -j "$THREADS" -g "$F" -G "$R" --discard-untrimmed $extra \
      -o "$t/a1.fq.gz" -p "$t/a2.fq.gz" "$a1" "$a2" > "$log" 2>&1
    trim3 "$t/a1.fq.gz" "$t/a2.fq.gz" "$dir/${run}_1.fastq.gz" "$dir/${run}_2.fastq.gz"
    rm -rf "$t"
    write_row main "cutadapt$([[ "$orient" == "R1_reverse" ]] && echo "_swapped")" \
      "$dir/${run}_1.fastq.gz" "$dir/${run}_2.fastq.gz"
  fi
  echo "$run ($batch): done"
done

echo "Done: $MANIFEST"
