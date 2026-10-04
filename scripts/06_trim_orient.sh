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
# Safety checks: runs whose R1/R2 files differ in read number are skipped (action
# "unpaired_input"); if any cutadapt call fails, partial outputs are deleted and the run is
# recorded as "cutadapt_failed". Neither is written as a successful run.
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
  : > "$log"

  if [[ ! -f "$in1" || ! -f "$in2" ]]; then
    echo -e "$run\t$batch\tNA\tmissing_paired_input\t0\t0\tNA\tNA" >> "$MANIFEST"
    echo "$run: paired input not found, skipped"; continue
  fi
  n_in=$(nreads "$in1"); n_in2=$(nreads "$in2")
  if [[ "$n_in" != "$n_in2" ]]; then
    echo -e "$run\t$batch\tNA\tunpaired_input\t$n_in\t0\tNA\tNA" >> "$MANIFEST"
    echo "$run: R1/R2 read numbers differ ($n_in vs $n_in2), skipped"; continue
  fi

  # R1_reverse: the reverse primer side is in R1, so swap
  if [[ "$orient" == "R1_reverse" ]]; then a1="$in2"; a2="$in1"; else a1="$in1"; a2="$in2"; fi

  ok=1
  cut5() {   # 5' primers, pairs without both primers discarded: in1 in2 out1 out2
    cutadapt -j "$THREADS" -g "$F" -G "$R" --discard-untrimmed $extra \
      -o "$3" -p "$4" "$1" "$2" >> "$log" 2>&1 || ok=0
  }
  trim3() {  # read-through into the opposite primer: in1 in2 out1 out2
    cutadapt -j "$THREADS" -a "$(rc "$R")" -A "$(rc "$F")" --minimum-length 50 \
      -o "$3" -p "$4" "$1" "$2" >> "$log" 2>&1 || ok=0
  }
  write_row() {  # set action r1 r2
    echo -e "$run\t$batch\t$1\t$2\t$n_in\t$(nreads "$3")\t$3\t$4" >> "$MANIFEST"
  }
  fail() {   # remove partial outputs and record the failure
    rm -rf "$dir/tmp_${run}" "$dir/${run}"_*.fastq.gz
    echo -e "$run\t$batch\tNA\tcutadapt_failed\t$n_in\t0\tNA\tNA" >> "$MANIFEST"
    echo "$run ($batch): cutadapt FAILED, see $log"
  }

  if [[ "$fstat" == "removed" && "$rstat" == "removed" ]]; then
    o1="$dir/${run}_1.fastq.gz"; o2="$dir/${run}_2.fastq.gz"
    cp "$a1" "$o1"; cp "$a2" "$o2"
    write_row main "copy$([[ "$orient" == "R1_reverse" ]] && echo "_swapped")" "$o1" "$o2"

  elif [[ "$orient" == "mixed" ]]; then
    t="$dir/tmp_${run}"; mkdir -p "$t"
    # set A: forward primer in R1; pairs not found are kept for set B
    cutadapt -j "$THREADS" -g "$F" -G "$R" $extra \
      --untrimmed-output "$t/u1.fq.gz" --untrimmed-paired-output "$t/u2.fq.gz" \
      -o "$t/a1.fq.gz" -p "$t/a2.fq.gz" "$in1" "$in2" >> "$log" 2>&1 || ok=0
    # set B: forward primer in R2 (swap the untrimmed pairs)
    [[ $ok == 1 ]] && cut5 "$t/u2.fq.gz" "$t/u1.fq.gz" "$t/b1.fq.gz" "$t/b2.fq.gz"
    [[ $ok == 1 ]] && trim3 "$t/a1.fq.gz" "$t/a2.fq.gz" "$dir/${run}_A_1.fastq.gz" "$dir/${run}_A_2.fastq.gz"
    [[ $ok == 1 ]] && trim3 "$t/b1.fq.gz" "$t/b2.fq.gz" "$dir/${run}_B_1.fastq.gz" "$dir/${run}_B_2.fastq.gz"
    if [[ $ok == 1 ]]; then
      rm -rf "$t"
      write_row A cutadapt_mixed "$dir/${run}_A_1.fastq.gz" "$dir/${run}_A_2.fastq.gz"
      write_row B cutadapt_mixed "$dir/${run}_B_1.fastq.gz" "$dir/${run}_B_2.fastq.gz"
    else
      fail; continue
    fi

  else
    t="$dir/tmp_${run}"; mkdir -p "$t"
    cut5 "$a1" "$a2" "$t/a1.fq.gz" "$t/a2.fq.gz"
    [[ $ok == 1 ]] && trim3 "$t/a1.fq.gz" "$t/a2.fq.gz" "$dir/${run}_1.fastq.gz" "$dir/${run}_2.fastq.gz"
    if [[ $ok == 1 ]]; then
      rm -rf "$t"
      write_row main "cutadapt$([[ "$orient" == "R1_reverse" ]] && echo "_swapped")" \
        "$dir/${run}_1.fastq.gz" "$dir/${run}_2.fastq.gz"
    else
      fail; continue
    fi
  fi
  echo "$run ($batch): done"
done

echo "Done: $MANIFEST"
awk -F'\t' 'NR > 1 && ($4 == "cutadapt_failed" || $4 == "unpaired_input" || $4 == "missing_paired_input")' "$MANIFEST" \
  | cut -f1,2,4 | sed 's/^/  NOT PROCESSED: /'
