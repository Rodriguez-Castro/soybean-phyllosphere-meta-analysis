#!/usr/bin/env bash
# Locate the reads of each SRA run on the E. coli 16S rRNA gene (J01859, standard numbering)
# with BLAST. Unlike step 01, this works whether primers were kept or removed before
# submission, with spacers/barcodes before the primers, and in any read orientation.
#
# Usage:
#   bash scripts/01b_align_reads.sh <run_list.tsv> <output.tsv> [n_reads]
#
# Output, one row per run and read file:
#   frac_16S     fraction of reads aligning to 16S (near 0: ITS, RNA-Seq, shotgun...)
#   frac_plus    fraction of aligned reads on the 16S forward strand (~0.5: mixed orientation)
#   start_plus   median E. coli position where forward-strand alignments start
#                (e.g. 515 = 515F still present; ~534 = 515F removed)
#   start_minus  median E. coli position where reverse-strand alignments start
#                (e.g. 806 = 806R still present; ~786 = 806R removed)
#   qstart_*     median read position where the alignment starts (>1: spacer/barcode before it)
#
# Resumable: runs already in <output.tsv> are skipped.
# Requires: sra-toolkit (fastq-dump), blast+ (blastn, makeblastdb), run from the repository root.

set -uo pipefail

# refuse to run on the login node (lab rule: no computation on iv12; use salloc --partition bio)
if [[ "$(hostname -s)" == iv* ]]; then
  echo "ERROR: you are on the login node ($(hostname -s)). Request a compute node first:" >&2
  echo "  salloc --time=24:00:00 --cpus-per-task=4 --mem=8G --partition bio" >&2
  exit 1
fi

RUNLIST=${1:?run list required}
OUT=${2:?output file required}
NREADS=${3:-200}
REF=ref/ecoli_16S_J01859.fasta

mkdir -p ref "$(dirname "$OUT")"
if [[ ! -s "$REF" ]]; then
  curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J01859&rettype=fasta&retmode=text" -o "$REF"
fi
grep -q '^>' "$REF" || { echo "ERROR: could not download the E. coli 16S reference"; exit 1; }
[[ -s "$REF.nsq" ]] || makeblastdb -in "$REF" -dbtype nucl > /dev/null

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

med() {  # median of numbers on stdin, NA if none
  sort -n | awk '{ a[NR] = $1 }
    END { if (NR == 0) print "NA"; else if (NR % 2) print a[(NR + 1) / 2]; else print int((a[NR / 2] + a[NR / 2 + 1]) / 2) }'
}

[[ -s "$OUT" ]] || echo -e "run\tbioproject\tread\tn_reads\tn_16S\tfrac_16S\tfrac_plus\tstart_plus\tstart_minus\tqstart_plus\tqstart_minus" > "$OUT"

total=$(($(wc -l < "$RUNLIST") - 1)); i=0

tail -n +2 "$RUNLIST" | while IFS=$'\t' read -r run bioproject _; do
  i=$((i + 1))
  if grep -q -P "^${run}\t" "$OUT"; then continue; fi
  echo "[$i/$total] $run ($bioproject)"

  if ! fastq-dump -X "$NREADS" --split-files -O "$TMP" "$run" < /dev/null > /dev/null 2>&1; then
    echo -e "$run\t$bioproject\tNA\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA" >> "$OUT"
    echo "  download failed"
    continue
  fi

  for f in "$TMP/${run}"*.fastq; do
    case "$f" in
      *_1.fastq) rd=R1 ;;
      *_2.fastq) rd=R2 ;;
      *)         rd=SE ;;
    esac
    awk 'NR % 4 == 1 { print ">r" NR } NR % 4 == 2 { print }' "$f" > "$TMP/q.fa"
    n=$(grep -c '^>' "$TMP/q.fa")

    # strand, E. coli start of the alignment, read position where it starts
    blastn -task blastn -query "$TMP/q.fa" -db "$REF" -evalue 1e-10 \
           -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid qstart qend sstart send" 2> /dev/null \
      | awk '!seen[$1]++ { if ($4 <= $5) print "+", $4, $2; else print "-", $4, $2 }' > "$TMP/hits.txt"

    nh=$(wc -l < "$TMP/hits.txt")
    np=$(grep -c '^+' "$TMP/hits.txt")
    frac=$(awk -v a="$nh" -v b="$n" 'BEGIN { printf "%.3f", (b ? a / b : 0) }')
    fplus=$(awk -v a="$np" -v b="$nh" 'BEGIN { if (b) printf "%.3f", a / b; else print "NA" }')
    sp=$(awk '$1 == "+" { print $2 }' "$TMP/hits.txt" | med)
    sm=$(awk '$1 == "-" { print $2 }' "$TMP/hits.txt" | med)
    qp=$(awk '$1 == "+" { print $3 }' "$TMP/hits.txt" | med)
    qm=$(awk '$1 == "-" { print $3 }' "$TMP/hits.txt" | med)

    echo -e "$run\t$bioproject\t$rd\t$n\t$nh\t$frac\t$fplus\t$sp\t$sm\t$qp\t$qm" >> "$OUT"
  done

  rm -f "$TMP/${run}"*.fastq
done

echo "Done: $OUT"
