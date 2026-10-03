# Soybean phyllosphere bacterial meta-analysis

Reproducible pipeline to reprocess public 16S rRNA amplicon data from the phyllosphere of *Glycine* spp. and integrate studies at genus level.

## Strategy

1. **Verify every run from its reads**, not from the published methods: detect which primers each run actually contains, and which sequencing run (flow cell) it comes from.
2. **One DADA2 batch per BioProject** (see `config/batches.tsv`). If a BioProject contains several sequencing runs (flow cells), error models are learned per sequencing run within that batch.
3. **Taxonomy assigned per batch** with identical settings (same SILVA release, `minBoot`, seed).
4. **Integration at genus level** across all batches and regions, with study/sequencing run kept as covariates.

## Repository structure

```
soybean-phyllosphere-meta-analysis/
├── README.md
├── config/
│   ├── batches.tsv              # one DADA2 batch per BioProject, reported vs verified primers
│   ├── primers.tsv              # primer catalogue (IUPAC) used for verification and trimming
│   └── selection_overrides.tsv  # documented exceptions to the original selection
├── metadata/
│   ├── runinfo/                 # SRA RunInfo CSVs, one per BioProject (input)
│   ├── run_list.tsv             # combined run list (step 00)
│   ├── phyllosphere_runs.csv    # runs selected for the meta-analysis (input of step 03)
│   └── sample_sheet.tsv         # final runs for DADA2 (step 04)
├── scripts/
│   ├── 00_build_run_list.py
│   ├── 01_verify_runs.sh
│   ├── 01b_align_reads.sh
│   ├── 02_summarise_verification.R
│   ├── 02b_summarise_alignment.R
│   ├── 03_cross_selection.R
│   └── 04_build_sample_sheet.R
├── results/
│   ├── 01_verification/
│   ├── 01b_alignment/
│   ├── 02_summary/
│   ├── 02b_summary/
│   ├── 03_selection/
│   └── 04_sample_sheet/
└── logs/
```

## Requirements (hpc-bio, Compute Canada software stack)

No computation is allowed on the login node (`iv12`). All steps run on the compute node `cv3401`, requested with `salloc` from inside a `tmux` session:

```bash
ssh rodl4348@hpc-bio.ccs.usherbrooke.ca
newgrp def-ilafores
tmux new -s verify
salloc --time=24:00:00 --cpus-per-task=1 --mem=2G --partition bio   # prompt changes to cv3401
module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9                  # steps 00-01b
module load blast+                                                  # step 01b (check version with: module spider blast+)
module load StdEnv/2023 r/4.4.0                                     # step 02
```

`cv3401` has internet access (needed by `fastq-dump`) and sees the home directory, but not `/jbod2`; keep the repository in `$HOME`.

Python 3 (standard library only) for step 00; base R only for step 02.

## Steps

### 00 — Build the run list

Put one SRA RunInfo CSV per BioProject in `metadata/runinfo/`, then:

```bash
python3 scripts/00_build_run_list.py metadata/runinfo metadata/run_list.tsv
# optional: keep only selected runs (one SRR per line)
python3 scripts/00_build_run_list.py metadata/runinfo metadata/run_list.tsv selected_runs.txt
```

### 01 — Verify primers and sequencing runs from the reads

Downloads the first 1,000 reads of each run, records the instrument / run / flow cell from the read headers, and counts how many reads start (first 30 bp) with each primer in `config/primers.tsv`. The script can be stopped and relaunched: runs already in the output are skipped.

```bash
tmux new -s verify
bash scripts/01_verify_runs.sh metadata/run_list.tsv results/01_verification/primer_counts.tsv 1000
```

Output columns: `run, bioproject, instrument, seq_run, flowcell, read, n_reads, primer, count`.

A run with no primer detected may have had primers removed before submission; check it manually.

### 01b — Locate reads on the E. coli 16S gene with BLAST

Step 01 only detects primers within the first 30 bp of the reads. It misses runs whose primers were removed before submission or are preceded by long spacers. Step 01b aligns the first 200 reads of each run to the *E. coli* 16S rRNA gene (J01859, standard numbering) and reports where the alignments start:

```bash
bash scripts/01b_align_reads.sh metadata/run_list.tsv results/01b_alignment/alignment.tsv 200
```

| Column | Meaning |
|---|---|
| `frac_16S` | fraction of reads aligning to 16S. Fungal ITS reads align partially (the end of the fungal 18S gene resembles the end of 16S), always near *E. coli* positions 1491/1533 |
| `frac_plus` | fraction of aligned reads on the forward strand; ~0.5 means mixed orientation |
| `start_plus` / `start_minus` | median *E. coli* position where forward / reverse alignments start, i.e. the amplicon ends (e.g. 515 = 515F kept, ~534 = 515F removed) |
| `qstart_plus` / `qstart_minus` | median read position where the alignment starts; > 1 means a spacer or barcode before the primer |

### 02 — Classify runs and summarise per BioProject

```bash
Rscript scripts/02_summarise_verification.R results/01_verification/primer_counts.tsv metadata/run_list.tsv results/02_summary
```

Each run is classified from the fraction of reads starting with each primer:

| Class | Rule |
|---|---|
| `bacterial_16S` | ≥ 80% of reads start with the assigned 16S primer pair |
| `16S_partial` | 20–80%; check manually |
| `fungal_ITS` | ≥ 50% of reads start with an ITS primer; excluded |
| `no_primer_detected` | < 20%; primers possibly removed before submission |
| `download_failed` / `no_reads` | relaunch step 01 for these runs |

Nested primers that match the same reads are resolved to the outer primer (338F over 341F, 806R over 785R). Orientation is `mixed` when both primers start ≥ 10% of R1 reads. The sequencing run is `instrument:run:flowcell` from the original read names.

Outputs: `run_summary.tsv` (one row per run), `bioproject_summary.tsv` (classes, primer pairs and sequencing runs per BioProject) and `bacterial_runs.tsv`.

### 02b — Interpret the alignment positions (final run classification)

```bash
Rscript scripts/02b_summarise_alignment.R results/01b_alignment/alignment.tsv metadata/run_list.tsv results/02b_summary results/02_summary/run_summary.tsv
```

Each run is classified as `16S`, `fungal_ITS_like` (all read ends beyond position 1450), `inconsistent` (one read end beyond 1450 and the other within 16S; excluded), `non_16S` (< 50% of reads align to 16S), `16S_not_amplicon` (library strategy is not AMPLICON, e.g. RNA-Seq) or `download_failed`. For 16S runs, the read ends are matched (±3 bp) to the primer binding sites on *E. coli* 16S:

- **present**: reads start at the primer's first base (primers must be trimmed);
- **removed**: reads start right after the primer (nothing to trim).

Primers with the same 3' end (338F/341F, 785R/806R) cannot be told apart once removed; they amplify the same region. The amplicon span gives the hypervariable regions (V1–V9, *E. coli* numbering) covered at least 80%; this tolerates reads trimmed a few bases beyond the primer (e.g. PRJNA1280517, reverse reads ending at ~1170 instead of 1175). `bacterial_runs.tsv` from this step is the input for DADA2.

### 03 — Cross the phyllosphere selection with the verified classification

`metadata/phyllosphere_runs.csv` lists the runs selected for the meta-analysis (column `run`; optional `bio_project`, `dada2_group`). It is copied from the lab storage, which is only mounted on the login node:

```bash
cp /jbod2/def-ilafores/analysis/meta_analysis_soybean/data/filtered/phyllosphere_runs.csv metadata/   # on iv12
Rscript scripts/03_cross_selection.R metadata/phyllosphere_runs.csv results/02b_summary/run_alignment_summary.tsv results/03_selection   # on cv3401
```

Outputs: `selection_check.tsv` (every selected run with its verified class), `selection_summary.tsv` (kept and lost runs per BioProject and previous group, by reason; `not_verified` = run outside the 10 verified BioProjects) and `selected_bacterial_runs.tsv` (input for DADA2).

### 04 — Build the final sample sheet

```bash
Rscript scripts/04_build_sample_sheet.R
```

Includes every selected run verified as bacterial 16S whose previous group is not `excluded_*`, plus the exceptions listed in `config/selection_overrides.tsv`. For each run, `metadata/sample_sheet.tsv` gives the DADA2 batch, sequencing run (error models are learned per sequencing run), region, leaf compartment, orientation, primer names and sequences, and whether each primer is `present` (to trim) or `removed` (nothing to trim). Runs whose read ends did not match a known primer site take the majority status of their batch and are flagged in `status_flag`. `results/04_sample_sheet/batch_summary.tsv` summarises each batch.

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-10-01 | DADA2 and taxonomy run separately per BioProject (10 batches); studies merged only at genus level. Hypervariable region is kept as a descriptive variable and covariate, not as a processing group | Each BioProject is, in practice, a separate sequencing run and error models are run-specific; ASVs from different primers or truncation lengths cannot be merged (e.g., 515F/806R vs 520F/799R in V4) |
| 2026-10-02 | PRJNA601979 relabelled from V2–V3 to V4 (verified) | The paper states that the "V2–V3" region was amplified with 520F/799R. These primers bind around *E. coli* positions 520 and 799, downstream of V3 (≈433–497) and upstream of V5 (≈822–879), so the amplicon can only span V4 (≈576–682). The literature consistently describes 520F/799R as a V4 primer pair. Confirmed: reads start with 520F (AGCAGCCGCGGTAAT) and 799R (CMGGGTATCTAATCCKGTT, reverse complement of 799F). |
| 2026-10-01 | Primers assigned from the reads, not from the papers | PRJNA1092852: paper reports 338F/806R, but bacterial libraries (suffix `.1b`, 134 runs) contain 799F/1193R (V5–V7) in ~99% of reads, in mixed orientation (~50% each). Libraries with suffix `.1` (129 runs) are fungal ITS1F/ITS2 and are excluded. Authors contacted for confirmation. |
| 2026-10-02 | Primer verification extended with BLAST against *E. coli* 16S (step 01b) | Step 01 left 2,174 of 3,700 runs as `no_primer_detected` and 1,076 as `16S_partial`: several BioProjects removed primers before submission (PRJNA1280517, PRJNA544311) or have spacers longer than the 30 bp window (PRJNA603147) |
| 2026-10-02 | PRJNA661376 flagged for exclusion | Its 27 runs are RNA-Seq, not 16S amplicons; the BioProject holding the 93 amplicon samples must be identified |
| 2026-10-02 | PRJNA603199 and PRJNA987554 flagged for review | Reads of sampled runs contain fungal ITS sequences (ITS1F/ITS4 sites; 18S end and 5.8S start) |
| 2026-10-02 | Final run classification based on BLAST positions (step 02b), not on primer matches alone | Within a single BioProject, runs can differ: PRJNA987554 contains runs with 338F/806R present, runs with primers removed, and fungal ITS runs |
| 2026-10-02 | PRJNA603199 excluded | 504 of 613 runs are fungal ITS and 95 are not 16S; the 14 remaining runs are inconsistent (R1 at the 18S end, R2 at 806R). It is the fungal dataset of the study, not bacterial |
| 2026-10-02 | Group_E_V5V6 BioProject corrected from PRJNA661376 to PRJNA662376; PRJNA390118 added | Step 03 showed the selection uses PRJNA662376 (the summary table had a typo; PRJNA661376 is RNA-Seq) and includes PRJNA390118, absent from the summary table. Both added to the verification |
| 2026-10-02 | Selected ITS runs removed from Groups A, C and D | Step 03: 36 of 150 selected runs of PRJNA603147, 125 of 126 of PRJNA603199, 37 of 78 of PRJNA601979 and 37 of 74 of PRJNA862265 are fungal ITS. They did not reach the previous ASV tables (no 16S primers), but inflated the reported sample numbers |
| 2026-10-02 | 86 runs of PRJNA1092852 recovered | All 172 selected runs had been labelled `excluded_ITS`, but 86 are bacterial 16S (V5–V7, 799F/1193R) from leaves: 42 `L` and 44 `LE` (leaf endosphere) |
| 2026-10-02 | Phyllosphere includes leaf surface and endosphere; each study labelled `epi`, `endo` or `endo+epi` (column `leaf_compartment` in `config/batches.tsv`) | Several studies do not separate epiphytes from endophytes; the label allows covariate and sensitivity analyses |
| 2026-10-02 | PRJNA390118 relabelled from V4–V5 to V4; its 72 selected ITS runs removed | Reads start right after 515F and end at ~781–785 (515F/806R removed); no 907R signal. Libraries are named `S-16S-*` and `S-ITS-*`; all 72 selected 16S runs are leaf (`L`) |
| 2026-10-02 | Verified selection: 983 bacterial phyllosphere runs | 897 of the 1,219 previously selected runs in included groups are 16S, plus 86 recovered runs of PRJNA1092852 |
