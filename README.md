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
│   ├── dada2_params.tsv         # per-batch processing parameters (cutadapt options, error model, truncLen)
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
│   ├── 04_build_sample_sheet.R
│   ├── 05_download_reads.sh
│   ├── 05b_check_pairs.sh
│   ├── 06_trim_orient.sh
│   └── 07_quality_profiles.R
├── results/
│   ├── 01_verification/
│   ├── 01b_alignment/
│   ├── 02_summary/
│   ├── 02b_summary/
│   ├── 03_selection/
│   ├── 04_sample_sheet/
│   ├── 06_trim/
│   └── 07_quality/
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

### 05 — Download the complete reads

```bash
bash scripts/05_download_reads.sh metadata/sample_sheet.tsv <raw_dir> 4
```

`fasterq-dump --split-3` for every run of the sample sheet, compressed to `<raw_dir>/<run>_1.fastq.gz` / `_2.fastq.gz` (reads without mate go to `<run>.fastq.gz` and are not used); runs whose `_1` and `_2` differ in read number are not marked as done. Resumable (`.done` markers); failures listed in `logs/05_download_failed.tsv`.

### 05b — Check read pairing

```bash
bash scripts/05b_check_pairs.sh metadata/sample_sheet.tsv <raw_dir> --reset
bash scripts/05_download_reads.sh metadata/sample_sheet.tsv <raw_dir> 4
```

Counts the reads of `_1` and `_2` for every run (`results/05_download/pair_check.tsv`). With `--reset`, mismatched runs are deleted so that step 05 downloads them again.

### 06 — Remove primers and orient reads

```bash
module load cutadapt   # check the version with: module spider cutadapt
bash scripts/06_trim_orient.sh metadata/sample_sheet.tsv <raw_dir> <trimmed_dir> 4
```

After this step, in every run file `_1` starts at the forward primer site and `_2` at the reverse primer site:

| Case (from the sample sheet) | Action |
|---|---|
| Primers present | cutadapt with `--discard-untrimmed` (non-anchored 5' primers, so spacers are removed too), then removal of reverse-complemented primers from 3' ends (read-through) |
| Primers removed | reads copied unchanged |
| `R1_reverse` | R1 and R2 swapped |
| `mixed` (B05) | two cutadapt passes; pairs found in each orientation are written as sets `<run>_A` and `<run>_B`, which get separate error models in DADA2 |

Extra cutadapt options per batch come from `config/dada2_params.tsv` (e.g. `--pair-filter=first` for B07, as in the predoc). `results/06_trim/trim_manifest.tsv` records reads in/out per run. Runs with unpaired input files or a failed cutadapt call are not written as successful (`unpaired_input`, `cutadapt_failed`) and are listed at the end of the run.

cutadapt v5.2 is installed in a dedicated virtual environment in `$HOME` (visible from `cv3401`):

```bash
module load StdEnv/2023 python/3.11
python -m venv ~/envs/cutadapt && source ~/envs/cutadapt/bin/activate && pip install cutadapt==5.2
```

### 07 — Quality profiles and truncation lengths

```bash
Rscript scripts/07_quality_profiles.R results/06_trim/trim_manifest.tsv results/07_quality 10 2000
```

For each batch (and each orientation set of B05), samples 2,000 reads from up to 10 runs and proposes `truncLen` values that (1) do not exceed the 5th percentile of trimmed read length, (2) stop where the smoothed median quality falls below Q25, and (3) keep `truncLen_F + truncLen_R ≥ insert length + 20` so that pairs can merge (insert length from the primer positions on *E. coli* 16S). When (2) and (3) conflict, overlap wins and R1 is extended first; the batch is flagged. Proposals (`quality_summary.tsv`) and plots (`quality_<batch>.png`) are reviewed before copying the values into `config/dada2_params.tsv`.

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
| 2026-10-03 | DADA2 parameters inherited from the predoc (Chapter 1 methods) | cutadapt v5.2 with `--discard-untrimmed`; DADA2 v1.28 in R 4.4.0; `filterAndTrim(maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE)` with `truncLen` set per batch; `mergePairs(maxMismatch = 0)` with `minOverlap` set per batch; `removeBimeraDenovo(method = "consensus")`; `drop_rare_asvs(at_least_n = 2)`; SILVA v138.2, `assignTaxonomy(minBoot = 80)` + `addSpecies()`; chloroplast and mitochondria removed; genus as taxonomic unit |
| 2026-10-03 | Changes to the predoc pipeline | (1) one batch per BioProject instead of 7 region groups; (2) reads oriented before DADA2 instead of reverse-complementing ASVs of Groups E and F before taxonomy; (3) `learnErrors(nbases = 3e8)` in every batch, not only Group_F; (4) `pool = "pseudo"` in every batch, as stated in the predoc (the Group_F rerun had used `pool = FALSE`); (5) error model fitted with enforced monotonicity for NovaSeq batches (B06, B10), whose binned quality scores break the default fit |
| 2026-10-03 | Error model unit: sequencing run when known (B05, B10, B11), otherwise BioProject; B06 pooled into one model | SRA did not keep original read names in 7 BioProjects; B06 has only 2–4 runs per flow cell, too few to train a model each |
| 2026-10-03 | Core definition pending | The predoc defined core genera as prevalence ≥ 70% in ≥ 3 of 7 region groups, which no longer exist. Proposed: same threshold with BioProject as the unit; to be decided with the supervisor |
| 2026-10-03 | Reads re-downloaded with `fasterq-dump --split-3`; pairing checked before trimming | The first download used `--split-files`, which put unmated reads into `_1`/`_2` and broke the pairing (e.g. SRR10966917: 23,675 vs 21,244 reads). cutadapt stopped at the first mismatch, so B01 and B03 kept only 63–66% of reads. Step 06 now refuses unpaired inputs and records failed cutadapt calls |
| 2026-10-03 | B06 (PRJNA987554) forward primer set to the 341F variant `CCTAYGGGRBGCASCAG` | With 338F only 7.7% of R1 reads matched; the predoc had trimmed the same runs with this primer (99.6% of pairs kept). 338F and 341F are indistinguishable by alignment position |
