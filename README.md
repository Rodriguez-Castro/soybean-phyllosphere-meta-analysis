# Soybean phyllosphere bacterial meta-analysis

Reproducible pipeline to reprocess public 16S rRNA amplicon data from the phyllosphere of *Glycine* spp. and integrate studies at genus level.

## Strategy

1. **Verify every run from its reads**, not from the published methods: detect which primers each run actually contains, and which sequencing run (flow cell) it comes from.
2. **One DADA2 batch per BioProject** (see `config/batches.tsv`). If a BioProject contains several sequencing runs (flow cells), error models are learned per sequencing run within that batch.
3. **Taxonomy assigned per batch** with identical settings (same SILVA release, `minBoot`, seed).
4. **Integration at genus level** across all batches and regions, with study/sequencing run kept as covariates.

## Repository structure

```
soy-phyllosphere-meta/
├── README.md
├── config/
│   ├── batches.tsv              # one DADA2 batch per BioProject, reported vs verified primers
│   └── primers.tsv              # primer catalogue (IUPAC) used for verification and trimming
├── metadata/
│   ├── runinfo/                 # SRA RunInfo CSVs, one per BioProject (input)
│   └── run_list.tsv             # combined run list (step 00)
├── scripts/
│   ├── 00_build_run_list.py
│   └── 01_verify_runs.sh
├── results/
│   └── 01_verification/
└── logs/
```

## Requirements (hpc-bio, Compute Canada software stack)

```bash
newgrp def-ilafores
module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9
```

Python 3 (standard library only) for step 00.

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

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-10-01 | DADA2 and taxonomy run separately per BioProject (10 batches); studies merged only at genus level. Hypervariable region is kept as a descriptive variable and covariate, not as a processing group | Each BioProject is, in practice, a separate sequencing run and error models are run-specific; ASVs from different primers or truncation lengths cannot be merged (e.g., 515F/806R vs 520F/799R in V4) |
| 2026-10-01 | PRJNA601979 relabelled from V2–V3 to V4 (pending read verification) | The paper states that the "V2–V3" region was amplified with 520F/799R. These primers bind around *E. coli* positions 520 and 799, downstream of V3 (≈433–497) and upstream of V5 (≈822–879), so the amplicon can only span V4 (≈576–682). The literature consistently describes 520F/799R as a V4 primer pair. To be confirmed in step 01 (primer match at read start) and by aligning reads to *E. coli* 16S (J01695). |
| 2026-10-01 | Primers assigned from the reads, not from the papers | PRJNA1092852: paper reports 338F/806R, but bacterial libraries (suffix `.1b`, 134 runs) contain 799F/1193R (V5–V7) in ~99% of reads, in mixed orientation (~50% each). Libraries with suffix `.1` (129 runs) are fungal ITS1F/ITS2 and are excluded. Authors contacted for confirmation. |
