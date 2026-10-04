#!/bin/bash
# Step 08: run DADA2 for every batch as a Slurm job array (one task per batch, 3 at a time).
#
# Usage (from the repository root, on the login node; Slurm runs it on the compute node):
#   newgrp def-ilafores
#   sbatch scripts/08_dada2.sbatch              # all batches
#   sbatch --array=5 scripts/08_dada2.sbatch    # only one batch (5 = B06), e.g. as a test
#
# Array index -> batch:  1 B01  2 B03  3 B04  4 B05  5 B06  6 B07  7 B08  8 B09  9 B10  10 B11
# Logs: logs/08_dada2_<jobid>_<index>.log ; results: results/08_dada2/<batch>/

#SBATCH --job-name=dada2
#SBATCH --partition=bio
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=1-10%3
#SBATCH --output=logs/08_dada2_%A_%a.log

BATCHES=(B01 B03 B04 B05 B06 B07 B08 B09 B10 B11)
BATCH=${BATCHES[$((SLURM_ARRAY_TASK_ID - 1))]}

module load StdEnv/2023 r/4.4.0
cd "$SLURM_SUBMIT_DIR" || exit 1
echo "Batch $BATCH on $(hostname -s), $SLURM_CPUS_PER_TASK CPUs"
Rscript scripts/08_dada2_batch.R "$BATCH" "$SLURM_CPUS_PER_TASK"
