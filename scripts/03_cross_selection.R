#!/usr/bin/env Rscript
# Step 03: cross the phyllosphere sample selection with the verified run classification (step 02b).
#
# Usage (from the repository root):
#   Rscript scripts/03_cross_selection.R [selection.csv] [run_alignment_summary.tsv] [outdir]
#
# selection.csv: the runs selected for the meta-analysis (phyllosphere), with at least a
# column "run"; optional columns "bio_project" and "dada2_group" (previous grouping).
#
# Outputs in <outdir>:
#   selection_check.tsv         every selected run with its verified class, region and primers
#   selection_summary.tsv       counts per BioProject (and previous group): kept and lost, by reason
#   selected_bacterial_runs.tsv selected runs verified as bacterial 16S (input for DADA2)
# Base R only.

args <- commandArgs(trailingOnly = TRUE)
sel_file <- if (length(args) >= 1) args[1] else "metadata/phyllosphere_runs.csv"
cls_file <- if (length(args) >= 2) args[2] else "results/02b_summary/run_alignment_summary.tsv"
outdir   <- if (length(args) >= 3) args[3] else "results/03_selection"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sel <- read.csv(sel_file, colClasses = "character", check.names = FALSE)
if (!"run" %in% names(sel)) stop("selection file needs a column named 'run'; columns found: ",
                                 paste(names(sel), collapse = ", "))
sel <- sel[!duplicated(sel$run), ]
cls <- read.delim(cls_file, colClasses = "character")

keep_cls <- intersect(c("run", "bioproject", "class", "region", "fwd", "rev", "spacer_fwd",
                        "spacer_rev", "orientation", "seq_run_id"), names(cls))
x <- merge(sel, cls[, keep_cls], by = "run", all.x = TRUE)
x$class[is.na(x$class)] <- "not_verified"   # selected run absent from the 10 verified BioProjects

if (!"bioproject" %in% names(x) || all(is.na(x$bioproject))) x$bioproject <- NA
if ("bio_project" %in% names(x)) x$bioproject[is.na(x$bioproject)] <- x$bio_project[is.na(x$bioproject)]
group_col <- if ("dada2_group" %in% names(x)) "dada2_group" else NULL

# ---- summary ----------------------------------------------------------------
key <- if (is.null(group_col)) x["bioproject"] else x[c("bioproject", group_col)]
summ <- aggregate(list(selected = x$run), key, length)
for (cl in c("16S", "fungal_ITS_like", "inconsistent", "non_16S", "16S_not_amplicon",
             "download_failed", "not_verified")) {
  cnt <- aggregate(list(n = x$class == cl), key, sum)
  summ[[cl]] <- cnt$n[match(do.call(paste, summ[names(key)]), do.call(paste, cnt[names(key)]))]
}
names(summ)[names(summ) == "16S"] <- "kept_16S"
summ <- summ[order(summ$bioproject), ]

# ---- write ------------------------------------------------------------------
write.table(x, file.path(outdir, "selection_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(summ, file.path(outdir, "selection_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
out_cols <- intersect(c("run", "bioproject", group_col, "seq_run_id", "region", "fwd", "rev",
                        "spacer_fwd", "spacer_rev", "orientation"), names(x))
write.table(x[x$class == "16S", out_cols], file.path(outdir, "selected_bacterial_runs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\nSelected runs:", nrow(x), "| kept as bacterial 16S:", sum(x$class == "16S"), "\n\n")
cat("Selected runs per verified class:\n"); print(table(x$class))
cat("\nPer BioProject:\n"); print(summ, row.names = FALSE)
cat("\nOutputs written to", outdir, "\n")
