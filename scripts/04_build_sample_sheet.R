#!/usr/bin/env Rscript
# Step 04: build the final sample sheet for DADA2.
#
# Usage (from the repository root):
#   Rscript scripts/04_build_sample_sheet.R
#
# Inputs:
#   results/03_selection/selection_check.tsv   selected runs with their verified class (step 03)
#   results/02b_summary/run_alignment_summary.tsv
#   metadata/run_list.tsv                       library layout and names (step 00)
#   config/batches.tsv                          one DADA2 batch per BioProject, verified primers
#   config/primers.tsv                          primer sequences
#   config/selection_overrides.tsv              documented exceptions to the original selection
#
# Included runs: verified 16S AND (previous group not "excluded_*" OR listed in selection_overrides).
#
# Outputs:
#   metadata/sample_sheet.tsv        one row per run: batch, sequencing run, primers to trim, etc.
#   results/04_sample_sheet/batch_summary.tsv
# Base R only.

sel_file  <- "results/03_selection/selection_check.tsv"
aln_file  <- "results/02b_summary/run_alignment_summary.tsv"
runs_file <- "metadata/run_list.tsv"
outdir    <- "results/04_sample_sheet"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sel     <- read.delim(sel_file, colClasses = "character")
aln     <- read.delim(aln_file, colClasses = "character")
runs    <- read.delim(runs_file, colClasses = "character")
batches <- read.delim("config/batches.tsv", colClasses = "character")
primers <- read.delim("config/primers.tsv", colClasses = "character", comment.char = "#")
over    <- read.delim("config/selection_overrides.tsv", colClasses = "character")

# ---- which runs are included ------------------------------------------------
group <- if ("dada2_group" %in% names(sel)) sel$dada2_group else rep("selected", nrow(sel))
override <- paste(sel$bioproject, group) %in% paste(over$bioproject, over$previous_group)
included <- sel$class == "16S" & (!grepl("^excluded", group) | override)
x <- sel[included, ]
x$previous_group <- group[included]
x$included_by    <- ifelse(override[included], "override", "selection")

# ---- add layout, library and sample names -----------------------------------
x <- merge(x, runs[, intersect(c("run", "library_layout", "library_name", "sample_name"), names(runs))],
           by = "run", all.x = TRUE)

# ---- batch-level information --------------------------------------------------
b <- batches[match(x$bioproject, batches$bioproject), ]
x$batch            <- b$batch
x$leaf_compartment <- b$leaf_compartment
x$region           <- b$verified_region
pair <- sub(" .*", "", b$verified_primers)            # e.g. "515F/806R"
x$fwd_primer <- sub("/.*", "", pair)
x$rev_primer <- sub(".*/", "", pair)
x$fwd_seq <- primers$sequence[match(x$fwd_primer, primers$name)]
x$rev_seq <- primers$sequence[match(x$rev_primer, primers$name)]

# ---- run-level primer status (present: trim; removed: nothing to trim) --------
status <- function(s) ifelse(is.na(s), "unknown",
                      ifelse(grepl("present", s), "present",
                      ifelse(grepl("removed", s), "removed", "unknown")))
x$fwd_status <- status(x$fwd)
x$rev_status <- status(x$rev)

# runs whose read ends did not match a known site take the majority status of their batch
fill_unknown <- function(st, bt) {
  for (k in unique(bt)) {
    i <- bt == k
    known <- st[i & st != "unknown"]
    if (length(known)) st[i & st == "unknown"] <- names(sort(table(known), decreasing = TRUE))[1]
  }
  st
}
x$status_flag <- ifelse(x$fwd_status == "unknown" | x$rev_status == "unknown",
                        "filled_from_batch_majority", "")
x$fwd_status <- fill_unknown(x$fwd_status, x$batch)
x$rev_status <- fill_unknown(x$rev_status, x$batch)

# ---- write sample sheet -------------------------------------------------------
cols <- c("run", "batch", "bioproject", "seq_run_id", "region", "leaf_compartment",
          "library_layout", "orientation", "fwd_primer", "fwd_seq", "fwd_status",
          "rev_primer", "rev_seq", "rev_status", "spacer_fwd", "spacer_rev",
          "status_flag", "previous_group", "included_by", "library_name", "sample_name")
x <- x[order(x$batch, x$seq_run_id, x$run), intersect(cols, names(x))]
write.table(x, "metadata/sample_sheet.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- summary per batch ----------------------------------------------------------
tab_string <- function(v) {
  v <- v[!is.na(v) & v != ""]
  if (!length(v)) return("")
  t <- sort(table(v), decreasing = TRUE)
  paste0(names(t), " (", as.integer(t), ")", collapse = "; ")
}
summ <- do.call(rbind, lapply(split(x, x$batch), function(d) data.frame(
  batch        = d$batch[1],
  bioproject   = d$bioproject[1],
  region       = d$region[1],
  primers      = paste0(d$fwd_primer[1], "/", d$rev_primer[1]),
  n_runs       = nrow(d),
  n_seq_runs   = length(unique(d$seq_run_id)),
  seq_runs     = tab_string(d$seq_run_id),
  primer_trim  = tab_string(paste0(d$fwd_status, "/", d$rev_status)),
  orientation  = tab_string(d$orientation),
  layout       = tab_string(d$library_layout),
  compartment  = d$leaf_compartment[1],
  flagged      = sum(d$status_flag != ""),
  stringsAsFactors = FALSE)))
write.table(summ, file.path(outdir, "batch_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nIncluded runs:", nrow(x), "(", sum(x$included_by == "override"), "by override )\n\n")
print(summ[, c("batch", "bioproject", "region", "primers", "n_runs", "n_seq_runs", "flagged")], row.names = FALSE)
for (i in seq_len(nrow(summ))) {
  cat("\n", summ$batch[i], " ", summ$bioproject[i], "\n",
      "  sequencing runs: ", summ$seq_runs[i], "\n",
      "  primer trimming: ", summ$primer_trim[i], "\n",
      "  orientation:     ", summ$orientation[i], "\n",
      "  layout:          ", summ$layout[i], "\n", sep = "")
}
cat("\nSample sheet written to metadata/sample_sheet.tsv\n")
