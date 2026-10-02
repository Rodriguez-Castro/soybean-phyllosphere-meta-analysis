#!/usr/bin/env Rscript
# Step 02: classify every run from the primer counts of step 01 and summarise per BioProject.
#
# Usage (from the repository root):
#   Rscript scripts/02_summarise_verification.R [counts.tsv] [run_list.tsv] [outdir]
#
# Outputs in <outdir>:
#   run_summary.tsv         one row per run: class, primer pair, region, orientation, sequencing run
#   bioproject_summary.tsv  one row per BioProject: run classes, primer pairs, sequencing runs
#   bacterial_runs.tsv      runs classified as bacterial 16S (input for DADA2 batches)
#
# Base R only, no packages required.

args <- commandArgs(trailingOnly = TRUE)
counts_file  <- if (length(args) >= 1) args[1] else "results/01_verification/primer_counts.tsv"
runlist_file <- if (length(args) >= 2) args[2] else "metadata/run_list.tsv"
outdir       <- if (length(args) >= 3) args[3] else "results/02_summary"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- thresholds -------------------------------------------------------------
MIN_PAIR    <- 0.80  # fraction of reads starting with the assigned primer pair -> bacterial_16S
MIN_PARTIAL <- 0.20  # below MIN_PAIR but above this -> 16S_partial (check manually)
MIN_ITS     <- 0.50  # fraction of reads starting with an ITS primer -> fungal_ITS
MIN_MIX     <- 0.10  # both primers present in R1 above this -> mixed orientation
TIE         <- 0.90  # primers within 90% of the top count are considered tied

# Nested primers: when both match the same reads, the outer one is the real primer.
# 338F reads also contain the 341F site, 515F reads the 520F site,
# and 806R/785R reads the 785R/799R sites.
NESTED <- list("338F" = "341F", "515F" = "520F", "806R" = c("785R", "799R"), "785R" = "799R")

# Primer pair -> hypervariable region
REGIONS <- data.frame(
  fwd    = c("27F",   "338F",  "341F",  "341F",  "515F", "515F",  "520F", "799F",  "799F"),
  rev    = c("534R",  "806R",  "806R",  "785R",  "806R", "907R",  "799R", "1115R", "1193R"),
  region = c("V1-V3", "V3-V4", "V3-V4", "V3-V4", "V4",   "V4-V5", "V4",   "V5-V6", "V5-V7"),
  stringsAsFactors = FALSE
)

# ---- input ------------------------------------------------------------------
counts  <- read.delim(counts_file, colClasses = "character")
counts$count   <- as.numeric(counts$count)
counts$n_reads <- as.numeric(counts$n_reads)
runs    <- read.delim(runlist_file, colClasses = "character")
primers <- read.delim("config/primers.tsv", colClasses = "character", comment.char = "#")

# ---- helpers ----------------------------------------------------------------
pick_primer <- function(tot, suffix) {
  d <- tot[tot$target == "16S" & grepl(paste0(suffix, "$"), tot$primer), ]
  if (nrow(d) == 0 || max(d$count) == 0) return(list(name = NA_character_, count = 0))
  tied <- d$primer[d$count >= TIE * max(d$count)]
  for (outer in names(NESTED)) {
    if (outer %in% tied) tied <- setdiff(tied, NESTED[[outer]])
  }
  list(name = paste(sort(tied), collapse = "|"), count = max(d$count))
}

region_of <- function(fwd, rev) {
  if (is.na(fwd) || is.na(rev)) return(NA_character_)
  hit <- REGIONS$region[REGIONS$fwd == fwd & REGIONS$rev == rev]
  if (length(hit)) return(hit[1])
  if (grepl("|", paste(fwd, rev), fixed = TRUE)) "ambiguous" else "unknown"
}

max_count <- function(d, names) {
  x <- d$count[d$primer %in% strsplit(names, "|", fixed = TRUE)[[1]]]
  if (length(x)) max(x) else 0
}

empty_row <- function(base, class) {
  cbind(base, class = class, fwd = NA, rev = NA, region = NA,
        pair_frac = NA, its_frac = NA, orientation = NA)
}

# ---- classify each run ------------------------------------------------------
classify_run <- function(d) {
  base <- data.frame(
    run        = d$run[1],
    bioproject = d$bioproject[1],
    seq_run_id = paste(d$instrument[1], d$seq_run[1], d$flowcell[1], sep = ":"),
    stringsAsFactors = FALSE
  )
  if (any(d$primer == "DOWNLOAD_FAILED")) return(empty_row(base, "download_failed"))

  d$target <- primers$target[match(d$primer, primers$name)]
  reads <- unique(d[, c("read", "n_reads")])
  total <- sum(reads$n_reads)
  if (total == 0) return(empty_row(base, "no_reads"))

  tot <- aggregate(count ~ primer + target, data = d, FUN = sum)
  f <- pick_primer(tot, "F")
  r <- pick_primer(tot, "R")
  pair_frac <- (f$count + r$count) / total
  its_frac  <- sum(tot$count[tot$target == "ITS"]) / total

  class <- if (its_frac >= MIN_ITS) "fungal_ITS"
           else if (pair_frac >= MIN_PAIR) "bacterial_16S"
           else if (pair_frac >= MIN_PARTIAL) "16S_partial"
           else "no_primer_detected"

  orientation <- NA_character_
  r1 <- d[d$read %in% c("R1", "SE"), ]
  if (nrow(r1) && !is.na(f$name) && !is.na(r$name)) {
    n1 <- r1$n_reads[1]
    f1 <- max_count(r1, f$name) / n1
    rr <- max_count(r1, r$name) / n1
    orientation <- if (f1 >= MIN_MIX && rr >= MIN_MIX) "mixed"
                   else if (f1 >= rr) "R1_forward"
                   else "R1_reverse"
  }

  cbind(base, class = class, fwd = f$name, rev = r$name,
        region = region_of(f$name, r$name),
        pair_frac = round(pair_frac, 3), its_frac = round(its_frac, 3),
        orientation = orientation)
}

per_run <- do.call(rbind, lapply(split(counts, counts$run), classify_run))
meta_cols <- intersect(c("run", "library_name", "library_strategy", "library_layout",
                         "model", "sample_name"), names(runs))
per_run <- merge(per_run, runs[, meta_cols], by = "run", all.x = TRUE)
per_run <- per_run[order(per_run$bioproject, per_run$run), ]

# ---- summarise per BioProject -----------------------------------------------
tab_string <- function(x) {
  if (!length(x)) return("")
  t <- sort(table(x), decreasing = TRUE)
  paste0(names(t), " (", as.integer(t), ")", collapse = "; ")
}

summarise_bp <- function(d) {
  b <- d[d$class == "bacterial_16S", ]
  data.frame(
    bioproject         = d$bioproject[1],
    n_runs             = nrow(d),
    bacterial_16S      = nrow(b),
    fungal_ITS         = sum(d$class == "fungal_ITS"),
    partial_16S        = sum(d$class == "16S_partial"),
    no_primer_detected = sum(d$class == "no_primer_detected"),
    failed             = sum(d$class %in% c("download_failed", "no_reads")),
    mixed_orientation  = sum(b$orientation == "mixed", na.rm = TRUE),
    primer_pairs       = if (nrow(b)) tab_string(paste0(b$fwd, "/", b$rev, " ", b$region)) else "",
    n_seq_runs         = length(unique(b$seq_run_id)),
    seq_runs           = tab_string(b$seq_run_id),
    stringsAsFactors   = FALSE
  )
}

bp <- do.call(rbind, lapply(split(per_run, per_run$bioproject), summarise_bp))

# ---- write ------------------------------------------------------------------
write.table(per_run, file.path(outdir, "run_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(bp, file.path(outdir, "bioproject_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(per_run[per_run$class == "bacterial_16S",
                    c("run", "bioproject", "seq_run_id", "fwd", "rev", "region", "orientation")],
            file.path(outdir, "bacterial_runs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\nRuns per class:\n"); print(table(per_run$class))
cat("\nPer BioProject:\n")
print(bp[, c("bioproject", "n_runs", "bacterial_16S", "fungal_ITS", "partial_16S",
             "no_primer_detected", "failed", "n_seq_runs")], row.names = FALSE)
cat("\nPrimer pairs (bacterial runs):\n")
for (i in seq_len(nrow(bp))) cat(" ", bp$bioproject[i], ":", bp$primer_pairs[i], "\n")
cat("\nOutputs written to", outdir, "\n")
