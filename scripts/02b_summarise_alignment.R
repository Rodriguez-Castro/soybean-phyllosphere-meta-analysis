#!/usr/bin/env Rscript
# Step 02b: interpret the BLAST positions of step 01b for every run and summarise per BioProject.
#
# Usage (from the repository root):
#   Rscript scripts/02b_summarise_alignment.R [alignment.tsv] [run_list.tsv] [outdir] [step02_run_summary.tsv]
#
# For each run:
#   class          16S | fungal_ITS_like | inconsistent | non_16S | 16S_not_amplicon | download_failed
#   amplicon       E. coli span covered by the reads (start_plus .. start_minus)
#   region         hypervariable regions covered at least MIN_COVER by that span
#   fwd / rev      primers whose position matches the read ends, with status:
#                  "present" (reads start at the primer) or "removed" (reads start right after it)
#   spacer_fwd/rev bases before the alignment start (spacer/barcode/heterogeneity spacer)
#   orientation    R1_forward | R1_reverse | mixed
#
# Outputs in <outdir>: run_alignment_summary.tsv, bioproject_alignment_summary.tsv, bacterial_runs.tsv
# Base R only.

args <- commandArgs(trailingOnly = TRUE)
aln_file     <- if (length(args) >= 1) args[1] else "results/01b_alignment/alignment.tsv"
runlist_file <- if (length(args) >= 2) args[2] else "metadata/run_list.tsv"
outdir       <- if (length(args) >= 3) args[3] else "results/02b_summary"
step02_file  <- if (length(args) >= 4) args[4] else "results/02_summary/run_summary.tsv"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- parameters -------------------------------------------------------------
MIN_16S  <- 0.50  # minimum fraction of reads aligning to 16S
ITS_POS  <- 1450  # all read ends beyond this position -> fungal 18S end (ITS amplicon)
TOL      <- 3     # tolerance (bp) when matching read ends to primer positions
MIX_LOW  <- 0.10  # R1 forward-strand fraction between MIX_LOW and 1 - MIX_LOW -> mixed
MIN_COVER <- 0.80 # a hypervariable region counts if the amplicon covers at least 80% of it

# Primer binding sites on E. coli 16S (J01859 numbering): first and last base covered.
PRIMERS <- data.frame(
  name  = c("27F", "338F", "341F", "515F", "520F", "799F",
            "534R", "785R", "799R", "806R", "907R", "1115R", "1193R"),
  dir   = c(rep("F", 6), rep("R", 7)),
  start = c(8,     338,    341,    515,    520,    781,
            518,   785,    781,    787,    907,    1100,   1176),
  end   = c(27,    357,    357,    533,    534,    799,
            534,   805,    799,    806,    926,    1115,   1193),
  stringsAsFactors = FALSE
)

VREGIONS <- data.frame(
  name  = paste0("V", 1:9),
  start = c(69, 137, 433, 576, 822, 986, 1117, 1243, 1435),
  end   = c(99, 242, 497, 682, 879, 1043, 1173, 1294, 1465)
)

# ---- input ------------------------------------------------------------------
aln <- read.delim(aln_file, colClasses = "character", na.strings = "NA")
num_cols <- c("n_reads", "n_16S", "frac_16S", "frac_plus", "start_plus",
              "start_minus", "qstart_plus", "qstart_minus")
for (cc in num_cols) aln[[cc]] <- suppressWarnings(as.numeric(aln[[cc]]))
runs <- read.delim(runlist_file, colClasses = "character")

# ---- helpers ----------------------------------------------------------------
match_primer <- function(pos, direction) {
  if (is.na(pos)) return(NA_character_)
  p <- PRIMERS[PRIMERS$dir == direction, ]
  if (direction == "F") {
    present <- p$name[abs(pos - p$start) <= TOL]
    removed <- p$name[abs(pos - (p$end + 1)) <= TOL]
  } else {
    present <- p$name[abs(pos - p$end) <= TOL]
    removed <- p$name[abs(pos - (p$start - 1)) <= TOL]
  }
  out <- c(if (length(present)) paste(paste(present, collapse = "|"), "present"),
           if (length(removed)) paste(paste(removed, collapse = "|"), "removed"))
  if (length(out)) paste(out, collapse = " / ") else paste0("unknown (", pos, ")")
}

region_of <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_character_)
  overlap <- pmax(0, pmin(VREGIONS$end, b) - pmax(VREGIONS$start, a) + 1)
  v <- VREGIONS$name[overlap / (VREGIONS$end - VREGIONS$start + 1) >= MIN_COVER]
  if (!length(v)) return("none")
  if (length(v) == 1) v else paste0(v[1], "-", v[length(v)])
}

mean_na <- function(x) if (all(is.na(x))) NA_real_ else round(mean(x, na.rm = TRUE))

# ---- classify each run ------------------------------------------------------
classify <- function(d) {
  base <- data.frame(run = d$run[1], bioproject = d$bioproject[1], stringsAsFactors = FALSE)
  if (all(is.na(d$read)) || all(d$n_reads == 0, na.rm = TRUE)) {
    return(cbind(base, class = "download_failed", frac_16S = NA, start_plus = NA,
                 start_minus = NA, amplicon = NA, region = NA, fwd = NA, rev = NA,
                 spacer_fwd = NA, spacer_rev = NA, orientation = NA))
  }
  frac <- max(d$frac_16S, na.rm = TRUE)
  sp <- mean_na(d$start_plus)
  sm <- mean_na(d$start_minus)
  ends <- c(sp, sm)
  ends <- ends[!is.na(ends)]

  class <- if (!length(ends)) "non_16S"
           else if (all(ends >= ITS_POS)) "fungal_ITS_like"
           else if (any(ends >= ITS_POS)) "inconsistent"
           else if (frac < MIN_16S) "non_16S"
           else "16S"

  r1 <- d[d$read %in% c("R1", "SE"), ]
  fp <- if (nrow(r1)) r1$frac_plus[1] else NA
  orientation <- if (is.na(fp)) NA_character_
                 else if (fp >= 1 - MIX_LOW) "R1_forward"
                 else if (fp <= MIX_LOW) "R1_reverse"
                 else "mixed"

  is16 <- class == "16S"
  cbind(base,
        class       = class,
        frac_16S    = round(frac, 3),
        start_plus  = sp,
        start_minus = sm,
        amplicon    = if (is16 && !is.na(sp) && !is.na(sm)) paste0(sp, "-", sm) else NA,
        region      = if (is16) region_of(sp, sm) else NA,
        fwd         = if (is16) match_primer(sp, "F") else NA,
        rev         = if (is16) match_primer(sm, "R") else NA,
        spacer_fwd  = if (is16) mean_na(d$qstart_plus) - 1 else NA,
        spacer_rev  = if (is16) mean_na(d$qstart_minus) - 1 else NA,
        orientation = if (is16) orientation else NA)
}

per_run <- do.call(rbind, lapply(split(aln, aln$run), classify))

# add metadata and, if available, the sequencing run from step 02
meta_cols <- intersect(c("run", "library_name", "library_strategy", "model", "sample_name"), names(runs))
per_run <- merge(per_run, runs[, meta_cols], by = "run", all.x = TRUE)
if (file.exists(step02_file)) {
  s2 <- read.delim(step02_file, colClasses = "character")
  per_run <- merge(per_run, s2[, c("run", "seq_run_id")], by = "run", all.x = TRUE)
}

# RNA-Seq, WGS and other non-amplicon libraries are never bacterial amplicons
if ("library_strategy" %in% names(per_run)) {
  not_amp <- !is.na(per_run$library_strategy) & per_run$library_strategy != "AMPLICON" &
             per_run$class == "16S"
  per_run$class[not_amp] <- "16S_not_amplicon"
}
per_run <- per_run[order(per_run$bioproject, per_run$run), ]

# ---- summarise per BioProject -----------------------------------------------
tab_string <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  t <- sort(table(x), decreasing = TRUE)
  paste0(names(t), " (", as.integer(t), ")", collapse = "; ")
}

bp <- do.call(rbind, lapply(split(per_run, per_run$bioproject), function(d) {
  b <- d[d$class == "16S", ]
  data.frame(
    bioproject       = d$bioproject[1],
    n_runs           = nrow(d),
    bacterial_16S    = nrow(b),
    fungal_ITS_like  = sum(d$class == "fungal_ITS_like"),
    inconsistent     = sum(d$class == "inconsistent"),
    non_16S          = sum(d$class %in% c("non_16S", "16S_not_amplicon")),
    failed           = sum(d$class == "download_failed"),
    region           = tab_string(b$region),
    primers          = if (nrow(b)) tab_string(paste(b$fwd, "+", b$rev)) else "",
    spacer_fwd       = tab_string(b$spacer_fwd),
    orientation      = tab_string(b$orientation),
    n_seq_runs       = if ("seq_run_id" %in% names(b)) length(unique(b$seq_run_id)) else NA,
    stringsAsFactors = FALSE
  )
}))

# ---- write ------------------------------------------------------------------
write.table(per_run, file.path(outdir, "run_alignment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(bp, file.path(outdir, "bioproject_alignment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
keep <- intersect(c("run", "bioproject", "seq_run_id", "region", "amplicon", "fwd", "rev",
                    "spacer_fwd", "spacer_rev", "orientation"), names(per_run))
write.table(per_run[per_run$class == "16S", keep], file.path(outdir, "bacterial_runs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\nRuns per class:\n"); print(table(per_run$class))
cat("\nPer BioProject:\n")
print(bp[, c("bioproject", "n_runs", "bacterial_16S", "fungal_ITS_like", "inconsistent", "non_16S", "failed", "n_seq_runs")],
      row.names = FALSE)
for (i in seq_len(nrow(bp))) {
  cat("\n", bp$bioproject[i], "\n",
      "  region:      ", bp$region[i], "\n",
      "  primers:     ", bp$primers[i], "\n",
      "  spacer_fwd:  ", bp$spacer_fwd[i], "\n",
      "  orientation: ", bp$orientation[i], "\n", sep = "")
}
cat("\nOutputs written to", outdir, "\n")
