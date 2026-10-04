#!/usr/bin/env Rscript
# Step 07: quality profiles per batch and proposed truncLen / minOverlap for DADA2.
#
# Usage (from the repository root):
#   Rscript scripts/07_quality_profiles.R [trim_manifest.tsv] [outdir] [runs_per_batch] [reads_per_run]
#
# For each batch (and each orientation set A/B in B05), reads a sample of trimmed reads from up to
# <runs_per_batch> runs and computes, for R1 and R2: median quality per position and read length
# distribution. It then proposes truncLen values that
#   (1) do not exceed the 5th percentile of read length (filterAndTrim drops shorter reads),
#   (2) stop where the median quality falls below Q_MIN, and
#   (3) keep enough overlap to merge: truncLen_F + truncLen_R >= insert length + MIN_OVERLAP + MARGIN.
# If (2) and (3) conflict, overlap wins (the quality cut is relaxed) and the batch is flagged.
#
# Outputs in <outdir>: quality_<batch>.png, quality_summary.tsv (proposals to review before
# copying them into config/dada2_params.tsv). Base R only.

# refuse to run on the login node (lab rule: no computation on iv12; use salloc --partition bio)
if (grepl("^iv", Sys.info()[["nodename"]])) {
  stop("you are on the login node (", Sys.info()[["nodename"]], "). Request a compute node first:\n",
       "  salloc --time=24:00:00 --cpus-per-task=4 --mem=8G --partition bio", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
manifest_file <- if (length(args) >= 1) args[1] else "results/06_trim/trim_manifest.tsv"
outdir        <- if (length(args) >= 2) args[2] else "results/07_quality"
RUNS          <- if (length(args) >= 3) as.integer(args[3]) else 10
READS         <- if (length(args) >= 4) as.integer(args[4]) else 2000
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

Q_MIN       <- 25   # median quality threshold for truncation
MIN_OVERLAP <- 12   # minimum overlap used in the predoc
MARGIN      <- 8    # extra overlap for length variation between taxa
set.seed(1)

# Insert length between primers on E. coli 16S (bases between forward primer end and
# reverse primer start), i.e. the length of a merged read once primers are removed.
COORD <- data.frame(
  primer = c("515F", "520F", "341Fv", "338F", "799F", "806R", "799R", "907R", "1115R", "1193R"),
  pos    = c(533,    534,    357,     357,    799,    787,    781,    907,    1100,    1176))
insert_len <- function(f, r) COORD$pos[COORD$primer == r] - COORD$pos[COORD$primer == f] - 1

man   <- read.delim(manifest_file, colClasses = "character")
man   <- man[!is.na(man$r1) & man$r1 != "NA" & as.numeric(man$reads_out) > 0, ]
sheet <- read.delim("metadata/sample_sheet.tsv", colClasses = "character")
pairs <- unique(sheet[, c("batch", "fwd_primer", "rev_primer")])
man$key <- ifelse(man$set %in% c("A", "B"), paste0(man$batch, "_", man$set), man$batch)

read_fastq <- function(path, n) {
  con <- gzfile(path, "r"); on.exit(close(con))
  x <- readLines(con, n = 4 * n)
  list(len = nchar(x[seq(2, length(x), 4)]), qual = x[seq(4, length(x), 4)])
}

profile <- function(paths) {
  q <- lapply(paths, read_fastq, n = READS)
  len  <- unlist(lapply(q, `[[`, "len"))
  qual <- unlist(lapply(q, `[[`, "qual"))
  L <- max(len)
  m <- matrix(NA_integer_, nrow = length(qual), ncol = L)
  for (i in seq_along(qual)) {
    v <- utf8ToInt(qual[i]) - 33L
    m[i, seq_along(v)] <- v
  }
  qs <- apply(m, 2, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  list(len = len, q25 = qs[1, ], med = qs[2, ], q75 = qs[3, ],
       cover = colMeans(!is.na(m)), p5 = unname(quantile(len, 0.05)))
}

last_good <- function(med, qmin) {   # last position before the smoothed median drops below qmin
  s <- stats::filter(med, rep(1 / 10, 10), sides = 2)
  s[is.na(s)] <- med[is.na(s)]
  bad <- which(s < qmin)
  if (length(bad)) max(1, min(bad) - 1) else length(med)
}

rows <- list()
for (key in sort(unique(man$key))) {
  d <- man[man$key == key, ]
  d <- d[sample(nrow(d), min(RUNS, nrow(d))), ]
  b <- d$batch[1]
  pr <- pairs[pairs$batch == b, ][1, ]
  L  <- insert_len(pr$fwd_primer, pr$rev_primer)
  p1 <- profile(d$r1); p2 <- profile(d$r2)

  tF <- min(floor(p1$p5), last_good(p1$med, Q_MIN))
  tR <- min(floor(p2$p5), last_good(p2$med, Q_MIN))
  need <- L + MIN_OVERLAP + MARGIN
  flag <- ""
  if (tF + tR < need) {      # relax the quality cut, never beyond the read length limit
    extra <- need - (tF + tR)
    # R1 first: its quality is usually better than R2's
    addF <- min(extra, floor(p1$p5) - tF); tF <- tF + addF; extra <- extra - addF
    addR <- min(extra, floor(p2$p5) - tR); tR <- tR + addR; extra <- extra - addR
    flag <- if (extra > 0) "CANNOT_MERGE: reads too short for this insert"
            else "quality cut relaxed to keep overlap"
  }
  overlap <- tF + tR - L
  rows[[key]] <- data.frame(
    key = key, batch = b, primers = paste0(pr$fwd_primer, "/", pr$rev_primer),
    insert_len = L, runs_sampled = nrow(d),
    len_p5_F = floor(p1$p5), len_p5_R = floor(p2$p5),
    q25_pos_F = last_good(p1$med, Q_MIN), q25_pos_R = last_good(p2$med, Q_MIN),
    truncLen_F = tF, truncLen_R = tR, expected_overlap = overlap,
    minOverlap = if (overlap >= 20 + MARGIN) 20 else MIN_OVERLAP, flag = flag)

  png(file.path(outdir, paste0("quality_", key, ".png")), width = 1400, height = 600, res = 110)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 4))
  for (k in 1:2) {
    p <- if (k == 1) p1 else p2; t <- if (k == 1) tF else tR
    x <- seq_along(p$med)
    plot(x, p$med, type = "l", lwd = 2, ylim = c(0, 42), xlab = "Position", ylab = "Quality",
         main = paste0(key, " ", c("R1 (forward)", "R2 (reverse)")[k], " - ", nrow(d), " runs"))
    lines(x, p$q25, lty = 2); lines(x, p$q75, lty = 2)
    abline(h = Q_MIN, col = "grey50", lty = 3)
    abline(v = t, col = "red", lwd = 2)
    par(new = TRUE)
    plot(x, p$cover, type = "l", col = "steelblue", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
    axis(4, col.axis = "steelblue"); mtext("Fraction of reads this long", 4, 2.5, col = "steelblue")
    legend("bottomleft", bty = "n", cex = 0.8,
           legend = c("median", "quartiles", paste("truncLen =", t), "read length coverage"),
           lty = c(1, 2, 1, 1), col = c("black", "black", "red", "steelblue"))
  }
  dev.off()
  cat(key, ": truncLen =", tF, "/", tR, "| insert", L, "| overlap", overlap,
      if (nzchar(flag)) paste("|", flag) else "", "\n")
}

summ <- do.call(rbind, rows)
write.table(summ, file.path(outdir, "quality_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nPlots and quality_summary.tsv written to", outdir, "\n")
