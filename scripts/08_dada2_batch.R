#!/usr/bin/env Rscript
# Step 08: DADA2 for one batch (one BioProject), from trimmed and oriented reads to a
# chimera-free ASV table.
#
# Usage (from the repository root, normally through scripts/08_dada2.sbatch):
#   Rscript scripts/08_dada2_batch.R <batch> [threads] [filtered_dir]
#
# Parameters (predoc, Chapter 1, unless noted):
#   filterAndTrim(maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE), truncLen per batch
#   learnErrors(nbases = 3e8) per error group: sequencing run (error_unit = seq_run) or the whole
#     BioProject (error_unit = bioproject); orientation sets A/B of B05 always separate.
#     Runs with unknown sequencing run join the largest sequencing run of their batch.
#     NovaSeq batches (error_model = novaseq): loess fit with monotonicity enforced.
#   dada(pool = "pseudo") per error group
#   mergePairs(maxMismatch = 0, minOverlap per batch, trimOverhang = TRUE)
#   removeBimeraDenovo(method = "consensus") on the whole batch
#
# Outputs in results/08_dada2/<batch>/:
#   seqtab.rds, seqtab_nochim.rds (runs x ASVs; sets A/B of a run summed), track.tsv (reads kept
#   at each step per run), asv_lengths.tsv, errors_<group>.pdf, summary.txt

# refuse to run on the login node (lab rule: no computation on iv12; use salloc --partition bio)
if (grepl("^iv", Sys.info()[["nodename"]])) {
  stop("you are on the login node (", Sys.info()[["nodename"]], "). Request a compute node first:\n",
       "  salloc --time=24:00:00 --cpus-per-task=4 --mem=8G --partition bio", call. = FALSE)
}

suppressPackageStartupMessages(library(dada2))

args     <- commandArgs(trailingOnly = TRUE)
BATCH    <- args[1]
THREADS  <- if (length(args) >= 2) as.integer(args[2]) else as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
FILT_DIR <- if (length(args) >= 3) args[3] else file.path(Sys.getenv("HOME"), "soy_meta_data", "filtered")
if (is.na(BATCH)) stop("usage: Rscript scripts/08_dada2_batch.R <batch> [threads] [filtered_dir]")

OUT <- file.path("results/08_dada2", BATCH)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FILT_DIR, BATCH), recursive = TRUE, showWarnings = FALSE)
set.seed(100)
t0 <- Sys.time()
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), BATCH, "-", ..., "\n")
msg("DADA2", as.character(packageVersion("dada2")), "|", THREADS, "threads")

# ---- inputs -------------------------------------------------------------------
man   <- read.delim("results/06_trim/trim_manifest.tsv", colClasses = "character")
man   <- man[man$batch == BATCH & man$r1 != "NA" & as.numeric(man$reads_out) > 0, ]
sheet <- read.delim("metadata/sample_sheet.tsv", colClasses = "character")
par   <- read.delim("config/dada2_params.tsv", colClasses = "character")
par   <- par[par$batch == BATCH, ]
if (nrow(par) != 1) stop("batch ", BATCH, " not found in config/dada2_params.tsv")
if (!nrow(man))     stop("no trimmed runs for batch ", BATCH)

# truncLen may be "274" or, for sets, "A:274;B:234"
pick <- function(v, set) {
  if (!grepl(":", v)) return(as.integer(v))
  kv <- do.call(rbind, strsplit(strsplit(v, ";")[[1]], ":"))
  as.integer(kv[kv[, 1] == set, 2])
}

man$set_key <- ifelse(man$set %in% c("A", "B"), man$set, "main")
man$sample  <- ifelse(man$set_key == "main", man$run, paste0(man$run, "_", man$set))
man$truncF  <- sapply(man$set_key, pick, v = par$truncLen_F)
man$truncR  <- sapply(man$set_key, pick, v = par$truncLen_R)
MINOV       <- as.integer(par$minOverlap)

# error groups
sr <- sheet$seq_run_id[match(man$run, sheet$run)]
sr[is.na(sr) | grepl("^NA", sr)] <- NA
if (par$error_unit == "seq_run" && any(!is.na(sr))) {
  sr[is.na(sr)] <- names(which.max(table(sr)))
} else {
  sr[] <- BATCH
}
man$err_group <- ifelse(man$set_key == "main", sr, paste0("set", man$set_key, "_", sr))
msg(nrow(man), "samples |", length(unique(man$err_group)), "error group(s):",
    paste(names(table(man$err_group)), table(man$err_group), sep = "=", collapse = ", "))

# ---- 1. filter and trim ---------------------------------------------------------
man$f1 <- file.path(FILT_DIR, BATCH, paste0(man$sample, "_F_filt.fastq.gz"))
man$f2 <- file.path(FILT_DIR, BATCH, paste0(man$sample, "_R_filt.fastq.gz"))
man$input <- NA_real_; man$filtered <- NA_real_
for (g in unique(paste(man$truncF, man$truncR))) {
  i <- which(paste(man$truncF, man$truncR) == g)
  msg("filterAndTrim truncLen =", man$truncF[i[1]], "/", man$truncR[i[1]], "|", length(i), "samples")
  r <- filterAndTrim(man$r1[i], man$f1[i], man$r2[i], man$f2[i],
                     truncLen = c(man$truncF[i[1]], man$truncR[i[1]]),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = THREADS)
  man$input[i] <- r[, 1]; man$filtered[i] <- r[, 2]
}
man$denoisedF <- 0; man$denoisedR <- 0; man$merged <- 0

# ---- error function for NovaSeq (binned quality scores) --------------------------
loess_monotone <- function(trans) {
  est <- loessErrfun(trans)
  nts <- c("A", "C", "G", "T")
  if (is.null(rownames(est))) rownames(est) <- paste0(rep(nts, each = 4), "2", nts)
  self <- rownames(est) %in% paste0(nts, "2", nts)
  # substitution rates must not increase with quality
  est[!self, ] <- t(apply(est[!self, , drop = FALSE], 1, function(x) rev(cummax(rev(x)))))
  for (nt in nts) {   # rows of each original nucleotide sum to 1
    sub <- startsWith(rownames(est), paste0(nt, "2")) & !self
    est[paste0(nt, "2", nt), ] <- 1 - colSums(est[sub, , drop = FALSE])
  }
  est
}
EFUN <- if (par$error_model == "novaseq") loess_monotone else loessErrfun
msg("error model:", par$error_model)

# ---- 2-4. errors, denoising and merging per error group --------------------------
getN <- function(x) sum(getUniques(x))
tables <- list()
for (g in unique(man$err_group)) {
  i <- which(man$err_group == g & man$filtered > 0 & file.exists(man$f1) & file.exists(man$f2))
  if (!length(i)) next
  msg("group", g, "|", length(i), "samples | learnErrors")
  errF <- learnErrors(man$f1[i], nbases = 3e8, errorEstimationFunction = EFUN,
                      multithread = THREADS, randomize = TRUE)
  errR <- learnErrors(man$f2[i], nbases = 3e8, errorEstimationFunction = EFUN,
                      multithread = THREADS, randomize = TRUE)
  pdf(file.path(OUT, paste0("errors_", gsub("[^A-Za-z0-9_.-]", "_", g), ".pdf")), width = 10, height = 8)
  print(plotErrors(errF, nominalQ = TRUE) + ggplot2::ggtitle(paste(BATCH, g, "forward")))
  print(plotErrors(errR, nominalQ = TRUE) + ggplot2::ggtitle(paste(BATCH, g, "reverse")))
  dev.off()

  msg("group", g, "| dada (pool = pseudo)")
  ddF <- dada(man$f1[i], err = errF, pool = "pseudo", multithread = THREADS)
  ddR <- dada(man$f2[i], err = errR, pool = "pseudo", multithread = THREADS)
  if (length(i) == 1) { ddF <- list(ddF); ddR <- list(ddR) }
  mg <- mergePairs(ddF, man$f1[i], ddR, man$f2[i],
                   maxMismatch = 0, minOverlap = MINOV, trimOverhang = TRUE)
  if (length(i) == 1) mg <- list(mg)
  names(mg) <- man$sample[i]

  man$denoisedF[i] <- sapply(ddF, getN)
  man$denoisedR[i] <- sapply(ddR, getN)
  man$merged[i]    <- sapply(mg, getN)
  tables[[g]] <- makeSequenceTable(mg)
}

# ---- 5. combine groups, sum sets A/B per run, remove chimeras ---------------------
seqtab <- if (length(tables) > 1) mergeSequenceTables(tables = tables) else tables[[1]]
seqtab <- rowsum(seqtab, man$run[match(rownames(seqtab), man$sample)])
msg("removeBimeraDenovo |", ncol(seqtab), "ASVs before")
nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = THREADS)

# ---- outputs ---------------------------------------------------------------------
saveRDS(seqtab, file.path(OUT, "seqtab.rds"))
saveRDS(nochim, file.path(OUT, "seqtab_nochim.rds"))

track <- aggregate(man[, c("input", "filtered", "denoisedF", "denoisedR", "merged")],
                   by = list(run = man$run), FUN = sum)
track$nonchim <- rowSums(nochim)[match(track$run, rownames(nochim))]
track$nonchim[is.na(track$nonchim)] <- 0
track$pct_kept <- round(100 * track$nonchim / track$input, 1)
write.table(track, file.path(OUT, "track.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

lens <- as.data.frame(table(length = nchar(colnames(nochim))))
lens$reads <- tapply(colSums(nochim), nchar(colnames(nochim)), sum)[as.character(lens$length)]
write.table(lens, file.path(OUT, "asv_lengths.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

top_len <- lens[order(-lens$reads), ][1:min(3, nrow(lens)), ]
summ <- c(
  paste("batch:", BATCH, "| bioproject:", par$bioproject),
  paste("runs:", nrow(track), "| error groups:", length(tables)),
  paste("reads: input", sum(track$input), "| filtered", sum(track$filtered),
        "| merged", sum(track$merged), "| non-chimeric", sum(track$nonchim)),
  paste0("kept from input: ", round(100 * sum(track$nonchim) / sum(track$input), 1), "% ",
         "(median per run ", median(track$pct_kept), "%)"),
  paste("ASVs:", ncol(seqtab), "before chimera removal |", ncol(nochim), "after",
        paste0("(", round(100 * sum(nochim) / sum(seqtab), 1), "% of merged reads non-chimeric)")),
  paste("most abundant ASV lengths:", paste0(top_len$length, " bp (",
        round(100 * top_len$reads / sum(nochim), 1), "% of reads)", collapse = ", ")),
  paste("runs with < 1000 non-chimeric reads:", sum(track$nonchim < 1000)),
  paste("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min"))
writeLines(summ, file.path(OUT, "summary.txt"))
cat("\n", paste(summ, collapse = "\n"), "\n")
