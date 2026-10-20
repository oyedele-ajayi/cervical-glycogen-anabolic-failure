# =============================================================================
# GDI Cervical Cancer Study - Re-Run from scratch
# Phase 0 + 1: environment, seed, directory tree, cohort load
# Runs top to bottom in a clean R session.
# =============================================================================

rm(list = ls())

# The project root is wherever this is run from, which is the repository
# root when the .Rproj file is open.
PROJ <- normalizePath(".", winslash = "/")

DOWNLOAD_IF_MISSING <- TRUE   # set FALSE to force local-only

# --- 0.1 reproducibility -----------------------------------------------------
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42)
GDI_SEED <- 42L

# --- 0.2 directories ---------------------------------------------------------
for (d in c("data/raw", "data/processed", "results/tables",
            "results/figures/qc", "results/figures/deg",
            "results/figures/gdi", "results/figures/validation",
            "logs", "scripts")) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- 0.3 packages ------------------------------------------------------------
cran <- c("ggplot2","dplyr","patchwork","pheatmap","pROC","tibble","tidyr","readr","RColorBrewer")
bioc <- c("GEOquery","limma","GSVA","fgsea","Biobase","AnnotationDbi","org.Hs.eg.db")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
invisible(lapply(c(cran, bioc), function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))

env_tab <- data.frame(
  package = c("GEOquery","limma","GSVA","fgsea","pROC"),
  version = vapply(c("GEOquery","limma","GSVA","fgsea","pROC"),
                   function(p) as.character(packageVersion(p)), character(1)),
  row.names = NULL, stringsAsFactors = FALSE)

writeLines(capture.output(sessionInfo()), "logs/session_info.txt")

# --- 0.4 decision record -----------------------------------------------------
decisions <- data.frame(id = c("D1","D2","D3","D4","D5","D6"), decision = c(
  "GSE27678 GPL570: exclude cell lines GSM685543 (MS751) and GSM685549 (CaSki). n = 31 (Normal 3 / Cancer 28).",
  "QC PCA on top 2000 most-variable genes of the collapsed gene-level matrix. All-probe PCA logged as concordance check only.",
  "Bootstrap 2000 iterations, stratified. Methods 2.9 amended from 1000.",
  "set.seed(42), RNGkind declared, re-set before every stochastic call.",
  "Primary index GDI_diff. GDI_ratio secondary, sign pathology documented. Pre-specified.",
  "Protocol section 7 values are a tripwire, not a contract. Tolerances: AUC 0.02, rho 0.05, module median 0.02, DEG counts 5%, any fgsea sign flip or padj crossing 0.05. Exception: the 7 limma DEG counts are a hard stop at Checkpoint 4."),
  stringsAsFactors = FALSE)
write.csv(decisions, "logs/decisions_record.csv", row.names = FALSE)

# =============================================================================
# PHASE 1 - load cohorts from data/raw
# =============================================================================

cat("\n--- Contents of data/raw ---\n")
print(list.files("data/raw"))
cat("\n")

# Loader: .rds cache -> series matrix .txt.gz -> download. Always returns a list.
load_gse <- function(gse) {
  rds <- file.path("data/raw", paste0(gse, ".rds"))
  if (file.exists(rds)) {
    cat("  [rds]  ", gse, "\n"); obj <- readRDS(rds)
    return(if (is.list(obj)) obj else list(obj))
  }
  sm <- list.files("data/raw", pattern = paste0("^", gse, ".*series_matrix\\.txt\\.gz$"), full.names = TRUE)
  if (length(sm)) {
    cat("  [matrix]", gse, "-", length(sm), "file(s)\n")
    obj <- lapply(sm, function(f) getGEO(filename = f, getGPL = TRUE, AnnotGPL = TRUE))
    saveRDS(obj, rds); return(obj)
  }
  if (!DOWNLOAD_IF_MISSING) stop(gse, " not found in data/raw and download disabled.")
  cat("  [dl]   ", gse, "\n")
  obj <- getGEO(gse, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/raw")
  saveRDS(obj, rds); obj
}

cat("--- Loading cohorts ---\n")
raw63514  <- load_gse("GSE63514")
raw138080 <- load_gse("GSE138080")
raw27678  <- load_gse("GSE27678")
cat("\n")

one <- function(x, g) { if (length(x) != 1L) stop(g, " has ", length(x), " elements, expected 1."); x[[1]] }

# GSE27678 carries two platforms. Resolve by annotation(), never by index.
if (length(raw27678) != 2L) stop("GSE27678 has ", length(raw27678), " elements, expected 2.")
ann <- vapply(raw27678, annotation, character(1))
cat("GSE27678 annotations:", paste(sprintf("[[%d]]=%s", seq_along(ann), ann), collapse = "  "), "\n")
i570 <- which(grepl("GPL570", ann, fixed = TRUE)); i571 <- which(grepl("GPL571", ann, fixed = TRUE))
if (length(i570) != 1L || length(i571) != 1L) stop("Cannot uniquely resolve GPL570/GPL571. Inspect by hand.")

esets <- list(
  GSE63514        = one(raw63514,  "GSE63514"),
  GSE138080       = one(raw138080, "GSE138080"),
  GSE27678_GPL571 = raw27678[[i571]],
  GSE27678_GPL570 = raw27678[[i570]]
)

# --- 1.2 D1: drop the two in vitro cell lines --------------------------------
CELL_LINES <- c("GSM685543", "GSM685549")
ttl <- as.character(pData(esets$GSE27678_GPL570)$title)
names(ttl) <- sampleNames(esets$GSE27678_GPL570)
hits <- grep("cell.?line", ttl, ignore.case = TRUE, value = TRUE)

cat("\n--- D1 cell-line check (GPL570) ---\n")
if (!length(hits)) stop("No cell-line titles matched. STOP: inspect pData titles by hand before excluding anything.")
for (i in seq_along(hits)) cat(sprintf("  %s  %s\n", names(hits)[i], hits[i]))
if (!identical(sort(names(hits)), sort(CELL_LINES)))
  stop("Matched accessions differ from the expected pair. STOP and reconcile.")

esets$GSE27678_GPL570 <- esets$GSE27678_GPL570[, setdiff(sampleNames(esets$GSE27678_GPL570), CELL_LINES)]
write.csv(data.frame(sample = names(ttl), title = unname(ttl), excluded = names(ttl) %in% CELL_LINES),
          "logs/GPL570_titles_D1.csv", row.names = FALSE)

# --- 1.3 scale diagnostic ----------------------------------------------------
# A linear-scale series matrix passes every structural check, then breaks
# limma's variance model and skews ssGSEA ranks without ever erroring.
scale_tab <- do.call(rbind, Map(function(e, nm) {
  m <- exprs(e); q <- quantile(m, c(0, .25, .5, .75, 1), na.rm = TRUE)
  data.frame(cohort = nm, features = nrow(m), samples = ncol(m),
             min = round(q[1],2), median = round(q[3],2), max = round(q[5],2),
             NAs = sum(is.na(m)),
             verdict = if (q[5] > 100) "LINEAR - STOP" else if (q[5] > 25) "CHECK" else "log OK",
             stringsAsFactors = FALSE)
}, esets, names(esets)))

# --- 1.4 CHECKPOINT ----------------------------------------------------------
expected <- c(GSE63514 = 128, GSE138080 = 35, GSE27678_GPL571 = 44, GSE27678_GPL570 = 31)
n_tab <- data.frame(cohort = names(expected), expected = unname(expected),
                    observed = vapply(esets[names(expected)], ncol, integer(1)),
                    row.names = NULL)
n_tab$match <- ifelse(n_tab$expected == n_tab$observed, "MATCH", "MISMATCH")

cat("\n================ CHECKPOINT 0 + 1 ================\n")
cat(R.version.string, "| Bioconductor", as.character(BiocManager::version()), "\n")
print(env_tab, row.names = FALSE)
cat("\nPlatforms:\n")
for (nm in names(esets)) cat(sprintf("  %-16s %s\n", nm, annotation(esets[[nm]])))
cat("\nDimensions and expression scale:\n"); print(scale_tab, row.names = FALSE)
cat("\nSample counts vs protocol section 4 (post-D1):\n"); print(n_tab, row.names = FALSE)
cat("\npData fields for group assignment:\n")
for (nm in names(esets)) {
  cl <- grep("title|characteristic|source", colnames(pData(esets[[nm]])), ignore.case = TRUE, value = TRUE)
  cat(sprintf("  %-16s %s\n", nm, paste(head(cl, 6), collapse = " | ")))
}
cat("==================================================\n\n")

if (any(n_tab$match == "MISMATCH") || any(scale_tab$verdict != "log OK"))
  cat(">> CHECKPOINT FAILED. Do not proceed to QC.\n\n")

saveRDS(esets, "data/processed/esets_raw.rds")
write.csv(scale_tab, "logs/checkpoint1_scale_dims.csv", row.names = FALSE)
cat("Done. Written: data/processed/esets_raw.rds, logs/\n")


# Phase 1b - fix the GSE138080 scale
# GSE138080 came in on a linear scale (median 110, max 650658).
# Everything else is already log. Transform it before any analysis.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(Biobase)
library(limma)

esets <- readRDS("data/processed/esets_raw.rds")
e <- esets$GSE138080


# 1. Check what platform this actually is.
# Table 1 calls it SurePrint G3 8x60K, but 27460 features looks like the 44K.

cat("--- GPL4133 identity ---\n")
cat("annotation :", annotation(e), "\n")
cat("features   :", nrow(e), "\n")
cat("samples    :", ncol(e), "\n")
cat("source     :", as.character(pData(e)$source_name_ch1[1]), "\n")

ch2 <- grep("_ch2$", colnames(pData(e)))
cat("ch2 fields :", length(ch2), " (0 means single channel)\n\n")


# 2. Log transform. Minimum is 50 so no offset is needed.

raw <- exprs(e)
logged <- log2(raw)
quantiled <- normalizeBetweenArrays(logged, method = "quantile")


# 3. Compare the two versions across the 35 arrays.

med1 <- apply(logged, 2, median)
med2 <- apply(quantiled, 2, median)
iqr1 <- apply(logged, 2, IQR)
iqr2 <- apply(quantiled, 2, IQR)

cat("--- Per-sample distribution ---\n")
cat("log2 only       median range:", round(max(med1) - min(med1), 3),
    " IQR", round(min(iqr1), 2), "to", round(max(iqr1), 2),
    " overall", round(min(logged), 2), "to", round(max(logged), 2), "\n")
cat("log2 + quantile median range:", round(max(med2) - min(med2), 3),
    " IQR", round(min(iqr2), 2), "to", round(max(iqr2), 2),
    " overall", round(min(quantiled), 2), "to", round(max(quantiled), 2), "\n\n")

cat("If the log2-only median range is under 0.5 the submitter already\n")
cat("normalised between arrays and quantile adds nothing.\n\n")


# 4. Save both versions. We pick one at the checkpoint.

esets_log <- esets
exprs(esets_log$GSE138080) <- logged
saveRDS(esets_log, "data/processed/esets_log2only.rds")

esets_qn <- esets
exprs(esets_qn$GSE138080) <- quantiled
saveRDS(esets_qn, "data/processed/esets_log2_quantile.rds")


# 5. Re-run the scale check on all four cohorts.

cat("--- Scale check after transform ---\n")
for (nm in names(esets_log)) {
  m <- exprs(esets_log[[nm]])
  top <- max(m, na.rm = TRUE)
  if (top > 100) {
    verdict <- "LINEAR - STOP"
  } else if (top > 25) {
    verdict <- "CHECK"
  } else {
    verdict <- "log OK"
  }
  cat(sprintf("  %-16s min %7.2f  median %6.2f  max %7.2f   %s\n",
              nm, min(m, na.rm = TRUE), median(m, na.rm = TRUE), top, verdict))
}
cat("\n")


# 6. Record the decision.

d <- read.csv("logs/decisions_record.csv", stringsAsFactors = FALSE)
new <- data.frame(id = "D7",
                  decision = "GSE138080 was deposited on a linear scale. log2 transformed before analysis. Between-array normalisation decided at the Phase 1b checkpoint. Not in the original protocol.")
d <- rbind(d, new)
write.csv(d, "logs/decisions_record.csv", row.names = FALSE)

cat("Saved esets_log2only.rds and esets_log2_quantile.rds\n")
