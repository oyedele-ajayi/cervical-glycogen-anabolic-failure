# Phase 2 - probe annotation, collapse to gene symbols, group assignment
#
# The protocol put the collapse in Phase 3, after QC. That does not work:
# the QC PCA runs on the top 2000 most-variable genes, so genes have to
# exist first. Collapse happens here.

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

esets <- readRDS("data/processed/esets_log2only.rds")


# ---------------------------------------------------------------------------
# 1. What annotation columns do we actually have?
# ---------------------------------------------------------------------------

cat("=== Feature annotation columns per cohort ===\n")
for (nm in names(esets)) {
  cols <- colnames(fData(esets[[nm]]))
  cat("\n", nm, " (", length(cols), " columns)\n", sep = "")
  cat("  ", paste(cols, collapse = " | "), "\n", sep = "")
}
cat("\n\n")


# ---------------------------------------------------------------------------
# 2. Collapse probes to gene symbols
# ---------------------------------------------------------------------------
# Affymetrix annot files use "Gene symbol", Agilent uses "GENE_SYMBOL".
# Find whichever is present rather than hard-coding it.

find_symbol_col <- function(eset) {
  cols <- colnames(fData(eset))
  hit <- cols[tolower(cols) %in% c("gene symbol", "gene_symbol", "symbol")]
  if (length(hit) == 0) {
    hit <- grep("symbol", cols, ignore.case = TRUE, value = TRUE)
  }
  if (length(hit) == 0) stop("No gene symbol column found. Columns: ", paste(cols, collapse = ", "))
  hit[1]
}

gene_mats <- list()

cat("=== Probe to gene collapse ===\n")
for (nm in names(esets)) {
  eset <- esets[[nm]]
  sym_col <- find_symbol_col(eset)
  symbols <- as.character(fData(eset)[[sym_col]])
  
  # Affymetrix multi-mapping probes look like "GYS1///GYS2". Take the first.
  symbols <- sub("///.*$", "", symbols)
  symbols <- trimws(symbols)
  symbols[symbols %in% c("", "NA", "---")] <- NA
  
  m <- exprs(eset)
  keep <- !is.na(symbols)
  m <- m[keep, , drop = FALSE]
  rownames(m) <- symbols[keep]
  
  # average probes that map to the same symbol
  collapsed <- avereps(m, ID = rownames(m))
  gene_mats[[nm]] <- collapsed
  
  cat(sprintf("  %-16s col '%s'  probes %6d -> annotated %6d -> genes %6d\n",
              nm, sym_col, nrow(exprs(eset)), nrow(m), nrow(collapsed)))
}
cat("\n")


# ---------------------------------------------------------------------------
# 3. Inventory of the sample annotation, so group regexes are not guesswork
# ---------------------------------------------------------------------------

cat("=== Sample annotation inventory ===\n")
for (nm in names(esets)) {
  pd <- pData(esets[[nm]])
  cat("\n---", nm, "---\n")
  for (field in c("title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1")) {
    if (field %in% colnames(pd)) {
      vals <- table(as.character(pd[[field]]))
      cat("  ", field, ":\n", sep = "")
      for (i in seq_along(vals)) {
        cat(sprintf("     %-60s n=%d\n", substr(names(vals)[i], 1, 60), vals[i]))
      }
    }
  }
}
cat("\n\n")


# ---------------------------------------------------------------------------
# 4. Assign groups
# ---------------------------------------------------------------------------
# Each cohort gets its own block. Regexes are applied to the whole pData row
# pasted together, so it does not matter which field carries the label.

pasted_annot <- function(eset) {
  pd <- pData(eset)
  fields <- intersect(c("title", "source_name_ch1", "characteristics_ch1",
                        "characteristics_ch1.1", "characteristics_ch1.2"),
                      colnames(pd))
  out <- rep("", nrow(pd))
  for (f in fields) out <- paste(out, as.character(pd[[f]]), sep = " | ")
  tolower(out)
}

groups <- list()

# --- GSE63514: Normal 24, CIN1 14, CIN2 22, CIN3 40, Cancer 28 -------------
a <- pasted_annot(esets$GSE63514)
g <- rep(NA_character_, length(a))
g[grepl("normal", a)] <- "Normal"
g[grepl("cin ?1|cin1", a)] <- "CIN1"
g[grepl("cin ?2|cin2", a)] <- "CIN2"
g[grepl("cin ?3|cin3", a)] <- "CIN3"
g[grepl("cancer|carcinoma|scc|tumor|tumour", a)] <- "Cancer"
groups$GSE63514 <- factor(g, levels = c("Normal", "CIN1", "CIN2", "CIN3", "Cancer"))

# --- GSE138080: Normal 10, HGCIN 15, SCC 10 --------------------------------
# These are histology on bulk tissue, so the series labels them Normal, HGCIN
# (high-grade CIN, meaning CIN2/3) and CxSCC. HSIL is a Bethesda cytology term
# and does not belong on histological specimens; it is kept only for GPL571.
# Five samples carry "cervical intraepithelial neoplasia, grade 2-3" in their
# characteristics rather than HGCIN in the title, so both patterns are matched.
a <- pasted_annot(esets$GSE138080)
g <- rep(NA_character_, length(a))
g[grepl("normal|healthy", a)] <- "Normal"
g[grepl("hgcin|high.grade|grade 2|grade 3|cin ?2|cin ?3", a)] <- "HGCIN"
g[grepl("scc|squamous cell carcinoma|cancer|carcinoma", a)] <- "SCC"
groups$GSE138080 <- factor(g, levels = c("Normal", "HGCIN", "SCC"))

# --- GSE27678 GPL571: Normal 12, LSIL 11, HSIL 21 --------------------------
a <- pasted_annot(esets$GSE27678_GPL571)
g <- rep(NA_character_, length(a))
g[grepl("normal", a)] <- "Normal"
g[grepl("lsil|low.grade|cin ?1", a)] <- "LSIL"
g[grepl("hsil|high.grade|cin ?2|cin ?3", a)] <- "HSIL"
groups$GSE27678_GPL571 <- factor(g, levels = c("Normal", "LSIL", "HSIL"))

# --- GSE27678 GPL570: Normal 3, Cancer 28 (cell lines already dropped) -----
a <- pasted_annot(esets$GSE27678_GPL570)
g <- rep(NA_character_, length(a))
g[grepl("normal", a)] <- "Normal"
g[grepl("scc|cancer|carcinoma|tumor|tumour", a)] <- "Cancer"
groups$GSE27678_GPL570 <- factor(g, levels = c("Normal", "Cancer"))


# ---------------------------------------------------------------------------
# 5. CHECKPOINT 2a - group counts against protocol section 4
# ---------------------------------------------------------------------------

expected <- list(
  GSE63514        = c(Normal = 24, CIN1 = 14, CIN2 = 22, CIN3 = 40, Cancer = 28),
  GSE138080       = c(Normal = 10, HGCIN = 15, SCC = 10),
  GSE27678_GPL571 = c(Normal = 12, LSIL = 11, HSIL = 21),
  GSE27678_GPL570 = c(Normal = 3,  Cancer = 28)
)

cat("=== CHECKPOINT 2a: group assignment ===\n")
all_ok <- TRUE

for (nm in names(groups)) {
  obs <- table(groups[[nm]])
  exp <- expected[[nm]]
  cat("\n", nm, "\n", sep = "")
  for (lv in names(exp)) {
    o <- as.integer(obs[lv])
    if (is.na(o)) o <- 0L
    flag <- if (o == exp[lv]) "ok" else "MISMATCH"
    cat(sprintf("   %-8s expected %3d  observed %3d   %s\n", lv, exp[lv], o, flag))
    if (o != exp[lv]) all_ok <- FALSE
  }
  n_na <- sum(is.na(groups[[nm]]))
  cat(sprintf("   unassigned: %d\n", n_na))
  if (n_na > 0) {
    all_ok <- FALSE
    cat("   -> unassigned samples:\n")
    bad <- sampleNames(esets[[nm]])[is.na(groups[[nm]])]
    pa <- pasted_annot(esets[[nm]])[is.na(groups[[nm]])]
    for (i in seq_along(bad)) cat("      ", bad[i], " :: ", substr(pa[i], 1, 90), "\n", sep = "")
  }
}

cat("\nGene-level matrix sizes:\n")
for (nm in names(gene_mats)) {
  cat(sprintf("   %-16s %6d genes x %3d samples\n",
              nm, nrow(gene_mats[[nm]]), ncol(gene_mats[[nm]])))
}

if (!all_ok) {
  cat("\n>> CHECKPOINT 2a FAILED. Fix the regexes above before continuing.\n")
  cat(">> The annotation inventory printed in section 3 shows the real strings.\n\n")
} else {
  cat("\n>> CHECKPOINT 2a PASSED.\n\n")
}


# ---------------------------------------------------------------------------
# 6. Save
# ---------------------------------------------------------------------------

# These matrices use the simple average collapse and are an intermediate only.
# Script 02c recollapses by MaxMean and overwrites this same file, so 02c must
# always follow 02. The stamp lets the analysis scripts refuse to run on these.
attr(gene_mats, "collapse") <- "avereps_stage02"

saveRDS(gene_mats, "data/processed/gene_matrices.rds")
saveRDS(groups, "data/processed/groups.rds")

cat("\n-------------------------------------------------------------\n")
cat("gene_matrices.rds now holds the STAGE 02 average collapse.\n")
cat("Run 02b, then 02c, before any analysis script. 02c replaces this\n")
cat("file with the MaxMean collapse the analysis is defined on.\n")
cat("-------------------------------------------------------------\n")

# one tidy file per cohort so nothing is ever retyped from the console
for (nm in names(groups)) {
  out <- data.frame(sample = colnames(gene_mats[[nm]]),
                    group = as.character(groups[[nm]]),
                    stringsAsFactors = FALSE)
  write.csv(out, paste0("results/tables/sample_groups_", nm, ".csv"), row.names = FALSE)
}

cat("Saved gene_matrices.rds, groups.rds and per-cohort group tables.\n")
