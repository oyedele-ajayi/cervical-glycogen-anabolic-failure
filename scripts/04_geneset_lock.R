# Phase 4 - lock the 33-gene panel and map it per cohort
#
# Module order is Synthesis, Degradation, Warburg, Differentiation and is used
# everywhere downstream: heatmap rows, table rows, figure panels.
#
# Coverage expectations corrected at Checkpoint 3. Protocol section 5 claimed
# 26 of 33 for GSE138080. The true figure is 28. G6PD (probe 32535, Entrez
# 2539) and PKM (probes 6623, 9906, 33586, Entrez 5315) are present as sole
# annotations, not multi-mapping, so the earlier count was wrong. Five genes
# are genuinely absent from the GPL4133 annotation: GYS2, PPP1R3A, PHKG1,
# SLC2A3, TGM3.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(Biobase)
library(AnnotationDbi)
library(org.Hs.eg.db)

esets <- readRDS("data/processed/esets_log2only.rds")
gene_mats <- readRDS("data/processed/gene_matrices.rds")
# The analysis is defined on the MaxMean collapse from 02c. Running 02 alone
# leaves the stage 02 average collapse in this file, which silently changes
# every score downstream, so refuse rather than compute on the wrong matrices.
tag <- attr(gene_mats, "collapse")
if (!identical(tag, "MaxMean")) {
  stop("gene_matrices.rds holds the '", if (is.null(tag)) "unstamped" else tag,
       "' collapse, not MaxMean. Run 02c_collapse_maxmean.R first.")
}



# ---------------------------------------------------------------------------
# 1. The locked panel
# ---------------------------------------------------------------------------

gene_set <- list(
  Synthesis = c("GYS1", "GYS2", "GBE1", "UGP2", "PGM1",
                "PPP1R3A", "PPP1R3B", "PPP1R3C", "PPP1R3D"),
  Degradation = c("PYGB", "PYGL", "PYGM", "AGL", "PHKA1",
                  "PHKA2", "PHKB", "PHKG1", "PHKG2"),
  Warburg = c("HK2", "PFKL", "PKM", "LDHA", "ENO1",
              "SLC2A1", "SLC2A3", "G6PD"),
  Differentiation = c("KRT4", "KRT13", "IVL", "SPRR1A", "TGM3", "CRNN", "ESR1")
)

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
gene_set <- gene_set[module_order]

all_genes <- unlist(gene_set, use.names = FALSE)

cat("=== Locked gene set ===\n")
for (mod in module_order) {
  cat(sprintf("  %-16s n=%d  %s\n", mod, length(gene_set[[mod]]),
              paste(gene_set[[mod]], collapse = ", ")))
}
cat("  total:", length(all_genes), "genes\n")
if (any(duplicated(all_genes))) stop("A gene appears in more than one module.")
cat("  no duplicates across modules\n\n")

saveRDS(gene_set, "data/processed/gene_set_locked.rds")

panel_tab <- data.frame(gene = all_genes,
                        module = rep(module_order, times = sapply(gene_set, length)),
                        stringsAsFactors = FALSE)
write.csv(panel_tab, "results/tables/gene_set_locked.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 2. Entrez rescue attempt for GSE138080
# ---------------------------------------------------------------------------
# Probes carrying a Gene ID but no usable symbol, mapped through org.Hs.eg.db.
# A zero recovery is the expected result and is reported as a negative.

cat("=== Entrez rescue, GSE138080 ===\n")

fd <- fData(esets$GSE138080)
sym <- trimws(sub("///.*$", "", as.character(fd[["Gene symbol"]])))
entrez <- trimws(sub("///.*$", "", as.character(fd[["Gene ID"]])))

no_symbol <- (is.na(sym) | sym %in% c("", "NA", "---")) &
  !(is.na(entrez) | entrez %in% c("", "NA", "---"))

cat("  probes with Entrez ID but no symbol:", sum(no_symbol), "\n")

rescued_symbols <- character(0)
if (sum(no_symbol) > 0) {
  ids <- unique(entrez[no_symbol])
  mapped <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = ids,
                          keytype = "ENTREZID", columns = "SYMBOL"))
  rescued_symbols <- unique(na.omit(mapped$SYMBOL))
  cat("  unique Entrez IDs:", length(ids),
      "-> symbols recovered:", length(rescued_symbols), "\n")
}

rescued_panel <- intersect(rescued_symbols, all_genes)
cat("  of the 33 panel genes, recovered:", length(rescued_panel))
if (length(rescued_panel) > 0) cat("  (", paste(rescued_panel, collapse = ", "), ")")
cat("\n\n")


# ---------------------------------------------------------------------------
# 3. Per-cohort coverage
# ---------------------------------------------------------------------------

cohort_genesets <- list()
coverage <- NULL

cat("=== Per-cohort coverage ===\n")

for (nm in names(gene_mats)) {
  
  available <- rownames(gene_mats[[nm]])
  cs <- list()
  for (mod in module_order) {
    cs[[mod]] <- intersect(gene_set[[mod]], available)
  }
  cohort_genesets[[nm]] <- cs
  
  found <- unlist(cs, use.names = FALSE)
  missing <- setdiff(all_genes, found)
  
  cat("\n", nm, "\n", sep = "")
  for (mod in module_order) {
    cat(sprintf("   %-16s %d of %d\n", mod, length(cs[[mod]]), length(gene_set[[mod]])))
  }
  cat(sprintf("   total            %d of %d\n", length(found), length(all_genes)))
  if (length(missing) > 0) {
    cat("   missing: ", paste(sort(missing), collapse = ", "), "\n", sep = "")
  } else {
    cat("   missing: none\n")
  }
  
  coverage <- rbind(coverage,
                    data.frame(cohort = nm,
                               gene = panel_tab$gene,
                               module = panel_tab$module,
                               present = panel_tab$gene %in% found,
                               stringsAsFactors = FALSE))
}

write.csv(coverage, "results/tables/gene_set_coverage.csv", row.names = FALSE)
saveRDS(cohort_genesets, "data/processed/cohort_genesets.rds")


# ---------------------------------------------------------------------------
# 4. CHECKPOINT 3
# ---------------------------------------------------------------------------
# Expectations below are the corrected ones, verified probe by probe in 04b.

cat("\n\n=== CHECKPOINT 3 ===\n\n")

expected_total <- c(GSE63514 = 33, GSE138080 = 28,
                    GSE27678_GPL571 = 32, GSE27678_GPL570 = 33)

expected_missing <- list(
  GSE63514 = character(0),
  GSE138080 = c("GYS2", "PPP1R3A", "PHKG1", "SLC2A3", "TGM3"),
  GSE27678_GPL571 = c("PPP1R3B"),
  GSE27678_GPL570 = character(0))

ok <- TRUE
for (nm in names(expected_total)) {
  
  found <- unlist(cohort_genesets[[nm]], use.names = FALSE)
  missing <- sort(setdiff(all_genes, found))
  exp_miss <- sort(expected_missing[[nm]])
  
  n_flag <- if (length(found) == expected_total[nm]) "ok" else "MISMATCH"
  m_flag <- if (identical(missing, exp_miss)) "ok" else "MISMATCH"
  if (n_flag != "ok" || m_flag != "ok") ok <- FALSE
  
  cat(sprintf("%-16s total expected %2d observed %2d  %s\n",
              nm, expected_total[nm], length(found), n_flag))
  cat("                 missing expected: ",
      ifelse(length(exp_miss) == 0, "none", paste(exp_miss, collapse = ", ")),
      "\n", sep = "")
  cat("                 missing observed: ",
      ifelse(length(missing) == 0, "none", paste(missing, collapse = ", ")),
      "   ", m_flag, "\n", sep = "")
}

cat("\nEntrez rescue for GSE138080 recovered", length(rescued_panel), "of the 33.\n")

if (ok) {
  cat("\n>> CHECKPOINT 3 PASSED.\n\n")
} else {
  cat("\n>> CHECKPOINT 3 MISMATCH. Diagnose before Phase 5.\n\n")
}

cat("Written: gene_set_locked.rds/.csv, cohort_genesets.rds,\n")
cat("  results/tables/gene_set_coverage.csv\n")
