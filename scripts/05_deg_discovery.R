# Phase 5 - differential expression, discovery cohort (GSE63514)
#
# Seven contrasts. Run three times:
#   full       all 128 samples
#   criteria4  excluding the four samples the pre-specified criteria flagged
#   legacy5    excluding the five samples named in the current manuscript
#
# The exclusion sets stay separate and are never merged. The criteria set is
# the pre-specified sensitivity analysis. The legacy set is reported only for
# continuity with the earlier version.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(limma)

gene_mats <- readRDS("data/processed/gene_matrices.rds")
# The analysis is defined on the MaxMean collapse from 02c. Running 02 alone
# leaves the stage 02 average collapse in this file, which silently changes
# every score downstream, so refuse rather than compute on the wrong matrices.
tag <- attr(gene_mats, "collapse")
if (!identical(tag, "MaxMean")) {
  stop("gene_matrices.rds holds the '", if (is.null(tag)) "unstamped" else tag,
       "' collapse, not MaxMean. Run 02c_collapse_maxmean.R first.")
}

groups <- readRDS("data/processed/groups.rds")
gene_set <- readRDS("data/processed/gene_set_locked.rds")

m_all <- gene_mats$GSE63514
grp_all <- groups$GSE63514

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
panel <- data.frame(gene = unlist(gene_set, use.names = FALSE),
                    module = rep(module_order, times = sapply(gene_set, length)),
                    stringsAsFactors = FALSE)

criteria4 <- c("GSM1551339", "GSM1551365", "GSM1551395", "GSM1551419")
legacy5 <- c("GSM1551334", "GSM1551316", "GSM1551349", "GSM1551350", "GSM1551351")

contrast_defs <- c("CIN1 - Normal", "CIN2 - Normal", "CIN3 - Normal",
                   "Cancer - Normal", "CIN2 - CIN1", "CIN3 - CIN2",
                   "Cancer - CIN3")
contrast_names <- c("CIN1vNormal", "CIN2vNormal", "CIN3vNormal", "CancervNormal",
                    "CIN2vCIN1", "CIN3vCIN2", "CancervCIN3")


# ---------------------------------------------------------------------------
# Fit one version of the cohort
# ---------------------------------------------------------------------------

run_deg <- function(m, grp, label) {
  
  grp <- factor(as.character(grp),
                levels = c("Normal", "CIN1", "CIN2", "CIN3", "Cancer"))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  
  cont <- makeContrasts(contrasts = contrast_defs, levels = design)
  colnames(cont) <- contrast_names
  
  fit <- lmFit(m, design)
  fit <- contrasts.fit(fit, cont)
  fit <- eBayes(fit)
  
  out <- NULL
  for (cn in contrast_names) {
    tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none", adjust.method = "BH")
    tt$gene <- rownames(tt)
    tt$contrast <- cn
    tt$version <- label
    out <- rbind(out, tt)
  }
  out
}

deg_full <- run_deg(m_all, grp_all, "full")

keep <- !(colnames(m_all) %in% criteria4)
deg_c4 <- run_deg(m_all[, keep], grp_all[keep], "criteria4")

keep <- !(colnames(m_all) %in% legacy5)
deg_l5 <- run_deg(m_all[, keep], grp_all[keep], "legacy5")

cols_out <- c("gene", "contrast", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
write.csv(deg_full[, cols_out], "results/tables/DEG_discovery_full.csv", row.names = FALSE)
write.csv(deg_c4[, cols_out], "results/tables/DEG_discovery_criteria4.csv", row.names = FALSE)
write.csv(deg_l5[, cols_out], "results/tables/DEG_discovery_legacy5.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# CHECKPOINT 4a - transcriptome-wide counts
# ---------------------------------------------------------------------------

cat("=== CHECKPOINT 4a: transcriptome-wide significant counts ===\n")
cat("Genes tested:", nrow(m_all), "\n\n")

reference <- c(CIN1vNormal = 1, CIN2vNormal = 266, CIN3vNormal = 5492,
               CancervNormal = 8449, CIN2vCIN1 = 0, CIN3vCIN2 = 1049,
               CancervCIN3 = 1861)

counts <- data.frame(contrast = contrast_names,
                     n_FDR05 = 0, min_FDR = 0,
                     reference = as.integer(reference[contrast_names]),
                     stringsAsFactors = FALSE)

for (i in seq_along(contrast_names)) {
  sub <- deg_full[deg_full$contrast == contrast_names[i], ]
  counts$n_FDR05[i] <- sum(sub$adj.P.Val < 0.05)
  counts$min_FDR[i] <- signif(min(sub$adj.P.Val), 4)
}
counts$diff <- counts$n_FDR05 - counts$reference
counts$pct <- round(100 * counts$diff / pmax(counts$reference, 1), 1)

print(counts, row.names = FALSE)

fail <- abs(counts$pct) > 5 & abs(counts$diff) > 2
if (any(fail)) {
  cat("\n>> HARD STOP at Checkpoint 4. Out of tolerance:",
      paste(counts$contrast[fail], collapse = ", "), "\n")
  cat(">> limma is stable across these versions, so suspect the group\n")
  cat(">> assignment or the probe collapse, not the package.\n\n")
} else {
  cat("\n>> Counts within tolerance.\n\n")
}

cat("Counts across the three cohort versions:\n")
sens <- data.frame(contrast = contrast_names, full = 0, criteria4 = 0, legacy5 = 0)
for (i in seq_along(contrast_names)) {
  cn <- contrast_names[i]
  sens$full[i] <- sum(deg_full$adj.P.Val[deg_full$contrast == cn] < 0.05)
  sens$criteria4[i] <- sum(deg_c4$adj.P.Val[deg_c4$contrast == cn] < 0.05)
  sens$legacy5[i] <- sum(deg_l5$adj.P.Val[deg_l5$contrast == cn] < 0.05)
}
print(sens, row.names = FALSE)
write.csv(sens, "results/tables/DEG_counts_by_version.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# CHECKPOINT 4b - panel genes, Cancer vs Normal
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 4b: panel genes, Cancer vs Normal ===\n\n")

cvn <- deg_full[deg_full$contrast == "CancervNormal", ]
pg <- merge(panel, cvn[, c("gene", "logFC", "P.Value", "adj.P.Val")],
            by = "gene", all.x = TRUE)
pg$module <- factor(pg$module, levels = module_order)
pg <- pg[order(pg$module, -abs(pg$logFC)), ]

show <- data.frame(gene = pg$gene, module = as.character(pg$module),
                   log2FC = round(pg$logFC, 3),
                   P = signif(pg$P.Value, 3),
                   FDR = signif(pg$adj.P.Val, 3))
print(show, row.names = FALSE)

both <- pg[!is.na(pg$adj.P.Val) & pg$adj.P.Val < 0.05 & abs(pg$logFC) > 1, ]
fdr_only <- pg[!is.na(pg$adj.P.Val) & pg$adj.P.Val < 0.05, ]

cat("\nFDR < 0.05 and |log2FC| > 1 :", nrow(both), "genes\n")
cat("  ", paste(both$gene, collapse = ", "), "\n", sep = "")
cat("FDR < 0.05 alone             :", nrow(fdr_only), "genes\n")

cat("\nModule means of log2FC, Cancer vs Normal:\n")
for (mod in module_order) {
  v <- pg$logFC[pg$module == mod]
  cat(sprintf("   %-16s mean %+.4f   n=%d\n", mod, mean(v, na.rm = TRUE), sum(!is.na(v))))
}

write.csv(pg, "results/tables/panel_DEG_CancervNormal.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# CHECKPOINT 4c - panel genes across the three grade-to-grade contrasts
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 4c: panel genes, grade to grade ===\n")

g2g_all <- NULL
for (cn in c("CIN2vCIN1", "CIN3vCIN2", "CancervCIN3")) {
  sub <- deg_full[deg_full$contrast == cn, ]
  tmp <- merge(panel, sub[, c("gene", "logFC", "adj.P.Val")], by = "gene", all.x = TRUE)
  tmp$contrast <- cn
  tmp$module <- factor(tmp$module, levels = module_order)
  tmp <- tmp[order(tmp$module, -abs(tmp$logFC)), ]
  g2g_all <- rbind(g2g_all, tmp)
  
  cat("\n--- ", cn, " ---\n", sep = "")
  print(data.frame(gene = tmp$gene, module = as.character(tmp$module),
                   log2FC = round(tmp$logFC, 3), FDR = signif(tmp$adj.P.Val, 3)),
        row.names = FALSE)
}
write.csv(g2g_all, "results/tables/panel_DEG_grade_to_grade.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# CHECKPOINT 4d - panel genes moving more than 0.8 at CIN2 vs CIN1
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 4d: |log2FC| > 0.8 at CIN2 vs CIN1 ===\n")

c21 <- g2g_all[g2g_all$contrast == "CIN2vCIN1", ]
big <- c21[!is.na(c21$logFC) & abs(c21$logFC) > 0.8, ]
big <- big[order(-abs(big$logFC)), ]

if (nrow(big) == 0) {
  cat("  none\n")
} else {
  print(data.frame(gene = big$gene, module = as.character(big$module),
                   log2FC = round(big$logFC, 3), FDR = signif(big$adj.P.Val, 3)),
        row.names = FALSE)
}

cat("\nMinimum FDR at CIN2vCIN1, transcriptome-wide:",
    signif(min(deg_full$adj.P.Val[deg_full$contrast == "CIN2vCIN1"]), 4), "\n")

saveRDS(list(full = deg_full, criteria4 = deg_c4, legacy5 = deg_l5),
        "data/processed/DEG_discovery.rds")

cat("\nWritten: DEG_discovery_*.csv, panel_DEG_*.csv, DEG_discovery.rds\n")
