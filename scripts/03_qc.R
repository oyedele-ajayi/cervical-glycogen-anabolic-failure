# Phase 3 - QC and outlier flagging
#
# Four diagnostics per cohort: inter-sample boxplots, density curves,
# sample-sample correlation heatmap, and PCA on the top 2000 most-variable
# genes (decision D2).
#
# Two pre-specified outlier criteria, applied to the discovery cohort:
#   C1  position outside the group cluster in PC1-PC2 space
#   C2  within-group correlation more than 3 SD below the group median
# Flagged samples are reported, not removed. The Phase 7 sensitivity analysis
# is what excludes them.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(pheatmap)

gene_mats <- readRDS("data/processed/gene_matrices.rds")
groups <- readRDS("data/processed/groups.rds")

N_VAR_GENES <- 2000
cohort_names <- names(gene_mats)

# colours, one per grade label across all cohorts
grade_colours <- c(Normal = "#4C72B0", CIN1 = "#55A868", CIN2 = "#C44E52",
                   CIN3 = "#8172B2", Cancer = "#CCB974",
                   LSIL = "#55A868", HSIL = "#C44E52", SCC = "#CCB974",
                   HGCIN = "#C44E52")


# ---------------------------------------------------------------------------
# 1. QC figures, one PNG per cohort
# ---------------------------------------------------------------------------

for (nm in cohort_names) {
  
  m <- gene_mats[[nm]]
  grp <- groups[[nm]]
  cols <- grade_colours[as.character(grp)]
  
  png(paste0("results/figures/qc/FigureS1_QC_", nm, ".png"),
      width = 2400, height = 2000, res = 180)
  par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
  
  # (A) boxplots. 128 boxes are unreadable, so subsample the discovery cohort
  # and say so in the legend.
  if (ncol(m) > 60) {
    set.seed(42)
    show <- sort(sample(seq_len(ncol(m)), 40))
    sub_title <- paste0("A. Expression by sample (", length(show), " of ", ncol(m), " shown)")
  } else {
    show <- seq_len(ncol(m))
    sub_title <- paste0("A. Expression by sample (all ", ncol(m), ")")
  }
  boxplot(m[, show], col = cols[show], outline = FALSE, las = 2,
          cex.axis = 0.5, main = sub_title, ylab = "log2 expression")
  
  # (B) density curves
  plot(density(m[, 1]), col = cols[1], lwd = 1,
       main = "B. Density", xlab = "log2 expression",
       ylim = c(0, max(apply(m, 2, function(x) max(density(x)$y)))))
  for (j in 2:ncol(m)) lines(density(m[, j]), col = cols[j], lwd = 1)
  legend("topright", legend = levels(grp),
         col = grade_colours[levels(grp)], lwd = 2, bty = "n", cex = 0.8)
  
  # (C) PCA on the top 2000 most-variable genes
  v <- apply(m, 1, var)
  top <- order(v, decreasing = TRUE)[1:min(N_VAR_GENES, nrow(m))]
  pca <- prcomp(t(m[top, ]), center = TRUE, scale. = FALSE)
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  plot(pca$x[, 1], pca$x[, 2], col = cols, pch = 19,
       xlab = paste0("PC1 (", pct[1], "%)"), ylab = paste0("PC2 (", pct[2], "%)"),
       main = paste0("C. PCA, top ", length(top), " variable genes"))
  legend("topright", legend = levels(grp),
         col = grade_colours[levels(grp)], pch = 19, bty = "n", cex = 0.8)
  
  # (D) sample-sample correlation
  cmat <- cor(m, method = "pearson")
  image(1:ncol(cmat), 1:ncol(cmat), cmat, col = hcl.colors(50, "YlOrRd", rev = TRUE),
        xlab = "sample", ylab = "sample", main = "D. Sample-sample correlation")
  
  dev.off()
  cat("wrote results/figures/qc/FigureS1_QC_", nm, ".png\n", sep = "")
}
cat("\n")


# ---------------------------------------------------------------------------
# 2. Outlier criteria on the discovery cohort
# ---------------------------------------------------------------------------

m <- gene_mats$GSE63514
grp <- groups$GSE63514
samples <- colnames(m)

# PCA, same input as the figure
v <- apply(m, 1, var)
top <- order(v, decreasing = TRUE)[1:N_VAR_GENES]
pca <- prcomp(t(m[top, ]), center = TRUE, scale. = FALSE)
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]

# C1: distance from the group centroid in PC1-PC2, z-scored within group
dist_z <- rep(NA, length(samples))
for (lv in levels(grp)) {
  idx <- which(grp == lv)
  cx <- mean(pc1[idx]); cy <- mean(pc2[idx])
  d <- sqrt((pc1[idx] - cx)^2 + (pc2[idx] - cy)^2)
  dist_z[idx] <- (d - mean(d)) / sd(d)
}

# C2: mean correlation to the rest of its own group, z-scored within group
cmat <- cor(m, method = "pearson")
corr_z <- rep(NA, length(samples))
for (lv in levels(grp)) {
  idx <- which(grp == lv)
  sub <- cmat[idx, idx]
  diag(sub) <- NA
  mc <- rowMeans(sub, na.rm = TRUE)
  corr_z[idx] <- (mc - median(mc)) / sd(mc)
}

qc <- data.frame(sample = samples,
                 group = as.character(grp),
                 PC1 = round(pc1, 3),
                 PC2 = round(pc2, 3),
                 PC1_rank = rank(pc1),
                 centroid_dist_z = round(dist_z, 3),
                 mean_corr_z = round(corr_z, 3),
                 stringsAsFactors = FALSE)

qc$flag_C1 <- qc$centroid_dist_z > 3
qc$flag_C2 <- qc$mean_corr_z < -3
qc$flagged <- qc$flag_C1 | qc$flag_C2

write.csv(qc, "results/tables/QC_discovery_metrics.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 3. CHECKPOINT 2b
# ---------------------------------------------------------------------------

cat("=== CHECKPOINT 2b: outlier flagging, discovery cohort ===\n\n")

cat("Samples flagged by the two pre-specified criteria:\n")
hits <- qc[qc$flagged, c("sample", "group", "PC1", "PC1_rank",
                         "centroid_dist_z", "mean_corr_z", "flag_C1", "flag_C2")]
if (nrow(hits) == 0) {
  cat("   none\n")
} else {
  print(hits, row.names = FALSE)
}

cat("\nThe five samples named in the current manuscript:\n")
manu <- c("GSM1551334", "GSM1551316", "GSM1551349", "GSM1551350", "GSM1551351")
print(qc[qc$sample %in% manu,
         c("sample", "group", "PC1", "PC1_rank", "centroid_dist_z",
           "mean_corr_z", "flagged")], row.names = FALSE)

cat("\nOpen question O1 - GSM1551331:\n")
print(qc[qc$sample == "GSM1551331",
         c("sample", "group", "PC1", "PC2", "PC1_rank", "centroid_dist_z",
           "mean_corr_z", "flagged")], row.names = FALSE)

cat("\nPC1 rank 1 is the most negative PC1. Lowest five PC1 overall:\n")
print(head(qc[order(qc$PC1), c("sample", "group", "PC1", "PC1_rank")], 5),
      row.names = FALSE)

cat("\nGroup counts, all cohorts:\n")
for (nm in cohort_names) {
  cat("  ", nm, ": ", paste(names(table(groups[[nm]])), table(groups[[nm]]),
                            sep = "=", collapse = "  "), "\n", sep = "")
}


# ---------------------------------------------------------------------------
# 4. All-probe PCA concordance check (logged, not a figure)
# ---------------------------------------------------------------------------

pca_all <- prcomp(t(m), center = TRUE, scale. = FALSE)
r <- cor(pca$x[, 1], pca_all$x[, 1], method = "spearman")

cat("\nAll-probe vs top-2000 PCA, PC1 Spearman correlation:", round(r, 4), "\n")
cat("(sign is arbitrary in PCA; magnitude near 1 means the two agree)\n")

writeLines(c(
  paste("top2000 PC1 vs all-gene PC1 Spearman:", round(r, 4)),
  paste("top2000 PC1 variance explained:", round(100 * pca$sdev[1]^2 / sum(pca$sdev^2), 2), "%"),
  paste("all-gene PC1 variance explained:", round(100 * pca_all$sdev[1]^2 / sum(pca_all$sdev^2), 2), "%")
), "logs/PCA_concordance_check.txt")

saveRDS(qc, "data/processed/QC_discovery.rds")

cat("\nWritten: results/tables/QC_discovery_metrics.csv,\n")
cat("  logs/PCA_concordance_check.txt, QC PNGs in results/figures/qc/\n")
