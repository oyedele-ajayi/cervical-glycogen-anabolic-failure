# Phase 2c - probe collapse by MaxMean selection  (decision D10)
#
# Replaces the avereps collapse. For each gene, keep the single probe with the
# highest mean expression across samples.
#
# Rationale (Miller et al. 2011, collapseRows; also the basis of jetset):
# probe failure from poor hybridisation kinetics, non-functional splice
# targets or 3' bias is common, while extreme cross-hybridisation producing a
# spuriously high signal is rare. So the highest-mean probe is the best single
# estimate of the dominant transcript. Averaging instead drags a gene toward
# background whenever its probes sit at very different expression levels
# (ENO1 carries probes at AveExpr 13.34 and 2.10), and that distortion feeds
# straight into the ssGSEA rank.
#
# Cost, stated in the limitations: MaxMean gives up the variance reduction of
# averaging genuine replicate probes. A systematic bias is worse than a modest
# precision loss, which is why we take this trade.
#
# Both matrices are written. The avereps version is kept so the supplement can
# report every headline statistic under each rule.

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

gene_set <- list(
  Synthesis = c("GYS1", "GYS2", "GBE1", "UGP2", "PGM1",
                "PPP1R3A", "PPP1R3B", "PPP1R3C", "PPP1R3D"),
  Degradation = c("PYGB", "PYGL", "PYGM", "AGL", "PHKA1",
                  "PHKA2", "PHKB", "PHKG1", "PHKG2"),
  Warburg = c("HK2", "PFKL", "PKM", "LDHA", "ENO1",
              "SLC2A1", "SLC2A3", "G6PD"),
  Differentiation = c("KRT4", "KRT13", "IVL", "SPRR1A", "TGM3", "CRNN", "ESR1"))
panel_genes <- unlist(gene_set, use.names = FALSE)


# ---------------------------------------------------------------------------
# 1. Collapse each cohort both ways
# ---------------------------------------------------------------------------

mats_maxmean <- list()
mats_avereps <- list()
probe_choices <- NULL

cat("=== Probe collapse ===\n\n")

for (nm in names(esets)) {
  
  eset <- esets[[nm]]
  m <- exprs(eset)
  
  sym <- trimws(sub("///.*$", "", as.character(fData(eset)[["Gene symbol"]])))
  sym[sym %in% c("", "NA", "---")] <- NA
  
  keep <- !is.na(sym)
  m <- m[keep, , drop = FALSE]
  sym <- sym[keep]
  probe_ids <- rownames(m)
  
  ave <- rowMeans(m)
  
  # MaxMean: for each gene, the probe with the highest mean expression.
  # order by gene then by descending mean, take the first row of each gene.
  ord <- order(sym, -ave)
  sym_ord <- sym[ord]
  first_of_gene <- !duplicated(sym_ord)
  pick <- ord[first_of_gene]
  
  m_max <- m[pick, , drop = FALSE]
  rownames(m_max) <- sym[pick]
  m_max <- m_max[order(rownames(m_max)), , drop = FALSE]
  
  # avereps for comparison
  m_avg <- m
  rownames(m_avg) <- sym
  m_avg <- avereps(m_avg, ID = rownames(m_avg))
  m_avg <- m_avg[order(rownames(m_avg)), , drop = FALSE]
  
  mats_maxmean[[nm]] <- m_max
  mats_avereps[[nm]] <- m_avg
  
  cat(sprintf("  %-16s probes %6d  genes %6d (MaxMean)  %6d (avereps)\n",
              nm, nrow(m), nrow(m_max), nrow(m_avg)))
  
  if (!identical(sort(rownames(m_max)), sort(rownames(m_avg))))
    stop("Gene sets differ between the two rules for ", nm, ". Investigate.")
  
  # record which probe was chosen for each panel gene
  for (g in intersect(panel_genes, sym)) {
    idx <- which(sym == g)
    chosen <- idx[which.max(ave[idx])]
    probe_choices <- rbind(probe_choices, data.frame(
      cohort = nm, gene = g,
      n_probes = length(idx),
      chosen_probe = probe_ids[chosen],
      chosen_AveExpr = round(ave[chosen], 2),
      lowest_AveExpr = round(min(ave[idx]), 2),
      AveExpr_spread = round(max(ave[idx]) - min(ave[idx]), 2),
      stringsAsFactors = FALSE))
  }
}

# stamped so the analysis scripts can tell these from the stage 02 matrices
attr(mats_maxmean, "collapse") <- "MaxMean"
attr(mats_avereps, "collapse") <- "avereps"

saveRDS(mats_maxmean, "data/processed/gene_matrices.rds")
saveRDS(mats_avereps, "data/processed/gene_matrices_avereps.rds")
write.csv(probe_choices, "results/tables/probe_selection_panel_genes.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 2. Which panel genes does the rule actually change?
# ---------------------------------------------------------------------------

cat("\n=== Panel genes with more than one probe, discovery cohort ===\n\n")

pc <- probe_choices[probe_choices$cohort == "GSE63514" & probe_choices$n_probes > 1, ]
pc <- pc[order(-pc$AveExpr_spread), ]
print(pc[, c("gene", "n_probes", "chosen_probe", "chosen_AveExpr",
             "lowest_AveExpr", "AveExpr_spread")], row.names = FALSE)

cat("\nGenes where the spread exceeds 4 log2 units are the ones averaging\n")
cat("would have diluted toward background:\n")
big <- pc$gene[pc$AveExpr_spread > 4]
cat("  ", ifelse(length(big) == 0, "none", paste(big, collapse = ", ")), "\n", sep = "")


# ---------------------------------------------------------------------------
# 3. Direct comparison of panel gene values under the two rules
# ---------------------------------------------------------------------------

cat("\n\n=== Panel gene expression under each rule, discovery cohort ===\n")
cat("(mean across all 128 samples)\n\n")

a <- mats_avereps$GSE63514
b <- mats_maxmean$GSE63514

cmp <- NULL
for (g in panel_genes) {
  if (g %in% rownames(a) && g %in% rownames(b)) {
    cmp <- rbind(cmp, data.frame(
      gene = g,
      module = names(gene_set)[sapply(gene_set, function(x) g %in% x)],
      avereps = round(mean(a[g, ]), 3),
      maxmean = round(mean(b[g, ]), 3),
      change = round(mean(b[g, ]) - mean(a[g, ]), 3),
      stringsAsFactors = FALSE))
  }
}
cmp <- cmp[order(-abs(cmp$change)), ]
print(head(cmp, 15), row.names = FALSE)

cat("\nMean absolute change by module:\n")
for (mod in names(gene_set)) {
  v <- abs(cmp$change[cmp$module == mod])
  cat(sprintf("   %-16s %.3f\n", mod, mean(v)))
}

write.csv(cmp, "results/tables/collapse_rule_gene_comparison.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 4. Record the decision
# ---------------------------------------------------------------------------

d <- read.csv("logs/decisions_record.csv", stringsAsFactors = FALSE)
d <- rbind(d, data.frame(id = "D10", decision = paste(
  "Probe collapse by MaxMean selection (highest mean expression per gene),",
  "replacing limma::avereps as specified in protocol section 5. Rationale:",
  "probe failure is far more common than extreme cross-hybridisation, so the",
  "highest-mean probe best represents the dominant transcript; averaging",
  "dilutes genes whose probes span wide expression ranges. Reference: Miller",
  "et al. 2011 (collapseRows). Caveat for limitations: MaxMean forgoes the",
  "variance reduction of averaging true replicate probes. Both rules are",
  "carried forward and every headline statistic is reported under each.")))
write.csv(d, "logs/decisions_record.csv", row.names = FALSE)

cat("\nWritten: gene_matrices.rds (MaxMean, primary),\n")
cat("  gene_matrices_avereps.rds (secondary),\n")
cat("  results/tables/probe_selection_panel_genes.csv,\n")
cat("  results/tables/collapse_rule_gene_comparison.csv\n")
cat("\nNow re-run 04, then 05, then 06.\n")
