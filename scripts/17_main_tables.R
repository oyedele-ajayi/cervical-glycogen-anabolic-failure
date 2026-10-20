# Phase 17 - the three main-text tables
#
# Table 1  cohort characteristics
# Table 2  the 33-gene panel, per-cohort availability, and the Cancer vs
#          Normal effect in the discovery cohort
# Table 3  trend statistics and discrimination, with group sizes attached
#
# Each is written to CSV and also printed as markdown, so the manuscript
# tables are produced from the analysis files rather than retyped.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(Biobase)

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")

deg <- readRDS("data/processed/DEG_discovery.rds")
gene_set <- readRDS("data/processed/gene_set_locked.rds")
cohort_genesets <- readRDS("data/processed/cohort_genesets.rds")
gene_mats <- readRDS("data/processed/gene_matrices.rds")
groups <- readRDS("data/processed/groups.rds")
ranks <- read.csv("results/tables/fgsea_member_ranks.csv", stringsAsFactors = FALSE)

panel_genes <- unlist(gene_set[module_order], use.names = FALSE)
panel_module <- rep(module_order, times = sapply(gene_set[module_order], length))


# helper: a data frame as a markdown pipe table
as_md <- function(df) {
  df[] <- lapply(df, as.character)
  df[is.na(df)] <- ""
  head <- paste("|", paste(colnames(df), collapse = " | "), "|")
  rule <- paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|")
  body <- apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))
  c(head, rule, body)
}


# ---------------------------------------------------------------------------
# Table 1 - cohorts
# ---------------------------------------------------------------------------

t1 <- data.frame(
  Role = c("Discovery", "Replication 1", "Replication 2", "Concordance check"),
  Accession = c("GSE63514", "GSE138080", "GSE27678", "GSE27678"),
  Platform = c("Affymetrix HG-U133 Plus 2.0 (GPL570)",
               "Agilent-014850 Whole Human Genome 4x44K (GPL4133)",
               "Affymetrix HG-U133A 2.0 (GPL571)",
               "Affymetrix HG-U133 Plus 2.0 (GPL570)"),
  Preparation = c("Laser-capture-microdissected epithelium",
                  "Bulk tissue", "Exfoliated cytology", "Exfoliated cytology"),
  stringsAsFactors = FALSE)

t1$Composition <- sapply(c("GSE63514", "GSE138080", "GSE27678_GPL571",
                           "GSE27678_GPL570"), function(nm) {
                             tb <- table(groups[[nm]])
                             paste(paste0(names(tb), " ", as.integer(tb)), collapse = ", ")
                           })
t1$n <- sapply(c("GSE63514", "GSE138080", "GSE27678_GPL571",
                 "GSE27678_GPL570"), function(nm) length(groups[[nm]]))
t1$Genes <- sapply(c("GSE63514", "GSE138080", "GSE27678_GPL571",
                     "GSE27678_GPL570"), function(nm) nrow(gene_mats[[nm]]))
t1$`Panel coverage` <- sapply(c("GSE63514", "GSE138080", "GSE27678_GPL571",
                                "GSE27678_GPL570"), function(nm)
                                  sprintf("%d of 33", length(unlist(cohort_genesets[[nm]], use.names = FALSE))))

write.csv(t1, "results/tables/Table1_cohorts.csv", row.names = FALSE)

cat("=== TABLE 1 ===\n\n")
writeLines(as_md(t1))
cat("\nTotal specimens:", sum(t1$n), "\n\n\n")


# ---------------------------------------------------------------------------
# Table 2 - the panel
# ---------------------------------------------------------------------------

cvn <- deg$full[deg$full$contrast == "CancervNormal", ]

t2 <- data.frame(Module = panel_module, Gene = panel_genes,
                 stringsAsFactors = FALSE)

short <- c(GSE63514 = "GSE63514", GSE138080 = "GSE138080",
           GSE27678_GPL571 = "GPL571", GSE27678_GPL570 = "GPL570")
for (nm in names(short)) {
  present <- unlist(cohort_genesets[[nm]], use.names = FALSE)
  t2[[short[[nm]]]] <- ifelse(t2$Gene %in% present, "yes", "-")
}

i <- match(t2$Gene, cvn$gene)
t2$log2FC <- sprintf("%+.3f", cvn$logFC[i])
t2$FDR <- ifelse(cvn$adj.P.Val[i] < 0.001,
                 formatC(cvn$adj.P.Val[i], format = "e", digits = 2),
                 sprintf("%.3f", cvn$adj.P.Val[i]))
t2$Rank <- ranks$rank[match(t2$Gene, ranks$gene)]
t2$Sig <- ifelse(cvn$adj.P.Val[i] < 0.05,
                 ifelse(abs(cvn$logFC[i]) > 1, "**", "*"), "")

# Detectability is a question about the tissue before disease, so this is the
# mean in the NORMAL group, not across all specimens. Averaging over all 128
# would drag GYS2 and TGM3 below any threshold simply because they are floored
# in cancer, and would wrongly mark them as undetected.
norm_idx <- groups$GSE63514 == "Normal"
t2$`MeanExpr (Normal)` <- sprintf("%.2f", sapply(t2$Gene, function(g)
  mean(gene_mats$GSE63514[g, norm_idx])))
t2$`MeanExpr (Cancer)` <- sprintf("%.2f", sapply(t2$Gene, function(g)
  mean(gene_mats$GSE63514[g, groups$GSE63514 == "Cancer"])))
t2$Detected <- ifelse(
  as.numeric(t2$`MeanExpr (Normal)`) < 5, "no", "yes")

write.csv(t2, "results/tables/Table2_panel.csv", row.names = FALSE)

cat("=== TABLE 2 ===\n\n")
writeLines(as_md(t2))
cat("\nFootnote: log2FC, FDR and Rank are Cancer vs Normal in the discovery\n")
cat("cohort across 21,755 genes; rank 1 is the most up-regulated. Mean log2\n")
cat("expression is given separately for the normal and cancer groups; a gene\n")
cat("is marked undetected where its mean in normal epithelium is below 5, the\n")
cat("array background, since detectability is a property of the tissue before\n")
cat("disease. ** FDR < 0.05 and |log2FC| > 1; * FDR < 0.05 only. Dash indicates\n")
cat("the gene is absent from that platform.\n\n")

cat("Counts: ", sum(t2$Sig == "**"), " genes at ** and ",
    sum(t2$Sig != ""), " at FDR < 0.05.\n", sep = "")
cat("Undetected in normal epithelium: ",
    paste(t2$Gene[t2$Detected == "no"], collapse = ", "), "\n", sep = "")
cat("(these should be PHKG1 and PPP1R3A only; GYS2 and TGM3 must NOT appear,\n")
cat(" since both are expressed in normal epithelium and lost in cancer)\n\n\n")


# ---------------------------------------------------------------------------
# Table 3 - trend and discrimination
# ---------------------------------------------------------------------------

trend <- read.csv("results/tables/GDI_trend_tests.csv", stringsAsFactors = FALSE)
roc <- read.csv("results/tables/GDI_ROC.csv", stringsAsFactors = FALSE)
frag <- read.csv("results/tables/ROC_fragility.csv", stringsAsFactors = FALSE)

trend <- trend[trend$GDI_form == "GDI_diff", ]
roc <- roc[roc$GDI_form == "GDI_diff", ]

t3 <- data.frame(Cohort = roc$cohort, Threshold = roc$threshold,
                 stringsAsFactors = FALSE)
t3$AUC <- sprintf("%.3f", roc$AUC)
t3$`95% CI` <- sprintf("%.3f to %.3f", roc$CI_lower, roc$CI_upper)
t3$`n pos` <- roc$n_pos
t3$`n neg` <- roc$n_neg

k <- match(t3$Cohort, trend$cohort)
t3$`Spearman rho` <- sprintf("%.3f", trend$Spearman_rho[k])
t3$`rho p` <- formatC(trend$Spearman_p[k], format = "e", digits = 2)
t3$`KW H (df)` <- sprintf("%.2f (%d)", trend$KW_H[k], trend$KW_df[k])
t3$`KW p` <- formatC(trend$KW_p[k], format = "e", digits = 2)

# trend statistics are per cohort, so show them once per cohort
dup <- duplicated(t3$Cohort)
t3$`Spearman rho`[dup] <- ""
t3$`rho p`[dup] <- ""
t3$`KW H (df)`[dup] <- ""
t3$`KW p`[dup] <- ""
t3$Cohort[dup] <- ""

write.csv(t3, "results/tables/Table3_trend_and_ROC.csv", row.names = FALSE)

cat("=== TABLE 3 ===\n\n")
writeLines(as_md(t3))
cat("\nFootnote: all values are for GDI_diff, the pre-specified primary form.\n")
cat("Confidence intervals are from 2,000 stratified bootstrap resamples. The\n")
cat("GPL570 set contains two groups, so its Spearman coefficient is a\n")
cat("two-group comparison rather than a trend test, and with three normal\n")
cat("specimens it is reported as a concordance check and not as a validation.\n")
cat("The ratio form is in Supplementary Table S9a.\n\n")

cat("=== Smaller group size per estimate (fragility) ===\n")
print(frag[, c("cohort", "threshold", "smaller_group", "discordant_pairs")],
      row.names = FALSE)

cat("\nWritten: Table1_cohorts.csv, Table2_panel.csv, Table3_trend_and_ROC.csv\n")
