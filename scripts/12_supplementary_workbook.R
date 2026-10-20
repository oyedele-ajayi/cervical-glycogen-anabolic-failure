# Phase 12 - supplementary workbook
#
# One workbook, Contents sheet first, carrying a legend for every table.
# Every number that appears in the manuscript must come from here, never from
# console output.
#
# S2 and S3 are transcriptome-wide across seven contrasts, so the workbook is
# large and takes a minute or two to write.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)
library(Biobase)

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")

deg <- readRDS("data/processed/DEG_discovery.rds")
gene_set <- readRDS("data/processed/gene_set_locked.rds")
cohort_genesets <- readRDS("data/processed/cohort_genesets.rds")
gdi_list <- readRDS("data/processed/GDI_results.rds")
gene_mats <- readRDS("data/processed/gene_matrices.rds")
fg <- readRDS("data/processed/fgsea_results.rds")

panel_genes <- unlist(gene_set, use.names = FALSE)
panel_module <- rep(module_order, times = sapply(gene_set[module_order], length))


# ---------------------------------------------------------------------------
# S1 - panel composition, per-cohort availability, detection status
# ---------------------------------------------------------------------------

S1 <- data.frame(gene = panel_genes, module = panel_module,
                 stringsAsFactors = FALSE)

for (nm in names(cohort_genesets)) {
  present <- unlist(cohort_genesets[[nm]], use.names = FALSE)
  S1[[paste0("in_", nm)]] <- ifelse(S1$gene %in% present, "yes", "no")
}

# detection status: highest mean expression for that gene in each cohort
for (nm in names(gene_mats)) {
  m <- gene_mats[[nm]]
  S1[[paste0("meanExpr_", nm)]] <- round(
    sapply(S1$gene, function(g) if (g %in% rownames(m)) mean(m[g, ]) else NA), 2)
}
S1$detection_note <- ifelse(
  !is.na(S1$meanExpr_GSE63514) & S1$meanExpr_GSE63514 < 5,
  "at or near background in discovery cohort", "")


# ---------------------------------------------------------------------------
# S2, S3 - transcriptome-wide DEG
# ---------------------------------------------------------------------------

cols <- c("gene", "contrast", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
S2 <- deg$full[, cols]
S3 <- deg$criteria4[, cols]
S3b <- deg$legacy5[, cols]


# ---------------------------------------------------------------------------
# S4 - trend statistics
# ---------------------------------------------------------------------------

S4 <- read.csv("results/tables/GDI_trend_tests.csv", stringsAsFactors = FALSE)


# ---------------------------------------------------------------------------
# S5 - panel genes across all seven contrasts
# ---------------------------------------------------------------------------

S5 <- NULL
for (cn in unique(deg$full$contrast)) {
  sub <- deg$full[deg$full$contrast == cn & deg$full$gene %in% panel_genes, ]
  sub$module <- panel_module[match(sub$gene, panel_genes)]
  S5 <- rbind(S5, sub[, c("gene", "module", "contrast", "logFC", "P.Value", "adj.P.Val")])
}
S5$module <- factor(S5$module, levels = module_order)
S5 <- S5[order(S5$contrast, S5$module, S5$gene), ]


# ---------------------------------------------------------------------------
# S6 - leave one gene out
# ---------------------------------------------------------------------------

S6 <- read.csv("results/tables/leave_one_gene_out.csv", stringsAsFactors = FALSE)


# ---------------------------------------------------------------------------
# S7 - per-sample module scores and both GDI forms
# ---------------------------------------------------------------------------

S7 <- NULL
for (nm in names(gdi_list)) {
  s <- gdi_list[[nm]]
  S7 <- rbind(S7, data.frame(
    cohort = nm, sample = s$sample, grade = as.character(s$group),
    Synthesis = round(s$Synthesis, 5), Degradation = round(s$Degradation, 5),
    Warburg = round(s$Warburg, 5), Differentiation = round(s$Differentiation, 5),
    GDI_diff = round(s$GDI_diff, 5), GDI_ratio = round(s$GDI_ratio, 5),
    stringsAsFactors = FALSE))
}


# ---------------------------------------------------------------------------
# S8 - per-grade quartiles, both forms
# ---------------------------------------------------------------------------

S8 <- NULL
for (nm in names(gdi_list)) {
  s <- gdi_list[[nm]]
  for (lv in levels(s$group)) {
    for (form in c("GDI_diff", "GDI_ratio")) {
      v <- s[[form]][s$group == lv]
      q <- quantile(v, c(0, 0.25, 0.5, 0.75, 1))
      S8 <- rbind(S8, data.frame(
        cohort = nm, grade = lv, GDI_form = form, n = length(v),
        min = round(q[1], 4), Q1 = round(q[2], 4), median = round(q[3], 4),
        Q3 = round(q[4], 4), max = round(q[5], 4),
        mean = round(mean(v), 4), SD = round(sd(v), 4),
        stringsAsFactors = FALSE))
    }
  }
}
rownames(S8) <- NULL


# ---------------------------------------------------------------------------
# S9 - ROC, fragility, module dropout, index variants
# ---------------------------------------------------------------------------

S9a <- read.csv("results/tables/GDI_ROC.csv", stringsAsFactors = FALSE)
S9b <- read.csv("results/tables/ROC_fragility.csv", stringsAsFactors = FALSE)
S9c <- read.csv("results/tables/module_dropout.csv", stringsAsFactors = FALSE)
S9d <- read.csv("results/tables/D14_index_variants_discovery.csv", stringsAsFactors = FALSE)
S9e <- read.csv("results/tables/D14_index_variants_all_cohorts.csv", stringsAsFactors = FALSE)
S9f <- read.csv("results/tables/outlier_sensitivity.csv", stringsAsFactors = FALSE)


# ---------------------------------------------------------------------------
# S10 - pathway enrichment
# ---------------------------------------------------------------------------

S10a <- read.csv("results/tables/fgsea_modules_CancervNormal.csv", stringsAsFactors = FALSE)
S10b <- read.csv("results/tables/fgsea_direction_split.csv", stringsAsFactors = FALSE)
S10c <- read.csv("results/tables/fgsea_member_ranks.csv", stringsAsFactors = FALSE)


# ---------------------------------------------------------------------------
# S11 - probe collapse comparison, and the decision record
# ---------------------------------------------------------------------------

S11a <- read.csv("results/tables/probe_selection_panel_genes.csv", stringsAsFactors = FALSE)
S11b <- read.csv("results/tables/collapse_rule_gene_comparison.csv", stringsAsFactors = FALSE)
S12 <- read.csv("logs/decisions_record.csv", stringsAsFactors = FALSE)


# ---------------------------------------------------------------------------
# Contents
# ---------------------------------------------------------------------------

contents <- data.frame(
  Sheet = c("S1", "S2", "S3", "S3b", "S4", "S5", "S6", "S7", "S8",
            "S9a", "S9b", "S9c", "S9d", "S9e", "S9f",
            "S10a", "S10b", "S10c", "S11a", "S11b", "S12"),
  Legend = c(
    "The 33-gene panel by module, with availability and mean expression in each cohort. Genes at or near background in the discovery cohort are flagged; they are retained because the panel was pre-specified, and the leave-one-gene-out analysis (S6) shows their removal changes the CIN3+ AUC by less than 0.003.",
    "Transcriptome-wide differential expression, discovery cohort, all 21,755 genes across the seven contrasts. Full cohort (n = 128).",
    "Transcriptome-wide differential expression, discovery cohort excluding the four samples flagged by the pre-specified outlier criteria (n = 124). This is the pre-specified sensitivity analysis.",
    "Transcriptome-wide differential expression, discovery cohort excluding the five samples named in the earlier version of the manuscript (n = 123). These five were not derived from the stated criteria; this analysis is reported for continuity only and is not merged with S3.",
    "Trend statistics for both GDI forms in all four cohorts: Spearman rho with p, and Kruskal-Wallis H with degrees of freedom and p. GSE27678 GPL570 has two groups, so its Spearman is a two-group comparison rather than a trend test.",
    "Panel-gene differential expression across all seven contrasts, with log2 fold change, nominal P and FDR.",
    "Leave-one-gene-out analysis. Each of the 33 genes removed in turn, ssGSEA and GDI_diff recomputed from the reduced set, and the CIN3+ AUC re-estimated in the discovery cohort.",
    "Per-sample ssGSEA module scores and both GDI forms for all four cohorts. ssGSEA normalises within a sample set, so scores are cohort-internal and absolute values are not comparable between cohorts.",
    "Per-grade quartiles, mean and SD for both GDI forms in all four cohorts.",
    "ROC estimates for both GDI forms at every pre-specified threshold, with 2000-iteration stratified bootstrap 95% confidence intervals and group sizes.",
    "Fragility of each ROC estimate: total and discordant pairs, and the size of the smaller group. Estimates whose smaller group is under 10 cannot support a validation claim.",
    "Each module alone as a CIN3+ classifier in the discovery cohort, against the composite index. Synthesis and Differentiation are inverted so low expression yields a high disease score.",
    "Index variants compared at all three discovery thresholds: the full four-module index, the index without Warburg, without both catabolic modules, and each of Synthesis and Differentiation alone.",
    "The same index variants compared across all four cohorts at the highest-grade threshold.",
    "Outlier sensitivity. Three AUCs under the full cohort and under each exclusion set, with confidence intervals. The exclusion sets are kept separate and are never merged.",
    "fgsea enrichment of the four modules on the Cancer vs Normal ranked list, with NES, nominal p, FDR and leading-edge genes.",
    "Significant module members split by direction of change. The leading edge reported by fgsea contains only members moving in the direction of the NES sign, so it is not a list of a module's moving genes; this table gives that information.",
    "Rank position and ranking metric of every panel gene in the Cancer vs Normal ranked list, with leading-edge membership.",
    "The probe selected for each panel gene under MaxMean collapse, with the number of probes available and the spread of their mean expression.",
    "Panel gene expression under MaxMean collapse against averaging across all probes, showing which genes the collapse rule changes.",
    "Record of every analytical decision taken during the re-run, with its rationale."),
  stringsAsFactors = FALSE)


# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------

cat("Building workbook. S2 and S3 are large; this takes a minute.\n")

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", halign = "left")

add <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  freezePane(wb, name, firstRow = TRUE)
}

addWorksheet(wb, "Contents")
writeData(wb, "Contents", contents, headerStyle = hdr)
setColWidths(wb, "Contents", cols = 1:2, widths = c(10, 130))

add("S1", S1);    add("S2", S2);    add("S3", S3);    add("S3b", S3b)
add("S4", S4);    add("S5", S5);    add("S6", S6);    add("S7", S7)
add("S8", S8);    add("S9a", S9a);  add("S9b", S9b);  add("S9c", S9c)
add("S9d", S9d);  add("S9e", S9e);  add("S9f", S9f)
add("S10a", S10a); add("S10b", S10b); add("S10c", S10c)
add("S11a", S11a); add("S11b", S11b); add("S12", S12)

saveWorkbook(wb, "results/tables/Supplementary_Tables.xlsx", overwrite = TRUE)


# ---------------------------------------------------------------------------
# CHECKPOINT 10
# ---------------------------------------------------------------------------

cat("\n=== CHECKPOINT 10: sheet names and row counts ===\n\n")
sheets <- list(Contents = contents, S1 = S1, S2 = S2, S3 = S3, S3b = S3b,
               S4 = S4, S5 = S5, S6 = S6, S7 = S7, S8 = S8,
               S9a = S9a, S9b = S9b, S9c = S9c, S9d = S9d, S9e = S9e, S9f = S9f,
               S10a = S10a, S10b = S10b, S10c = S10c,
               S11a = S11a, S11b = S11b, S12 = S12)

for (nm in names(sheets)) {
  cat(sprintf("  %-9s %7d rows x %2d cols\n", nm,
              nrow(sheets[[nm]]), ncol(sheets[[nm]])))
}

cat("\nEvery supplementary item must be cited in the body text. Five of eight\n")
cat("were orphaned in the previous version.\n")
cat("\nWritten: results/tables/Supplementary_Tables.xlsx\n")
