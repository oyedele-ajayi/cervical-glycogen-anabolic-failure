# Phase 8 - robustness battery
#
# (a) Outlier sensitivity. Four separate exclusion sets, never merged:
#       criteria4   the pre-specified statistical criteria
#       legacy5     the five named in the current manuscript (ad hoc)
#       posthoc1    GSM1551331 alone (signature-informed, post hoc)
#       posthoc2    GSM685524 in GPL570 (contradictory metadata, post hoc)
# (b) Leave-one-gene-out across all 33 genes.
# (c) Module dropout, each module alone as a CIN3+ classifier.
# (d) Reduced-index comparison, which settles whether Warburg earns its place.
#
# ssGSEA normalises against the sample set, so every subset recomputes from the
# subset matrix. Nothing here slices a previously computed score.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42)

library(GSVA)
library(pROC)

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
cohort_genesets <- readRDS("data/processed/cohort_genesets.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
N_BOOT <- 2000

m_disc <- gene_mats$GSE63514
grp_disc <- groups$GSE63514
gs_disc <- cohort_genesets$GSE63514[module_order]


# ---------------------------------------------------------------------------
# Helper: ssGSEA on a matrix, returning module scores plus GDI_diff
# ---------------------------------------------------------------------------

score_matrix <- function(m, gs) {
  set.seed(42)
  sc <- gsva(ssgseaParam(m, gs), verbose = FALSE)
  out <- as.data.frame(t(sc))
  out$GDI_diff <- (out$Warburg + out$Degradation) -
    (out$Synthesis + out$Differentiation)
  out
}

auc_for <- function(score, grp, positive) {
  outcome <- ifelse(as.character(grp) %in% positive, 1, 0)
  set.seed(42)
  r <- roc(outcome, score, quiet = TRUE, direction = "<")
  as.numeric(r$auc)
}

auc_ci_for <- function(score, grp, positive) {
  outcome <- ifelse(as.character(grp) %in% positive, 1, 0)
  set.seed(42)
  r <- roc(outcome, score, quiet = TRUE, direction = "<")
  set.seed(42)
  ci <- ci.auc(r, method = "bootstrap", boot.n = N_BOOT, boot.stratified = TRUE)
  c(auc = as.numeric(r$auc), lo = as.numeric(ci[1]), hi = as.numeric(ci[3]))
}

cin2plus <- c("CIN2", "CIN3", "Cancer")
cin3plus <- c("CIN3", "Cancer")
cancer_only <- c("Cancer")


# ---------------------------------------------------------------------------
# (a) Outlier sensitivity
# ---------------------------------------------------------------------------

exclusions <- list(
  full      = character(0),
  criteria4 = c("GSM1551339", "GSM1551365", "GSM1551395", "GSM1551419"),
  legacy5   = c("GSM1551334", "GSM1551316", "GSM1551349", "GSM1551350", "GSM1551351"),
  posthoc1  = c("GSM1551331")
)

cat("=== CHECKPOINT 7a: outlier sensitivity, discovery cohort ===\n\n")

sens <- NULL
for (lab in names(exclusions)) {
  drop <- exclusions[[lab]]
  keep <- !(colnames(m_disc) %in% drop)
  m <- m_disc[, keep]
  g <- grp_disc[keep]
  
  s <- score_matrix(m, gs_disc)
  
  a2 <- auc_ci_for(s$GDI_diff, g, cin2plus)
  a3 <- auc_ci_for(s$GDI_diff, g, cin3plus)
  ac <- auc_ci_for(s$GDI_diff, g, cancer_only)
  
  sens <- rbind(sens, data.frame(
    version = lab, n = ncol(m), n_excluded = length(drop),
    CIN2plus = round(a2["auc"], 4), CIN2_lo = round(a2["lo"], 4), CIN2_hi = round(a2["hi"], 4),
    CIN3plus = round(a3["auc"], 4), CIN3_lo = round(a3["lo"], 4), CIN3_hi = round(a3["hi"], 4),
    Cancer = round(ac["auc"], 4), Cancer_lo = round(ac["lo"], 4), Cancer_hi = round(ac["hi"], 4),
    stringsAsFactors = FALSE))
}
rownames(sens) <- NULL
print(sens[, c("version", "n", "CIN2plus", "CIN3plus", "Cancer")], row.names = FALSE)

cat("\nWith confidence intervals:\n")
print(sens, row.names = FALSE)

cat("\nCI half-widths for the full cohort, as the comparison yardstick:\n")
cat(sprintf("  CIN2+  %.4f\n", (sens$CIN2_hi[1] - sens$CIN2_lo[1]) / 2))
cat(sprintf("  CIN3+  %.4f\n", (sens$CIN3_hi[1] - sens$CIN3_lo[1]) / 2))
cat(sprintf("  Cancer %.4f\n", (sens$Cancer_hi[1] - sens$Cancer_lo[1]) / 2))
cat("Any shift smaller than its half-width is inside sampling noise.\n")

write.csv(sens, "results/tables/outlier_sensitivity.csv", row.names = FALSE)


# GPL570 post hoc: drop the metadata-ambiguous normal
cat("\n--- Post hoc, GPL570: excluding GSM685524 ---\n")
m570 <- gene_mats$GSE27678_GPL570
g570 <- groups$GSE27678_GPL570
gs570 <- cohort_genesets$GSE27678_GPL570[module_order]

s_all <- score_matrix(m570, gs570)
a_all <- auc_for(s_all$GDI_diff, g570, cancer_only)

keep <- colnames(m570) != "GSM685524"
s_cut <- score_matrix(m570[, keep], gs570)
a_cut <- auc_for(s_cut$GDI_diff, g570[keep], cancer_only)

cat(sprintf("  all 31 samples (3 normals): AUC %.4f\n", a_all))
cat(sprintf("  excluding GSM685524 (2 normals): AUC %.4f\n", a_cut))
cat("  Labelled post hoc. Never merged with the pre-specified analysis.\n")


# ---------------------------------------------------------------------------
# (b) Leave-one-gene-out
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 7b: leave-one-gene-out, CIN3+ AUC ===\n\n")

baseline_s <- score_matrix(m_disc, gs_disc)
baseline_auc <- auc_for(baseline_s$GDI_diff, grp_disc, cin3plus)
cat(sprintf("Baseline CIN3+ AUC: %.4f\n\n", baseline_auc))

loo <- NULL
for (mod in module_order) {
  for (g in gs_disc[[mod]]) {
    
    gs_mod <- gs_disc
    gs_mod[[mod]] <- setdiff(gs_mod[[mod]], g)
    
    if (length(gs_mod[[mod]]) == 0) {
      cat("  skipping", g, "- removal would empty", mod, "\n")
      next
    }
    
    s <- score_matrix(m_disc, gs_mod)
    a <- auc_for(s$GDI_diff, grp_disc, cin3plus)
    
    loo <- rbind(loo, data.frame(gene = g, module = mod,
                                 AUC = round(a, 4),
                                 delta = round(a - baseline_auc, 4),
                                 stringsAsFactors = FALSE))
  }
}

cat(sprintf("Genes tested: %d\n", nrow(loo)))
cat(sprintf("  mean AUC  %.4f\n", mean(loo$AUC)))
cat(sprintf("  SD        %.4f\n", sd(loo$AUC)))
cat(sprintf("  range     %.4f to %.4f\n", min(loo$AUC), max(loo$AUC)))

cat("\nFive largest drops:\n")
print(head(loo[order(loo$delta), ], 5), row.names = FALSE)

cat("\nFive largest gains (removal improves the classifier):\n")
print(head(loo[order(-loo$delta), ], 5), row.names = FALSE)

cat("\nPHKG1 and PPP1R3A, the two genes at background expression (D11):\n")
print(loo[loo$gene %in% c("PHKG1", "PPP1R3A"), ], row.names = FALSE)

write.csv(loo, "results/tables/leave_one_gene_out.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# (c) Module dropout
# ---------------------------------------------------------------------------
# Each module alone as a CIN3+ classifier. Synthesis and Differentiation are
# inverted, because for those low expression means high disease score.

cat("\n\n=== CHECKPOINT 7c: module dropout, CIN3+ AUC ===\n\n")

drop_tab <- NULL
for (mod in module_order) {
  v <- baseline_s[[mod]]
  if (mod %in% c("Synthesis", "Differentiation")) v <- -v
  a <- auc_ci_for(v, grp_disc, cin3plus)
  drop_tab <- rbind(drop_tab, data.frame(
    predictor = mod,
    inverted = mod %in% c("Synthesis", "Differentiation"),
    n_genes = length(gs_disc[[mod]]),
    AUC = round(a["auc"], 4),
    CI_lower = round(a["lo"], 4), CI_upper = round(a["hi"], 4),
    stringsAsFactors = FALSE))
}

a <- auc_ci_for(baseline_s$GDI_diff, grp_disc, cin3plus)
drop_tab <- rbind(drop_tab, data.frame(
  predictor = "GDI_diff (composite)", inverted = NA, n_genes = 33,
  AUC = round(a["auc"], 4), CI_lower = round(a["lo"], 4), CI_upper = round(a["hi"], 4),
  stringsAsFactors = FALSE))

rownames(drop_tab) <- NULL
print(drop_tab, row.names = FALSE)

write.csv(drop_tab, "results/tables/module_dropout.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# (d) Does Warburg earn its place?  (decision D14)
# ---------------------------------------------------------------------------
# O2 showed Warburg never separates from Normal at any grade. Test whether
# dropping it from the index changes discrimination, at all three thresholds
# and in all four cohorts.

cat("\n\n=== D14: reduced index comparison ===\n\n")

variants <- function(s) {
  list(
    full        = (s$Warburg + s$Degradation) - (s$Synthesis + s$Differentiation),
    no_warburg  = s$Degradation - (s$Synthesis + s$Differentiation),
    no_catab    = -(s$Synthesis + s$Differentiation),
    synth_only  = -s$Synthesis,
    diff_only   = -s$Differentiation
  )
}

cat("Discovery cohort, all three thresholds:\n\n")
vd <- variants(baseline_s)
d14 <- NULL
for (vn in names(vd)) {
  d14 <- rbind(d14, data.frame(
    index = vn,
    CIN2plus = round(auc_for(vd[[vn]], grp_disc, cin2plus), 4),
    CIN3plus = round(auc_for(vd[[vn]], grp_disc, cin3plus), 4),
    Cancer = round(auc_for(vd[[vn]], grp_disc, cancer_only), 4),
    stringsAsFactors = FALSE))
}
print(d14, row.names = FALSE)

cat("\nAll four cohorts, highest-grade threshold:\n\n")
top_pos <- list(GSE63514 = cancer_only, GSE138080 = c("SCC"),
                GSE27678_GPL571 = c("HSIL"), GSE27678_GPL570 = cancer_only)

cross <- NULL
for (nm in names(gene_mats)) {
  s <- score_matrix(gene_mats[[nm]], cohort_genesets[[nm]][module_order])
  vv <- variants(s)
  row <- data.frame(cohort = nm, stringsAsFactors = FALSE)
  for (vn in names(vv)) {
    row[[vn]] <- round(auc_for(vv[[vn]], groups[[nm]], top_pos[[nm]]), 4)
  }
  cross <- rbind(cross, row)
}
print(cross, row.names = FALSE)

write.csv(d14, "results/tables/D14_index_variants_discovery.csv", row.names = FALSE)
write.csv(cross, "results/tables/D14_index_variants_all_cohorts.csv", row.names = FALSE)

cat("\nIf no_warburg matches or beats full across the board, Warburg is a\n")
cat("passenger term and the index should be reduced. If full wins anywhere\n")
cat("by more than the CI half-width, it stays.\n")

cat("\nWritten: outlier_sensitivity.csv, leave_one_gene_out.csv,\n")
cat("  module_dropout.csv, D14_index_variants_*.csv\n")
