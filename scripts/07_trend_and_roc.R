# Phase 7 - trend tests and ROC
#
# Trend: Spearman (exact = FALSE) and Kruskal-Wallis of each GDI form against
# ordinal grade, all four cohorts. H and df are written out, not just p.
#
# ROC: pROC::roc, stratified bootstrap with 2000 iterations (decision D3).
#
# GSE27678 GPL570 has two groups, so its Spearman is a two-group comparison
# and not a trend test. It is reported but must be described that way.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42)

library(pROC)

gdi_list <- readRDS("data/processed/GDI_results.rds")

N_BOOT <- 2000
cohorts <- names(gdi_list)


# ---------------------------------------------------------------------------
# 1. Trend tests
# ---------------------------------------------------------------------------

trend <- NULL

for (nm in cohorts) {
  s <- gdi_list[[nm]]
  ord <- as.integer(s$group)          # ordinal grade, 1..k
  
  for (form in c("GDI_diff", "GDI_ratio")) {
    v <- s[[form]]
    
    sp <- cor.test(v, ord, method = "spearman", exact = FALSE)
    kw <- kruskal.test(v ~ s$group)
    
    trend <- rbind(trend, data.frame(
      cohort = nm,
      GDI_form = form,
      n = length(v),
      n_groups = length(levels(s$group)),
      Spearman_rho = round(as.numeric(sp$estimate), 4),
      Spearman_p = signif(sp$p.value, 4),
      KW_H = round(as.numeric(kw$statistic), 4),
      KW_df = as.integer(kw$parameter),
      KW_p = signif(kw$p.value, 4),
      stringsAsFactors = FALSE))
  }
}

cat("=== CHECKPOINT 6a: trend tests ===\n\n")
print(trend, row.names = FALSE)
cat("\nNote: GSE27678_GPL570 has 2 groups. Its Spearman is a two-group\n")
cat("comparison, not a trend test, and must be described as such in the text.\n\n")

write.csv(trend, "results/tables/GDI_trend_tests.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 2. ROC thresholds
# ---------------------------------------------------------------------------
# Each entry: which grades count as positive for that threshold.

thresholds <- list(
  list(cohort = "GSE63514",        label = "CIN2+",       pos = c("CIN2", "CIN3", "Cancer")),
  list(cohort = "GSE63514",        label = "CIN3+",       pos = c("CIN3", "Cancer")),
  list(cohort = "GSE63514",        label = "Cancer",      pos = c("Cancer")),
  list(cohort = "GSE138080",       label = "HGCIN+",      pos = c("HGCIN", "SCC")),
  list(cohort = "GSE138080",       label = "SCC",         pos = c("SCC")),
  list(cohort = "GSE27678_GPL571", label = "LSIL+",       pos = c("LSIL", "HSIL")),
  list(cohort = "GSE27678_GPL571", label = "HSIL",        pos = c("HSIL")),
  list(cohort = "GSE27678_GPL570", label = "Cancer",      pos = c("Cancer"))
)


# ---------------------------------------------------------------------------
# 3. ROC estimates
# ---------------------------------------------------------------------------

roc_tab <- NULL
roc_objects <- list()

for (th in thresholds) {
  s <- gdi_list[[th$cohort]]
  outcome <- ifelse(as.character(s$group) %in% th$pos, 1, 0)
  
  for (form in c("GDI_diff", "GDI_ratio")) {
    
    set.seed(42)
    r <- roc(outcome, s[[form]], quiet = TRUE, direction = "<")
    
    set.seed(42)
    ci <- ci.auc(r, method = "bootstrap", boot.n = N_BOOT, boot.stratified = TRUE)
    
    key <- paste(th$cohort, th$label, form, sep = "|")
    roc_objects[[key]] <- r
    
    roc_tab <- rbind(roc_tab, data.frame(
      cohort = th$cohort,
      threshold = th$label,
      positive_grades = paste(th$pos, collapse = "+"),
      GDI_form = form,
      AUC = round(as.numeric(r$auc), 4),
      CI_lower = round(as.numeric(ci[1]), 4),
      CI_upper = round(as.numeric(ci[3]), 4),
      n_pos = sum(outcome == 1),
      n_neg = sum(outcome == 0),
      stringsAsFactors = FALSE))
  }
}

cat("=== CHECKPOINT 6b: ROC ===\n\n")
print(roc_tab, row.names = FALSE)

write.csv(roc_tab, "results/tables/GDI_ROC.csv", row.names = FALSE)
saveRDS(roc_objects, "data/processed/ROC_objects.rds")


# ---------------------------------------------------------------------------
# 4. How thin are the small-group estimates?
# ---------------------------------------------------------------------------
# An AUC resting on three negatives is not a validation. Report the count of
# discordant pairs so the text can say how fragile each estimate is.

cat("\n\n=== Fragility of each estimate (GDI_diff only) ===\n\n")

frag <- NULL
for (th in thresholds) {
  s <- gdi_list[[th$cohort]]
  outcome <- ifelse(as.character(s$group) %in% th$pos, 1, 0)
  v <- s$GDI_diff
  
  pos <- v[outcome == 1]
  neg <- v[outcome == 0]
  
  discordant <- 0
  for (a in pos) discordant <- discordant + sum(a < neg)
  
  frag <- rbind(frag, data.frame(
    cohort = th$cohort, threshold = th$label,
    n_pos = length(pos), n_neg = length(neg),
    total_pairs = length(pos) * length(neg),
    discordant_pairs = discordant,
    smaller_group = min(length(pos), length(neg)),
    stringsAsFactors = FALSE))
}
print(frag, row.names = FALSE)

cat("\nAny estimate whose smaller group is under 10 cannot carry a validation\n")
cat("claim. State the group size in the same sentence as the AUC.\n")

write.csv(frag, "results/tables/ROC_fragility.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 5. Primary table for the manuscript, GDI_diff only
# ---------------------------------------------------------------------------

cat("\n\n=== Primary ROC table (GDI_diff) ===\n\n")
primary <- roc_tab[roc_tab$GDI_form == "GDI_diff",
                   c("cohort", "threshold", "AUC", "CI_lower", "CI_upper",
                     "n_pos", "n_neg")]
print(primary, row.names = FALSE)

cat("\n\nSecondary (GDI_ratio), reported once in the supplement:\n\n")
secondary <- roc_tab[roc_tab$GDI_form == "GDI_ratio",
                     c("cohort", "threshold", "AUC", "CI_lower", "CI_upper")]
print(secondary, row.names = FALSE)

cat("\nWritten: GDI_trend_tests.csv, GDI_ROC.csv, ROC_fragility.csv,\n")
cat("  data/processed/ROC_objects.rds\n")
