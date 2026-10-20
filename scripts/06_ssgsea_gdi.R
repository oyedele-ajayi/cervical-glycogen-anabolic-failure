# Phase 6 - ssGSEA module scores and the two GDI forms
#
# GDI_ratio = (Warburg + Degradation) / (Synthesis + Differentiation)
# GDI_diff  = (Warburg + Degradation) - (Synthesis + Differentiation)
#
# GDI_diff is the primary index (decision D5). The ratio is secondary because
# its denominator crosses zero in some samples, which inverts the sign and
# makes the quantity undefined as a scale.
#
# ssGSEA normalises against the whole sample set, so scores are cohort-internal
# and not comparable across cohorts on an absolute scale. Every subset analysis
# downstream recomputes from the subset matrix rather than slicing these scores.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(GSVA)

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
cohort_genesets <- readRDS("data/processed/cohort_genesets.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
cohorts <- names(gene_mats)

grade_levels <- list(
  GSE63514        = c("Normal", "CIN1", "CIN2", "CIN3", "Cancer"),
  GSE138080       = c("Normal", "HGCIN", "SCC"),
  GSE27678_GPL571 = c("Normal", "LSIL", "HSIL"),
  GSE27678_GPL570 = c("Normal", "Cancer")
)


# ---------------------------------------------------------------------------
# 1. Run ssGSEA per cohort
# ---------------------------------------------------------------------------

gdi_list <- list()

cat("=== ssGSEA ===\n")

for (nm in cohorts) {
  
  m <- gene_mats[[nm]]
  gs <- cohort_genesets[[nm]]
  gs <- gs[module_order]
  
  set.seed(42)
  param <- ssgseaParam(m, gs)
  scores <- gsva(param, verbose = FALSE)
  
  # samples as rows
  s <- as.data.frame(t(scores))
  s <- s[, module_order]
  
  denom <- s$Synthesis + s$Differentiation
  numer <- s$Warburg + s$Degradation
  
  s$GDI_ratio <- numer / denom
  s$GDI_diff <- numer - denom
  s$sample <- rownames(s)
  s$group <- factor(as.character(groups[[nm]]), levels = grade_levels[[nm]])
  s$cohort <- nm
  
  gdi_list[[nm]] <- s
  
  write.csv(s, paste0("results/tables/GDI_scores_", nm, ".csv"), row.names = FALSE)
  cat(sprintf("  %-16s %d samples, modules %s\n", nm, nrow(s),
              paste(sapply(gs, length), collapse = "/")))
}

saveRDS(gdi_list, "data/processed/GDI_results.rds")
cat("\n")


# ---------------------------------------------------------------------------
# 2. CHECKPOINT 5a - per grade summaries
# ---------------------------------------------------------------------------

cat("=== CHECKPOINT 5a: module and GDI medians by grade ===\n")

summary_all <- NULL

for (nm in cohorts) {
  s <- gdi_list[[nm]]
  cat("\n--- ", nm, " ---\n", sep = "")
  
  rows <- NULL
  for (lv in levels(s$group)) {
    sub <- s[s$group == lv, ]
    rows <- rbind(rows, data.frame(
      cohort = nm, grade = lv, n = nrow(sub),
      Synthesis = round(median(sub$Synthesis), 3),
      Degradation = round(median(sub$Degradation), 3),
      Warburg = round(median(sub$Warburg), 3),
      Differentiation = round(median(sub$Differentiation), 3),
      GDI_diff = round(median(sub$GDI_diff), 3),
      GDI_ratio = round(median(sub$GDI_ratio), 3),
      SD_diff = round(sd(sub$GDI_diff), 4),
      SD_ratio = round(sd(sub$GDI_ratio), 4),
      stringsAsFactors = FALSE))
  }
  print(rows, row.names = FALSE)
  summary_all <- rbind(summary_all, rows)
}

write.csv(summary_all, "results/tables/GDI_summary_by_grade.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 3. CHECKPOINT 5b - monotonicity
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 5b: consecutive deltas ===\n")

for (nm in cohorts) {
  sub <- summary_all[summary_all$cohort == nm, ]
  cat("\n--- ", nm, " ---\n", sep = "")
  
  for (form in c("GDI_diff", "GDI_ratio")) {
    v <- sub[[form]]
    d <- diff(v)
    rev <- which(d < 0)
    cat(sprintf("  %-10s medians: %s\n", form, paste(v, collapse = " -> ")))
    cat(sprintf("  %-10s deltas : %s", form, paste(round(d, 3), collapse = "  ")))
    if (length(rev) > 0) {
      cat("   REVERSAL at ",
          paste(paste0(sub$grade[rev], "->", sub$grade[rev + 1]), collapse = ", "), "\n", sep = "")
    } else {
      cat("   monotonic\n")
    }
  }
}


# ---------------------------------------------------------------------------
# 4. CHECKPOINT 5c - sign pathology in the ratio denominator
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 5c: ratio denominator sign ===\n")

for (nm in cohorts) {
  s <- gdi_list[[nm]]
  denom <- s$Synthesis + s$Differentiation
  bad <- which(denom < 0)
  
  cat("\n--- ", nm, " ---\n", sep = "")
  cat("  samples with (Synthesis + Differentiation) < 0 : ", length(bad), "\n", sep = "")
  if (length(bad) > 0) {
    print(data.frame(sample = s$sample[bad], group = as.character(s$group[bad]),
                     denominator = round(denom[bad], 4),
                     GDI_ratio = round(s$GDI_ratio[bad], 3),
                     GDI_diff = round(s$GDI_diff[bad], 3)),
          row.names = FALSE)
  }
  cat(sprintf("  GDI_ratio range: %.3f to %.3f\n", min(s$GDI_ratio), max(s$GDI_ratio)))
  cat(sprintf("  GDI_diff  range: %.3f to %.3f\n", min(s$GDI_diff), max(s$GDI_diff)))
  cat(sprintf("  denominator min: %.4f\n", min(denom)))
}


# ---------------------------------------------------------------------------
# 5. CHECKPOINT 5d - dispersion in the highest grade
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 5d: dispersion, highest grade only ===\n\n")

disp <- NULL
for (nm in cohorts) {
  s <- gdi_list[[nm]]
  top <- levels(s$group)[length(levels(s$group))]
  sub <- s[s$group == top, ]
  sd_d <- sd(sub$GDI_diff)
  sd_r <- sd(sub$GDI_ratio)
  disp <- rbind(disp, data.frame(cohort = nm, grade = top, n = nrow(sub),
                                 SD_GDI_diff = round(sd_d, 4),
                                 SD_GDI_ratio = round(sd_r, 4),
                                 fold = round(sd_r / sd_d, 2),
                                 stringsAsFactors = FALSE))
}
print(disp, row.names = FALSE)
write.csv(disp, "results/tables/GDI_dispersion_highest_grade.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 6. Open question O2 - module-level comparisons against Normal
# ---------------------------------------------------------------------------
# The gene-level CIN1 to CIN2 claim does not survive (zero genes significant).
# Section 3.3 has to be rebuilt on module-level separation instead. This tests
# each module against Normal at every grade, in the discovery cohort.

cat("\n\n=== O2: module scores vs Normal, discovery cohort ===\n")
cat("Mann-Whitney, two-sided, each grade against Normal\n\n")

s <- gdi_list$GSE63514
normal <- s[s$group == "Normal", ]

o2 <- NULL
for (mod in c(module_order, "GDI_diff")) {
  for (lv in c("CIN1", "CIN2", "CIN3", "Cancer")) {
    sub <- s[s$group == lv, ]
    w <- wilcox.test(sub[[mod]], normal[[mod]], exact = FALSE)
    o2 <- rbind(o2, data.frame(
      module = mod, grade = lv,
      median_grade = round(median(sub[[mod]]), 4),
      median_normal = round(median(normal[[mod]]), 4),
      W = as.numeric(w$statistic),
      p = signif(w$p.value, 4),
      stringsAsFactors = FALSE))
  }
}

for (mod in c(module_order, "GDI_diff")) {
  cat("--- ", mod, " ---\n", sep = "")
  print(o2[o2$module == mod, c("grade", "median_grade", "median_normal", "p")],
        row.names = FALSE)
  first_sig <- o2$grade[o2$module == mod & o2$p < 0.05][1]
  cat("  first grade separating from Normal at p < 0.05: ",
      ifelse(is.na(first_sig), "none", first_sig), "\n\n", sep = "")
}

write.csv(o2, "results/tables/O2_module_vs_normal.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 7. O4 - is the Normal to CIN1 reversal noise?
# ---------------------------------------------------------------------------

cat("=== O4: Normal vs CIN1, discovery cohort ===\n")

n_diff <- s$GDI_diff[s$group == "Normal"]
c1_diff <- s$GDI_diff[s$group == "CIN1"]

w <- wilcox.test(c1_diff, n_diff, exact = FALSE)
cat(sprintf("  GDI_diff  Normal median %.4f  CIN1 median %.4f  Mann-Whitney p = %.4f\n",
            median(n_diff), median(c1_diff), w$p.value))

set.seed(42)
boot <- numeric(2000)
for (i in 1:2000) {
  boot[i] <- median(sample(c1_diff, replace = TRUE)) - median(sample(n_diff, replace = TRUE))
}
ci <- quantile(boot, c(0.025, 0.975))
cat(sprintf("  bootstrap 95%% CI on the median difference (CIN1 - Normal): [%.3f, %.3f]\n",
            ci[1], ci[2]))
cat(sprintf("  CIN1 exceeds Normal in %.1f%% of 2000 resamples\n", 100 * mean(boot > 0)))

cat("\nWritten: GDI_scores_*.csv, GDI_summary_by_grade.csv,\n")
cat("  GDI_dispersion_highest_grade.csv, O2_module_vs_normal.csv\n")
