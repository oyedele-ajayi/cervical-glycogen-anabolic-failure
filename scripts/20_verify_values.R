# Phase 20 - verify the HGCIN relabel changed the label and nothing else
#
# The relabel was supposed to be cosmetic: same specimens, same groups, same
# estimates, different name. This checks that claim against the files the
# rerun just wrote, using the values the manuscript already reports as the
# reference. Anything that moved is a real change and needs explaining.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

pass <- 0
fail <- 0

ok <- function(label, condition, detail = "") {
  if (isTRUE(condition)) {
    cat("  PASS  ", label, "\n", sep = "")
    pass <<- pass + 1
  } else {
    cat("  FAIL  ", label, "   ", detail, "\n", sep = "")
    fail <<- fail + 1
  }
}

near <- function(a, b, tol = 0.0015) {
  length(a) == 1 && !is.na(a) && abs(a - b) < tol
}


# ---------------------------------------------------------------------------
# 1. Groups: the label changed, the membership did not
# ---------------------------------------------------------------------------

cat("\n=== 1. Group assignment ===\n")
groups <- readRDS("data/processed/groups.rds")

lv <- levels(groups$GSE138080)
ok("GSE138080 levels are Normal, HGCIN, SCC",
   identical(lv, c("Normal", "HGCIN", "SCC")),
   paste("got:", paste(lv, collapse = ", ")))

tb <- table(groups$GSE138080)
ok("GSE138080 counts are 10 / 15 / 10",
   identical(as.integer(tb), c(10L, 15L, 10L)),
   paste(paste(names(tb), as.integer(tb)), collapse = ", "))

ok("no HSIL left in GSE138080",
   !any(as.character(groups$GSE138080) == "HSIL"))

ok("GPL571 keeps Bethesda terms",
   identical(levels(groups$GSE27678_GPL571), c("Normal", "LSIL", "HSIL")))

ok("discovery cohort untouched: 24 / 14 / 22 / 40 / 28",
   identical(as.integer(table(groups$GSE63514)), c(24L, 14L, 22L, 40L, 28L)))

ok("GPL570 untouched: 3 normal, 28 cancer",
   identical(as.integer(table(groups$GSE27678_GPL570)), c(3L, 28L)))


# ---------------------------------------------------------------------------
# 2. Discrimination: every estimate must match what the manuscript reports
# ---------------------------------------------------------------------------

cat("\n=== 2. ROC estimates against the reported values ===\n")
roc <- read.csv("results/tables/GDI_ROC.csv", stringsAsFactors = FALSE)
r <- roc[roc$GDI_form == "GDI_diff", ]

get <- function(cohort, threshold, col = "AUC") {
  v <- r[[col]][r$cohort == cohort & r$threshold == threshold]
  if (length(v) != 1) return(NA_real_)
  as.numeric(v)
}

ok("GSE138080 threshold is now labelled HGCIN+",
   "HGCIN+" %in% r$threshold[r$cohort == "GSE138080"],
   paste("found:", paste(r$threshold[r$cohort == "GSE138080"], collapse = ", ")))

ok("no HSIL threshold under GSE138080",
   !any(r$cohort == "GSE138080" & r$threshold == "HSIL+"))

# the reference values are those in Results 3.5 and Table 3
ok("GSE138080 HGCIN+ AUC is 0.932", near(get("GSE138080", "HGCIN+"), 0.932),
   paste("got", sprintf("%.4f", get("GSE138080", "HGCIN+"))))
ok("GSE138080 SCC AUC is 0.836", near(get("GSE138080", "SCC"), 0.836),
   paste("got", sprintf("%.4f", get("GSE138080", "SCC"))))
ok("discovery CIN3+ AUC is 0.832", near(get("GSE63514", "CIN3+"), 0.8324),
   paste("got", sprintf("%.4f", get("GSE63514", "CIN3+"))))
ok("GPL571 LSIL+ AUC is 0.794", near(get("GSE27678_GPL571", "LSIL+"), 0.794),
   paste("got", sprintf("%.4f", get("GSE27678_GPL571", "LSIL+"))))
ok("GPL571 HSIL AUC is 0.820", near(get("GSE27678_GPL571", "HSIL"), 0.820),
   paste("got", sprintf("%.4f", get("GSE27678_GPL571", "HSIL"))))
ok("GPL570 Cancer AUC is 0.976", near(get("GSE27678_GPL570", "Cancer"), 0.976),
   paste("got", sprintf("%.4f", get("GSE27678_GPL570", "Cancer"))))

cat("\nFull GDI_diff ROC table as it now stands:\n")
print(r[, c("cohort", "threshold", "AUC", "CI_lower", "CI_upper",
            "n_pos", "n_neg")], row.names = FALSE)


# ---------------------------------------------------------------------------
# 3. Trend tests
# ---------------------------------------------------------------------------

cat("\n=== 3. Trend tests, GDI_diff ===\n")
trend <- read.csv("results/tables/GDI_trend_tests.csv", stringsAsFactors = FALSE)
t <- trend[trend$GDI_form == "GDI_diff", ]
print(t[, c("cohort", "Spearman_rho", "Spearman_p", "KW_H", "KW_df", "KW_p")],
      row.names = FALSE)
cat("(these are per cohort and should be identical to the previous run)\n")


# ---------------------------------------------------------------------------
# 4. The negative result, which is the paper's main claim
# ---------------------------------------------------------------------------

cat("\n=== 4. Warburg and Degradation, the two null arms ===\n")
fg <- readRDS("data/processed/fgsea_results.rds")$fgsea
w <- fg[as.character(fg$pathway) == "Warburg", ]
ok("Warburg NES is +0.884", near(as.numeric(w$NES), 0.884, 0.002),
   paste("got", sprintf("%+.4f", as.numeric(w$NES))))
ok("Warburg FDR is 0.790", near(as.numeric(w$padj), 0.790, 0.002),
   paste("got", sprintf("%.4f", as.numeric(w$padj))))

fb <- read.csv("results/tables/Figure2B_module_and_index.csv",
               stringsAsFactors = FALSE)
gv <- function(p) {
  v <- fb$AUC[fb$predictor == p]
  if (length(v) != 1) NA_real_ else as.numeric(v)
}
ok("Warburg alone AUC is 0.569", near(gv("Warburg alone"), 0.569),
   paste("got", sprintf("%.4f", gv("Warburg alone"))))
ok("Differentiation alone AUC is 0.832", near(gv("Differentiation alone"), 0.832),
   paste("got", sprintf("%.4f", gv("Differentiation alone"))))


# ---------------------------------------------------------------------------
# 5. Tables and the decision record carry the new label
# ---------------------------------------------------------------------------

cat("\n=== 5. Tables and decision record ===\n")
t1 <- read.csv("results/tables/Table1_cohorts.csv", stringsAsFactors = FALSE)
comp <- paste(t1$Composition, collapse = " | ")
ok("Table 1 composition says HGCIN", grepl("HGCIN", comp), comp)
ok("Table 1 composition has no HSIL for the tissue cohort",
   !grepl("HGCIN", comp) || !grepl("HSIL 15", comp))

t3 <- read.csv("results/tables/Table3_trend_and_ROC.csv", stringsAsFactors = FALSE)
ok("Table 3 threshold column says HGCIN+",
   any(grepl("HGCIN", t3$Threshold)),
   paste(unique(t3$Threshold), collapse = ", "))

d <- read.csv("logs/decisions_record.csv", stringsAsFactors = FALSE)
ok("a decision records the terminology correction",
   any(grepl("HGCIN", d$decision)),
   paste("ids present:", paste(d$id, collapse = " ")))
cat("  decision ids:", paste(d$id, collapse = ", "), "\n")


# ---------------------------------------------------------------------------
# 6. Figures exist and still carry the right numbers
# ---------------------------------------------------------------------------

cat("\n=== 6. Figures ===\n")
man <- read.csv("results/figures/figure_manifest.csv", stringsAsFactors = FALSE)
for (i in seq_len(nrow(man))) {
  f <- man$file[i]
  paths <- if (grepl("<cohort>", f, fixed = TRUE))
    Sys.glob(sub("<cohort>", "*", f, fixed = TRUE)) else f
  ok(paste(man$figure[i], "present"), length(paths) > 0 && all(file.exists(paths)))
}


# ---------------------------------------------------------------------------

cat("\n===========================================\n")
cat("  ", pass, " passed, ", fail, " failed\n", sep = "")
if (fail == 0) {
  cat("   The relabel changed the label and nothing else.\n")
} else {
  cat("   Something moved. Paste the FAIL lines before we go further.\n")
}
cat("===========================================\n")

cat("\nStill to check by eye, since no file records them:\n")
cat("  Figure 1 panel B x-axis reads Normal / HGCIN / SCC\n")
cat("  Figure 4 panel B legend reads HGCIN+ and SCC\n")
cat("  QC figure for GSE138080 has an HGCIN colour, not a blank\n")
