# Phase 2b - fix the GPL570 group assignment
#
# One Normal was pulled into Cancer. The rules overwrite in sequence with the
# cancer pattern last, so a normal specimen whose annotation mentions carcinoma
# anywhere gets reclassified. GPL570 titles are explicit (normal_* vs SCC_*),
# so this cohort is assigned from the title alone.
#
# The same edit goes into 02_annotate_and_groups.R for the archival record.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(Biobase)

esets <- readRDS("data/processed/esets_log2only.rds")
gene_mats <- readRDS("data/processed/gene_matrices.rds")
groups <- readRDS("data/processed/groups.rds")

eset <- esets$GSE27678_GPL570
pd <- pData(eset)


# ---------------------------------------------------------------------------
# 1. Diagnose. Show everything GEO carries for the three normal-titled samples.
# ---------------------------------------------------------------------------

titles <- as.character(pd$title)
is_normal_title <- grepl("^normal", titles, ignore.case = TRUE)

cat("=== Full annotation for the normal-titled GPL570 samples ===\n")
fields <- intersect(c("title", "source_name_ch1", "characteristics_ch1",
                      "characteristics_ch1.1", "description"), colnames(pd))

for (i in which(is_normal_title)) {
  cat("\n", sampleNames(eset)[i], "\n", sep = "")
  for (f in fields) {
    cat("   ", f, ": ", as.character(pd[[f]][i]), "\n", sep = "")
  }
  cat("    old assignment: ", as.character(groups$GSE27678_GPL570[i]), "\n", sep = "")
}
cat("\n")


# ---------------------------------------------------------------------------
# 2. Reassign from the title
# ---------------------------------------------------------------------------

g <- rep(NA_character_, length(titles))
g[grepl("^normal", titles, ignore.case = TRUE)] <- "Normal"
g[grepl("^scc", titles, ignore.case = TRUE)] <- "Cancer"

groups$GSE27678_GPL570 <- factor(g, levels = c("Normal", "Cancer"))


# ---------------------------------------------------------------------------
# 3. Re-check all four cohorts
# ---------------------------------------------------------------------------

expected <- list(
  GSE63514        = c(Normal = 24, CIN1 = 14, CIN2 = 22, CIN3 = 40, Cancer = 28),
  GSE138080       = c(Normal = 10, HGCIN = 15, SCC = 10),
  GSE27678_GPL571 = c(Normal = 12, LSIL = 11, HSIL = 21),
  GSE27678_GPL570 = c(Normal = 3,  Cancer = 28)
)

cat("=== CHECKPOINT 2a, rerun ===\n")
all_ok <- TRUE

for (nm in names(groups)) {
  obs <- table(groups[[nm]])
  exp <- expected[[nm]]
  cat("\n", nm, "\n", sep = "")
  for (lv in names(exp)) {
    o <- as.integer(obs[lv])
    if (is.na(o)) o <- 0L
    flag <- if (o == exp[lv]) "ok" else "MISMATCH"
    cat(sprintf("   %-8s expected %3d  observed %3d   %s\n", lv, exp[lv], o, flag))
    if (o != exp[lv]) all_ok <- FALSE
  }
  n_na <- sum(is.na(groups[[nm]]))
  cat(sprintf("   unassigned: %d\n", n_na))
  if (n_na > 0) all_ok <- FALSE
}

if (all_ok) {
  cat("\n>> CHECKPOINT 2a PASSED. All four cohorts match protocol section 4.\n\n")
} else {
  cat("\n>> STILL FAILING. Stop here.\n\n")
}


# ---------------------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------------------

saveRDS(groups, "data/processed/groups.rds")

for (nm in names(groups)) {
  out <- data.frame(sample = colnames(gene_mats[[nm]]),
                    group = as.character(groups[[nm]]),
                    stringsAsFactors = FALSE)
  write.csv(out, paste0("results/tables/sample_groups_", nm, ".csv"), row.names = FALSE)
}

cat("groups.rds and the per-cohort group tables rewritten.\n")
