# Phase 17b - correct D11 in the decision record
#
# D11 named two genes at background. Measuring detectability in the normal
# group rather than across all specimens identifies three, and all three are
# the skeletal-muscle isoforms.
#
# Revised: the expression values are now read from Table2_panel.csv rather
# than typed in. They were typed in before, which meant the decision text
# could disagree with the table it describes without anything noticing. Run
# this after 17, which writes that table.
#
# Assigning D11 rather than appending it, so rerunning is harmless.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

t2 <- read.csv("results/tables/Table2_panel.csv", stringsAsFactors = FALSE,
               check.names = FALSE)
groups <- readRDS("data/processed/groups.rds")

norm_col <- "MeanExpr (Normal)"
canc_col <- "MeanExpr (Cancer)"
stopifnot(all(c(norm_col, canc_col, "Detected", "Gene") %in% names(t2)))

n_discovery <- length(groups$GSE63514)


# --- which genes are actually undetected, per the table --------------------

und <- t2$Gene[t2$Detected == "no"]
cat("Undetected in normal epithelium, per Table 2:",
    paste(und, collapse = ", "), "\n")

expected <- c("PYGM", "PHKG1", "PPP1R3A")
if (!setequal(und, expected)) {
  stop("The undetected set is ", paste(und, collapse = ", "),
       ", not ", paste(expected, collapse = ", "),
       ". D11 describes the three muscle isoforms, so do not write it until\n",
       "  this difference is understood. Check that 02c ran before 17.")
}

# order them as D11 reads, lowest first
und <- expected[order(as.numeric(t2[[norm_col]][match(expected, t2$Gene)]))]

val <- function(gene, col) as.numeric(t2[[col]][match(gene, t2$Gene)])

phrase <- paste0(und, " (mean log2 ", sprintf("%.2f", val(und, norm_col)), ")")
undetected_text <- paste0(paste(phrase[-length(phrase)], collapse = ", "),
                          " and ", phrase[length(phrase)])

# the two genes that a whole-cohort average would wrongly flag
floored <- c("GYS2", "TGM3")
floored_text <- paste0(
  sprintf("%s (%.2f in normal, %.2f in cancer)", floored[1],
          val(floored[1], norm_col), val(floored[1], canc_col)),
  " and ",
  sprintf("%s (%.2f and %.2f)", floored[2],
          val(floored[2], norm_col), val(floored[2], canc_col)))

cat("Built from the table:\n  undetected:", undetected_text,
    "\n  floored:   ", floored_text, "\n\n")


# --- write D11 -------------------------------------------------------------

d <- read.csv("logs/decisions_record.csv", stringsAsFactors = FALSE)
if (!any(d$id == "D11")) stop("No D11 row in the decision record.")

cat("Replacing D11.\nOld:\n", d$decision[d$id == "D11"][1], "\n\n", sep = "")

d$decision[d$id == "D11"] <- paste0(
  "The 33-gene panel is not trimmed. Three members are not detected above ",
  "array background in normal cervical epithelium: ", undetected_text, ". ",
  "All three are the skeletal-muscle isoforms of their enzymes, while the ",
  "liver and brain isoforms are well expressed, so their absence is expected ",
  "behaviour of squamous epithelium and serves as an internal check that the ",
  "panel reports tissue-appropriate biology. Detectability is assessed in the ",
  "normal group rather than across all specimens, because averaging over all ",
  n_discovery, " would place ", floored_text, " below the same threshold ",
  "purely because they are floored in cancer, which is the finding rather ",
  "than a detection failure. Trimming after inspecting detection levels ",
  "would be post hoc panel surgery, so detection status is reported per gene ",
  "in Table 2 and the leave-one-gene-out analysis quantifies the ",
  "contribution of each.")

write.csv(d, "logs/decisions_record.csv", row.names = FALSE)
cat("New:\n", d$decision[d$id == "D11"][1], "\n\n", sep = "")

if (sum(d$id == "D11") > 1) {
  cat("NOTE: the record holds", sum(d$id == "D11"),
      "D11 rows. All were set to the same text, so 21_fix_decision_record.R\n")
  cat("will report them as identical and keep one.\n\n")
}


# --- the leave-one-gene-out values the Results cites ------------------------

loo <- read.csv("results/tables/leave_one_gene_out.csv", stringsAsFactors = FALSE)
sel <- loo$gene %in% expected

cat("=== Leave-one-gene-out, the three undetected genes ===\n")
print(loo[sel, c("gene", "module", "AUC", "delta")], row.names = FALSE)

cat("\nLargest absolute change among the three:",
    sprintf("%.4f", max(abs(loo$delta[sel]))), "\n")

cat("\nRerun 12_supplementary_tables.R so sheet S12 picks up the revised D11.\n")
