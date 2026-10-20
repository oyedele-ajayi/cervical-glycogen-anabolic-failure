# Phase 22 - renumber the figures to match order of first citation
#
# Figures were numbered by the order they were built, not the order the
# Results cites them, so the manuscript ran 2, 3, 1, 4. Journals number by
# first citation. The supplementary pair was inverted the same way, because
# fgsea is cited in 3.6 and robustness in 3.7.
#
#   panel genes        Figure 2 -> Figure 1    cited 3.2
#   signal location    Figure 3 -> Figure 2    cited 3.3 and 3.7
#   trajectories       Figure 1 -> Figure 3    cited 3.4
#   ROC                Figure 4 unchanged      cited 3.5
#   QC                 Figure S1 unchanged     cited 3.1
#   volcanoes          Figure S2 unchanged     cited 3.2
#   fgsea              Figure S4 -> Figure S3  cited 3.6
#   robustness         Figure S3 -> Figure S4  cited 3.7
#
# The main set is a three-cycle, so renaming in place would collide. Every
# move goes through a temporary name first.
#
# The producing scripts have been changed to write the new names, so a clean
# run from 00 needs none of this. Run it once to fix the files already on disk.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

fig <- "results/figures"
tab <- "results/tables"

moves <- rbind(
  # figures, both formats
  data.frame(from = "gdi/Figure1_module_trajectories.png", to = "gdi/Figure3_module_trajectories.png"),
  data.frame(from = "gdi/Figure1_module_trajectories.pdf", to = "gdi/Figure3_module_trajectories.pdf"),
  data.frame(from = "deg/Figure2_panel_genes.png",         to = "deg/Figure1_panel_genes.png"),
  data.frame(from = "deg/Figure2_panel_genes.pdf",         to = "deg/Figure1_panel_genes.pdf"),
  data.frame(from = "gdi/Figure3_signal_location.png",     to = "gdi/Figure2_signal_location.png"),
  data.frame(from = "gdi/Figure3_signal_location.pdf",     to = "gdi/Figure2_signal_location.pdf"),
  data.frame(from = "validation/FigureS3_robustness.png",  to = "validation/FigureS4_robustness.png"),
  data.frame(from = "validation/FigureS3_robustness.pdf",  to = "validation/FigureS4_robustness.pdf"),
  data.frame(from = "gdi/FigureS4_fgsea.png",              to = "gdi/FigureS3_fgsea.png"),
  data.frame(from = "gdi/FigureS4_fgsea.pdf",              to = "gdi/FigureS3_fgsea.pdf"),
  stringsAsFactors = FALSE)

table_moves <- data.frame(
  from = c("Figure1_module_trajectories.csv", "Figure3A_module_separation.csv",
           "Figure3B_module_and_index.csv"),
  to   = c("Figure3_module_trajectories.csv", "Figure2A_module_separation.csv",
           "Figure2B_module_and_index.csv"),
  stringsAsFactors = FALSE)


# two-phase move, so a file never lands on a name still in use
two_phase <- function(root, m) {
  tmp <- file.path(root, paste0(".renumber_tmp_", seq_len(nrow(m))))
  done <- 0
  
  for (i in seq_len(nrow(m))) {
    src <- file.path(root, m$from[i])
    dst <- file.path(root, m$to[i])
    if (!file.exists(src)) {
      if (file.exists(dst)) cat("  already done:", m$to[i], "\n")
      else                  cat("  not found:   ", m$from[i], "\n")
      next
    }
    if (!file.rename(src, tmp[i])) {
      cat("  FAILED to stage:", m$from[i], "- close whatever has it open\n")
      next
    }
    done <- done + 1
  }
  
  for (i in seq_len(nrow(m))) {
    if (!file.exists(tmp[i])) next
    dst <- file.path(root, m$to[i])
    if (file.exists(dst)) {
      cat("  target occupied, restoring:", m$from[i], "\n")
      file.rename(tmp[i], file.path(root, m$from[i]))
      next
    }
    ok <- file.rename(tmp[i], dst)
    cat(ifelse(ok, "  renamed:     ", "  FAILED:      "),
        m$from[i], " -> ", m$to[i], "\n", sep = "")
  }
  invisible(done)
}

cat("=== Figures ===\n")
two_phase(fig, moves)

cat("\n=== Result tables named after their figure ===\n")
two_phase(tab, table_moves)

# the trajectory PDF and PNG for the superseded five-panel ROC stay put
cat("\nUntouched: Figure4_ROC, FigureS1_QC_*, FigureS2_volcanoes, superseded/\n")


# ---------------------------------------------------------------------------
# Manifest, rewritten here so this script leaves the project consistent even
# if 15 is not rerun
# ---------------------------------------------------------------------------

manifest <- data.frame(
  figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4",
             "Figure S1", "Figure S2", "Figure S3", "Figure S4"),
  content = c("A: rank position of the 33 panel genes; B: panel gene heatmap, discovery",
              "A: grade at which each module separates from normal; B: discrimination by module and reduced index",
              "Module score trajectories across grade, four cohorts",
              "ROC curves for the index, four cohorts",
              "Quality control: boxplots, density, PCA",
              "Volcano plots, four grade contrasts against normal",
              "fgsea running enrichment, four modules",
              "Robustness: outlier exclusion sets and leave-one-gene-out"),
  file = c("results/figures/deg/Figure1_panel_genes.png",
           "results/figures/gdi/Figure2_signal_location.png",
           "results/figures/gdi/Figure3_module_trajectories.png",
           "results/figures/validation/Figure4_ROC.png",
           "results/figures/qc/FigureS1_QC_<cohort>.png",
           "results/figures/deg/FigureS2_volcanoes.png",
           "results/figures/gdi/FigureS3_fgsea.png",
           "results/figures/validation/FigureS4_robustness.png"),
  script = c("15_assemble_figures.R", "14_figure3.R", "10_figures_main.R",
             "15_assemble_figures.R", "03_QC.R", "11_figures_deg_heatmap.R",
             "16_supplementary_figures.R", "16_supplementary_figures.R"),
  first_cited = c("3.2", "3.3", "3.4", "3.5", "3.1", "3.2", "3.6", "3.7"),
  stringsAsFactors = FALSE)

write.csv(manifest, file.path(fig, "figure_manifest.csv"), row.names = FALSE)


cat("\n=== Check: present, correctly named, and in citation order ===\n")
fail <- 0
for (i in seq_len(nrow(manifest))) {
  f <- manifest$file[i]
  paths <- if (grepl("<cohort>", f, fixed = TRUE))
    Sys.glob(sub("<cohort>", "*", f, fixed = TRUE)) else f
  
  present <- length(paths) > 0 && all(file.exists(paths))
  num <- sub(" ", "", manifest$figure[i])
  named <- present && all(grepl(num, basename(paths), fixed = TRUE))
  if (!present || !named) fail <- fail + 1
  
  cat(sprintf("  %-10s %-8s %s\n", manifest$figure[i],
              ifelse(present, "present", "MISSING"),
              ifelse(!present, "", ifelse(named, "name matches the number",
                                          "NAME DOES NOT MATCH THE NUMBER"))))
}

for (set in list(c("Figure 1", "Figure 2", "Figure 3", "Figure 4"),
                 c("Figure S1", "Figure S2", "Figure S3", "Figure S4"))) {
  k <- match(set, manifest$figure)
  cites <- as.numeric(manifest$first_cited[k])
  lab <- if (grepl("S", set[1])) "supplementary" else "main"
  if (is.unsorted(cites)) {
    fail <- fail + 1
    cat("\n  CITATION ORDER WRONG in the ", lab, " set\n", sep = "")
  } else {
    cat("\n  ", lab, " figures run in order of first citation (",
        paste(manifest$first_cited[k], collapse = ", "), ")\n", sep = "")
  }
}

cat("\n")
if (fail == 0) {
  cat("All eight present, correctly named, and in citation order.\n")
} else {
  cat(fail, "problems. Paste the output before changing the manuscript.\n")
}

cat("\nStill to change by hand in the manuscript: every figure citation and\n")
cat("legend. Figure 3A and 3B become Figure 2A and 2B.\n")
