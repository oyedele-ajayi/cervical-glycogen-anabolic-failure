# Phase 18 - bring the figure filenames into line with the manuscript numbers
#
# The figures were built before the numbering was settled, so three of them sit
# on disk under their draft numbers. The trajectory panel is the one that
# matters: it is saved as Figure3 and cited as Figure 1. This renames the files
# that are already there, so the whole figure set does not have to be rebuilt.
#
# The ggsave paths in 03_QC.R, 10_figures_main.R and 11_figures_deg_heatmap.R
# have also been changed, so a clean re-run from script 00 now writes the
# correct names and this script is not needed again.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

fig <- "results/figures"
dir.create(file.path(fig, "superseded"), showWarnings = FALSE, recursive = TRUE)


# one rename, reported either way
move <- function(from, to) {
  from <- file.path(fig, from)
  to   <- file.path(fig, to)
  
  if (!file.exists(from)) {
    if (file.exists(to)) {
      cat("  already done:", basename(to), "\n")
    } else {
      cat("  not found:   ", basename(from), "\n")
    }
    return(invisible(FALSE))
  }
  if (file.exists(to)) {
    cat("  target exists, leaving both:", basename(to), "\n")
    return(invisible(FALSE))
  }
  ok <- file.rename(from, to)
  cat(ifelse(ok, "  renamed:     ", "  FAILED:      "),
      basename(from), " -> ", basename(to), "\n", sep = "")
  invisible(ok)
}


cat("=== Figure 1, module trajectories (was Figure3) ===\n")
move("gdi/Figure3_module_trajectories.png", "gdi/Figure1_module_trajectories.png")
move("gdi/Figure3_module_trajectories.pdf", "gdi/Figure1_module_trajectories.pdf")

cat("\n=== Figure S2, volcano panels (was Figure1_v2) ===\n")
move("deg/Figure1_volcanoes_v2.png", "deg/FigureS2_volcanoes.png")
move("deg/Figure1_volcanoes_v2.pdf", "deg/FigureS2_volcanoes.pdf")

cat("\n=== Figure S1, quality control, one file per cohort ===\n")
qc <- list.files(file.path(fig, "qc"), pattern = "^QC_.*\\.png$")
if (length(qc) == 0) cat("  nothing to rename\n")
for (f in qc) move(file.path("qc", f), file.path("qc", paste0("FigureS1_", f)))

cat("\n=== Superseded, so it cannot be deposited by mistake ===\n")
# the five-panel ROC: script 15 rebuilt it as the four-panel Figure 4, and its
# panel E became Figure 3B in script 14
move("validation/Figure4_ROC_and_dropout.png", "superseded/ROC_5panel_superseded.png")
move("validation/Figure4_ROC_and_dropout.pdf", "superseded/ROC_5panel_superseded.pdf")
# legend text files written by scripts 10 and 11 under the draft numbering; the
# manuscript holds the legends now
move("Figure3_Figure4_legends.txt", "superseded/legends_draft_numbering_10.txt")
move("Figure1_Figure2_legends.txt", "superseded/legends_draft_numbering_11.txt")


# ---------------------------------------------------------------------------
# Rewrite the manifest and check it against what is actually on disk
# ---------------------------------------------------------------------------

manifest <- data.frame(
  figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4",
             "Figure S1", "Figure S2", "Figure S3", "Figure S4"),
  content = c("Module score trajectories across grade, four cohorts",
              "A: rank position of the 33 panel genes; B: panel gene heatmap, discovery",
              "A: grade at which each module separates from normal; B: discrimination by module and reduced index",
              "ROC curves for the index, four cohorts",
              "Quality control: boxplots, density, PCA",
              "Volcano plots, four grade contrasts against normal",
              "Robustness: outlier exclusion sets and leave-one-gene-out",
              "fgsea running enrichment, four modules"),
  file = c("results/figures/gdi/Figure1_module_trajectories.png",
           "results/figures/deg/Figure2_panel_genes.png",
           "results/figures/gdi/Figure3_signal_location.png",
           "results/figures/validation/Figure4_ROC.png",
           "results/figures/qc/FigureS1_QC_<cohort>.png",
           "results/figures/deg/FigureS2_volcanoes.png",
           "results/figures/validation/FigureS3_robustness.png",
           "results/figures/gdi/FigureS4_fgsea.png"),
  script = c("10_figures_main.R", "15_assemble_figures.R", "14_figure3.R",
             "15_assemble_figures.R", "03_QC.R", "11_figures_deg_heatmap.R",
             "16_supplementary_figures.R", "16_supplementary_figures.R"),
  cited_in = c("3.4", "3.2", "3.3 and 3.7", "3.5", "3.1", "3.2", "3.7", "3.6"),
  stringsAsFactors = FALSE)

write.csv(manifest, file.path(fig, "figure_manifest.csv"), row.names = FALSE)

cat("\n=== Check: does every file exist, and does its name carry its number? ===\n")
fail <- 0
for (i in seq_len(nrow(manifest))) {
  f <- manifest$file[i]
  paths <- if (grepl("<cohort>", f, fixed = TRUE)) {
    Sys.glob(sub("<cohort>", "*", f, fixed = TRUE))
  } else f
  
  present <- length(paths) > 0 && all(file.exists(paths))
  num <- sub(" ", "", manifest$figure[i])
  named <- present && all(grepl(num, basename(paths), fixed = TRUE))
  if (!present || !named) fail <- fail + 1
  
  cat(sprintf("  %-10s %-8s %s\n", manifest$figure[i],
              ifelse(present, "present", "MISSING"),
              ifelse(!present, "", ifelse(named, "name matches the number",
                                          "NAME DOES NOT MATCH THE NUMBER"))))
}

cat("\n")
if (fail == 0) {
  cat("All eight figures present and correctly named.\n")
} else {
  cat(fail, "entries need attention before deposit.\n")
}

cat("\nFiles kept under their existing names because they were already correct:\n")
cat("  Figure2_panel_genes, Figure3_signal_location, Figure4_ROC,\n")
cat("  FigureS3_robustness, FigureS4_fgsea\n")
cat("\nIntermediates not in the manifest, consumed by script 15:\n")
cat("  deg/Figure2_heatmap_v2.png, deg/Figure2_rank_strip.png\n")
