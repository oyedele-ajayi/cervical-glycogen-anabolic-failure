# Phase 24 - give the scripts names a stranger can read
#
# The scripts are currently named things like "2c. collapse by MaxMean.R" and
# "Phase_22.R". Three problems for a public repository: spaces break
# command-line use, capitalisation is inconsistent and Linux treats
# phase_20.R and Phase_22.R as unrelated names, and four numbers are used
# twice with nothing to say which runs first.
#
# This renames them to lowercase, underscore-separated, two-digit names, and
# names them for what they do rather than which figure they currently make,
# because the figure numbers have already changed twice.
#
# It also retires the stale copy of the phase 10 script, which still writes
# the superseded five-panel ROC into the validation folder, and strips the
# no-op setwd(PROJ) line left in script 1 by the previous patch.
#
# Run from the project root.

rm(list = ls())

if (!dir.exists("scripts")) stop("Run this from the project root.")

map <- rbind(
  c("1.Phase 1.R",                         "01_setup_and_load.R"),
  c("2. Probe_annotation.R",               "02_annotate_and_groups.R"),
  c("2b_fix.R",                            "02b_fix_gpl570_groups.R"),
  c("2c. collapse by MaxMean.R",           "02c_collapse_maxmean.R"),
  c("3. QC.R",                             "03_qc.R"),
  c("4. Gene_Mapping.R",                   "04_geneset_lock.R"),
  c("5. DEG_Discovery.R",                  "05_deg_discovery.R"),
  c("6. ssGSEA_module.R",                  "06_ssgsea_gdi.R"),
  c("7. trend_test_ROC.R",                 "07_trend_and_roc.R"),
  c("8. Robustness_battery.R",             "08_robustness.R"),
  c("9. fgsea.R",                          "09_fgsea.R"),
  c("10. Main_figure_v2.R",                "10_module_trajectories.R"),
  c("11. Deg_heatmap.R",                   "11_volcanoes_and_heatmap.R"),
  c("11. Deg_heatmap_b.R",                 "11b_heatmap_colour.R"),
  c("12. Supplementary_workbook_tables.R", "12_supplementary_workbook.R"),
  c("12b. Supplementary_rerun_b4_12.R",    "12b_complete_decisions.R"),
  c("13. Strip_plot.R",                    "13_rank_strip.R"),
  c("14. figure_3.R",                      "14_signal_location.R"),
  c("15. Final_figure.R",                  "15_assemble_figures.R"),
  c("16. supplementary_figures.R",         "16_supplementary_figures.R"),
  c("17. Main_tables.R",                   "17_main_tables.R"),
  c("17. Correct_D11.R",                   "17b_update_d11.R"),
  c("18.R",                                "18_rename_figures.R"),
  c("19.R",                                "19_add_decision_terminology.R"),
  c("phase_20.R",                          "20_verify_values.R"),
  c("phase_21.R",                          "21_fix_decision_record.R"),
  c("Phase_22.R",                          "22_renumber_figures.R")
)
colnames(map) <- c("from", "to")

retire <- "10. Main_figure.R"   # superseded copy, writes the old Figure 4


# --- look before leaping ---------------------------------------------------

present <- list.files("scripts", pattern = "\\.R$")

cat("=== Before ===\n")
cat("  scripts present:", length(present), "\n")

missing <- setdiff(map[, "from"], present)
if (length(missing) > 0) {
  cat("\n  NOT FOUND, so not renamed:\n")
  for (m in missing) cat("    ", m, "\n")
}

extra <- setdiff(present, c(map[, "from"], retire, "23_prepare_for_repo.R",
                            "24_rename_scripts.R"))
if (length(extra) > 0) {
  cat("\n  present but not in the map, left alone:\n")
  for (e in extra) cat("    ", e, "\n")
}

clash <- intersect(map[, "to"], present)
if (length(clash) > 0) {
  stop("Target name already exists: ", paste(clash, collapse = ", "),
       ". Resolve before rerunning.")
}


# --- rename ----------------------------------------------------------------

cat("\n=== Renaming ===\n")
done <- 0
for (i in seq_len(nrow(map))) {
  from <- file.path("scripts", map[i, "from"])
  to   <- file.path("scripts", map[i, "to"])

  if (!file.exists(from)) next

  if (file.rename(from, to)) {
    cat(sprintf("  %-38s -> %s\n", map[i, "from"], map[i, "to"]))
    done <- done + 1
  } else {
    cat("  FAILED:", map[i, "from"], "- is it open in RStudio?\n")
  }
}
cat("\n  ", done, " renamed\n", sep = "")


# --- retire the stale phase 10 copy ----------------------------------------

if (file.exists(file.path("scripts", retire))) {
  dir.create("scripts_superseded", showWarnings = FALSE)
  ok <- file.rename(file.path("scripts", retire),
                    file.path("scripts_superseded", retire))
  cat("\n", if (ok) "  moved to scripts_superseded/: " else "  FAILED to move: ",
      retire, "\n", sep = "")
  cat("  It writes Figure4_ROC_and_dropout, which the four-panel Figure 4\n")
  cat("  replaced. Kept outside scripts/ rather than deleted.\n")
}


# --- strip the no-op setwd left by the previous patch ----------------------

f <- "scripts/01_setup_and_load.R"
if (file.exists(f)) {
  txt <- readLines(f, warn = FALSE)
  hit <- which(trimws(txt) == "setwd(PROJ)")
  if (length(hit) > 0) {
    txt <- txt[-hit]
    writeLines(txt, f)
    cat("\n  removed the no-op setwd(PROJ) from 01_setup_and_load.R\n")
  }
}


# --- after ------------------------------------------------------------------

cat("\n=== After ===\n")
after <- sort(list.files("scripts", pattern = "\\.R$"))
for (a in after) cat("  ", a, "\n")

bad <- after[grepl("[ A-Z]", after)]
if (length(bad) > 0) {
  cat("\n  Still carrying a space or a capital:\n")
  for (b in bad) cat("    ", b, "\n")
} else {
  cat("\n  All lowercase, no spaces. Safe on every filesystem.\n")
}

cat("\nNext: the corrected run_all.R, which names these files.\n")
