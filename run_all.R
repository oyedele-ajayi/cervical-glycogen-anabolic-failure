# Reproduce every reported value, in order
#
# Each script runs in its own R session. That matters: the scripts begin with
# rm(list = ls()), so sourcing them into one session would clear this one's
# state halfway through. Running them separately also means a failure stops the
# pipeline instead of leaving later scripts to compute on stale files, which is
# how the two worst errors in this project's history happened.
#
# Usage, from the repository root:
#   Rscript run_all.R            run everything
#   Rscript run_all.R 06 07      run only the scripts whose numbers you name
#
# Expect a long run. Script 08 recomputes enrichment 33 times and script 12
# writes a workbook with two transcriptome-wide sheets.

args <- commandArgs(trailingOnly = TRUE)

if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}


# Not alphabetical order, and four dependencies are invisible from the
# numbering alone. See RUN_ORDER.md.
order <- c(
  "01_setup_and_load.R",
  "02_annotate_and_groups.R",
  "02b_fix_gpl570_groups.R",
  "02c_collapse_maxmean.R",        # must follow 02: it overwrites 02's matrices
  "03_qc.R",
  "04_geneset_lock.R",
  "05_deg_discovery.R",
  "09_fgsea.R",
  "06_ssgsea_gdi.R",
  "07_trend_and_roc.R",
  "08_robustness.R",
  "10_module_trajectories.R",
  "11_volcanoes_and_heatmap.R",
  "11b_heatmap_colour.R",
  "13_rank_strip.R",
  "14_signal_location.R",
  "15_assemble_figures.R",
  "16_supplementary_figures.R",
  "17_main_tables.R",
  "12b_complete_decisions.R",      # must precede 17b, whose D11 it would undo
  "17b_update_d11.R",              # must follow 17: reads Table2_panel.csv
  "19_add_decision_terminology.R",
  "12_supplementary_workbook.R",   # last: collects everything above
  "20_verify_values.R"             # checks the output against the manuscript
)

if (length(args) > 0) {
  keep <- vapply(order, function(s) any(startsWith(s, args)), logical(1))
  order <- order[keep]
  if (length(order) == 0) stop("No script matched: ", paste(args, collapse = ", "))
}


rscript <- file.path(R.home("bin"), "Rscript")
started <- Sys.time()
failed <- character(0)
dir.create("logs", showWarnings = FALSE)

cat("Running", length(order), "scripts\n")
cat(strrep("-", 64), "\n", sep = "")

for (s in order) {
  path <- file.path("scripts", s)
  if (!file.exists(path)) {
    cat(sprintf("  %-34s MISSING\n", s))
    failed <- c(failed, s)
    next
  }

  t0 <- Sys.time()
  log <- file.path("logs", sub("\\.R$", ".log", s))

  status <- system2(rscript, c("--vanilla", shQuote(path)),
                    stdout = log, stderr = log)

  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %-34s %s %6.1fs\n", s,
              if (status == 0) "ok     " else "FAILED ", secs))

  if (status != 0) {
    failed <- c(failed, s)
    cat("\n  Stopped. The last 20 lines of ", log, ":\n\n", sep = "")
    cat(paste0("    ", tail(readLines(log, warn = FALSE), 20), collapse = "\n"), "\n\n")
    break
  }
}

cat(strrep("-", 64), "\n", sep = "")
cat(sprintf("Elapsed %.1f minutes\n",
            as.numeric(difftime(Sys.time(), started, units = "mins"))))

if (length(failed) == 0) {
  cat("\nAll scripts completed. 20_verify_values.R compares the output against\n")
  cat("the values reported in the manuscript; read its log for 30 passes.\n")
} else {
  cat("\nFailed: ", paste(failed, collapse = ", "), "\n", sep = "")
  cat("Per-script logs are in logs/. Fix the first failure and restart there:\n")
  cat("  Rscript run_all.R ", substr(failed[1], 1, 3), "\n", sep = "")
  quit(status = 1)
}
