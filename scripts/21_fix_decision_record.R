# Phase 21 - repair the decision record
#
# Script 20 printed the ids as
#   D1 D2 D3 D4 D5 D6 D7 D10 D8 D9 D11 D12 D13 D14 D15 D16 D8 D9 D11 D12 D13 D14 D15
# so D8, D9 and D11 to D15 each appear twice. Something appended the later
# decisions a second time instead of updating them in place. This keeps one
# row per id, puts them in numerical order, and reports anything it had to
# choose between rather than deciding quietly.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

f <- "logs/decisions_record.csv"
d <- read.csv(f, stringsAsFactors = FALSE)

cat("Rows in:", nrow(d), "\n")
cat("Ids in: ", paste(d$id, collapse = ", "), "\n\n")

backup <- sub("\\.csv$", "_before_dedup.csv", f)
write.csv(d, backup, row.names = FALSE)
cat("Backup written to", backup, "\n\n")


# --- do the duplicates actually agree? ------------------------------------

dups <- unique(d$id[duplicated(d$id)])

if (length(dups) == 0) {
  cat("No duplicate ids. Nothing to repair.\n")
} else {
  cat("=== Duplicated ids:", paste(dups, collapse = ", "), "===\n\n")
  
  keep_last <- character(0)
  for (id in dups) {
    rows <- d[d$id == id, ]
    same <- all(duplicated(rows)[-1])
    cat(id, if (same) ": copies are identical, keeping one\n"
        else ": copies DIFFER, see below\n", sep = "")
    if (!same) {
      for (i in seq_len(nrow(rows))) {
        cat("   copy ", i, ": ", substr(rows$decision[i], 1, 160), "...\n", sep = "")
      }
      cat("   keeping the LAST copy, on the assumption it is the corrected one.\n")
      keep_last <- c(keep_last, id)
    }
  }
  cat("\n")
  
  # keep the last occurrence of each id, which is the corrected one where
  # they differ, and identical either way where they do not
  d <- d[!duplicated(d$id, fromLast = TRUE), ]
  
  # numerical order, so D10 stops sitting between D7 and D8
  num <- suppressWarnings(as.integer(sub("^D", "", d$id)))
  d <- d[order(is.na(num), num), ]
  
  write.csv(d, f, row.names = FALSE)
  
  cat("Rows out:", nrow(d), "\n")
  cat("Ids out: ", paste(d$id, collapse = ", "), "\n\n")
  
  if (length(keep_last) > 0) {
    cat("CHECK BY HAND: the copies of ",
        paste(keep_last, collapse = ", "),
        " were not identical.\n", sep = "")
    cat("Open the backup and confirm the kept version is the one you want.\n\n")
  }
}


# --- the terminology decision should be there exactly once -----------------

hg <- grep("HGCIN", d$decision)
cat("Decisions mentioning HGCIN:", length(hg),
    if (length(hg) == 1) " (correct)\n" else " (expected 1)\n", sep = "")
if (length(hg) >= 1) cat("  id:", paste(d$id[hg], collapse = ", "), "\n")

cat("\nRerun 12_supplementary_tables.R so the decisions sheet matches.\n")
