# Phase 23 - make the scripts runnable by someone who is not you
#
# Every script begins setwd("C:/Users/USER/Documents/GDI_project"). That path
# exists on one computer in the world, so as the repository stands nobody can
# run a line of it. This replaces the setwd call with a check that the working
# directory is the project root, which is what opening the .Rproj file in
# RStudio already gives you.
#
# Nothing about the analysis changes. The scripts still read and write the same
# relative paths; they just stop insisting on an absolute one.
#
# Run this from the project root. It edits the scripts in place and keeps a
# backup of each under scripts_before_repo_patch/.

rm(list = ls())

scripts_dir <- "scripts"
if (!dir.exists(scripts_dir)) {
  scripts_dir <- "."
  cat("No scripts/ folder; patching the .R files in the current directory.\n\n")
}

backup_dir <- "scripts_before_repo_patch"
dir.create(backup_dir, showWarnings = FALSE)

files <- list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE)
files <- files[!grepl("23_prepare_for_repo\\.R$|run_all\\.R$", files)]

cat("Found", length(files), "scripts.\n\n")


old_setwd <- 'setwd("C:/Users/USER/Documents/GDI_project")'

new_block <- paste(
  '# Run from the repository root. Opening the .Rproj file in RStudio puts you',
  '# there. The check stops early rather than writing results into whatever',
  '# folder happens to be current.',
  'if (!dir.exists("scripts") || !dir.exists("data")) {',
  '  stop("Run this from the repository root, the folder holding scripts/ and data/.")',
  '}',
  sep = "\n")

patched <- 0
skipped <- 0
special <- character(0)

for (f in files) {
  txt <- readLines(f, warn = FALSE)
  hit <- grep(old_setwd, txt, fixed = TRUE)

  # script 00 builds the directory tree from a PROJ variable instead
  proj <- grep('^PROJ <- "C:/Users/USER/Documents/GDI_project"', txt)

  if (length(hit) == 0 && length(proj) == 0) {
    cat("  no change needed:", basename(f), "\n")
    skipped <- skipped + 1
    next
  }

  file.copy(f, file.path(backup_dir, basename(f)), overwrite = TRUE)

  if (length(hit) > 0) {
    txt[hit] <- new_block
  }
  if (length(proj) > 0) {
    txt[proj] <- paste(
      '# The project root is wherever this is run from, which is the repository',
      '# root when the .Rproj file is open.',
      'PROJ <- normalizePath(".", winslash = "/")',
      sep = "\n")
    special <- c(special, basename(f))
  }

  writeLines(txt, f)
  cat("  patched:", basename(f), "\n")
  patched <- patched + 1
}

cat("\n", patched, " patched, ", skipped, " unchanged.\n", sep = "")
cat("Backups in ", backup_dir, "/\n", sep = "")

if (length(special) > 0) {
  cat("\nPROJ rewritten in: ", paste(special, collapse = ", "), "\n", sep = "")
  cat("Read that script once and confirm the directory tree it builds still\n")
  cat("lands where you expect.\n")
}


# --- verify nothing absolute is left ---------------------------------------

cat("\n=== Check: any absolute paths remaining ===\n")
left <- 0
for (f in list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE)) {
  txt <- readLines(f, warn = FALSE)
  bad <- grep("C:/|/Users/|setwd\\(", txt)
  if (length(bad) > 0) {
    left <- left + length(bad)
    for (i in bad) cat("  ", basename(f), " line ", i, ": ", trimws(txt[i]), "\n", sep = "")
  }
}
if (left == 0) {
  cat("  none. The scripts will run from any clone of the repository.\n")
} else {
  cat("\n", left, " lines still carry an absolute path or a setwd call.\n", sep = "")
  cat("Fix those by hand before committing.\n")
}

cat("\nNext: run_all.R, which runs the scripts in the right order.\n")
