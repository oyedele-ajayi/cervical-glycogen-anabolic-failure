# Phase 19 - record the GSE138080 terminology correction as a decision
#
# The grades in GSE138080 were called HSIL. That is a Bethesda cytological
# term and the series is histology on bulk tissue. The grouping was right, the
# name was wrong. Since every other analytical choice sits in the decision
# record, this one belongs there too.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

d <- read.csv("logs/decisions_record.csv", stringsAsFactors = FALSE)

# next free D number, rather than a hard-coded one
used <- as.integer(sub("^D", "", d$id[grepl("^D[0-9]+$", d$id)]))
new_id <- paste0("D", max(used) + 1)

cat("Existing decisions:", length(used), "- adding", new_id, "\n\n")

decision <- paste(
  "Grades are named in the system that applies to the specimen. GSE138080 is",
  "histology on bulk tissue, so its grades are Normal, HGCIN (high-grade",
  "cervical intraepithelial neoplasia, CIN2/3) and SCC, following the labels",
  "the series itself uses. An earlier version of the pipeline called the",
  "middle group HSIL, which is a Bethesda cytological term and does not apply",
  "to histological specimens. The group assignment was unaffected: 10 Normal,",
  "15 high-grade and 10 carcinoma, matching the series exactly, so no",
  "estimate changes. HSIL and LSIL are retained for GSE27678 GPL571, which is",
  "exfoliated cytology and where those terms are correct. No cohort is",
  "relabelled to match another, because HGCIN and HSIL are different",
  "observations made on different material.")

row <- d[1, ]
row[] <- NA
row$id <- new_id
row$decision <- decision
if ("stage" %in% names(d))     row$stage <- "Phase 2, group assignment"
if ("rationale" %in% names(d)) row$rationale <- "Histological specimens require histological nomenclature."
if ("timing" %in% names(d))    row$timing <- "post hoc, after manuscript review"

d <- rbind(d, row)
write.csv(d, "logs/decisions_record.csv", row.names = FALSE)

cat(new_id, "added.\n\n")
print(d[d$id == new_id, ], row.names = FALSE)

cat("\nColumns in the record:", paste(names(d), collapse = ", "), "\n")
cat("Fill any column this script left blank by hand if your record has one\n")
cat("it does not know about.\n")
cat("\nRerun 12_supplementary_tables.R so the decisions sheet picks this up.\n")
