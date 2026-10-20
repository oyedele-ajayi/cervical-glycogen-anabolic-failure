# Results

Every table and figure the manuscript reports is here, except four large files
that the pipeline regenerates.

Not tracked in git:

| File | Size | Why |
| --- | --- | --- |
| `tables/Supplementary_Tables.xlsx` | 32 MB | the supplementary workbook; sheets S2, S3 and S3b carry the three tables below |
| `tables/DEG_discovery_full.csv` | 17 MB | transcriptome-wide differential expression, full cohort |
| `tables/DEG_discovery_criteria4.csv` | 17 MB | the same, pre-specified exclusion set |
| `tables/DEG_discovery_legacy5.csv` | 17 MB | the same, earlier version's exclusion set |

Git stores every version of every file permanently, so tracking these would
burden every clone forever with output that `Rscript run_all.R` reproduces
exactly. They are published as the supplementary file accompanying the paper, where the three tables appear as sheets S2, S3 and S3b of the workbook.

`figures/figure_manifest.csv` maps each figure to its file, the script that
makes it, and the Results subsection where it is first cited. The numbering
follows order of first citation and the manifest is checked on every run.
