# Execution order

`Rscript run_all.R` runs everything in the right order. This file explains the
order, because four of the dependencies cannot be inferred from the file names
and each of them has silently corrupted results at least once during this
project.

## The order

```
00_setup_and_load.R          environment, directory tree, decision record
01_data_acquisition.R        retrieve four GEO series
01b_scale_fix.R              log2 transform GSE138080, which was deposited linear
02_annotate_and_groups.R     probe annotation and disease group assignment
02b_fix_gpl570_groups.R      correct one GPL570 group from its title field
02c_collapse_maxmean.R       collapse probes to genes by MaxMean
03_QC.R                      quality control, Figure S1
04_geneset_lock.R            lock the 33-gene panel, per-cohort availability
05_DEG_discovery.R           transcriptome-wide differential expression
09_fgsea.R                   module enrichment on the ranked list
06_ssGSEA_GDI.R              module scores and both index forms
07_trend_and_ROC.R           trend statistics and discrimination
08_robustness.R              exclusion sets, leave-one-gene-out, module dropout
10_figures_main.R            Figure 3, module trajectories
11_figures_deg_heatmap.R     Figure S2 volcanoes, heatmap source
11b_figure2_heatmap.R        heatmap panel
13_figure2_rank_strip.R      rank strip panel
14_figure3.R                 Figure 2, location of the signal
15_assemble_figures.R        Figure 1, Figure 4, figure manifest
16_supplementary_figures.R   Figures S3 and S4
17_main_tables.R             Tables 1 to 3
12b_complete_decisions.R     decisions settled during review
17b_update_D11.R             rewrite D11 from Table 2
12_supplementary_tables.R    the supplementary workbook
20_verify_relabel.R          check every reported value against the manuscript
```

## The four constraints that are not obvious

**02c must follow 02.** Script 02 writes `data/processed/gene_matrices.rds`
using a simple average collapse, as an intermediate. Script 02c overwrites that
same file with the MaxMean collapse the analysis is defined on. Run 02 without
02c and every downstream score is computed on the wrong matrices. Scripts 04,
05, 06 and 08 now check the collapse attribute and refuse to run rather than
produce quiet nonsense.

**12b must precede 17b.** Script 12b writes the decisions settled during review,
including an older version of D11 that names two background genes. Script 17b
rewrites D11 from Table 2, where there are three. Run them the other way round
and the correction is undone.

**17b must follow 17.** It builds the D11 text from `Table2_panel.csv`, which 17
writes. It refuses to run if the undetected gene set is not the three expected
skeletal-muscle isoforms.

**12 runs last.** The supplementary workbook collects the output of every other
script, including the decision record. Run it earlier and the sheets are stale.

## Scripts not in the pipeline

These were diagnostic or one-off repairs. They are kept because the decision
record refers to them, but `run_all.R` does not call them and a clean run does
not need them.

| Script | What it did |
| --- | --- |
| `04b_diagnose_coverage.R` | established that GSE138080 covers 28 panel genes, not 26 |
| `05b_diagnose_warburg.R` | traced a Warburg gene mismatch to a sparser probe annotation |
| `12_supplementary_tables.R` is in the pipeline; `18`, `19`, `21`, `22` are not | |
| `18_rename_figures.R` | renamed figure files after the first renumbering |
| `19_add_decision_terminology.R` | recorded the HGCIN terminology correction |
| `21_fix_decision_record.R` | removed duplicate decision ids left by running 12b twice |
| `22_renumber_figures.R` | renumbered figures into order of first citation |
| `23_prepare_for_repo.R` | removed the hardcoded Windows paths |

Scripts 18 and 22 are obsolete by design: the figure scripts now write the
correct names directly, so a clean run never needs a rename pass.

## Runtime

Script 08 recomputes single-sample enrichment 33 times for the
leave-one-gene-out analysis and is the slowest step. Script 12 writes two
transcriptome-wide sheets and takes a minute or two. Everything else is quick.

## If something fails

`run_all.R` stops at the first failure, prints the tail of that script's log and
exits non-zero. Per-script logs are in `logs/`. Fix the cause and restart from
that script:

```
Rscript run_all.R 06
```

Do not skip a failed script and continue. Three of the four constraints above
exist because a later script ran happily on a file an earlier one should have
rewritten.
