# Anabolic failure, not catabolic activation, underlies glycogen depletion in cervical carcinogenesis

Analysis code and results for a multi-cohort transcriptomic study of why
glycogen disappears from cervical epithelium as disease grade rises.

Glycogen loss from dysplastic cervical epithelium has been demonstrable since
Schiller introduced iodine staining, and it is the tissue change that visual
inspection with Lugol's iodine reads directly. What has never been established
is which part of glycogen metabolism accounts for the loss. Histochemistry
reports that glycogen is absent; it cannot report why.

This study measures all four candidate mechanisms on a single scale in the same
specimens: glycogen synthesis, glycogen degradation, Warburg-axis glycolysis,
and squamous differentiation. The answer is that depletion is anabolic failure
inside a collapsing differentiation programme. The catabolic and glycolytic
arms, which the cancer-metabolism literature would have predicted, do not move
at any disease grade in any cohort.

## Main findings

- Synthesis and Differentiation module scores fall with grade in all four
  cohorts. Degradation and Warburg stay within 0.076 of their normal baseline
  everywhere.
- The synthetic failure is selective. Of nine Synthesis members only two change:
  GYS2, the hepatic-type synthase, and PPP1R3C, the PTG targeting subunit. The
  ubiquitous synthase GYS1 does not move. Glycogen storage in this tissue
  depends on a restricted axis, and it is that axis which fails.
- All seven Differentiation members fall and all seven sit in the top 5% of the
  ranked transcriptome, with CRNN the most down-regulated gene of 21,755. That
  is a programme being switched off rather than genes drifting.
- Three independent analyses agree on the negative. Warburg separates from
  normal at no grade (p 0.45 to 0.62), classifies CIN3 or worse at AUC 0.569
  with an interval spanning chance, and shows no enrichment in carcinoma
  (NES +0.884, FDR 0.790).
- The composite index is an instrument for decomposing the mechanism, not a
  classifier. Every form containing the Differentiation module falls within
  0.007 of every other, so the discriminating information sits in one module and
  no combination is advanced as a test.

## Cohorts

| Role | Accession | Platform | Preparation | n |
| --- | --- | --- | --- | --- |
| Discovery | GSE63514 | Affymetrix HG-U133 Plus 2.0 (GPL570) | Laser-capture-microdissected epithelium | 128 |
| Replication 1 | GSE138080 | Agilent 4x44K (GPL4133) | Bulk tissue | 35 |
| Replication 2 | GSE27678 | Affymetrix HG-U133A 2.0 (GPL571) | Exfoliated cytology | 44 |
| Concordance check | GSE27678 | Affymetrix HG-U133 Plus 2.0 (GPL570) | Exfoliated cytology | 31 |

238 specimens in total, after excluding two immortalised cell lines from the
GPL570 set.

Platform and specimen preparation vary independently. The same Affymetrix array
covers both microdissected epithelium and exfoliated cytology, and a single
cytology series spans two array chemistries, so a platform effect cannot be
mistaken for a preparation effect.

Grades are named in the system that belongs to the specimen: histological terms
for microdissected epithelium and bulk tissue, Bethesda terms for cytology.
HGCIN and HSIL are not interchangeable labels for one observation.

## Reproducing the analysis

```
git clone https://github.com/oyedele-ajayi/cervical-glycogen-anabolic-failure.git
cd cervical-glycogen-anabolic-failure
Rscript run_all.R
```

Open the `.Rproj` file first if you work in RStudio; the scripts expect the
repository root as the working directory and stop with a clear message if they
do not get it.

The raw GEO data is not in this repository. Script 01 downloads it, which takes
some minutes and about 2 GB. See `data/README.md`.

Execution order matters more than the numbering suggests. `RUN_ORDER.md` sets
out the four dependencies that are invisible from the file names, each of which
has corrupted results at least once during this project.

Script 20 is the last step and compares what the pipeline just produced against
every value reported in the manuscript. A clean run returns 30 passes.

## What is in here

```
scripts/        the analysis, numbered in dependency order
run_all.R       runs them in a fresh session each, stops at the first failure
RUN_ORDER.md    the order and the four non-obvious constraints
results/
  tables/       every CSV the manuscript reports, plus the supplementary workbook
  figures/      the four main and four supplementary figures, PNG and PDF
  figures/figure_manifest.csv  which file is which figure, and where it is cited
logs/
  decisions_record.csv   every analytical decision taken, with its rationale
  session_info.txt       R and package versions the results were produced under
data/           empty; see data/README.md for what script 01 puts here
```

## Figures

Figures are numbered by order of first citation. The manifest in
`results/figures/figure_manifest.csv` records the mapping and the scripts check
it on every run.

| Figure | Content | First cited |
| --- | --- | --- |
| 1 | Panel gene rank position and expression heatmap | Results 3.2 |
| 2 | Grade at which each module separates; discrimination by module | 3.3 and 3.7 |
| 3 | Module score trajectories across grade, four cohorts | 3.4 |
| 4 | ROC curves for the index, four cohorts | 3.5 |
| S1 | Quality control per cohort | 3.1 |
| S2 | Volcano plots, four grade contrasts | 3.2 |
| S3 | fgsea running enrichment, four modules | 3.6 |
| S4 | Robustness: exclusion sets and leave-one-gene-out | 3.7 |

## The decision record

`logs/decisions_record.csv` holds every analytical decision taken during the
analysis with its rationale, including the ones taken after inspecting the data
and labelled as such. It is in the repository because a negative result is worth
only as much as a reader's ability to check how it was reached.

It records, among others: why probes are collapsed by MaxMean rather than
averaged; why the 33-gene panel is not trimmed despite three members sitting at
background; why four outlier exclusion sets are reported separately and never
pooled; why one specimen with contradictory metadata is retained; and why the
composite index was demoted from deliverable to instrument.

## Reproducibility

Every value in the manuscript was written to a file by this code and read back
from that file. None was transcribed from console output.

The random number generator is fixed with `RNGkind("Mersenne-Twister",
"Inversion", "Rejection")` and the seed set to 42 at the start of every script
and again immediately before each stochastic call, so results do not depend on
execution order or on the R version's default sampling algorithm.

Produced under R 4.5.3 with Bioconductor 3.22. Exact package versions are in
`logs/session_info.txt`.

## Citation

If you use this code or these results, cite the paper. Citation metadata is in
`CITATION.cff`.

## Licence

Code is MIT. Figures, tables and the decision record are CC BY 4.0. See
`LICENSE` and `LICENSE-CC-BY-4.0.md`.

The underlying expression data belongs to the original depositors and is subject
to the terms of the Gene Expression Omnibus. This repository redistributes none
of it.
