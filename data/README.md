# Data

This folder is empty in the repository. Nothing here is redistributed.

The four expression series belong to their original depositors and are subject
to the terms of the Gene Expression Omnibus. Script 01 retrieves them with
GEOquery and writes them here.

```
data/
  raw/         series matrices and platform annotation as downloaded
  processed/   log2-scaled matrices, gene-level matrices, groups, scores
```

Both subfolders are in `.gitignore`. Expect roughly 2 GB after a full run,
mostly the two transcriptome-wide Affymetrix series.

## The series

| Accession | Platform | Preparation | n | Groups |
| --- | --- | --- | --- | --- |
| GSE63514 | GPL570 | Laser-capture-microdissected epithelium | 128 | Normal 24, CIN1 14, CIN2 22, CIN3 40, Cancer 28 |
| GSE138080 | GPL4133 | Bulk tissue | 35 | Normal 10, HGCIN 15, SCC 10 |
| GSE27678 | GPL571 | Exfoliated cytology | 44 | Normal 12, LSIL 11, HSIL 21 |
| GSE27678 | GPL570 | Exfoliated cytology | 31 | Normal 3, Cancer 28 |

GSE27678 contains two sub-platforms. Script 02 separates them by querying the
platform annotation of each returned expression set, not by position in the
list GEOquery returns, which is not guaranteed stable.

## Two things script 01 and 02 handle that are easy to miss

**GSE138080 is deposited on a linear scale**, median 110 and maximum 650,658,
while the other three are already log2. Script 01b transforms it. Running the
analysis without that step compares a linear matrix against three log ones.

**Two GPL570 specimens are immortalised cell lines**, MS751 (GSM685543) and
CaSki (GSM685549). Both are excluded, because a cell line has no stratified
epithelium and so cannot carry the differentiation gradient the panel measures.
That is why the GPL570 set is 31 and not 33.

## If GEO is unreachable

Script 01 falls back through three routes: a cached `.rds` if one exists, then
the series matrix, then a full download. If all three fail the script stops
rather than continuing on partial data.

The series can also be fetched by hand from
`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63514` and placed in
`data/raw/` under the names script 01 expects, which it prints on failure.
