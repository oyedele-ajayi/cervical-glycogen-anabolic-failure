# Phase 12b - complete the decision record
#
# D1-D7 and D10 were written to file as they were taken. D8, D9, D11, D12,
# D13, D14 and D15 were settled during review and never recorded. Added here
# so the file, not the conversation, is the authority.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

d <- read.csv("logs/decisions_record.csv", stringsAsFactors = FALSE)

new <- data.frame(
  id = c("D8", "D9", "D11", "D12", "D13", "D14", "D15"),
  decision = c(
    
    "GSE27678 GPL570 sample GSM685524 carries contradictory metadata: title and characteristics say normal uterine cervix and type Normal, while source_name_ch1 says Squamous cell carcinoma. Two fields against one, so it is retained as Normal in the primary analysis. This is also the conservative choice, since including it lowers the AUC. It is additionally the sample generating all three discordant pairs in that cohort, so the expression data agrees with the minority label. The conflict is reported, and a post hoc analysis excluding it is labelled post hoc and kept separate from every pre-specified analysis.",
    
    "Outlier exclusion sets. The four samples flagged by the pre-specified criteria (GSM1551339, GSM1551365, GSM1551395, GSM1551419) constitute the pre-specified sensitivity analysis. The five named in the earlier manuscript (GSM1551334, GSM1551316, GSM1551349, GSM1551350, GSM1551351) meet neither stated criterion and are reported separately as an ad hoc analysis for continuity. GSM1551331 and GSM685524 are each post hoc and separate again. The four sets are never merged. Note also that no sample met the PCA criterion; every flag came from the correlation criterion, and the threshold was not moved after inspection.",
    
    "The 33-gene panel is not trimmed. PHKG1 (max probe AveExpr 1.64) and PPP1R3A (2.96) sit at or near background in the discovery cohort and contribute noise rather than signal. Trimming after inspecting detection levels would be post hoc panel surgery. Instead, detection status is reported per gene per cohort in Supplementary Table S1, and the leave-one-gene-out analysis demonstrates the point empirically: removing PPP1R3A changes the CIN3+ AUC by -0.0002 and removing PHKG1 by +0.0022.",
    
    "The Checkpoint 4a hard stop on transcriptome-wide DEG counts is retired. It was justified on limma being version-stable, which holds, but it was conditional on an unchanged input matrix. D10 changed the input matrix deliberately. CIN2vNormal came in at 311 against a reference of 266, and it is the weakest contrast (minimum FDR 4.5e-6), so a modest absolute change lands hardest there proportionally. The two structurally diagnostic counts are unchanged: CIN1vNormal 1 and CIN2vCIN1 0.",
    
    "Monotonicity is not claimed for either GDI form, in either direction. Under MaxMean collapse the discovery medians increase monotonically; under averaging they reverse at Normal to CIN1. Same data, same index, different probe rule, so monotonicity is an artefact of an analytic choice. The Spearman trend is reported as the evidence. The Normal to CIN1 step is not significant under either rule (Mann-Whitney p = 0.44; bootstrap 95% CI on the median difference [-0.074, +0.141]), and the text states that median ordering at that step is sensitive to probe collapse.",
    
    "The composite index is demoted from deliverable to summary statistic. Differentiation alone gives a CIN3+ AUC of 0.8319 against the composite's 0.8324, and equals or beats the composite in three of four cohorts, so no classification advantage is claimed for combining modules. The 33-gene panel and all four modules are retained, because the negative results are load-bearing: one cannot report that catabolism does not move without having measured it. Degradation alone gives 0.5044 and Warburg 0.5689, both with CIs spanning 0.5, and Warburg separates from Normal at no grade (p 0.45 to 0.62) and shows null enrichment (NES +0.884, padj 0.790). Module dropout moves from supplementary into the main text as a primary result.",
    
    "Cohort framing. GSE63514 is the primary analysis and carries the headline estimates; it is the only cohort with adequate group sizes at every threshold. GSE138080 and GSE27678 GPL571 are replication cohorts in which the trend test is the primary evidence because it uses every sample, and each AUC is reported with its confidence interval and both group sizes in the same sentence. GSE27678 GPL570 is a binary concordance check, not a validation: it rests on three negatives and two discordant pairs out of eighty-four. It is reported once with n = 3 negatives stated inline and is excluded from the abstract."),
  
  stringsAsFactors = FALSE)

d <- rbind(d, new)
write.csv(d, "logs/decisions_record.csv", row.names = FALSE)

cat("Decision record now holds", nrow(d), "entries:\n")
cat("  ", paste(d$id, collapse = ", "), "\n\n")

cat("Rerun 12_supplementary_tables.R so sheet S12 picks these up.\n")
