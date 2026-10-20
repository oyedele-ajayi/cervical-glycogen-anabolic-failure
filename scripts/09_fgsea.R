# Phase 9 - pathway enrichment on the Cancer vs Normal contrast
#
# Ranking metric: sign(log2FC) * -log10(P + 1e-300), using the nominal P,
# not the FDR. Non-finite and exactly-zero values dropped, sorted decreasing.
#
# Open question O3. Because the Degradation NES is negative, fgsea's leading
# edge reports only the down-regulated members and structurally excludes the
# up-regulated ones. So the leading edge is not a list of the module's moving
# genes. We therefore also report, per module, how many members are
# significant in each direction, and the rank position of every member.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42)

library(fgsea)

deg <- readRDS("data/processed/DEG_discovery.rds")
gene_set <- readRDS("data/processed/gene_set_locked.rds")
cohort_genesets <- readRDS("data/processed/cohort_genesets.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
gs <- cohort_genesets$GSE63514[module_order]

cvn <- deg$full[deg$full$contrast == "CancervNormal", ]


# ---------------------------------------------------------------------------
# 1. Build the ranked list
# ---------------------------------------------------------------------------

metric <- sign(cvn$logFC) * -log10(cvn$P.Value + 1e-300)
names(metric) <- cvn$gene

keep <- is.finite(metric) & metric != 0
dropped <- sum(!keep)
metric <- metric[keep]
metric <- sort(metric, decreasing = TRUE)

cat("=== CHECKPOINT 8a: ranked list ===\n\n")
cat("Genes in contrast      :", nrow(cvn), "\n")
cat("Dropped (zero/non-finite):", dropped, "\n")
cat("Ranked list length     :", length(metric), "\n\n")

cat("Top 5 (most up in cancer):\n")
print(round(head(metric, 5), 4))
cat("\nBottom 5 (most down in cancer):\n")
print(round(tail(metric, 5), 4))

write.csv(data.frame(gene = names(metric), rank_metric = round(metric, 6),
                     rank = seq_along(metric)),
          "results/tables/fgsea_ranked_list.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 2. fgsea
# ---------------------------------------------------------------------------

set.seed(42)
res <- fgsea(pathways = gs, stats = metric,
             minSize = 5, maxSize = 50, nPermSimple = 10000)

res <- as.data.frame(res)
res$pathway <- factor(res$pathway, levels = module_order)
res <- res[order(res$pathway), ]

cat("\n\n=== CHECKPOINT 8b: module enrichment ===\n\n")
out <- data.frame(module = as.character(res$pathway),
                  size = res$size,
                  NES = round(res$NES, 4),
                  p_nominal = signif(res$pval, 4),
                  padj = signif(res$padj, 4),
                  stringsAsFactors = FALSE)
print(out, row.names = FALSE)

cat("\nLeading-edge genes per module:\n")
for (i in seq_len(nrow(res))) {
  le <- res$leadingEdge[[i]]
  cat(sprintf("  %-16s (%d genes) %s\n", as.character(res$pathway[i]),
              length(le), paste(le, collapse = ", ")))
}

res_out <- out
res_out$leading_edge <- sapply(res$leadingEdge, paste, collapse = ", ")
res_out$leading_edge_n <- sapply(res$leadingEdge, length)
write.csv(res_out, "results/tables/fgsea_modules_CancervNormal.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 3. O3 - significant members by direction
# ---------------------------------------------------------------------------
# This is what the leading edge cannot tell you.

cat("\n\n=== CHECKPOINT 8c: significant members by direction (O3) ===\n\n")

dir_tab <- NULL
for (mod in module_order) {
  members <- gs[[mod]]
  sub <- cvn[cvn$gene %in% members, ]
  
  sig <- sub[sub$adj.P.Val < 0.05, ]
  up <- sig$gene[sig$logFC > 0]
  dn <- sig$gene[sig$logFC < 0]
  
  dir_tab <- rbind(dir_tab, data.frame(
    module = mod, n_members = length(members),
    n_sig = nrow(sig),
    n_sig_up = length(up), n_sig_down = length(dn),
    genes_up = paste(up, collapse = ", "),
    genes_down = paste(dn, collapse = ", "),
    NES = round(res$NES[res$pathway == mod], 3),
    stringsAsFactors = FALSE))
}

print(dir_tab[, c("module", "n_members", "n_sig", "n_sig_up", "n_sig_down", "NES")],
      row.names = FALSE)

cat("\nDetail:\n")
for (i in seq_len(nrow(dir_tab))) {
  cat("  ", dir_tab$module[i], "\n", sep = "")
  cat("    up  : ", ifelse(dir_tab$genes_up[i] == "", "none", dir_tab$genes_up[i]), "\n", sep = "")
  cat("    down: ", ifelse(dir_tab$genes_down[i] == "", "none", dir_tab$genes_down[i]), "\n", sep = "")
}

cat("\nWhere a module has significant members in both directions, the module\n")
cat("null is a null about coordination, not about activity. The leading edge\n")
cat("for a negative NES lists only down-regulated members, so Table 3 needs a\n")
cat("footnote defining it and this direction split as a column.\n")

write.csv(dir_tab, "results/tables/fgsea_direction_split.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# 4. Rank position of every module member
# ---------------------------------------------------------------------------

cat("\n\n=== CHECKPOINT 8d: rank positions ===\n\n")

ranks <- data.frame(gene = names(metric), rank = seq_along(metric),
                    metric = as.numeric(metric), stringsAsFactors = FALSE)
N <- nrow(ranks)

pos <- NULL
for (mod in module_order) {
  for (g in gs[[mod]]) {
    r <- ranks[ranks$gene == g, ]
    if (nrow(r) == 0) {
      pos <- rbind(pos, data.frame(module = mod, gene = g, rank = NA,
                                   percentile = NA, metric = NA,
                                   in_leading_edge = NA, stringsAsFactors = FALSE))
    } else {
      le <- res$leadingEdge[[which(res$pathway == mod)]]
      pos <- rbind(pos, data.frame(
        module = mod, gene = g, rank = r$rank[1],
        percentile = round(100 * r$rank[1] / N, 1),
        metric = round(r$metric[1], 3),
        in_leading_edge = g %in% le,
        stringsAsFactors = FALSE))
    }
  }
}

for (mod in module_order) {
  cat("--- ", mod, " ---\n", sep = "")
  sub <- pos[pos$module == mod, ]
  sub <- sub[order(sub$rank), ]
  print(sub[, c("gene", "rank", "percentile", "metric", "in_leading_edge")],
        row.names = FALSE)
  cat("\n")
}

cat("Ranked list runs 1 (most up in cancer) to", N, "(most down).\n")

write.csv(pos, "results/tables/fgsea_member_ranks.csv", row.names = FALSE)

saveRDS(list(metric = metric, fgsea = res, direction = dir_tab, ranks = pos),
        "data/processed/fgsea_results.rds")

cat("\nWritten: fgsea_ranked_list.csv, fgsea_modules_CancervNormal.csv,\n")
cat("  fgsea_direction_split.csv, fgsea_member_ranks.csv\n")
