# Phase 13 - Figure 2A, rank position of the panel genes  (revision 2)
#
# Fixes from the first draft: label only the five genes the Results text
# names, so labels stop colliding across rows; jitter the points slightly so
# the Differentiation clump resolves into seven; carry the count of extreme
# members in the margin beside the enrichment statistic.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(ggplot2)
library(ggrepel)

ranks <- read.csv("results/tables/fgsea_member_ranks.csv", stringsAsFactors = FALSE)
fg <- read.csv("results/tables/fgsea_modules_CancervNormal.csv", stringsAsFactors = FALSE)

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
module_cols <- c(Synthesis = "#0072B2", Degradation = "#D55E00",
                 Warburg = "#009E73", Differentiation = "#CC79A7")

N <- 21755
ranks$module <- factor(ranks$module, levels = rev(module_order))
ranks$y <- as.numeric(ranks$module)

# label only the genes carried in the Results text
named <- c("GYS2", "PPP1R3C", "PYGB", "SLC2A3", "CRNN")
ranks$lab <- ifelse(ranks$gene %in% named, ranks$gene, "")

# small deterministic vertical offset so overlapping points resolve
set.seed(42)
ranks$yj <- ranks$y + runif(nrow(ranks), -0.17, 0.17)

# right margin: enrichment statistic, and how many members sit in a tail
ann <- data.frame(module = factor(module_order, levels = rev(module_order)),
                  stringsAsFactors = FALSE)
ann$y <- as.numeric(ann$module)
ann$NES <- fg$NES[match(as.character(ann$module), fg$module)]
ann$padj <- fg$padj[match(as.character(ann$module), fg$module)]
ann$stat <- sprintf("NES %+.2f, FDR %s", ann$NES,
                    ifelse(ann$padj < 0.001,
                           formatC(ann$padj, format = "e", digits = 1),
                           sprintf("%.3f", ann$padj)))
ann$count <- sapply(as.character(ann$module), function(m) {
  sub <- ranks[as.character(ranks$module) == m, ]
  sprintf("%d of %d beyond the 90th percentile",
          sum(sub$percentile > 90 | sub$percentile < 10), nrow(sub))
})

XMAX <- N * 1.34

p <- ggplot(ranks, aes(x = rank, y = yj, colour = module)) +
  geom_hline(yintercept = 1:4, linewidth = 0.25, colour = "grey93") +
  geom_vline(xintercept = N / 2, linewidth = 0.3, linetype = "dashed",
             colour = "grey75") +
  geom_point(size = 2.1, alpha = 0.9) +
  geom_text_repel(aes(label = lab), size = 2.4, colour = "grey15",
                  segment.size = 0.2, segment.colour = "grey55",
                  min.segment.length = 0, box.padding = 0.4,
                  nudge_y = 0.28, direction = "both",
                  max.overlaps = 20, show.legend = FALSE) +
  geom_text(data = ann, aes(x = N * 1.03, y = y + 0.10, label = stat),
            hjust = 0, size = 2.5, colour = "grey25", inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = N * 1.03, y = y - 0.16, label = count),
            hjust = 0, size = 2.2, colour = "grey50", inherit.aes = FALSE) +
  scale_colour_manual(values = module_cols, guide = "none") +
  scale_x_continuous(limits = c(-400, XMAX),
                     breaks = c(1, 5000, 10000, 15000, 20000, N),
                     labels = c("1", "5,000", "10,000", "15,000", "20,000", "21,755"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:4, labels = levels(ranks$module),
                     limits = c(0.45, 4.6), expand = c(0, 0)) +
  labs(x = "rank in the Cancer vs Normal ranked transcriptome", y = NULL,
       title = "A   Position of the 33 panel genes among 21,755 ranked genes",
       subtitle = "rank 1 = most up-regulated in cancer; 21,755 = most down-regulated") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0),
        plot.subtitle = element_text(size = 7.5, colour = "grey45",
                                     margin = margin(b = 8)),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks.x = element_line(linewidth = 0.3, colour = "grey30"),
        axis.text = element_text(colour = "grey20"),
        axis.text.y = element_text(size = 8.5, colour = "grey10"),
        panel.grid = element_blank(),
        plot.margin = margin(6, 6, 6, 6))

ggsave("results/figures/deg/Figure2A_rank_strip_v2.png", p,
       width = 9.6, height = 3.6, dpi = 400)
ggsave("results/figures/deg/Figure2A_rank_strip_v2.pdf", p,
       width = 9.6, height = 3.6)

cat("wrote Figure2A_rank_strip_v2.png and .pdf\n\n")

cat("=== Rank summary by module ===\n\n")
for (m in module_order) {
  sub <- ranks[as.character(ranks$module) == m, ]
  cat(sprintf("%-16s n=%d  median percentile %.1f  range %.1f to %.1f\n",
              m, nrow(sub), median(sub$percentile),
              min(sub$percentile), max(sub$percentile)))
  hi <- sub$gene[sub$percentile > 90]
  lo <- sub$gene[sub$percentile < 10]
  cat("   above 90th: ", ifelse(length(hi) == 0, "none", paste(hi, collapse = ", ")), "\n", sep = "")
  cat("   below 10th: ", ifelse(length(lo) == 0, "none", paste(lo, collapse = ", ")), "\n", sep = "")
}

write.csv(ranks[order(ranks$rank), c("module", "gene", "rank", "percentile", "metric")],
          "results/tables/Figure2A_rank_positions.csv", row.names = FALSE)
cat("\nWritten: results/tables/Figure2A_rank_positions.csv\n")
