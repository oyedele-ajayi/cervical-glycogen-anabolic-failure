# Phase 11 - Figure 1 (volcanoes) and Figure 2 (panel heatmap)
#
# Figure 1  four volcano panels, CIN1/CIN2/CIN3/Cancer each against Normal.
#           Grey = not significant, blue = FDR<0.05 and |log2FC|>1,
#           vermillion = panel gene. Panel genes labelled where significant.
#
# Figure 2  33-gene heatmap, discovery cohort. Columns ordered
#           Normal -> CIN1 -> CIN2 -> CIN3 -> Cancer, and within each grade by
#           ascending GDI_diff. Rows in module order. Z-scored within gene,
#           clipped to -3..+3.
#
# Grade order is taken from the saved factor levels, never rebuilt after a
# read.csv. That is where the previous run's two ordering bugs came from.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(ggplot2)
library(patchwork)
library(ggrepel)
library(pheatmap)

deg <- readRDS("data/processed/DEG_discovery.rds")
gene_set <- readRDS("data/processed/gene_set_locked.rds")
gene_mats <- readRDS("data/processed/gene_matrices.rds")
groups <- readRDS("data/processed/groups.rds")
gdi_list <- readRDS("data/processed/GDI_results.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
panel_genes <- unlist(gene_set, use.names = FALSE)

COL_NS <- "grey78"
COL_SIG <- "#0072B2"      # blue
COL_PANEL <- "#D55E00"    # vermillion

base_theme <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 4)),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        axis.text = element_text(colour = "grey20"),
        axis.title = element_text(colour = "grey20"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = "grey20"),
        plot.margin = margin(5, 8, 5, 5))


# ---------------------------------------------------------------------------
# Figure 1 - volcanoes
# ---------------------------------------------------------------------------

contrasts <- c("CIN1vNormal", "CIN2vNormal", "CIN3vNormal", "CancervNormal")
titles <- c("A   CIN1 vs Normal", "B   CIN2 vs Normal",
            "C   CIN3 vs Normal", "D   Cancer vs Normal")

full <- deg$full
xlim_v <- max(abs(full$logFC[full$contrast %in% contrasts]))
ylim_v <- max(-log10(full$adj.P.Val[full$contrast %in% contrasts]))

panels1 <- list()
for (i in seq_along(contrasts)) {
  
  d <- full[full$contrast == contrasts[i], ]
  d$y <- -log10(d$adj.P.Val)
  d$is_panel <- d$gene %in% panel_genes
  d$sig <- d$adj.P.Val < 0.05 & abs(d$logFC) > 1
  
  d$class <- "not significant"
  d$class[d$sig] <- "FDR<0.05 & |log2FC|>1"
  d$class[d$is_panel] <- "panel gene"
  d$class <- factor(d$class,
                    levels = c("not significant", "FDR<0.05 & |log2FC|>1", "panel gene"))
  
  lab <- d[d$is_panel & d$sig, ]
  
  p <- ggplot(d[!d$is_panel, ], aes(x = logFC, y = y, colour = class)) +
    geom_point(size = 0.45, alpha = 0.55) +
    geom_point(data = d[d$is_panel, ], size = 1.5, alpha = 0.95) +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.25,
               linetype = "dashed", colour = "grey55") +
    geom_vline(xintercept = c(-1, 1), linewidth = 0.25,
               linetype = "dashed", colour = "grey55") +
    scale_colour_manual(values = c("not significant" = COL_NS,
                                   "FDR<0.05 & |log2FC|>1" = COL_SIG,
                                   "panel gene" = COL_PANEL),
                        drop = FALSE) +
    coord_cartesian(xlim = c(-xlim_v, xlim_v), ylim = c(0, ylim_v)) +
    labs(title = titles[i], x = "log2 fold change", y = "-log10 FDR") +
    base_theme
  
  if (nrow(lab) > 0) {
    p <- p + geom_text_repel(data = lab, aes(label = gene), size = 2.2,
                             colour = "grey15", segment.size = 0.2,
                             min.segment.length = 0, max.overlaps = 30,
                             show.legend = FALSE)
  }
  panels1[[i]] <- p
}

fig1 <- (panels1[[1]] | panels1[[2]]) / (panels1[[3]] | panels1[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("results/figures/deg/Figure1_volcanoes_v2.png", fig1,
       width = 9.0, height = 7.2, dpi = 400)
ggsave("results/figures/deg/Figure1_volcanoes_v2.pdf", fig1,
       width = 9.0, height = 7.2)
cat("wrote Figure 1\n")


# ---------------------------------------------------------------------------
# Figure 2 - 33-gene heatmap
# ---------------------------------------------------------------------------

m <- gene_mats$GSE63514
grp <- groups$GSE63514                       # factor, levels already correct
gdi <- gdi_list$GSE63514

# column order: by grade, then by ascending GDI_diff within grade
ord <- order(as.integer(grp), gdi$GDI_diff[match(colnames(m), gdi$sample)])
m <- m[, ord]
grp <- grp[ord]

# row order: module order, genes within module as locked
row_genes <- character(0)
row_module <- character(0)
for (mod in module_order) {
  present <- gene_set[[mod]][gene_set[[mod]] %in% rownames(m)]
  row_genes <- c(row_genes, present)
  row_module <- c(row_module, rep(mod, length(present)))
}

mat <- m[row_genes, ]

# z-score within gene, clipped
z <- t(scale(t(mat)))
z[z > 3] <- 3
z[z < -3] <- -3

ann_col <- data.frame(Grade = grp)
rownames(ann_col) <- colnames(z)

ann_row <- data.frame(Module = factor(row_module, levels = module_order))
rownames(ann_row) <- row_genes

grade_cols <- c(Normal = "#DEEBF7", CIN1 = "#9ECAE1", CIN2 = "#4292C6",
                CIN3 = "#2171B5", Cancer = "#08306B")
module_cols <- c(Synthesis = "#0072B2", Degradation = "#D55E00",
                 Warburg = "#009E73", Differentiation = "#CC79A7")

png("results/figures/gdi/Figure2_heatmap_v2.png",
    width = 2800, height = 1700, res = 300)
pheatmap(z,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, fontsize_row = 7,
         annotation_col = ann_col, annotation_row = ann_row,
         annotation_colors = list(Grade = grade_cols, Module = module_cols),
         gaps_row = cumsum(table(factor(row_module, levels = module_order)))[1:3],
         gaps_col = cumsum(table(grp))[1:4],
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
         breaks = seq(-3, 3, length.out = 102),
         border_color = NA,
         main = "")
dev.off()

pdf("results/figures/gdi/Figure2_heatmap_v2.pdf", width = 9.3, height = 5.7)
pheatmap(z,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_colnames = FALSE, fontsize_row = 7,
         annotation_col = ann_col, annotation_row = ann_row,
         annotation_colors = list(Grade = grade_cols, Module = module_cols),
         gaps_row = cumsum(table(factor(row_module, levels = module_order)))[1:3],
         gaps_col = cumsum(table(grp))[1:4],
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
         breaks = seq(-3, 3, length.out = 102),
         border_color = NA, main = "")
dev.off()
cat("wrote Figure 2\n")

# module counts for the legend, computed rather than typed
counts <- table(factor(row_module, levels = module_order))
cat("\nModule counts for the Figure 2 legend: ",
    paste(names(counts), counts, sep = " ", collapse = " / "), "\n", sep = "")


# ---------------------------------------------------------------------------
# Legends
# ---------------------------------------------------------------------------

legends <- c(
  "Figure 1. Transcriptome-wide differential expression across disease grade,",
  "discovery cohort (GSE63514, n = 128: Normal 24, CIN1 14, CIN2 22, CIN3 40,",
  "Cancer 28). Each panel contrasts one grade against Normal across 21,755",
  "genes. Grey, not significant; blue, FDR < 0.05 and |log2FC| > 1; vermillion,",
  "a member of the 33-gene panel. Panel genes meeting both thresholds are",
  "labelled. Axes are shared across panels. Dashed lines mark FDR 0.05 and",
  "|log2FC| = 1.",
  "",
  paste0("Figure 2. Expression of the 33 panel genes across disease grade,",
         " discovery cohort."),
  "Columns are the 128 samples ordered by grade and, within each grade, by",
  "ascending GDI_diff. Rows follow the locked module order: Synthesis (9),",
  paste0("Degradation (", counts[["Degradation"]], "), Warburg (",
         counts[["Warburg"]], "), Differentiation (", counts[["Differentiation"]], ")."),
  "Values are z-scores computed within gene and clipped to -3 to +3.")

writeLines(legends, "results/figures/Figure1_Figure2_legends.txt")

cat("\nWritten: Figure1_volcanoes_v2.png/.pdf, Figure2_heatmap_v2.png/.pdf,\n")
cat("  results/figures/Figure1_Figure2_legends.txt\n")
