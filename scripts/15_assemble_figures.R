# Phase 15 - assemble the final figure set under the new numbering
#
# Figure 1  module trajectories, four cohorts            (was Figure 3)
# Figure 2  A rank strip + B 33-gene heatmap             (new assembly)
# Figure 3  A module separation + B forest               (built in script 14)
# Figure 4  ROC curves, four cohorts, A to D             (panel E removed;
#           single-module discrimination now lives in Figure 3B)
#
# S1  quality control
# S2  volcanoes, four panels                             (was Figure 1)
# S3  robustness: outlier comparison and LOO distribution
# S4  fgsea running enrichment, four modules
#
# Figure 2 combines a ggplot panel with a pheatmap grob, so it is assembled
# with cowplot rather than patchwork.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42)

for (p in c("cowplot", "ggplotify")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(ggplot2)
library(patchwork)
library(pROC)
library(pheatmap)
library(cowplot)
library(ggplotify)
library(ggrepel)

# refuse to run if an output is locked by a viewer
check_writable <- function(path) {
  if (file.exists(path) && file.access(path, 2) != 0)
    stop("Cannot write ", path, " - close whatever has it open.")
}

gdi_list <- readRDS("data/processed/GDI_results.rds")
roc_objects <- readRDS("data/processed/ROC_objects.rds")
gene_set <- readRDS("data/processed/gene_set_locked.rds")
gene_mats <- readRDS("data/processed/gene_matrices.rds")
groups <- readRDS("data/processed/groups.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
module_cols <- c(Synthesis = "#0072B2", Degradation = "#D55E00",
                 Warburg = "#009E73", Differentiation = "#CC79A7")

base_theme <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 4)),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        axis.text = element_text(colour = "grey20"),
        axis.title = element_text(colour = "grey20"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = "grey20"),
        plot.margin = margin(5, 8, 5, 5))


# ===========================================================================
# FIGURE 2 - A rank strip, B heatmap
# ===========================================================================

# --- panel A ---------------------------------------------------------------

ranks <- read.csv("results/tables/fgsea_member_ranks.csv", stringsAsFactors = FALSE)
fg <- read.csv("results/tables/fgsea_modules_CancervNormal.csv", stringsAsFactors = FALSE)

N <- 21755
ranks$module <- factor(ranks$module, levels = rev(module_order))
ranks$y <- as.numeric(ranks$module)
named <- c("GYS2", "PPP1R3C", "PYGB", "SLC2A3", "CRNN")
ranks$lab <- ifelse(ranks$gene %in% named, ranks$gene, "")
set.seed(42)
ranks$yj <- ranks$y + runif(nrow(ranks), -0.17, 0.17)

ann <- data.frame(module = factor(module_order, levels = rev(module_order)))
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

p2a <- ggplot(ranks, aes(x = rank, y = yj, colour = module)) +
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
            hjust = 0, size = 2.4, colour = "grey25", inherit.aes = FALSE) +
  geom_text(data = ann, aes(x = N * 1.03, y = y - 0.16, label = count),
            hjust = 0, size = 2.1, colour = "grey50", inherit.aes = FALSE) +
  scale_colour_manual(values = module_cols, guide = "none") +
  scale_x_continuous(limits = c(-400, N * 1.34),
                     breaks = c(1, 5000, 10000, 15000, 20000, N),
                     labels = c("1", "5,000", "10,000", "15,000", "20,000", "21,755"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:4, labels = levels(ranks$module),
                     limits = c(0.45, 4.6), expand = c(0, 0)) +
  labs(x = "rank in the Cancer vs Normal ranked transcriptome", y = NULL,
       title = "A   Position of the 33 panel genes among 21,755 ranked genes") +
  base_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 8.5, colour = "grey10"))

# --- panel B ---------------------------------------------------------------

m <- gene_mats$GSE63514
grp <- groups$GSE63514
gdi <- gdi_list$GSE63514

ord <- order(as.integer(grp), gdi$GDI_diff[match(colnames(m), gdi$sample)])
m <- m[, ord]; grp <- grp[ord]

row_genes <- character(0); row_module <- character(0)
for (mod in module_order) {
  present <- gene_set[[mod]][gene_set[[mod]] %in% rownames(m)]
  row_genes <- c(row_genes, present)
  row_module <- c(row_module, rep(mod, length(present)))
}

z <- t(scale(t(m[row_genes, ])))
LIM <- as.numeric(round(quantile(abs(z), 0.98, na.rm = TRUE), 1))
z[z > LIM] <- LIM; z[z < -LIM] <- -LIM

ann_col <- data.frame(Grade = grp); rownames(ann_col) <- colnames(z)
ann_row <- data.frame(Module = factor(row_module, levels = module_order))
rownames(ann_row) <- row_genes

grade_cols <- c(Normal = "#DEEBF7", CIN1 = "#9ECAE1", CIN2 = "#4292C6",
                CIN3 = "#2171B5", Cancer = "#08306B")
ramp <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))(255)

hm <- pheatmap(z, cluster_rows = FALSE, cluster_cols = FALSE,
               show_colnames = FALSE, fontsize_row = 7,
               annotation_col = ann_col, annotation_row = ann_row,
               annotation_colors = list(Grade = grade_cols, Module = module_cols),
               annotation_names_col = FALSE, annotation_names_row = FALSE,
               gaps_row = cumsum(table(factor(row_module, levels = module_order)))[1:3],
               gaps_col = cumsum(table(grp))[1:4],
               color = ramp, breaks = seq(-LIM, LIM, length.out = 256),
               border_color = NA,
               legend_breaks = c(-LIM, 0, LIM),
               silent = TRUE, main = "")

p2b <- as.ggplot(hm) +
  labs(title = "B   Panel gene expression across grade, discovery cohort") +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 3)))

fig2 <- plot_grid(p2a, p2b, ncol = 1, rel_heights = c(0.78, 1))

check_writable("results/figures/deg/Figure1_panel_genes.png")
ggsave("results/figures/deg/Figure1_panel_genes.png", fig2,
       width = 9.6, height = 9.0, dpi = 400)
ggsave("results/figures/deg/Figure1_panel_genes.pdf", fig2,
       width = 9.6, height = 9.0)
cat("wrote Figure 1 (panel genes)\n")


# ===========================================================================
# FIGURE 4 - ROC curves only, panels A to D
# ===========================================================================

roc_panel <- function(cohort, labels, title) {
  d <- NULL; lv <- character(0)
  for (lab in labels) {
    r <- roc_objects[[paste(cohort, lab, "GDI_diff", sep = "|")]]
    s <- sprintf("%s  (AUC %.3f)", lab, as.numeric(r$auc))
    lv <- c(lv, s)
    d <- rbind(d, data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities,
                             series = s, stringsAsFactors = FALSE))
  }
  d$series <- factor(d$series, levels = lv)
  d <- d[order(d$series, d$fpr), ]
  cols <- c("#0072B2", "#D55E00", "#CC79A7")[seq_along(lv)]
  names(cols) <- lv
  
  ggplot(d, aes(x = fpr, y = tpr, colour = series)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3,
                linetype = "dashed", colour = "grey70") +
    geom_path(linewidth = 0.7) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "1 - specificity", y = "sensitivity") +
    base_theme +
    theme(legend.position = c(0.97, 0.06), legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key.height = unit(3.4, "mm"),
          legend.text = element_text(size = 6.8))
}

fig4 <- (roc_panel("GSE63514", c("CIN2+", "CIN3+", "Cancer"), "A   GSE63514 (n=128)") |
           roc_panel("GSE138080", c("HGCIN+", "SCC"), "B   GSE138080 (n=35)")) /
  (roc_panel("GSE27678_GPL571", c("LSIL+", "HSIL"), "C   GSE27678 GPL571 (n=44)") |
     roc_panel("GSE27678_GPL570", c("Cancer"), "D   GSE27678 GPL570 (n=31; only 3 normals)"))

check_writable("results/figures/validation/Figure4_ROC.png")
ggsave("results/figures/validation/Figure4_ROC.png", fig4,
       width = 8.2, height = 7.4, dpi = 400)
ggsave("results/figures/validation/Figure4_ROC.pdf", fig4,
       width = 8.2, height = 7.4)
cat("wrote Figure 4\n")


# ===========================================================================
# Manifest, so the numbering is recorded rather than remembered
# ===========================================================================

manifest <- data.frame(
  figure = c("Figure 1", "Figure 2", "Figure 3", "Figure 4",
             "Figure S1", "Figure S2", "Figure S3", "Figure S4"),
  content = c("A: rank position of the 33 panel genes; B: panel gene heatmap, discovery",
              "A: grade at which each module separates from normal; B: discrimination by module and reduced index",
              "Module score trajectories across grade, four cohorts",
              "ROC curves for the index, four cohorts",
              "Quality control: boxplots, density, PCA",
              "Volcano plots, four grade contrasts against normal",
              "fgsea running enrichment, four modules",
              "Robustness: outlier exclusion sets and leave-one-gene-out"),
  file = c("results/figures/deg/Figure1_panel_genes.png",
           "results/figures/gdi/Figure2_signal_location.png",
           "results/figures/gdi/Figure3_module_trajectories.png",
           "results/figures/validation/Figure4_ROC.png",
           "results/figures/qc/FigureS1_QC_<cohort>.png",
           "results/figures/deg/FigureS2_volcanoes.png",
           "results/figures/gdi/FigureS3_fgsea.png",
           "results/figures/validation/FigureS4_robustness.png"),
  script = c("15_assemble_figures.R", "14_figure3.R", "10_figures_main.R",
             "15_assemble_figures.R", "03_QC.R", "11_figures_deg_heatmap.R",
             "16_supplementary_figures.R", "16_supplementary_figures.R"),
  first_cited = c("3.2", "3.3", "3.4", "3.5", "3.1", "3.2", "3.6", "3.7"),
  stringsAsFactors = FALSE)

write.csv(manifest, "results/figures/figure_manifest.csv", row.names = FALSE)

cat("\n=== Figure manifest ===\n")
print(manifest[, c("figure", "file", "first_cited")], row.names = FALSE)

# Three things have to agree, and all three have gone wrong at some point:
# the file must exist, its name must carry its number, and the numbers must
# run in order of first citation. The last one is the rule journals apply and
# the one nothing here used to check.
cat("\n=== Manifest check ===\n")
for (i in seq_len(nrow(manifest))) {
  f <- manifest$file[i]
  paths <- if (grepl("<cohort>", f, fixed = TRUE)) {
    Sys.glob(sub("<cohort>", "*", f, fixed = TRUE))
  } else f
  ok <- length(paths) > 0 && all(file.exists(paths))
  num <- sub(" ", "", manifest$figure[i])
  named <- ok && all(grepl(num, basename(paths), fixed = TRUE))
  cat(sprintf("  %-10s %-8s %s\n", manifest$figure[i],
              ifelse(ok, "present", "MISSING"),
              ifelse(!ok, "", ifelse(named, "name matches the number",
                                     "NAME DOES NOT MATCH THE NUMBER"))))
}

for (set in list(c("Figure 1", "Figure 2", "Figure 3", "Figure 4"),
                 c("Figure S1", "Figure S2", "Figure S3", "Figure S4"))) {
  k <- match(set, manifest$figure)
  cites <- as.numeric(manifest$first_cited[k])
  lab <- if (grepl("S", set[1])) "supplementary" else "main"
  if (is.unsorted(cites)) {
    cat("\n  CITATION ORDER WRONG in the ", lab, " set: ",
        paste(set, "cited", manifest$first_cited[k], collapse = "; "), "\n", sep = "")
  } else {
    cat("\n  ", lab, " figures are numbered in order of first citation\n", sep = "")
  }
}

cat("\nIf a figure is renumbered, change its ggsave path, this manifest and\n")
cat("the manuscript citations together.\n")
