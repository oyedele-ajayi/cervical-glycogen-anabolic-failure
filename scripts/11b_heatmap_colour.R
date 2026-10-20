# Phase 11b - Figure 2 heatmap, colour scale fixed
#
# The first version clipped at +/-3. Within-gene z-scores across 128 samples
# sit almost entirely inside +/-1, so most of the ramp went unused and every
# cell rendered pale. Here the limit is taken from the data (98th percentile
# of |z|), the ramp has 255 steps across the full eleven-stop RdBu, and the
# clipping is stated in the legend.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}
set.seed(42)

library(pheatmap)

gene_set <- readRDS("data/processed/gene_set_locked.rds")
gene_mats <- readRDS("data/processed/gene_matrices.rds")
groups <- readRDS("data/processed/groups.rds")
gdi_list <- readRDS("data/processed/GDI_results.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")

m <- gene_mats$GSE63514
grp <- groups$GSE63514
gdi <- gdi_list$GSE63514

# columns: by grade, then ascending GDI_diff within grade
ord <- order(as.integer(grp), gdi$GDI_diff[match(colnames(m), gdi$sample)])
m <- m[, ord]
grp <- grp[ord]

# rows: locked module order
row_genes <- character(0)
row_module <- character(0)
for (mod in module_order) {
  present <- gene_set[[mod]][gene_set[[mod]] %in% rownames(m)]
  row_genes <- c(row_genes, present)
  row_module <- c(row_module, rep(mod, length(present)))
}

z <- t(scale(t(m[row_genes, ])))


# ---------------------------------------------------------------------------
# Where does the data actually sit?
# ---------------------------------------------------------------------------

cat("=== z-score distribution, 33 genes x 128 samples ===\n")
qs <- quantile(abs(z), c(0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 1), na.rm = TRUE)
print(round(qs, 3))
cat("\nproportion of cells inside |z| < 1 :",
    round(mean(abs(z) < 1, na.rm = TRUE), 3), "\n")
cat("proportion inside |z| < 2 :",
    round(mean(abs(z) < 2, na.rm = TRUE), 3), "\n")
cat("proportion inside |z| < 3 :",
    round(mean(abs(z) < 3, na.rm = TRUE), 3), "\n\n")

LIM <- as.numeric(round(quantile(abs(z), 0.98, na.rm = TRUE), 1))
cat("Colour limit set to +/-", LIM, " (98th percentile of |z|)\n", sep = "")
cat("Cells beyond the limit are clipped and shown at the extreme colour.\n")
cat("Clipped cells:", sum(abs(z) > LIM, na.rm = TRUE), "of", length(z),
    sprintf("(%.1f%%)\n\n", 100 * mean(abs(z) > LIM, na.rm = TRUE)))

z[z > LIM] <- LIM
z[z < -LIM] <- -LIM


# ---------------------------------------------------------------------------
# Draw
# ---------------------------------------------------------------------------

ann_col <- data.frame(Grade = grp)
rownames(ann_col) <- colnames(z)

ann_row <- data.frame(Module = factor(row_module, levels = module_order))
rownames(ann_row) <- row_genes

grade_cols <- c(Normal = "#DEEBF7", CIN1 = "#9ECAE1", CIN2 = "#4292C6",
                CIN3 = "#2171B5", Cancer = "#08306B")
module_cols <- c(Synthesis = "#0072B2", Degradation = "#D55E00",
                 Warburg = "#009E73", Differentiation = "#CC79A7")

# full eleven-stop diverging ramp, neutral at the midpoint
ramp <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                           "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582",
                           "#D6604D", "#B2182B", "#67001F"))(255)

draw_it <- function() {
  pheatmap(z,
           cluster_rows = FALSE, cluster_cols = FALSE,
           show_colnames = FALSE, fontsize_row = 7.5,
           annotation_col = ann_col, annotation_row = ann_row,
           annotation_colors = list(Grade = grade_cols, Module = module_cols),
           annotation_names_col = FALSE, annotation_names_row = FALSE,
           gaps_row = cumsum(table(factor(row_module, levels = module_order)))[1:3],
           gaps_col = cumsum(table(grp))[1:4],
           color = ramp,
           breaks = seq(-LIM, LIM, length.out = 256),
           border_color = NA,
           legend_breaks = c(-LIM, -LIM/2, 0, LIM/2, LIM),
           main = "")
}

png("results/figures/gdi/Figure2_heatmap_v3.png",
    width = 3000, height = 1800, res = 300)
draw_it()
dev.off()

pdf("results/figures/gdi/Figure2_heatmap_v3.pdf", width = 10, height = 6)
draw_it()
dev.off()

cat("wrote Figure2_heatmap_v3.png and .pdf\n")

counts <- table(factor(row_module, levels = module_order))
cat("\nModule counts for the legend: ",
    paste(names(counts), counts, sep = " ", collapse = " / "), "\n", sep = "")

legend_txt <- c(
  "Figure 2. Expression of the 33 panel genes across disease grade, discovery",
  "cohort (GSE63514, n = 128). Columns are samples ordered by grade and, within",
  "each grade, by ascending GDI_diff. Rows follow the locked module order:",
  paste0("Synthesis (", counts[["Synthesis"]], "), Degradation (", counts[["Degradation"]],
         "), Warburg (", counts[["Warburg"]], "), Differentiation (",
         counts[["Differentiation"]], ")."),
  paste0("Values are z-scores computed within each gene across all samples, with the",
         " colour scale clipped at +/-", LIM, ", the 98th percentile of |z|;"),
  paste0(sprintf("%.1f%%", 100 * mean(abs(t(scale(t(gene_mats$GSE63514[row_genes, ])))) > LIM, na.rm = TRUE)),
         " of cells exceed this limit and are drawn at the extreme colour."))

writeLines(legend_txt, "results/figures/Figure2_legend.txt")
cat("legend written to results/figures/Figure2_legend.txt\n")
