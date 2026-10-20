# Phase 16 - supplementary figures S3 and S4  (revision 2)
#
# Robustness (manuscript Figure S4) panel B rebuilt: geom_dotplot sized its dots to the bin width, which on a
#    0.019-wide axis overflowed the panel. Now a jittered point strip.
# fgsea (manuscript Figure S3) rebuilt: the running enrichment is computed here rather than taken from
#    plotEnrichment, whose internal layer order is not stable across versions
#    and whose colours could not be addressed reliably.

rm(list = ls())
# Run from the repository root. Opening the .Rproj file in RStudio puts you
# there. The check stops early rather than writing results into whatever
# folder happens to be current.
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("Run this from the repository root, the folder holding scripts/ and data/.")
}

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(42)

library(ggplot2)
library(patchwork)
library(ggrepel)

check_writable <- function(path) {
  if (file.exists(path) && file.access(path, 2) != 0)
    stop("Cannot write ", path, " - close whatever has it open.")
}

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
module_cols <- c(Synthesis = "#0072B2", Degradation = "#D55E00",
                 Warburg = "#009E73", Differentiation = "#CC79A7")

base_theme <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 7.5, colour = "grey45"),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        axis.text = element_text(colour = "grey20"),
        axis.title = element_text(colour = "grey20"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = "grey20"),
        plot.margin = margin(5, 8, 5, 5))


# ===========================================================================
# Manuscript Figure S4 - robustness, first cited in Results 3.7
# ===========================================================================

sens <- read.csv("results/tables/outlier_sensitivity.csv", stringsAsFactors = FALSE)
loo <- read.csv("results/tables/leave_one_gene_out.csv", stringsAsFactors = FALSE)

lab <- c(full = "Full cohort (n=128)",
         criteria4 = "Pre-specified criteria (n=124)",
         legacy5 = "Earlier version's five (n=123)",
         posthoc1 = "Post hoc, GSM1551331 (n=127)")

sA <- rbind(
  data.frame(version = sens$version, threshold = "CIN2+",
             AUC = sens$CIN2plus, lo = sens$CIN2_lo, hi = sens$CIN2_hi),
  data.frame(version = sens$version, threshold = "CIN3+",
             AUC = sens$CIN3plus, lo = sens$CIN3_lo, hi = sens$CIN3_hi),
  data.frame(version = sens$version, threshold = "Cancer",
             AUC = sens$Cancer, lo = sens$Cancer_lo, hi = sens$Cancer_hi))

sA$version <- factor(lab[sA$version], levels = rev(lab))
sA$threshold <- factor(sA$threshold, levels = c("CIN2+", "CIN3+", "Cancer"))

pS3a <- ggplot(sA, aes(x = AUC, y = version, colour = threshold)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 linewidth = 0.45, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c("CIN2+" = "#0072B2", "CIN3+" = "#D55E00",
                                 "Cancer" = "#CC79A7")) +
  scale_x_continuous(limits = c(0.6, 1.0), breaks = seq(0.6, 1.0, 0.1)) +
  labs(title = "A   Discrimination under each exclusion set, discovery cohort",
       subtitle = "the four sets are defined by different procedures and are never pooled",
       x = "AUC (95% bootstrap CI)", y = NULL) +
  base_theme +
  theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey93"),
        legend.position = "bottom")


# --- B: leave one gene out, as a jittered strip ----------------------------

# read the baseline rather than typing it, so it cannot go stale when the
# discovery estimate is recomputed
roc_all <- read.csv("results/tables/GDI_ROC.csv", stringsAsFactors = FALSE)
baseline <- roc_all$AUC[roc_all$GDI_form == "GDI_diff" &
                          roc_all$cohort == "GSE63514" &
                          roc_all$threshold == "CIN3+"]
stopifnot(length(baseline) == 1)
loo$module <- factor(loo$module, levels = module_order)

set.seed(42)
loo$yj <- runif(nrow(loo), -0.33, 0.33)

# label the five largest drops and the largest gain
loo <- loo[order(loo$delta), ]
loo$lab <- ""
loo$lab[1:5] <- loo$gene[1:5]
loo$lab[nrow(loo)] <- loo$gene[nrow(loo)]

pS3b <- ggplot(loo, aes(x = AUC, y = yj, colour = module)) +
  geom_vline(xintercept = baseline, linewidth = 0.45, colour = "grey25") +
  annotate("text", x = baseline, y = 0.52,
           label = sprintf("baseline %.4f", baseline),
           hjust = -0.06, size = 2.4, colour = "grey25") +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_text_repel(aes(label = lab), size = 2.3, colour = "grey15",
                  segment.size = 0.2, segment.colour = "grey60",
                  min.segment.length = 0, box.padding = 0.35,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = module_cols) +
  scale_y_continuous(limits = c(-0.62, 0.62), breaks = NULL) +
  labs(title = "B   Leave-one-gene-out, CIN3+ discrimination",
       subtitle = sprintf("each of the 33 genes removed in turn; mean %.4f, SD %.4f, range %.4f to %.4f",
                          mean(loo$AUC), sd(loo$AUC), min(loo$AUC), max(loo$AUC)),
       x = "AUC after removing one gene", y = NULL) +
  base_theme +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey93"),
        legend.position = "bottom")

figS3 <- pS3a / pS3b + plot_layout(heights = c(1, 0.8))

check_writable("results/figures/validation/FigureS4_robustness.png")
ggsave("results/figures/validation/FigureS4_robustness.png", figS3,
       width = 8.2, height = 6.6, dpi = 400)
ggsave("results/figures/validation/FigureS4_robustness.pdf", figS3,
       width = 8.2, height = 6.6)
cat("wrote Figure S4 (robustness)\n")


# ===========================================================================
# Manuscript Figure S3 - running enrichment, first cited in Results 3.6
# ===========================================================================
# Standard weighted GSEA running score with p = 1, which is fgsea's default.
# Computing it here removes any dependence on plotEnrichment's layer layout.

fgres <- readRDS("data/processed/fgsea_results.rds")
metric <- fgres$metric                 # already sorted decreasing
res <- fgres$fgsea
cohort_genesets <- readRDS("data/processed/cohort_genesets.rds")
gs <- cohort_genesets$GSE63514[module_order]

running_es <- function(stats, set) {
  genes <- names(stats)
  hit <- genes %in% set
  nr <- sum(abs(stats[hit]))
  n_miss <- length(genes) - sum(hit)
  
  step <- ifelse(hit, abs(stats) / nr, -1 / n_miss)
  es <- cumsum(step)
  data.frame(rank = seq_along(es), es = es, hit = hit)
}

curves <- list()
for (mod in module_order) {
  d <- running_es(metric, gs[[mod]])
  d$module <- mod
  curves[[mod]] <- d
}

ylo <- min(sapply(curves, function(d) min(d$es)))
yhi <- max(sapply(curves, function(d) max(d$es)))
pad <- 0.08 * (yhi - ylo)
ylims <- c(ylo - pad * 2.2, yhi + pad)

panels <- list()
for (i in seq_along(module_order)) {
  mod <- module_order[i]
  d <- curves[[mod]]
  r <- res[as.character(res$pathway) == mod, ]
  hits <- d[d$hit, ]
  tick_y <- ylims[1] + 0.06 * diff(ylims)
  
  panels[[mod]] <- ggplot(d, aes(x = rank, y = es)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey75") +
    geom_line(linewidth = 0.6, colour = module_cols[[mod]]) +
    geom_segment(data = hits,
                 aes(x = rank, xend = rank,
                     y = tick_y - 0.035 * diff(ylims),
                     yend = tick_y + 0.035 * diff(ylims)),
                 linewidth = 0.3, colour = "grey25", inherit.aes = FALSE) +
    scale_x_continuous(breaks = c(1, 5000, 10000, 15000, 20000),
                       labels = c("1", "5,000", "10,000", "15,000", "20,000")) +
    coord_cartesian(ylim = ylims) +
    labs(title = sprintf("%s   %s", LETTERS[i], mod),
         subtitle = sprintf("NES %+.2f, FDR %s, %d genes", r$NES,
                            ifelse(r$padj < 0.001,
                                   formatC(r$padj, format = "e", digits = 1),
                                   sprintf("%.3f", r$padj)),
                            r$size),
         x = "rank in the ranked list", y = "running enrichment score") +
    base_theme
}

figS4 <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]])

check_writable("results/figures/gdi/FigureS3_fgsea.png")
ggsave("results/figures/gdi/FigureS3_fgsea.png", figS4,
       width = 8.4, height = 6.0, dpi = 400)
ggsave("results/figures/gdi/FigureS3_fgsea.pdf", figS4,
       width = 8.4, height = 6.0)
cat("wrote Figure S3 (fgsea)\n")


# --- check the computed curves against the reported scores -----------------

cat("\n=== S4 check: peak running score against reported ES ===\n")
for (mod in module_order) {
  d <- curves[[mod]]
  peak <- d$es[which.max(abs(d$es))]
  cat(sprintf("  %-16s computed peak %+.4f   reported NES %+.3f\n",
              mod, peak, res$NES[as.character(res$pathway) == mod]))
}
cat("\n(peak is the raw enrichment score, NES is normalised, so they differ\n")
cat(" in magnitude but must agree in sign)\n")

cat("\nShared y limits:", sprintf("%.3f to %.3f", ylims[1], ylims[2]), "\n")
cat("\nFigure set complete: 1, 2, 3, 4, S1, S2, S3, S4.\n")
