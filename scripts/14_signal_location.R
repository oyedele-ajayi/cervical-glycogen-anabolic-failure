# Phase 14 - Figure 2, where the signal lives
#
# A  At which grade does each module separate from normal epithelium?
#    Mann-Whitney against the normal group, discovery cohort, plotted as
#    -log10(p) so the crossing point is readable. Synthesis crosses at CIN2,
#    Differentiation at CIN3, Degradation spikes at CIN2 and returns, Warburg
#    never crosses.
#
# B  Single-module and reduced-index discrimination of CIN3+, with bootstrap
#    intervals computed here rather than carried over, so the index variants
#    and the single modules are on the same footing.

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
library(pROC)
library(ggrepel)

o2 <- read.csv("results/tables/O2_module_vs_normal.csv", stringsAsFactors = FALSE)
gdi_list <- readRDS("data/processed/GDI_results.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")
module_cols <- c(Synthesis = "#0072B2", Degradation = "#D55E00",
                 Warburg = "#009E73", Differentiation = "#CC79A7")
N_BOOT <- 2000

base_theme <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 5)),
        axis.line = element_line(linewidth = 0.3, colour = "grey30"),
        axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
        axis.text = element_text(colour = "grey20"),
        axis.title = element_text(colour = "grey20"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8, colour = "grey20"),
        plot.margin = margin(6, 8, 6, 6))


# ---------------------------------------------------------------------------
# Panel A - grade at which each module separates from normal
# ---------------------------------------------------------------------------

d <- o2[o2$module %in% module_order, ]
d$module <- factor(d$module, levels = module_order)
d$grade <- factor(d$grade, levels = c("CIN1", "CIN2", "CIN3", "Cancer"))
d$grade_i <- as.integer(d$grade)
d$y <- -log10(d$p)

last <- d[d$grade_i == max(d$grade_i), ]

pA <- ggplot(d, aes(x = grade_i, y = y, colour = module, group = module)) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.3,
             linetype = "dashed", colour = "grey55") +
  annotate("text", x = 0.75, y = -log10(0.05) + 0.25, label = "p = 0.05",
           size = 2.3, colour = "grey45", hjust = 0) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.1) +
  ggrepel::geom_text_repel(data = last, aes(label = module), hjust = 0,
                           nudge_x = 0.10, direction = "y", size = 2.4,
                           segment.size = 0.2, segment.colour = "grey65",
                           min.segment.length = 0, box.padding = 0.15,
                           show.legend = FALSE) +
  scale_colour_manual(values = module_cols, guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = levels(d$grade),
                     limits = c(0.7, 5.25)) +
  labs(title = "A   Grade at which each module separates from normal epithelium",
       x = NULL, y = "-log10 p (Mann-Whitney vs Normal)") +
  base_theme

write.csv(d[, c("module", "grade", "median_grade", "median_normal", "p")],
          "results/tables/Figure2A_module_separation.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# Panel B - single modules and reduced index forms
# ---------------------------------------------------------------------------

s <- gdi_list$GSE63514
grp <- s$group
outcome <- ifelse(as.character(grp) %in% c("CIN3", "Cancer"), 1, 0)

score_list <- list(
  "Synthesis alone"            = -s$Synthesis,
  "Degradation alone"          = s$Degradation,
  "Warburg alone"              = s$Warburg,
  "Differentiation alone"      = -s$Differentiation,
  "Synthesis + Differentiation" = -(s$Synthesis + s$Differentiation),
  "Index without Warburg"      = s$Degradation - (s$Synthesis + s$Differentiation),
  "Full four-module index"     = (s$Warburg + s$Degradation) - (s$Synthesis + s$Differentiation)
)

fb <- NULL
for (nm in names(score_list)) {
  set.seed(42)
  r <- roc(outcome, score_list[[nm]], quiet = TRUE, direction = "<")
  set.seed(42)
  ci <- ci.auc(r, method = "bootstrap", boot.n = N_BOOT, boot.stratified = TRUE)
  fb <- rbind(fb, data.frame(predictor = nm,
                             AUC = as.numeric(r$auc),
                             lo = as.numeric(ci[1]),
                             hi = as.numeric(ci[3]),
                             stringsAsFactors = FALSE))
}

fb$kind <- ifelse(grepl("alone", fb$predictor), "single module", "composite")
fb$predictor <- factor(fb$predictor, levels = rev(names(score_list)))

point_cols <- c("Synthesis alone" = module_cols[["Synthesis"]],
                "Degradation alone" = module_cols[["Degradation"]],
                "Warburg alone" = module_cols[["Warburg"]],
                "Differentiation alone" = module_cols[["Differentiation"]],
                "Synthesis + Differentiation" = "grey45",
                "Index without Warburg" = "grey45",
                "Full four-module index" = "grey25")
fb$col <- point_cols[as.character(fb$predictor)]

pB <- ggplot(fb, aes(x = AUC, y = predictor)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "dashed",
             colour = "grey60") +
  geom_vline(xintercept = fb$AUC[fb$predictor == "Full four-module index"],
             linewidth = 0.3, linetype = "dotted", colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 linewidth = 0.5, colour = "grey40") +
  geom_point(size = 2.5, colour = fb$col) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -1.15,
            size = 2.4, colour = "grey20") +
  scale_x_continuous(limits = c(0.35, 1.0),
                     breaks = seq(0.4, 1.0, 0.1)) +
  labs(title = "B   Discrimination of CIN3+ by module and by reduced index (discovery cohort)",
       x = "AUC (95% bootstrap CI)", y = NULL) +
  base_theme +
  theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey93"))

write.csv(fb[, c("predictor", "AUC", "lo", "hi", "kind")],
          "results/tables/Figure2B_module_and_index.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------

fig3 <- pA / pB + plot_layout(heights = c(1, 1.05))

ggsave("results/figures/gdi/Figure2_signal_location.png", fig3,
       width = 8.4, height = 7.0, dpi = 400)
ggsave("results/figures/gdi/Figure2_signal_location.pdf", fig3,
       width = 8.4, height = 7.0)

cat("wrote Figure2_signal_location.png and .pdf\n\n")

cat("=== Panel A: first grade separating from normal at p < 0.05 ===\n")
for (m in module_order) {
  sub <- d[d$module == m, ]
  first <- as.character(sub$grade[sub$p < 0.05])[1]
  cat(sprintf("  %-16s %s\n", m, ifelse(is.na(first), "none", first)))
}

cat("\n=== Panel B ===\n")
print(data.frame(predictor = as.character(fb$predictor),
                 AUC = round(fb$AUC, 4),
                 CI = sprintf("%.3f to %.3f", fb$lo, fb$hi)),
      row.names = FALSE)

cat("\nWritten: Figure2A_module_separation.csv, Figure2B_module_and_index.csv\n")
