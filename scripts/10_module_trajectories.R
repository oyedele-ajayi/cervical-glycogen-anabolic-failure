# Phase 10 - main figures, reframed version  (revision 2)
#
# Figure 3  module trajectories across grade, four cohorts. The thesis figure.
#           Plotted as change from each cohort's Normal median, on a shared
#           y-scale so the four panels are directly comparable.
#
# Figure 4  ROC curves (A-D) plus single-module discrimination (E). Panel E is
#           a forest plot, not bars: an AUC with a CI is a point estimate with
#           uncertainty, and bars on a truncated axis would misrepresent it.
#
# Okabe-Ito colourblind-safe palette. Grade order comes from saved factor
# levels, never from a read.csv.

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
library(pROC)

gdi_list <- readRDS("data/processed/GDI_results.rds")
roc_objects <- readRDS("data/processed/ROC_objects.rds")

module_order <- c("Synthesis", "Degradation", "Warburg", "Differentiation")

module_cols <- c(Synthesis = "#0072B2", Degradation = "#D55E00",
                 Warburg = "#009E73", Differentiation = "#CC79A7")

cohort_labels <- c(
  GSE63514 = "A   GSE63514 - microdissected epithelium (n=128)",
  GSE138080 = "B   GSE138080 - bulk tissue (n=35)",
  GSE27678_GPL571 = "C   GSE27678 GPL571 - exfoliated cytology (n=44)",
  GSE27678_GPL570 = "D   GSE27678 GPL570 - exfoliated cytology (n=31)")

base_theme <- theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 4)),
    axis.line = element_line(linewidth = 0.3, colour = "grey30"),
    axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
    axis.text = element_text(colour = "grey20"),
    axis.title = element_text(colour = "grey20"),
    legend.title = element_blank(),
    legend.text = element_text(colour = "grey20", size = 8),
    legend.key.width = unit(7, "mm"),
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
    plot.margin = margin(5, 8, 5, 5))


# ---------------------------------------------------------------------------
# Figure 3 - module trajectories
# ---------------------------------------------------------------------------

traj <- NULL
for (nm in names(gdi_list)) {
  s <- gdi_list[[nm]]
  for (mod in module_order) {
    base <- median(s[[mod]][s$group == "Normal"])
    for (lv in levels(s$group)) {
      v <- s[[mod]][s$group == lv]
      traj <- rbind(traj, data.frame(
        cohort = nm, module = mod, grade = lv,
        grade_i = which(levels(s$group) == lv),
        delta = median(v) - base, n = length(v),
        stringsAsFactors = FALSE))
    }
  }
}
traj$module <- factor(traj$module, levels = module_order)
write.csv(traj, "results/tables/Figure3_module_trajectories.csv", row.names = FALSE)

# one shared y-scale, so the four panels can be compared directly
ylo <- floor(min(traj$delta) * 20) / 20
yhi <- ceiling(max(traj$delta) * 20) / 20

panels3 <- list()
for (nm in names(gdi_list)) {
  d <- traj[traj$cohort == nm, ]
  lv <- levels(gdi_list[[nm]]$group)
  last <- d[d$grade_i == max(d$grade_i), ]
  
  p <- ggplot(d, aes(x = grade_i, y = delta, colour = module, group = module)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.9) +
    geom_text(data = last, aes(label = module), hjust = 0, nudge_x = 0.08,
              size = 2.3, show.legend = FALSE) +
    scale_colour_manual(values = module_cols, drop = FALSE) +
    scale_x_continuous(breaks = seq_along(lv), labels = lv,
                       limits = c(0.9, length(lv) + 1.25)) +
    scale_y_continuous(limits = c(ylo, yhi)) +
    labs(title = cohort_labels[[nm]], x = NULL,
         y = "change from Normal median") +
    base_theme
  
  panels3[[nm]] <- p
}

fig3 <- (panels3$GSE63514 | panels3$GSE138080) /
  (panels3$GSE27678_GPL571 | panels3$GSE27678_GPL570) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Manuscript Figure 3, first cited in Results 3.4. Figures are numbered by
# order of first citation, which is what the filename carries.
ggsave("results/figures/gdi/Figure3_module_trajectories.png", fig3,
       width = 9.4, height = 6.6, dpi = 400)
ggsave("results/figures/gdi/Figure3_module_trajectories.pdf", fig3,
       width = 9.4, height = 6.6)
cat("wrote Figure 3 (module trajectories)\n")


# ---------------------------------------------------------------------------
# Figure 4 - ROC panels plus single-module forest plot
# ---------------------------------------------------------------------------

roc_panel <- function(cohort, labels, title) {
  d <- NULL
  series_levels <- character(0)
  for (lab in labels) {
    key <- paste(cohort, lab, "GDI_diff", sep = "|")
    r <- roc_objects[[key]]
    s <- sprintf("%s  (AUC %.3f)", lab, as.numeric(r$auc))
    series_levels <- c(series_levels, s)
    d <- rbind(d, data.frame(fpr = 1 - r$specificities,
                             tpr = r$sensitivities,
                             series = s, stringsAsFactors = FALSE))
  }
  # legend follows severity order, not colour order
  d$series <- factor(d$series, levels = series_levels)
  d <- d[order(d$series, d$fpr), ]
  
  cols <- c("#0072B2", "#D55E00", "#CC79A7")[seq_along(series_levels)]
  names(cols) <- series_levels
  
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
    theme(legend.position = c(0.97, 0.06),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key.height = unit(3.4, "mm"),
          legend.text = element_text(size = 6.8),
          panel.grid.major.y = element_blank())
}

p4a <- roc_panel("GSE63514", c("CIN2+", "CIN3+", "Cancer"), "A   GSE63514 (n=128)")
p4b <- roc_panel("GSE138080", c("HGCIN+", "SCC"), "B   GSE138080 (n=35)")
p4c <- roc_panel("GSE27678_GPL571", c("LSIL+", "HSIL"), "C   GSE27678 GPL571 (n=44)")
p4d <- roc_panel("GSE27678_GPL570", c("Cancer"),
                 "D   GSE27678 GPL570 (n=31; only 3 normals)")


# Panel E as a forest plot. An AUC with a bootstrap CI is a point estimate
# with uncertainty; bars on a truncated axis would misrepresent it.

dropout <- read.csv("results/tables/module_dropout.csv", stringsAsFactors = FALSE)
dropout$predictor <- factor(dropout$predictor, levels = rev(dropout$predictor))
dropout$col <- ifelse(grepl("composite", dropout$predictor), "grey35",
                      module_cols[as.character(dropout$predictor)])

p4e <- ggplot(dropout, aes(x = AUC, y = predictor)) +
  geom_vline(xintercept = 0.5, linewidth = 0.3, linetype = "dashed",
             colour = "grey60") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0,
                 linewidth = 0.5, colour = "grey35") +
  geom_point(size = 2.6, colour = dropout$col) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -1.1,
            size = 2.5, colour = "grey20") +
  scale_x_continuous(limits = c(0.35, 1.0),
                     breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(title = "E   Single-module discrimination of CIN3+ (discovery cohort)",
       x = "AUC (95% bootstrap CI)", y = NULL) +
  base_theme +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey92"))

fig4 <- (p4a | p4b) / (p4c | p4d) / p4e +
  plot_layout(heights = c(1, 1, 0.62))

# Superseded. Script 15 rebuilds this as the four-panel manuscript Figure 4,
# and the single-module forest plot here became Figure 3B in script 14. Kept
# because the five-panel version is what the first draft discussed, but written
# to a superseded folder so it cannot be deposited as Figure 4 by mistake.
dir.create("results/figures/superseded", showWarnings = FALSE, recursive = TRUE)
ggsave("results/figures/superseded/ROC_5panel_superseded.png", fig4,
       width = 8.2, height = 9.4, dpi = 400)
ggsave("results/figures/superseded/ROC_5panel_superseded.pdf", fig4,
       width = 8.2, height = 9.4)
cat("wrote the superseded five-panel ROC figure\n")


# ---------------------------------------------------------------------------
# Legends, written to file
# ---------------------------------------------------------------------------

# Legends are not written here any more. They drifted out of step with the
# manuscript once the figures were renumbered, and two sources of truth is what
# caused that. The manuscript document holds the legends; this script holds only
# the figures and the trajectory table.

cat("\nWritten: Figure3_module_trajectories.png/.pdf,\n")
cat("  superseded/ROC_5panel_superseded.png/.pdf, trajectory table\n")
