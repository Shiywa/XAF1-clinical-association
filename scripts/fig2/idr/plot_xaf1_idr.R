#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: plot_xaf1_idr.R <position_summary.tsv> <output_directory>")
}

input_file <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(input_file, check.names = FALSE)
n_species <- unique(d$species_total)
stopifnot(length(n_species) == 1)

cross_label <- sprintf(
  "%d non-human mammalian orthologues\n(median; 10th–90th percentiles)",
  n_species
)

main_plot <- ggplot(d, aes(x = human_position)) +
  annotate(
    "rect", xmin = 188.5, xmax = 257.5, ymin = -Inf, ymax = Inf,
    fill = "#E64B35", alpha = 0.065
  ) +
  geom_ribbon(
    aes(ymin = cross_species_q10, ymax = cross_species_q90, fill = cross_label),
    alpha = 0.34, colour = NA
  ) +
  geom_line(
    aes(y = cross_species_median, colour = cross_label),
    linewidth = 0.82
  ) +
  geom_line(
    aes(y = human_iupred2_long, colour = "Human XAF1"),
    linewidth = 0.72
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#3C65E5", linewidth = 0.55) +
  annotate(
    "segment", x = 189, xend = 257, y = 0.035, yend = 0.035,
    colour = "#E64B35", linewidth = 1.6, lineend = "butt"
  ) +
  annotate(
    "text", x = 223, y = 0.072, label = "IDR 189–257",
    colour = "#E64B35", size = 3.0
  ) +
  scale_colour_manual(
    values = setNames(c("#E64B35", "#171717"), c("Human XAF1", cross_label)),
    breaks = c("Human XAF1", cross_label)
  ) +
  scale_fill_manual(
    values = setNames("#9B9B9B", cross_label),
    breaks = cross_label
  ) +
  scale_x_continuous(
    breaks = c(1, 50, 100, 150, 200, 250, 301),
    limits = c(1, 301), expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2), limits = c(0, 1),
    expand = expansion(mult = c(0, 0.015))
  ) +
  labs(x = "Human XAF1 residue position", y = "IUPred2 long disorder score") +
  guides(
    colour = guide_legend(order = 1, override.aes = list(linewidth = c(0.8, 1.1))),
    fill = "none"
  ) +
  theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(1.6, "mm"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.title = element_blank(),
    legend.key.width = unit(8, "mm"),
    legend.key.height = unit(4, "mm"),
    legend.spacing.x = unit(2.0, "mm"),
    plot.margin = margin(4, 5, 3, 4, unit = "mm")
  )

coverage_plot <- ggplot(d, aes(x = human_position)) +
  annotate(
    "rect", xmin = 188.5, xmax = 257.5, ymin = -Inf, ymax = Inf,
    fill = "#E64B35", alpha = 0.065
  ) +
  geom_line(aes(y = fraction_score_gt_0.5), colour = "#4B6F91", linewidth = 0.65) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "#777777", linewidth = 0.45) +
  scale_x_continuous(
    breaks = c(1, 50, 100, 150, 200, 250, 301),
    limits = c(1, 301), expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25), limits = c(0, 1),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = "Human XAF1 residue position",
    y = "Orthologues with\nscore > 0.5"
  ) +
  theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = unit(1.6, "mm"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.margin = margin(4, 5, 3, 4, unit = "mm")
  )

prefix <- file.path(output_dir, "XAF1_IDR_maximal_nonhuman_mammals")
ggsave(paste0(prefix, ".pdf"), main_plot, width = 7.1, height = 3.45, units = "in", device = cairo_pdf)
ggsave(paste0(prefix, ".png"), main_plot, width = 7.1, height = 3.45, units = "in", dpi = 600, bg = "white")
svg(paste0(prefix, ".svg"), width = 7.1, height = 3.45, family = "sans")
print(main_plot)
dev.off()

qc_prefix <- file.path(output_dir, "XAF1_IDR_fraction_disordered_QC")
ggsave(paste0(qc_prefix, ".pdf"), coverage_plot, width = 7.1, height = 2.45, units = "in", device = cairo_pdf)
ggsave(paste0(qc_prefix, ".png"), coverage_plot, width = 7.1, height = 2.45, units = "in", dpi = 600, bg = "white")

message(sprintf("Created figures using %d non-human mammalian orthologues", n_species))
