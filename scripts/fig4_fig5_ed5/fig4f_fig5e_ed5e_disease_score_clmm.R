#!/usr/bin/env Rscript

# Fig. 4f, Fig. 5e and Extended Data Fig. 5e
# Longitudinal ordinal disease-score analysis using cumulative-link mixed models.

suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
    library(ordinal)
})

project_root <- Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd())
source(file.path(project_root, "scripts", "common", "config.R"))

input_file <- xaf1_env_path(
    "XAF1_DISEASE_SCORE_FILE",
    xaf1_path("data", "mouse", "disease_score", "disease_score.xlsx"),
    must_exist = TRUE
)
output_dir <- xaf1_output_dir("fig4_fig5_ed5", "disease_score")

# Workbook order follows the original analysis file.
analysis_spec <- tibble::tribble(
    ~sheet, ~panel,   ~model_day_min, ~output_stem,
    1L,     "Fig4f", 3,              "Fig4f_EV71_disease_score",
    2L,     "Fig5e", -Inf,           "Fig5e_HSV1_brain_disease_score",
    3L,     "ED5e",  -Inf,           "ED5e_HSV1_systemic_disease_score"
)

analyse_sheet <- function(sheet, panel, model_day_min, output_stem) {
    raw <- read_excel(input_file, sheet = sheet)
    if (ncol(raw) < 3L) {
        stop("Sheet ", sheet, " must contain a group column and at least two day columns.")
    }
    names(raw)[1] <- "Group_label"

    long <- raw |>
        mutate(MouseID = paste0(panel, "_M", row_number())) |>
        pivot_longer(
            cols = -c(Group_label, MouseID),
            names_to = "Day",
            values_to = "Score"
        ) |>
        mutate(
            Day = parse_number(as.character(Day)),
            Score = as.numeric(Score),
            Group = if_else(grepl("^WT", Group_label, ignore.case = TRUE), "WT", "Xaf1_KO"),
            Group = factor(Group, levels = c("WT", "Xaf1_KO")),
            MouseID = factor(MouseID)
        )

    model_data <- long |>
        filter(!is.na(Score), !is.na(Day), Day >= model_day_min) |>
        mutate(Score_ordered = ordered(Score, levels = sort(unique(Score))))

    if (nlevels(model_data$Score_ordered) < 2L) {
        stop("Sheet ", sheet, " has fewer than two observed disease-score levels.")
    }

    model <- clmm(
        Score_ordered ~ Group * Day + (1 | MouseID),
        data = model_data,
        link = "logit",
        Hess = TRUE
    )

    coefficient_table <- as.data.frame(coef(summary(model))) |>
        tibble::rownames_to_column("term")
    write_tsv(
        coefficient_table,
        file.path(output_dir, paste0(output_stem, "_CLMM_coefficients.tsv"))
    )

    summary_table <- long |>
        filter(!is.na(Score), !is.na(Day)) |>
        group_by(Group, Group_label, Day) |>
        summarise(
            mean_score = mean(Score),
            sd = sd(Score),
            n_observed = n(),
            se = sd / sqrt(n_observed),
            ci_lower = mean_score - 1.96 * se,
            ci_upper = mean_score + 1.96 * se,
            .groups = "drop"
        )
    write_tsv(
        summary_table,
        file.path(output_dir, paste0(output_stem, "_daily_summary.tsv"))
    )

    mouse_counts <- raw |>
        count(Group_label, name = "n_mice") |>
        mutate(label = paste0(Group_label, " (n=", n_mice, ")"))
    label_lookup <- setNames(mouse_counts$label, mouse_counts$Group_label)

    plot_data <- summary_table |>
        mutate(legend_label = factor(label_lookup[Group_label], levels = label_lookup))

    p <- ggplot(plot_data, aes(x = Day, y = mean_score, group = Group, colour = Group, fill = Group)) +
        geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, colour = NA) +
        geom_line(linewidth = 0.85) +
        geom_point(size = 1.8) +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 0.35) +
        scale_colour_manual(values = c(WT = "#222222", Xaf1_KO = "#E64B35")) +
        scale_fill_manual(values = c(WT = "#222222", Xaf1_KO = "#E64B35")) +
        labs(x = "Day post infection", y = "Mean disease score", colour = NULL, fill = NULL) +
        theme_classic(base_size = 11) +
        theme(legend.position = "top", axis.text = element_text(colour = "black"))

    output_pdf <- file.path(output_dir, paste0(output_stem, ".pdf"))
    ggsave(output_pdf, p, width = 5.2, height = 4.2, device = cairo_pdf)

    saveRDS(model, file.path(output_dir, paste0(output_stem, "_CLMM.rds")))
    message(panel, ": created ", output_pdf)
    invisible(model)
}

models <- Map(
    analyse_sheet,
    analysis_spec$sheet,
    analysis_spec$panel,
    analysis_spec$model_day_min,
    analysis_spec$output_stem
)
names(models) <- analysis_spec$panel
