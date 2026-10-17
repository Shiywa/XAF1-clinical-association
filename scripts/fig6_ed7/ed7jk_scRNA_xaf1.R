#!/usr/bin/env Rscript

# Extended Data Fig. 7j-k: XAF1 expression in GSE271284 PBMC single-cell data.
# Cell-level distributions are shown in the figure; donor-level means are used
# for the accompanying statistical table so that cells from one donor are not
# treated as independent biological replicates.

suppressPackageStartupMessages({
    library(Seurat)
    library(harmony)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
    library(patchwork)
    library(viridis)
})

project_root <- Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd())
source(file.path(project_root, "scripts", "common", "config.R"))

python_path <- Sys.getenv("XAF1_RETICULATE_PYTHON", unset = "")
conda_binary <- Sys.getenv("XAF1_CONDA_BINARY", unset = "")
if (nzchar(python_path)) options(RETICULATE_PYTHON = python_path)
if (nzchar(conda_binary)) options(reticulate.conda_binary = conda_binary)

input_rds <- xaf1_env_path(
    "XAF1_GSE271284_RDS",
    xaf1_path("data", "public", "GSE271284", "230217_recast_after_qc.rds"),
    must_exist = TRUE
)
output_dir <- xaf1_output_dir("ed7", "GSE271284")

sce <- readRDS(input_rds)
required_metadata <- c("Patient.Id", "condition_unified")
missing_metadata <- setdiff(required_metadata, colnames(sce@meta.data))
if (length(missing_metadata)) {
    stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}
if (!"XAF1" %in% rownames(sce)) stop("XAF1 was not found in the expression matrix.")

DefaultAssay(sce) <- "RNA"
sce <- NormalizeData(sce, verbose = FALSE)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sce <- ScaleData(sce, features = VariableFeatures(sce), verbose = FALSE)
sce <- RunPCA(sce, features = VariableFeatures(sce), verbose = FALSE)
sce <- RunHarmony(
    sce,
    group.by.vars = "Patient.Id",
    sigma = 0.4,
    theta = 2,
    lambda = 1,
    plot_convergence = FALSE,
    verbose = FALSE
)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30, seed.use = 1234, verbose = FALSE)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:30, verbose = FALSE)
sce <- FindClusters(sce, resolution = 0.6, random.seed = 1234, verbose = FALSE)

cluster_id <- as.character(sce$seurat_clusters)
sce$cluster_anno <- case_when(
    cluster_id %in% c("0", "3", "10", "18", "21", "30", "31", "32", "36") ~ "B cell",
    cluster_id %in% c("2", "6", "13", "17", "24", "28") ~ "CD4+ T",
    cluster_id %in% c("4", "5", "7", "12", "16", "22", "25", "27") ~ "CD8+ T",
    cluster_id %in% c("1", "8", "9", "15", "19", "20", "26", "34") ~ "Monocyte/macrophage",
    TRUE ~ "Other"
)

condition_levels <- c(
    "healthy child", "asymptomatic child", "symptomatic child",
    "healthy adult", "mild adult", "severe adult"
)
celltype_levels <- c("CD4+ T", "CD8+ T", "B cell", "Monocyte/macrophage")
sce$condition_unified <- factor(sce$condition_unified, levels = condition_levels)
sce$cluster_anno <- factor(sce$cluster_anno, levels = c(celltype_levels, "Other"))

keep_cells <- colnames(sce)[sce$cluster_anno %in% celltype_levels & !is.na(sce$condition_unified)]
sce_use <- subset(sce, cells = keep_cells)

embedding <- Embeddings(sce_use, "umap")
xaf1_expression <- FetchData(sce_use, vars = "XAF1")[[1]]
plot_data <- tibble(
    cell = colnames(sce_use),
    UMAP_1 = embedding[, 1],
    UMAP_2 = embedding[, 2],
    XAF1 = xaf1_expression,
    cell_type = factor(sce_use$cluster_anno, levels = celltype_levels),
    condition = factor(sce_use$condition_unified, levels = condition_levels),
    donor = as.character(sce_use$Patient.Id)
)
write_tsv(plot_data, file.path(output_dir, "ED7jk_cell_level_plot_data.tsv.gz"))

# ED7j: density background plus normalized XAF1 expression.
label_positions <- plot_data |>
    group_by(cell_type) |>
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

p_umap <- ggplot(plot_data, aes(UMAP_1, UMAP_2)) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE, n = 250) +
    scale_fill_viridis(option = "magma", guide = "none") +
    geom_point(
        data = filter(plot_data, XAF1 > 0),
        aes(colour = XAF1),
        size = 0.08,
        alpha = 0.55
    ) +
    scale_colour_gradient(low = "#D9F0D3", high = "#008837", name = "XAF1") +
    geom_label(
        data = label_positions,
        aes(label = cell_type),
        size = 2.5,
        colour = "white",
        fill = "#222222CC",
        label.size = 0
    ) +
    coord_equal(expand = FALSE) +
    theme_void(base_size = 9) +
    theme(legend.position = "right")

ggsave(
    file.path(output_dir, "ED7j_XAF1_UMAP.pdf"),
    p_umap,
    width = 6.2,
    height = 5.2,
    device = cairo_pdf
)

# ED7k: cell-level display.
p_box <- ggplot(plot_data, aes(condition, XAF1, fill = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.68, linewidth = 0.35) +
    facet_wrap(~cell_type, ncol = 2, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(x = NULL, y = "Normalized XAF1 expression") +
    theme_classic(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 40, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
    )

ggsave(
    file.path(output_dir, "ED7k_XAF1_by_celltype_and_severity.pdf"),
    p_box,
    width = 8.2,
    height = 6.2,
    device = cairo_pdf
)

# Donor-level means and prespecified within-age comparisons.
donor_summary <- plot_data |>
    group_by(donor, condition, cell_type) |>
    summarise(mean_XAF1 = mean(XAF1), n_cells = n(), .groups = "drop")
write_tsv(donor_summary, file.path(output_dir, "ED7k_donor_level_XAF1.tsv"))

comparisons <- tibble::tribble(
    ~group_1, ~group_2,
    "healthy child", "asymptomatic child",
    "asymptomatic child", "symptomatic child",
    "healthy adult", "mild adult",
    "mild adult", "severe adult"
)

test_donor_groups <- function(cell_type_value, group_1_value, group_2_value) {
    x <- donor_summary$mean_XAF1[
        donor_summary$cell_type == cell_type_value & donor_summary$condition == group_1_value
    ]
    y <- donor_summary$mean_XAF1[
        donor_summary$cell_type == cell_type_value & donor_summary$condition == group_2_value
    ]
    tibble(
        cell_type = cell_type_value,
        group_1 = group_1_value,
        group_2 = group_2_value,
        n_1 = length(x),
        n_2 = length(y),
        p_value = if (length(x) && length(y)) wilcox.test(x, y, exact = FALSE)$p.value else NA_real_
    )
}

donor_tests <- bind_rows(lapply(celltype_levels, function(cell_type_value) {
    bind_rows(Map(
        test_donor_groups,
        cell_type_value = cell_type_value,
        group_1_value = comparisons$group_1,
        group_2_value = comparisons$group_2
    ))
})) |>
    mutate(p_adjusted_BH = p.adjust(p_value, method = "BH"))
write_tsv(donor_tests, file.path(output_dir, "ED7k_donor_level_Wilcoxon_tests.tsv"))

message("Created Extended Data Fig. 7j-k outputs in ", output_dir)
