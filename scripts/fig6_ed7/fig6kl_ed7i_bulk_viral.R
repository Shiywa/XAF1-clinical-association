#!/usr/bin/env Rscript

# Fig. 6k-l and Extended Data Fig. 7i
# XAF1 expression and pathway associations in four human viral-infection cohorts.

suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(GSVA)
    library(BiocParallel)
    library(SCP)
})

project_root <- Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd())
source(file.path(project_root, "scripts", "common", "config.R"))

data_root <- xaf1_env_path(
    "XAF1_VIRAL_DATA_DIR",
    xaf1_path("data", "public", "viral_bulk")
)
output_dir <- xaf1_output_dir("fig6_ed7", "viral_bulk")
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)

normalise_microarray <- function(matrix) {
    qx <- quantile(matrix, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
    log_transform <- (qx[[5]] > 100) || (qx[[6]] - qx[[1]] > 50 && qx[[2]] > 0)
    if (log_transform) {
        matrix[matrix <= 0] <- NA
        matrix <- log2(matrix + 1)
    }
    normalizeBetweenArrays(matrix)
}

collapse_probes <- function(matrix, mapping) {
    mapping <- mapping |>
        transmute(probe_id = as.character(probe_id), symbol = as.character(symbol)) |>
        filter(!is.na(symbol), nzchar(symbol)) |>
        distinct(probe_id, .keep_all = TRUE)

    joined <- as.data.frame(matrix) |>
        rownames_to_column("probe_id") |>
        inner_join(mapping, by = "probe_id")
    sample_columns <- setdiff(colnames(joined), c("probe_id", "symbol"))
    joined$mean_expression <- rowMeans(joined[, sample_columns, drop = FALSE], na.rm = TRUE)
    selected <- joined |>
        group_by(symbol) |>
        slice_max(mean_expression, n = 1, with_ties = FALSE) |>
        ungroup()
    result <- as.matrix(selected[, sample_columns, drop = FALSE])
    rownames(result) <- selected$symbol
    storage.mode(result) <- "numeric"
    result
}

collapse_gene_ids <- function(matrix, annotation) {
    annotation <- annotation |>
        transmute(GeneID = as.character(GeneID), symbol = as.character(Symbol)) |>
        filter(!is.na(symbol), nzchar(symbol)) |>
        distinct(GeneID, .keep_all = TRUE)
    joined <- as.data.frame(matrix) |>
        rownames_to_column("GeneID") |>
        inner_join(annotation, by = "GeneID")
    sample_columns <- setdiff(colnames(joined), c("GeneID", "symbol"))
    joined$mean_expression <- rowMeans(joined[, sample_columns, drop = FALSE], na.rm = TRUE)
    selected <- joined |>
        group_by(symbol) |>
        slice_max(mean_expression, n = 1, with_ties = FALSE) |>
        ungroup()
    result <- as.matrix(selected[, sample_columns, drop = FALSE])
    rownames(result) <- selected$symbol
    storage.mode(result) <- "numeric"
    result
}

read_geo <- function(accession) {
    accession_dir <- file.path(data_root, accession)
    dir.create(accession_dir, recursive = TRUE, showWarnings = FALSE)
    getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE, destdir = accession_dir)[[1]]
}

load_gse101702 <- function() {
    gset <- read_geo("GSE101702")
    matrix <- normalise_microarray(exprs(gset))
    platform <- read_delim(
        xaf1_require_file(file.path(data_root, "GSE101702", "GPL21185-21174.txt"), "GSE101702 platform annotation"),
        delim = "\t", skip = 18, show_col_types = FALSE
    )
    expression <- collapse_probes(matrix, setNames(platform[, c(1, 6)], c("probe_id", "symbol")))
    metadata <- pData(gset)
    metadata <- metadata[match(colnames(expression), metadata$geo_accession), , drop = FALSE]
    group <- recode(metadata$`severity:ch1`, hlty_ctrl = "Healthy", flu_mod = "Moderate", flu_svre = "Severe")
    plot_data <- tibble(sample = metadata$geo_accession, XAF1 = expression["XAF1", ], group = group) |>
        filter(!is.na(group)) |>
        mutate(group = factor(group, levels = c("Healthy", "Moderate", "Severe")))
    list(accession = "GSE101702", virus = "IAV", expression = expression[, plot_data$sample, drop = FALSE], plot_data = plot_data)
}

load_gse111368 <- function() {
    gset <- read_geo("GSE111368")
    matrix <- normalise_microarray(exprs(gset))
    ids <- idmap2::getIDs("gpl10558") |>
        transmute(probe_id = as.character(probe_id), symbol = as.character(symbol))
    expression <- collapse_probes(matrix, ids)
    metadata <- pData(gset)
    metadata <- metadata[match(colnames(expression), metadata$geo_accession), , drop = FALSE]
    group <- recode(as.character(metadata$`t1severity:ch1`), HC = "Healthy", `1` = "Mild", `2` = "Mild", `3` = "Severe")
    keep <- metadata$`timepoint:ch1` %in% c("HC", "T1") & !is.na(group)
    plot_data <- tibble(sample = metadata$geo_accession[keep], XAF1 = expression["XAF1", metadata$geo_accession[keep]], group = group[keep]) |>
        mutate(group = factor(group, levels = c("Healthy", "Mild", "Severe")))
    list(accession = "GSE111368", virus = "IAV", expression = expression[, plot_data$sample, drop = FALSE], plot_data = plot_data)
}

load_gse51808 <- function() {
    gset <- read_geo("GSE51808")
    matrix <- normalise_microarray(exprs(gset))
    platform <- read_delim(
        xaf1_require_file(file.path(data_root, "GSE51808", "GPL13158.annot.gz"), "GSE51808 platform annotation"),
        delim = "\t", skip = 27, show_col_types = FALSE
    )
    expression <- collapse_probes(matrix, setNames(platform[, c(1, 3)], c("probe_id", "symbol")))
    metadata <- pData(gset)
    metadata <- metadata[match(colnames(expression), metadata$geo_accession), , drop = FALSE]
    group <- recode(metadata$`status:ch1`, control = "Healthy", convalescent = "Convalescent", DF = "DF", DHF = "DHF")
    keep <- group %in% c("Healthy", "DF", "DHF")
    plot_data <- tibble(sample = metadata$geo_accession[keep], XAF1 = expression["XAF1", metadata$geo_accession[keep]], group = group[keep]) |>
        mutate(group = factor(group, levels = c("Healthy", "DF", "DHF")))
    list(accession = "GSE51808", virus = "DENV", expression = expression[, plot_data$sample, drop = FALSE], plot_data = plot_data)
}

load_gse94892 <- function() {
    gset <- read_geo("GSE94892")
    metadata <- pData(gset)
    counts <- read_tsv(
        xaf1_require_file(file.path(data_root, "GSE94892", "GSE94892_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"), "GSE94892 TPM matrix"),
        show_col_types = FALSE
    ) |>
        as.data.frame()
    rownames(counts) <- as.character(counts$GeneID)
    counts$GeneID <- NULL
    common_samples <- intersect(colnames(counts), metadata$geo_accession)
    counts <- as.matrix(counts[, common_samples, drop = FALSE])
    metadata <- metadata[match(common_samples, metadata$geo_accession), , drop = FALSE]
    annotation <- read_tsv(
        xaf1_require_file(file.path(data_root, "annotation", "Human.GRCh38.p13.annot.tsv.gz"), "NCBI human annotation"),
        show_col_types = FALSE
    )
    expression <- collapse_gene_ids(counts, annotation)
    group <- recode(
        metadata$`disease:ch1`,
        `Other Febrile Illness` = "OFI",
        `Dengue fever` = "DF",
        `Dengue hemorrhagic fever` = "DHF",
        `Dengue shock syndrome` = "DSS"
    )
    keep <- group %in% c("OFI", "DF", "DHF")
    plot_data <- tibble(sample = metadata$geo_accession[keep], XAF1 = expression["XAF1", metadata$geo_accession[keep]], group = group[keep]) |>
        mutate(group = factor(group, levels = c("OFI", "DF", "DHF")))
    list(accession = "GSE94892", virus = "DENV", expression = expression[, plot_data$sample, drop = FALSE], plot_data = plot_data)
}

cohorts <- list(load_gse101702(), load_gse111368(), load_gse51808(), load_gse94892())
names(cohorts) <- vapply(cohorts, `[[`, character(1), "accession")

expression_plot_data <- bind_rows(lapply(cohorts, function(x) {
    x$plot_data |>
        mutate(accession = x$accession, virus = x$virus)
}))
write_tsv(expression_plot_data, file.path(output_dir, "Fig6kl_XAF1_expression_plot_data.tsv"))

p_expression <- ggplot(expression_plot_data, aes(group, XAF1, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.65) +
    facet_wrap(~accession, scales = "free", nrow = 1) +
    labs(x = NULL, y = "XAF1 normalized expression") +
    theme_classic(base_size = 9) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
    )
ggsave(file.path(output_dir, "Fig6kl_XAF1_viral_severity.pdf"), p_expression, width = 10, height = 3.4, device = cairo_pdf)

comparison_map <- list(
    GSE101702 = list(c("Healthy", "Moderate"), c("Moderate", "Severe")),
    GSE111368 = list(c("Healthy", "Mild"), c("Mild", "Severe")),
    GSE51808 = list(c("Healthy", "DF"), c("DF", "DHF")),
    GSE94892 = list(c("OFI", "DF"), c("DF", "DHF"))
)
expression_tests <- bind_rows(lapply(names(comparison_map), function(accession) {
    data <- filter(expression_plot_data, .data$accession == accession)
    bind_rows(lapply(comparison_map[[accession]], function(groups) {
        x <- data$XAF1[data$group == groups[[1]]]
        y <- data$XAF1[data$group == groups[[2]]]
        tibble(
            accession = accession,
            group_1 = groups[[1]],
            group_2 = groups[[2]],
            n_1 = length(x),
            n_2 = length(y),
            p_value = wilcox.test(x, y, exact = FALSE)$p.value
        )
    }))
})) |>
    mutate(p_adjusted_BH = p.adjust(p_value, method = "BH"))
write_tsv(expression_tests, file.path(output_dir, "Fig6kl_Wilcoxon_tests.tsv"))

# Extended Data Fig. 7i: ssGSEA followed by within-cohort Spearman correlation.
kegg_db <- SCP::PrepareDB(species = "Homo_sapiens", db = "KEGG")$Homo_sapiens$KEGG
kegg_names <- kegg_db$TERM2NAME
gene_sets <- split(kegg_db$TERM2GENE$symbol, kegg_db$TERM2GENE$Term)
gene_sets <- gene_sets[lengths(gene_sets) >= 10]
workers <- min(xaf1_threads(), 8L)
bpparam <- if (.Platform$OS.type == "windows" || workers == 1L) {
    SerialParam()
} else {
    MulticoreParam(workers = workers)
}

correlation_results <- bind_rows(lapply(cohorts, function(cohort) {
    matrix <- cohort$expression
    parameter <- ssgseaParam(as.matrix(matrix), gene_sets)
    scores <- gsva(parameter, BPPARAM = bpparam, verbose = FALSE)
    xaf1 <- as.numeric(matrix["XAF1", colnames(scores)])
    tests <- lapply(seq_len(nrow(scores)), function(i) {
        score <- as.numeric(scores[i, ])
        keep <- is.finite(xaf1) & is.finite(score)
        if (sum(keep) < 4L) return(c(rho = NA_real_, p_value = NA_real_, n = sum(keep)))
        test <- suppressWarnings(cor.test(xaf1[keep], score[keep], method = "spearman", exact = FALSE))
        c(rho = unname(test$estimate), p_value = test$p.value, n = sum(keep))
    })
    tests <- as.data.frame(do.call(rbind, tests)) |>
        mutate(
            pathway_ID = rownames(scores),
            adjusted_p = p.adjust(p_value, method = "BH"),
            accession = cohort$accession,
            virus = cohort$virus
        ) |>
        left_join(kegg_names, by = c("pathway_ID" = "Term"))
    tests
}))
write_tsv(correlation_results, file.path(output_dir, "ED7i_XAF1_KEGG_correlations.tsv"))

pathway_rank <- correlation_results |>
    filter(!is.na(adjusted_p)) |>
    group_by(Name) |>
    summarise(n_significant = sum(adjusted_p < 0.05), minimum_adjusted_p = min(adjusted_p), .groups = "drop") |>
    filter(n_significant >= 3L) |>
    arrange(minimum_adjusted_p)
top_pathways <- head(pathway_rank$Name, 30L)
if (!length(top_pathways)) stop("No KEGG pathways passed the ED7i selection rule.")

heatmap_data <- correlation_results |>
    filter(Name %in% top_pathways) |>
    mutate(
        accession = factor(accession, levels = names(cohorts)),
        Name = factor(Name, levels = rev(top_pathways))
    )

p_heatmap <- ggplot(heatmap_data, aes(accession, Name)) +
    geom_tile(fill = "#D9D9D9", colour = "white", linewidth = 0.3) +
    geom_tile(
        data = filter(heatmap_data, adjusted_p < 0.05),
        aes(fill = rho),
        colour = "white",
        linewidth = 0.3
    ) +
    geom_text(data = filter(heatmap_data, adjusted_p < 0.05), label = "*", colour = "#C51B7D", size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Spearman rho") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 8.5) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 0, colour = "black"),
        axis.text.y = element_text(colour = "black")
    )
ggsave(file.path(output_dir, "ED7i_XAF1_KEGG_heatmap.pdf"), p_heatmap, width = 7.2, height = 8.5, device = cairo_pdf)

message("Created Fig. 6k-l and Extended Data Fig. 7i outputs in ", output_dir)
