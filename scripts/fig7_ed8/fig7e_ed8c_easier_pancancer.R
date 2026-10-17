#!/usr/bin/env Rscript

# Fig. 7e and Extended Data Fig. 8c: pan-cancer XAF1-stratified
# transcriptomic immune features and EaSIeR-predicted ICB response scores.

library(easier)
library(biomaRt)
library(readr)
library(dplyr)

project_root <- Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd())
source(file.path(project_root, "scripts", "common", "config.R"))
easier_root <- xaf1_env_path("XAF1_EASIER_ROOT", must_exist = TRUE)
gdc_root <- xaf1_env_path("XAF1_GDC_ROOT", must_exist = TRUE)
output_dir <- xaf1_output_dir("fig7_ed8", "easier")
all_easier_dir <- file.path(output_dir, "per_cancer")
dir.create(all_easier_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

pan_cancer_set <- c(
    "BLCA", "BRCA", "CESC", "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD",
    "LUSC", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC", "COAD", "READ", "ESCA"
)
stopifnot(length(pan_cancer_set) == 19L, length(unique(pan_cancer_set)) == 19L)

# 3. 预测 ICB 响应可能性
predict_immune_response2 <- function(pathways = NULL, immunecells = NULL, tfs = NULL, lrpairs = NULL,
                                     ccpairs = NULL, cancer_type, verbose = TRUE) {
    if (missing(cancer_type)) {
        stop("cancer type needs to be specified")
    }
    available_cancer_types <- c(
        "BLCA", "BRCA", "CESC", "CRC",
        "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC",
        "NSCLC", "PAAD", "PRAD", "SKCM", "STAD", "THCA",
        "UCEC"
    )
    if (!cancer_type %in% available_cancer_types) {
        stop("cancer-type specific model not available for this cancer type")
    }
    if (all(
        is.null(pathways), is.null(immunecells), is.null(tfs),
        is.null(lrpairs), is.null(ccpairs)
    )) {
        stop("none signature specified")
    }
    views <- c(
        pathways = "gaussian", immunecells = "gaussian",
        tfs = "gaussian", lrpairs = "gaussian", ccpairs = "gaussian"
    )
    miss_views <- c(ifelse(missing(pathways), NA, 1), ifelse(missing(immunecells),
        NA, 2
    ), ifelse(missing(tfs), NA, 3), ifelse(missing(lrpairs),
        NA, 4
    ), ifelse(missing(ccpairs), NA, 5))
    view_simples <- lapply(miss_views[!is.na(miss_views)], function(X) {
        tmp <- views[X]
        return(tmp)
    })
    view_combinations <- view_simples
    opt_models <- readRDS(file.path(easier_root, "opt_models.rds"))
    opt_xtrain_stats <- readRDS(file.path(easier_root, "opt_xtrain_stats.rds"))
    compute_prediction <- function(view, verbose, opt_models,
                                   opt_xtrain_stats, cancer_type) {
        view_info <- view_combinations[[view]]
        view_name <- paste(names(view_info), collapse = "_")
        view_data <- lapply(tolower(names(view_info)), function(x) as.data.frame(get(x)))
        names(view_data) <- names(view_info)
        if (verbose) {
            message(
                "Computing predictions using ", view_name,
                "...\n"
            )
        }
        opt_model_cancer_view_spec <- lapply(view_name, function(X) {
            return(opt_models[[cancer_type]][[X]])
        })
        names(opt_model_cancer_view_spec) <- view_name
        opt_xtrain_stats_cancer_view_spec <- lapply(
            view_name,
            function(X) {
                return(opt_xtrain_stats[[cancer_type]][[X]])
            }
        )
        names(opt_xtrain_stats_cancer_view_spec) <- view_name
        prediction_view <- predict_with_rmtlr(
            view_name = view_name,
            view_info = view_info, view_data = view_data, opt_model_cancer_view_spec = opt_model_cancer_view_spec,
            opt_xtrain_stats_cancer_view_spec = opt_xtrain_stats_cancer_view_spec
        )
        return(prediction_view)
    }
    BiocParallel::register(BiocParallel::MulticoreParam(workers = 2))
    all_predictions <- BiocParallel::bplapply(seq_len(length(view_combinations)),
        FUN = compute_prediction, verbose = verbose, opt_models,
        opt_xtrain_stats, cancer_type
    )
    names(all_predictions) <- vapply(seq_len(length(view_combinations)),
        function(X) {
            paste(names(view_combinations[[X]]), collapse = "_")
        },
        FUN.VALUE = character(1)
    )
    return(all_predictions)
}
predict_with_rmtlr <- function(view_name,
                               view_info,
                               view_data,
                               opt_model_cancer_view_spec,
                               opt_xtrain_stats_cancer_view_spec,
                               verbose = TRUE) {
    # initialize variables
    P <- length(view_info)
    K <- ncol(opt_model_cancer_view_spec[[1]][[1]])
    standardize_any <- TRUE
    tasks <- names(opt_model_cancer_view_spec[[1]])

    # Algorithm do not deal with NA values: here we removed features with NA values,
    # patients with all NA values are not removed
    view_data_new <- lapply(names(view_data), function(x) {
        tmp_data <- view_data[[x]]
        if (all(is.na(apply(tmp_data, 1, sum)))) {
            NA_sum <- apply(tmp_data, 2, sum)
            tmp_data <- tmp_data[, !is.na(NA_sum)]
        }
        return(tmp_data)
    })
    names(view_data_new) <- names(view_data)

    state <- opt_model_cancer_view_spec
    prediction_X <- view_data_new

    # standardize
    if (standardize_any == TRUE) {
        for (m in seq_len(P)) {
            # Check features availability (we account for mismatches in
            # new releases of progeny, dorothea or quantiseqr)
            keep_pos <- stats::na.omit(match(
                colnames(prediction_X[[m]]),
                rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean)
            ))
            keep_names <- intersect(
                colnames(prediction_X[[m]]),
                rownames(opt_xtrain_stats_cancer_view_spec[[m]]$mean)
            )
            prediction_X[[m]] <- prediction_X[[m]][, keep_names]

            # Normalization should be done taking into account the train set
            opt_xtrain_stats_cancer_view_spec[[m]]$mean <- opt_xtrain_stats_cancer_view_spec[[m]]$mean[keep_pos, ]
            opt_xtrain_stats_cancer_view_spec[[m]]$sd <- opt_xtrain_stats_cancer_view_spec[[m]]$sd[keep_pos, ]

            prediction_X_norm <- lapply(seq_len(K), function(k) {
                prediction_X[[m]] <- calc_z_score(
                    X = prediction_X[[m]], mean = opt_xtrain_stats_cancer_view_spec[[m]]$mean[, k],
                    sd = opt_xtrain_stats_cancer_view_spec[[m]]$sd[, k]
                )
            })
        }
    }
    # perform prediction
    prediction_cv <- lapply(seq_len(K), function(k) {
        coef_matrix <- vapply(tasks, function(task) {
            state[[view_name]][[task]][, k]
        }, FUN.VALUE = numeric(nrow(state[[1]][[1]])))
        rmtlr_test(prediction_X_norm[[k]], coef_matrix)
    })

    predictions_all_tasks_cv <- lapply(tasks, function(task) {
        prediction_task_cv <- vapply(seq_len(K), function(k) {
            prediction_cv[[k]][, task]
        }, FUN.VALUE = numeric(nrow(prediction_X[[1]])))
        return(prediction_task_cv)
    })
    names(predictions_all_tasks_cv) <- tasks
    return(predictions_all_tasks_cv)
}

calc_z_score <- function(X,
                         mean,
                         sd) {
    X_scale <- matrix(0, nrow(X), ncol(X),
        dimnames = list(rownames(X), colnames(X))
    )

    if (missing(mean) & missing(sd)) {
        mean_X <- colMeans(X, na.rm = TRUE)
        sd_X <- matrixStats::colSds(as.matrix(X), na.rm = TRUE)
        X_scale <- sweep(X, 2, mean_X, FUN = "-")
        X_scale <- sweep(X_scale, 2, sd_X, FUN = "/")
    } else {
        mean <- mean[na.omit(match(colnames(X), names(mean)))]
        sd <- sd[na.omit(match(colnames(X), names(sd)))]

        X <- X[, na.omit(match(names(sd), colnames(X)))]
        X_scale <- sweep(X, 2, mean, FUN = "-")
        X_scale <- sweep(X_scale, 2, sd, FUN = "/")
    }
    return(as.matrix(X_scale))
}
rmtlr_test <- function(x_test,
                       coef_matrix) {
    # Keep intercept
    intercept <- as.matrix(coef_matrix)[1, ]

    # Remove intercept from coef matrix
    coef <- as.matrix(coef_matrix)[-1, , drop = FALSE]

    # match features properly
    pos <- stats::na.omit(match(colnames(x_test), rownames(coef)))
    coef <- coef[pos, , drop = FALSE]

    if (length(coef) > 1) {
        slope <- coef
        fit_pred <- t(matrix(as.matrix(intercept),
            nrow = ncol(slope),
            ncol = nrow(x_test)
        )
        + t(slope) %*% t(as.matrix(x_test[, rownames(slope)])))
    } else {
        slope <- 0
        fit_pred <- matrix(as.matrix(intercept),
            nrow = ncol(slope),
            ncol = nrow(x_test[[1]])
        )
    }
    return(fit_pred)
}

###### for CRC ##############
setwd(all_easier_dir)
library(tibble)

cancer_use <- setdiff(pan_cancer_set, c("COAD", "READ", "ESCA"))

for (i in c("COAD", "READ")) {
    setwd(all_easier_dir)
    dir.create(i, recursive = TRUE, showWarnings = FALSE)
    setwd(i)

    file_path <- file.path(gdc_root, i)

    all_file <- list.files(file_path, full.names = T)
    tpm_file <- all_file[grep("Symbol_tpm.tsv", all_file)]
    count_file <- all_file[grep("Symbol_count.tsv", all_file)]
    meta_file <- all_file[grep("survival.tsv.gz", all_file)]

    RNA_tpm <- read_table(tpm_file, col_names = T)
    RNA_tpm <- as.data.frame(RNA_tpm)
    rownames(RNA_tpm) <- RNA_tpm$hgnc_symbol
    RNA_tpm <- RNA_tpm[, -1]

    RNA_counts <- read_table(count_file, col_names = T)
    RNA_counts <- as.data.frame(RNA_counts)
    rownames(RNA_counts) <- RNA_counts$hgnc_symbol
    RNA_counts <- RNA_counts[, -1]

    RNA_meta <- read_table(gzfile(meta_file), col_names = T)
    RNA_tpm <- RNA_tpm[, intersect(colnames(RNA_tpm), RNA_meta$sample)]
    RNA_counts <- RNA_counts[, intersect(colnames(RNA_counts), RNA_meta$sample)]

    cancer_type <- "CRC"

    # 2. TME 特征计算
    cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
    pathway_activities <- compute_pathway_activity(RNA_counts = RNA_counts)
    tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)
    lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm, cancer_type = cancer_type)
    ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, cancer_type = cancer_type)


    # 3. 预测 ICB 响应可能性

    predictions <- predict_immune_response2(
        pathways = pathway_activities,
        immunecells = cell_fractions,
        tfs = tf_activities,
        lrpairs = lrpair_weights,
        ccpairs = ccpair_scores,
        cancer_type = cancer_type
    )

    # 4. 提取 easier 分数（可选融合 TMB）
    easier_scores <- retrieve_easier_score(predictions_immune_response = predictions)


    #
    XAF1_factor <- ifelse(RNA_tpm["XAF1", ] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "R", "NR")
    table(XAF1_factor)
    names(XAF1_factor) <- colnames(RNA_tpm)

    output_eval_with_resp <- assess_immune_response(
        predictions_immune_response = predictions,
        patient_response = XAF1_factor,
        RNA_tpm = RNA_tpm
    )

    pdf("2.AUC of different factor in easier score.pdf", height = 6, width = 7)
    print(output_eval_with_resp[[1]])
    dev.off()

    pdf("2.AUC of different factor in easier score2.pdf", height = 6, width = 10)
    print(output_eval_with_resp[[2]])
    dev.off()

    # output_eval_with_resp[[3]]

    output_biomarkers <- explore_biomarkers(
        pathways = pathway_activities,
        immunecells = cell_fractions,
        tfs = tf_activities,
        lrpairs = lrpair_weights,
        ccpairs = ccpair_scores,
        cancer_type = cancer_type,
        patient_label = XAF1_factor
    )

    pdf("3.pathway descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[1]])
    dev.off()
    pdf("3.immune cells descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[2]])
    dev.off()
    pdf("3.TFs descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[3]])
    dev.off()
    pdf("3.L-R pairs descriptor.pdf", height = 6, width = 12)
    print(output_biomarkers[[4]])
    dev.off()
    pdf("3.C-C pairs descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[5]])
    dev.off()
    pdf("3.All descriptor.pdf", height = 10, width = 12)
    print(output_biomarkers[[6]])
    dev.off()

    XAF1_factor2 <- ifelse(RNA_tpm["XAF1", ] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "XAF1_high", "XAF1_low")

    easier_df <- data.frame(easier_scores = easier_scores$easier_score, group = as.character(XAF1_factor2))

    library(ggplot2)
    library(ggpubr)

    # 自定义颜色，可根据喜好调整
    group_colors <- c("XAF1_high" = "#377eb8", "XAF1_low" = "#e41a1c")

    # 画图
    p <- ggplot(easier_df, aes(x = group, y = easier_scores, fill = group)) +
        geom_boxplot(width = 0.5, alpha = 0.7, color = "gray30", outlier.shape = NA) +
        geom_jitter(aes(color = group), width = 0.15, size = 2, alpha = 0.8) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.format", vjust = 0, hjust = -1,
            label.y = max(easier_df$easier_scores) * 1.05,
            size = 5
        ) +
        scale_fill_manual(values = group_colors) +
        scale_color_manual(values = group_colors) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        labs(
            y = "Relative easier immune response score",
            title = "Immune response score in XAF1 high vs low groups"
        )

    pdf("1.easier immune response score.pdf", height = 6, width = 7)
    print(p)
    dev.off()

    write.table(easier_df, file = "easier_score.txt", quote = F, row.names = F, sep = "\t")

    cell_fractions <- rownames_to_column(cell_fractions, var = "patient")
    pathway_activities <- rownames_to_column(pathway_activities, var = "patient")
    tf_activities <- rownames_to_column(tf_activities, var = "patient")
    lrpair_weights <- rownames_to_column(lrpair_weights, var = "patient")
    ccpair_scores <- rownames_to_column(ccpair_scores, var = "patient")

    write.table(cell_fractions, file = "cell_fraction.txt", quote = F, sep = "\t", row.names = F)
    write.table(pathway_activities, file = "pathway_activities.txt", quote = F, sep = "\t", row.names = F)
    write.table(tf_activities, file = "tf_activities.txt", quote = F, sep = "\t", row.names = F)
    write.table(lrpair_weights, file = "lrpair_weights.txt", quote = F, sep = "\t", row.names = F)
    write.table(ccpair_scores, file = "ccpair_scores.txt", quote = F, sep = "\t", row.names = F)
}

##### for all other cancer #######

for (i in cancer_use) {
    setwd(all_easier_dir)
    dir.create(i, recursive = TRUE, showWarnings = FALSE)
    setwd(i)

    file_path <- file.path(gdc_root, i)

    all_file <- list.files(file_path, full.names = T)
    tpm_file <- all_file[grep("Symbol_tpm.tsv", all_file)]
    count_file <- all_file[grep("Symbol_count.tsv", all_file)]
    meta_file <- all_file[grep("survival.tsv.gz", all_file)]

    RNA_tpm <- read_table(tpm_file, col_names = T)
    RNA_tpm <- as.data.frame(RNA_tpm)
    rownames(RNA_tpm) <- RNA_tpm$hgnc_symbol
    RNA_tpm <- RNA_tpm[, -1]

    RNA_counts <- read_table(count_file, col_names = T)
    RNA_counts <- as.data.frame(RNA_counts)
    rownames(RNA_counts) <- RNA_counts$hgnc_symbol
    RNA_counts <- RNA_counts[, -1]

    RNA_meta <- read_table(gzfile(meta_file), col_names = T)
    RNA_tpm <- RNA_tpm[, intersect(colnames(RNA_tpm), RNA_meta$sample)]
    RNA_counts <- RNA_counts[, intersect(colnames(RNA_counts), RNA_meta$sample)]

    cancer_type <- i

    # 2. TME 特征计算
    cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
    pathway_activities <- compute_pathway_activity(RNA_counts = RNA_counts)
    tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)
    lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm, cancer_type = "pancan")
    ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, cancer_type = "pancan")

    # 3. 预测 ICB 响应可能性

    predictions <- predict_immune_response2(
        pathways = pathway_activities,
        immunecells = cell_fractions,
        tfs = tf_activities,
        lrpairs = lrpair_weights,
        ccpairs = ccpair_scores,
        cancer_type = cancer_type
    )

    # 4. 提取 easier 分数（可选融合 TMB）
    easier_scores <- retrieve_easier_score(predictions_immune_response = predictions)

    #
    XAF1_factor <- ifelse(RNA_tpm["XAF1", ] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "R", "NR")
    table(XAF1_factor)
    names(XAF1_factor) <- colnames(RNA_tpm)

    output_eval_with_resp <- assess_immune_response(
        predictions_immune_response = predictions,
        patient_response = XAF1_factor,
        RNA_tpm = RNA_tpm
    )

    pdf("2.AUC of different factor in easier score.pdf", height = 6, width = 7)
    print(output_eval_with_resp[[1]])
    dev.off()

    pdf("2.AUC of different factor in easier score2.pdf", height = 6, width = 10)
    print(output_eval_with_resp[[2]])
    dev.off()

    # output_eval_with_resp[[3]]

    output_biomarkers <- explore_biomarkers(
        pathways = pathway_activities,
        immunecells = cell_fractions,
        tfs = tf_activities,
        lrpairs = lrpair_weights,
        ccpairs = ccpair_scores,
        cancer_type = cancer_type,
        patient_label = XAF1_factor
    )

    pdf("3.pathway descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[1]])
    dev.off()
    pdf("3.immune cells descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[2]])
    dev.off()
    pdf("3.TFs descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[3]])
    dev.off()
    pdf("3.L-R pairs descriptor.pdf", height = 6, width = 12)
    print(output_biomarkers[[4]])
    dev.off()
    pdf("3.C-C pairs descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[5]])
    dev.off()
    pdf("3.All descriptor.pdf", height = 10, width = 12)
    print(output_biomarkers[[6]])
    dev.off()

    XAF1_factor2 <- ifelse(RNA_tpm["XAF1", ] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "XAF1_high", "XAF1_low")

    easier_df <- data.frame(easier_scores = easier_scores$easier_score, group = as.character(XAF1_factor2))

    library(ggplot2)
    library(ggpubr)

    # 自定义颜色，可根据喜好调整
    group_colors <- c("XAF1_high" = "#377eb8", "XAF1_low" = "#e41a1c")

    # 画图
    p <- ggplot(easier_df, aes(x = group, y = easier_scores, fill = group)) +
        geom_boxplot(width = 0.5, alpha = 0.7, color = "gray30", outlier.shape = NA) +
        geom_jitter(aes(color = group), width = 0.15, size = 2, alpha = 0.8) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.format", vjust = 0, hjust = -1,
            label.y = max(easier_df$easier_scores) * 1.05,
            size = 5
        ) +
        scale_fill_manual(values = group_colors) +
        scale_color_manual(values = group_colors) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        labs(
            y = "Relative easier immune response score",
            title = "Immune response score in XAF1 high vs low groups"
        )

    pdf("1.easier immune response score.pdf", height = 6, width = 7)
    print(p)
    dev.off()

    write.table(easier_df, file = "easier_score.txt", quote = F, row.names = F, sep = "\t")

    cell_fractions <- rownames_to_column(cell_fractions, var = "patient")
    pathway_activities <- rownames_to_column(pathway_activities, var = "patient")
    tf_activities <- rownames_to_column(tf_activities, var = "patient")
    lrpair_weights <- rownames_to_column(lrpair_weights, var = "patient")
    ccpair_scores <- rownames_to_column(ccpair_scores, var = "patient")

    write.table(cell_fractions, file = "cell_fraction.txt", quote = F, sep = "\t", row.names = F)
    write.table(pathway_activities, file = "pathway_activities.txt", quote = F, sep = "\t", row.names = F)
    write.table(tf_activities, file = "tf_activities.txt", quote = F, sep = "\t", row.names = F)
    write.table(lrpair_weights, file = "lrpair_weights.txt", quote = F, sep = "\t", row.names = F)
    write.table(ccpair_scores, file = "ccpair_scores.txt", quote = F, sep = "\t", row.names = F)
}

##### add ESCA 等肿瘤类型 #####

cancer_use <- intersect("ESCA", pan_cancer_set)

for (i in cancer_use) {
    setwd(all_easier_dir)
    dir.create(i, recursive = TRUE, showWarnings = FALSE)
    setwd(i)

    file_path <- file.path(gdc_root, i)

    all_file <- list.files(file_path, full.names = T)
    tpm_file <- all_file[grep("Symbol_tpm.tsv", all_file)]
    count_file <- all_file[grep("Symbol_count.tsv", all_file)]
    meta_file <- all_file[grep("survival.tsv.gz", all_file)]

    RNA_tpm <- read_table(tpm_file, col_names = T)
    RNA_tpm <- as.data.frame(RNA_tpm)
    rownames(RNA_tpm) <- RNA_tpm$hgnc_symbol
    RNA_tpm <- RNA_tpm[, -1]

    RNA_counts <- read_table(count_file, col_names = T)
    RNA_counts <- as.data.frame(RNA_counts)
    rownames(RNA_counts) <- RNA_counts$hgnc_symbol
    RNA_counts <- RNA_counts[, -1]

    RNA_meta <- read_table(gzfile(meta_file), col_names = T)
    RNA_tpm <- RNA_tpm[, intersect(colnames(RNA_tpm), RNA_meta$sample)]
    RNA_counts <- RNA_counts[, intersect(colnames(RNA_counts), RNA_meta$sample)]

    cancer_type <- "STAD"

    # 2. TME 特征计算
    cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
    pathway_activities <- compute_pathway_activity(RNA_counts = RNA_counts)
    tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)
    lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm, cancer_type = "pancan")
    ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, cancer_type = "pancan")

    # 3. 预测 ICB 响应可能性

    predictions <- predict_immune_response2(
        pathways = pathway_activities,
        immunecells = cell_fractions,
        tfs = tf_activities,
        lrpairs = lrpair_weights,
        ccpairs = ccpair_scores,
        cancer_type = cancer_type
    )

    # 4. 提取 easier 分数（可选融合 TMB）
    easier_scores <- retrieve_easier_score(predictions_immune_response = predictions)

    #
    XAF1_factor <- ifelse(RNA_tpm["XAF1", ] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "R", "NR")
    table(XAF1_factor)
    names(XAF1_factor) <- colnames(RNA_tpm)

    output_eval_with_resp <- assess_immune_response(
        predictions_immune_response = predictions,
        patient_response = XAF1_factor,
        RNA_tpm = RNA_tpm
    )

    pdf("2.AUC of different factor in easier score.pdf", height = 6, width = 7)
    print(output_eval_with_resp[[1]])
    dev.off()

    pdf("2.AUC of different factor in easier score2.pdf", height = 6, width = 10)
    print(output_eval_with_resp[[2]])
    dev.off()

    # output_eval_with_resp[[3]]

    output_biomarkers <- explore_biomarkers(
        pathways = pathway_activities,
        immunecells = cell_fractions,
        tfs = tf_activities,
        lrpairs = lrpair_weights,
        ccpairs = ccpair_scores,
        cancer_type = cancer_type,
        patient_label = XAF1_factor
    )

    pdf("3.pathway descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[1]])
    dev.off()
    pdf("3.immune cells descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[2]])
    dev.off()
    pdf("3.TFs descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[3]])
    dev.off()
    pdf("3.L-R pairs descriptor.pdf", height = 6, width = 12)
    print(output_biomarkers[[4]])
    dev.off()
    pdf("3.C-C pairs descriptor.pdf", height = 6, width = 8)
    print(output_biomarkers[[5]])
    dev.off()
    pdf("3.All descriptor.pdf", height = 10, width = 12)
    print(output_biomarkers[[6]])
    dev.off()

    XAF1_factor2 <- ifelse(RNA_tpm["XAF1", ] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "XAF1_high", "XAF1_low")

    easier_df <- data.frame(easier_scores = easier_scores$easier_score, group = as.character(XAF1_factor2))

    library(ggplot2)
    library(ggpubr)

    # 自定义颜色，可根据喜好调整
    group_colors <- c("XAF1_high" = "#377eb8", "XAF1_low" = "#e41a1c")

    # 画图
    p <- ggplot(easier_df, aes(x = group, y = easier_scores, fill = group)) +
        geom_boxplot(width = 0.5, alpha = 0.7, color = "gray30", outlier.shape = NA) +
        geom_jitter(aes(color = group), width = 0.15, size = 2, alpha = 0.8) +
        stat_compare_means(
            method = "wilcox.test",
            label = "p.format", vjust = 0, hjust = -1,
            label.y = max(easier_df$easier_scores) * 1.05,
            size = 5
        ) +
        scale_fill_manual(values = group_colors) +
        scale_color_manual(values = group_colors) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        labs(
            y = "Relative easier immune response score",
            title = "Immune response score in XAF1 high vs low groups"
        )

    pdf("1.easier immune response score.pdf", height = 6, width = 7)
    print(p)
    dev.off()

    write.table(easier_df, file = "easier_score.txt", quote = F, row.names = F, sep = "\t")

    cell_fractions <- rownames_to_column(cell_fractions, var = "patient")
    pathway_activities <- rownames_to_column(pathway_activities, var = "patient")
    tf_activities <- rownames_to_column(tf_activities, var = "patient")
    lrpair_weights <- rownames_to_column(lrpair_weights, var = "patient")
    ccpair_scores <- rownames_to_column(ccpair_scores, var = "patient")

    write.table(cell_fractions, file = "cell_fraction.txt", quote = F, sep = "\t", row.names = F)
    write.table(pathway_activities, file = "pathway_activities.txt", quote = F, sep = "\t", row.names = F)
    write.table(tf_activities, file = "tf_activities.txt", quote = F, sep = "\t", row.names = F)
    write.table(lrpair_weights, file = "lrpair_weights.txt", quote = F, sep = "\t", row.names = F)
    write.table(ccpair_scores, file = "ccpair_scores.txt", quote = F, sep = "\t", row.names = F)
}


###### 整合所有的easier结果 ####
library(rstatix)
library(ggpubr)
setwd(output_dir)

keep_pan_cancer_files <- function(paths) {
    paths[basename(dirname(paths)) %in% pan_cancer_set]
}

validate_pan_cancer_files <- function(paths, label) {
    paths <- keep_pan_cancer_files(paths)
    cancers <- basename(dirname(paths))
    if (length(paths) != 19L || !setequal(cancers, pan_cancer_set)) {
        stop(label, " files do not match the required 19-cancer set.")
    }
    paths
}

all.score <- list.files(all_easier_dir, recursive = TRUE, pattern = "easier_score", full.names = TRUE) |>
    validate_pan_cancer_files("EaSIeR score")

score.list <- list()
for (i in seq(all.score)) {
    d1 <- read.table(all.score[i], header = T, sep = "\t")
    d1$cancer <- basename(dirname(all.score[i]))
    score.list[[i]] <- d1
}
all.df.score <- do.call(rbind, score.list)

stat_results <- all.df.score %>%
    group_by(cancer) %>%
    wilcox_test(easier_scores ~ group) %>% # from rstatix
    mutate(p.adj = p.adjust(p, method = "BH")) %>%
    add_significance("p.adj")

top_cancers <- stat_results %>%
    arrange(p.adj) %>%
    pull(cancer)

all.df.score$cancer <- factor(all.df.score$cancer, levels = top_cancers)

group_colors <- c("XAF1_high" = "#377eb8", "XAF1_low" = "#e41a1c")

p1 <- ggplot(all.df.score, aes(x = group, y = easier_scores, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    facet_wrap(~cancer, scales = "free_y") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.4, size = 6) +
    theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 16),
        strip.text = element_text(hjust = 0.5, face = "bold", size = 15)
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(
        y = "Relative easier immune response score"
    )

pdf("4.all easier.pdf", height = 8, width = 12)
print(p1)
dev.off()

all.pathway <- list.files(all_easier_dir, recursive = TRUE, pattern = "pathway_activities", full.names = TRUE) |>
    validate_pan_cancer_files("Pathway activity")
all.tf_activity <- list.files(all_easier_dir, recursive = TRUE, pattern = "tf_activities", full.names = TRUE) |>
    validate_pan_cancer_files("TF activity")
all.LR_weight <- list.files(all_easier_dir, recursive = TRUE, pattern = "lrpair_weights", full.names = TRUE) |>
    validate_pan_cancer_files("Ligand-receptor")
all.cell_fraction <- list.files(all_easier_dir, recursive = TRUE, pattern = "cell_fraction", full.names = TRUE) |>
    validate_pan_cancer_files("Cell fraction")

used_pathway <- c("Hypoxia", "JAK-STAT", "MAPK", "NFkB", "p53", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT")
used_tf <- c("STAT4", "IRF1", "STAT1", "RELB", "STAT5B", "RUNX3", "WT1", "STAT2", "RFX5", "TP63")
used_LR <- c("CXCL10_SDC4", "IFNG_IFNGR1_IFNGR2", "CCL5_SDC4", "B2M_HLA-F", "CXCL16_CXCR6", "HLA-B_HLA-E_KLRD1", "GZMB_IGF2R", "ICAM1_IL2RG", "CCL8_CCR1", "FASLG_FAS")
used_cell <- c("CD8+ T")

library(effsize)

all.score.list <- list()

file_for_cancer <- function(paths, cancer_type) {
    matched <- paths[basename(dirname(paths)) == cancer_type]
    if (length(matched) != 1L) {
        stop("Expected one file for ", cancer_type, "; found ", length(matched))
    }
    matched
}

for (i in seq(all.pathway)) {
    cancer_type <- basename(dirname(all.pathway[i]))

    ### calculate easier score
    easier_score <- read.table(file_for_cancer(all.score, cancer_type), header = T, sep = "\t")
    easier_score$group <- factor(easier_score$group, levels = c("XAF1_high", "XAF1_low"))
    cohen.d.value <- cohen.d(easier_scores ~ group, data = easier_score)
    cliff.delta.value <- cliff.delta(easier_scores ~ group, data = easier_score)
    wilcox.p <- wilcox.test(easier_scores ~ group, data = easier_score)

    all.score.list <- append(all.score.list, list(c(cohen.d.value$estimate, cliff.delta.value$estimate, wilcox.p$p.value, "ICB_response_score", "ICB response score", cancer_type)))

    ### pathway score
    pathway_activities <- read.table(all.pathway[i], header = T, sep = "\t", check.names = F)

    ### group

    file_path <- file.path(gdc_root, cancer_type)
    all_file <- list.files(file_path, full.names = T)
    tpm_file <- all_file[grep("Symbol_tpm.tsv", all_file)]
    RNA_tpm <- read_table(tpm_file, col_names = T)
    RNA_tpm <- as.data.frame(RNA_tpm)
    rownames(RNA_tpm) <- RNA_tpm$hgnc_symbol
    RNA_tpm <- RNA_tpm[, -1]
    # RNA_tpm <- RNA_tpm[,pathway_activities$patient]

    XAF1_factor <- ifelse(RNA_tpm["XAF1", pathway_activities$patient] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "XAF1_high", "XAF1_low")

    pathway_activities$group <- as.character(XAF1_factor)

    ### pathway score

    for (j in used_pathway) {
        pathway_df <- pathway_activities[, c(j, "group")]
        colnames(pathway_df)[1] <- "score"
        pathway_df$group <- factor(pathway_df$group, levels = c("XAF1_high", "XAF1_low"))
        cohen.d.value <- cohen.d(score ~ group, data = pathway_df)
        cliff.delta.value <- cliff.delta(score ~ group, data = pathway_df)
        wilcox.p <- wilcox.test(score ~ group, data = pathway_df)

        all.score.list <- append(all.score.list, list(c(cohen.d.value$estimate, cliff.delta.value$estimate, wilcox.p$p.value, j, "Pathway activity", cancer_type)))
    }

    tf_activities <- read.table(file_for_cancer(all.tf_activity, cancer_type), header = T, sep = "\t", check.names = F)

    XAF1_factor <- ifelse(RNA_tpm["XAF1", tf_activities$patient] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "XAF1_high", "XAF1_low")
    tf_activities$group <- as.character(XAF1_factor)

    for (j in used_tf) {
        pathway_df <- tf_activities[, c(j, "group")]
        colnames(pathway_df)[1] <- "score"
        pathway_df$group <- factor(pathway_df$group, levels = c("XAF1_high", "XAF1_low"))
        cohen.d.value <- cohen.d(score ~ group, data = pathway_df)
        cliff.delta.value <- cliff.delta(score ~ group, data = pathway_df)
        wilcox.p <- wilcox.test(score ~ group, data = pathway_df)

        all.score.list <- append(all.score.list, list(c(cohen.d.value$estimate, cliff.delta.value$estimate, wilcox.p$p.value, j, "TF activity", cancer_type)))
    }

    lrpair_weights <- read.table(file_for_cancer(all.LR_weight, cancer_type), header = T, sep = "\t", check.names = F)

    XAF1_factor <- ifelse(RNA_tpm["XAF1", lrpair_weights$patient] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "XAF1_high", "XAF1_low")
    lrpair_weights$group <- as.character(XAF1_factor)

    for (j in used_LR) {
        pathway_df <- lrpair_weights[, c(j, "group")]
        colnames(pathway_df)[1] <- "score"
        pathway_df$group <- factor(pathway_df$group, levels = c("XAF1_high", "XAF1_low"))
        cohen.d.value <- cohen.d(score ~ group, data = pathway_df)
        cliff.delta.value <- cliff.delta(score ~ group, data = pathway_df)
        wilcox.p <- wilcox.test(score ~ group, data = pathway_df)

        all.score.list <- append(all.score.list, list(c(cohen.d.value$estimate, cliff.delta.value$estimate, wilcox.p$p.value, j, "L-R expression", cancer_type)))
    }

    cell_fractions <- read.table(file_for_cancer(all.cell_fraction, cancer_type), header = T, sep = "\t", check.names = F)

    XAF1_factor <- ifelse(RNA_tpm["XAF1", cell_fractions$patient] > median(as.numeric(as.character(RNA_tpm["XAF1", ]))), "XAF1_high", "XAF1_low")
    cell_fractions$group <- as.character(XAF1_factor)

    for (j in used_cell) {
        pathway_df <- cell_fractions[, c(j, "group")]
        colnames(pathway_df)[1] <- "score"
        pathway_df$group <- factor(pathway_df$group, levels = c("XAF1_high", "XAF1_low"))
        cohen.d.value <- cohen.d(score ~ group, data = pathway_df)
        cliff.delta.value <- cliff.delta(score ~ group, data = pathway_df)
        wilcox.p <- wilcox.test(score ~ group, data = pathway_df)

        all.score.list <- append(all.score.list, list(c(cohen.d.value$estimate, cliff.delta.value$estimate, wilcox.p$p.value, j, "Cell fraction", cancer_type)))
    }
}

all.score.df <- do.call(rbind, all.score.list) %>% as.data.frame()
colnames(all.score.df) <- c("cohen.d", "cliff.delta", "wilcox.p", "subtype", "group", "cancer")
all.score.df$p.adj <- p.adjust(all.score.df$wilcox.p, method = "BH")

all.score.df <- all.score.df %>%
    mutate(
        cancer = factor(cancer, levels = unique(cancer)),
        cohen.d = as.numeric(cohen.d),
        cliff.delta = as.numeric(cliff.delta)
    ) %>%
    mutate(signif = ifelse(p.adj < 0.05, "Significant", ""))

all.score.df$group <- factor(all.score.df$group, levels = c("ICB response score", "Cell fraction", "Pathway activity", "TF activity", "L-R expression"))

subtype_order <- all.score.df %>%
    group_by(group, subtype) %>%
    summarise(mean_cliff = mean(cliff.delta, na.rm = TRUE)) %>%
    arrange(group, desc(mean_cliff)) %>%
    mutate(subtype = factor(subtype, levels = unique(subtype)))

all.score.df <- all.score.df %>%
    mutate(subtype = factor(subtype, levels = rev(subtype_order$subtype)))

p <- ggplot() +
    geom_tile(
        data = all.score.df, aes(x = cancer, y = subtype),
        fill = "#d9d9d9", color = "#ffffff", width = 0.65, height = 0.65
    ) +
    geom_tile(
        data = all.score.df %>% na.omit() %>% dplyr::filter(p.adj < 0.05),
        aes(x = cancer, y = subtype, fill = cliff.delta),
        color = NA, width = 0.65, height = 0.65
    ) +
    #  geom_tile(data = all.score.df %>%
    #              filter(!group %in% c("ICB response score","Cell fraction")) %>%
    #              group_by(cancer, group) %>%
    #              slice_max(order_by = cliff.delta, n = 1, with_ties = FALSE) %>%
    #              ungroup(),
    #            aes(x = cancer, y = subtype, fill = cliff.delta),
    #            color = "black", width = 0.65, height = 0.65,linewidth = 0.3)+
    geom_star(
        data = all.score.df %>% na.omit() %>% dplyr::filter(p.adj < 0.05),
        aes(x = cancer, y = subtype, starshape = signif),
        size = 2, fill = "#d7301f", color = "#d7301f"
    ) +
    scale_fill_gradient2(
        low = "#984EA3", mid = "white", high = "#FF7F00",
        midpoint = 0, name = bquote(~Cliff ~ "'" ~ s ~ "Delta" ~ " (" ~ italic(XAF1)^high / ~ italic(XAF1)^low ~ ")")
    ) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "") +
    coord_equal(clip = "off") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face = "bold", color = "#000000", size = 12, angle = 45, hjust = 0),
        legend.background = element_rect(color = "#737373")
    )
p

for (i in c(10.5, 20.5, 31.5, 32.5)) {
    p <- p + annotation_custom(
        grob = linesGrob(
            x = unit(c(0, 1), "npc"),
            y = unit(c(0, 1), "npc"),
            default.units = "npc",
            arrow = NULL, name = NULL,
            gp = gpar(), vp = NULL
        ),
        xmin = 0.5, xmax = 20.5, ymin = i, ymax = i
    )
}

p1 <- all.score.df %>%
    dplyr::select(subtype, group) %>%
    # dplyr::mutate(method = factor(method, levels = rev(unique(method)), ordered = T)) %>%
    ggplot() +
    geom_tile(aes(x = 1, y = subtype, fill = group), color = NA, width = 0.5, height = 0.5) +
    scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#A65628")) +
    coord_equal() +
    labs(x = "", y = "") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", color = "#000000", size = 12),
        legend.background = element_rect(color = "#737373")
    )
p1

p_combine <- p1 + p + plot_layout(guides = "collect")
p_combine

setwd(output_dir)

pdf("5.all easier detector.pdf", height = 12, width = 12)
print(p_combine)
dev.off()
