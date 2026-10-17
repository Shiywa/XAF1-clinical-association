#!/usr/bin/env Rscript

# Fig. 7f: pan-cancer XAF1 methylation, expression and EaSIeR associations.
# Run from the repository root. External data and model locations are supplied
# through config/environment.example rather than embedded workstation paths.

source(file.path(
    Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd()),
    "scripts", "common", "config.R"
))

library(readr)
library(dplyr)
library(patchwork)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

###

library(survival)

# --- Input and output locations
base_dir <- xaf1_env_path("XAF1_GDC_ROOT", must_exist = TRUE)
easier_root <- xaf1_env_path("XAF1_EASIER_ROOT", must_exist = TRUE)
result.path <- xaf1_output_dir("fig7_ed8", "methylation")

pan_cancer_set <- c(
    "BLCA", "BRCA", "CESC", "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD",
    "LUSC", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC", "COAD", "READ", "ESCA"
)
stopifnot(length(pan_cancer_set) == 19L, length(unique(pan_cancer_set)) == 19L)

# 你的注释表（已筛选出 XAF1 相关 probes，且包含 Name 列）
# 假设 anno_use 已存在于环境。如果没有，请先生成：
# anno_use <- anno[,c("chr","pos","strand","Name", ... )] |> as.data.frame() |> filter(grepl("XAF1", UCSC_RefGene_Name))

# 如果 anno_use 还不是 data.frame，请确保如下：
# anno_use <- as.data.frame(anno_use)

# 文件列表（按你之前方式获取）
mythe.file <- list.files(base_dir, pattern = "methylation450", full.names = TRUE, recursive = TRUE)
clinical.file <- list.files(base_dir, pattern = "survival", full.names = TRUE, recursive = TRUE)
tpm.file <- list.files(base_dir, pattern = "Symbol_tpm", full.names = TRUE, recursive = TRUE)
count.file <- list.files(base_dir, pattern = "Symbol_count", full.names = TRUE, recursive = TRUE)

# 将文件按癌种名配对：尝试从路径中提取癌种简称（目录名或文件名包含癌种）
get_cancer_name <- function(path) {
    # 先尝试父目录名
    pdir <- basename(dirname(path))
    if (nzchar(pdir)) {
        return(pdir)
    }
    # 否则用文件名的 TCGA_XXX 前缀
    fname <- basename(path)
    m <- regexpr("TCGA_([A-Z]+)", fname)
    if (m[1] > 0) {
        substr(fname, m[1] + 5, m[1] + attr(m, "match.length") - 1)
    } else {
        # fallback: 文件名中提取第一个全大写词
        mm <- regmatches(fname, regexpr("[A-Z]{2,5}", fname))
        if (length(mm) && nzchar(mm)) mm else fname
    }
}

# 建立 data.frame 列表并合并匹配（按癌种名）
methy_df <- data.frame(path = mythe.file, cancer = sapply(mythe.file, get_cancer_name), stringsAsFactors = FALSE)
tpm_df <- data.frame(path = tpm.file, cancer = sapply(tpm.file, get_cancer_name), stringsAsFactors = FALSE)
clin_df <- data.frame(path = clinical.file, cancer = sapply(clinical.file, get_cancer_name), stringsAsFactors = FALSE)
count_df <- data.frame(path = count.file, cancer = sapply(count.file, get_cancer_name), stringsAsFactors = FALSE)

# 找到三者都存在的癌种
common_cancers <- Reduce(intersect, list(methy_df$cancer, tpm_df$cancer, clin_df$cancer, count_df$cancer))
if (length(common_cancers) == 0) stop("没有找到同时包含 methylation/tpm/clinical 的相同癌种文件，请检查文件名与路径。")
missing_cancers <- setdiff(pan_cancer_set, common_cancers)
if (length(missing_cancers)) {
    stop("Missing methylation/expression/clinical inputs for: ", paste(missing_cancers, collapse = ", "))
}
common_cancers <- pan_cancer_set


# 获取甲基化芯片注释信息
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno[1:5, 1:33]

anno_use <- anno[, c("chr", "pos", "strand", "Name", "Probe_rs", "Probe_maf", "CpG_rs", "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Enhancer", "Regulatory_Feature_Group", "DHS")]
anno_use[1:5, ]

anno_use <- anno_use |>
    as.data.frame() |>
    filter(grepl("XAF1", UCSC_RefGene_Name))

# 如果 anno_use 未定义，报错提醒
if (!exists("anno_use")) stop("请先准备好 anno_use（包含 XAF1 相关 probes，且含 Name 列）。")

# 确保 anno_use 的 Name 列为字符
if (!"Name" %in% colnames(anno_use)) stop("anno_use 必须包含 'Name' 列（probe ID）。")

probes <- unique(as.character(anno_use$Name))

# 输出结果容器
surv_results <- list()
corr_results <- list()

### 分析每个癌种的甲基化的分析 ###########
# 循环每个癌种
for (c in common_cancers) {
    message("Processing cancer: ", c)
    mpath <- methy_df$path[methy_df$cancer == c][1]
    tpath <- tpm_df$path[tpm_df$cancer == c][1]
    clpath <- clin_df$path[clin_df$cancer == c][1]

    # 读取数据（兼容 gz）
    methy <- read_tsv(mpath, col_names = TRUE, show_col_types = FALSE)
    tpm <- read_tsv(tpath, col_names = TRUE, show_col_types = FALSE)
    clin <- read_tsv(clpath, col_names = TRUE, show_col_types = FALSE)

    clin[1:5, 1:4]
    # 确认列名
    # methy 第一列为探针名（Composite Element REF）
    probe_col <- colnames(methy)[1]
    # tpm 第一列为基因名（hgnc_symbol）
    gene_col <- colnames(tpm)[1]
    # clinical 的样本列名可能为 "sample"（你示例是 sample）
    # 确保使用 OS.time 和 OS 字段
    if (!("OS.time" %in% colnames(clin) && "OS" %in% colnames(clin) && "sample" %in% colnames(clin))) {
        stop("临床文件需要包含列：sample, OS.time, OS。癌种: ", c)
    }

    # 找到 XAF1 在 tpm 行（exact match），若没有则跳过相关性分析
    xaf1_row_idx <- which(tpm[[gene_col]] %in% c("XAF1", "xaf1"))
    if (length(xaf1_row_idx) == 0) {
        message("Warning: XAF1 not found in tpm for ", c, " — correlation results will be NA.")
    }

    # 只保留在三类数据中同时存在的样本
    methy_samples <- colnames(methy)[-1]
    tpm_samples <- colnames(tpm)[-1]
    clinic_samples <- clin$sample

    common_samples <- intersect(intersect(methy_samples, tpm_samples), clinic_samples)
    if (length(common_samples) < 10) {
        message("Warning: common samples < 10 for ", c, " (", length(common_samples), "). Skipping if too few.")
    }

    # 以 common_samples 的顺序构建向量
    # 把 methy 转为以 probe 为行、samples 为列快速查找（但不用转置全部矩阵）
    # 遍历 probes
    for (p in probes) {
        # 在 methy 中找到该 probe 行
        idx <- which(methy[[probe_col]] == p)
        if (length(idx) == 0) {
            # probe 不存在于该癌种甲基化数据，记录 NA（可选）
            # 记录在 surv_results/corr_results 中
            for (split_method in c("mean", "Q75", "Q25")) {
                surv_results[[length(surv_results) + 1]] <- data.frame(
                    cancer = c, probe = p, split = split_method,
                    HR = NA_real_, p_value = NA_real_, n_high = NA_integer_, n_low = NA_integer_,
                    stringsAsFactors = FALSE
                )
            }
            corr_results[[length(corr_results) + 1]] <- data.frame(
                cancer = c, probe = p, rho = NA_real_, p_value = NA_real_, n = NA_integer_, stringsAsFactors = FALSE
            )
            next
        }

        # 提取 beta 向量（按 common_samples 顺序）
        beta_vec <- as.numeric(methy[idx, common_samples])
        names(beta_vec) <- common_samples

        # 若大部分为 NA，跳过
        if (sum(!is.na(beta_vec)) < 5) {
            for (split_method in c("mean", "Q75", "Q25")) {
                surv_results[[length(surv_results) + 1]] <- data.frame(
                    cancer = c, probe = p, split = split_method,
                    HR = NA_real_, p_value = NA_real_, n_high = sum(!is.na(beta_vec) & beta_vec > -Inf), n_low = NA_integer_,
                    stringsAsFactors = FALSE
                )
            }
            corr_results[[length(corr_results) + 1]] <- data.frame(
                cancer = c, probe = p, rho = NA_real_, p_value = NA_real_, n = sum(!is.na(beta_vec)), stringsAsFactors = FALSE
            )
            next
        }

        # 提取 XAF1 表达（按 common_samples 顺序）
        expr_vec <- rep(NA_real_, length(common_samples))
        names(expr_vec) <- common_samples
        if (length(xaf1_row_idx) > 0) {
            expr_row <- tpm[xaf1_row_idx[1], common_samples, drop = TRUE]
            expr_vec <- as.numeric(expr_row)
            names(expr_vec) <- common_samples
        }

        # 临床数据按 sample
        clin_sub <- clin %>%
            filter(sample %in% common_samples) %>%
            arrange(match(sample, common_samples))

        if (nrow(clin_sub) != length(common_samples)) {
            # 如果 clinical 有重复或缺失，重新对齐
            clin_map <- clin %>% filter(sample %in% common_samples)
            clin_sub <- clin_map[match(common_samples, clin_map$sample), ]
        }

        # 确保 OS.time 和 OS 非空
        if (any(is.na(clin_sub$OS.time)) || any(is.na(clin_sub$OS))) {
            # 仍可继续但会跳过 Cox 若事件/时间不足
        }

        # 三种分割阈值
        thresholds <- c(
            mean = mean(beta_vec, na.rm = TRUE),
            Q75  = as.numeric(quantile(beta_vec, 0.75, na.rm = TRUE)),
            Q25  = as.numeric(quantile(beta_vec, 0.25, na.rm = TRUE))
        )

        for (split_method in names(thresholds)) {
            th <- thresholds[[split_method]]
            # 分组： high = 1 if beta > th, low = 0 if beta <= th
            group <- as.integer(beta_vec > th)
            # 只保留有观察值的样本（非 NA）
            ok <- !is.na(group) & !is.na(clin_sub$OS.time) & !is.na(clin_sub$OS)
            if (sum(ok) < 10) {
                # 样本太少，记录 NA
                surv_results[[length(surv_results) + 1]] <- data.frame(
                    cancer = c, probe = p, split = split_method,
                    HR = NA_real_, p_value = NA_real_, n_high = sum(group[ok] == 1, na.rm = TRUE),
                    n_low = sum(group[ok] == 0, na.rm = TRUE), stringsAsFactors = FALSE
                )
                next
            }

            # 组内至少要有事件或样本量，否则 Cox 会失败
            n_high <- sum(group[ok] == 1, na.rm = TRUE)
            n_low <- sum(group[ok] == 0, na.rm = TRUE)
            if (n_high < 3 || n_low < 3) {
                surv_results[[length(surv_results) + 1]] <- data.frame(
                    cancer = c, probe = p, split = split_method,
                    HR = NA_real_, p_value = NA_real_, n_high = n_high, n_low = n_low, stringsAsFactors = FALSE
                )
                next
            }

            # 构建生存对象并拟合 Cox（group as factor）
            try(
                {
                    df_cox <- data.frame(
                        time = clin_sub$OS.time,
                        status = clin_sub$OS,
                        group = group
                    )[ok, , drop = FALSE]

                    cox_model <- coxph(Surv(time, status) ~ group, data = df_cox)
                    s <- summary(cox_model)
                    beta_coef <- s$coefficients[, "coef"]
                    hr <- exp(beta_coef)
                    pval <- s$coefficients[, "Pr(>|z|)"] # Wald p-value

                    surv_results[[length(surv_results) + 1]] <- data.frame(
                        cancer = c, probe = p, split = split_method,
                        HR = as.numeric(hr), p_value = as.numeric(pval),
                        n_high = n_high, n_low = n_low, stringsAsFactors = FALSE
                    )
                },
                silent = TRUE
            )
        } # end split methods

        # ---- 相关性分析（Spearman）: beta vs XAF1 expression
        ok_corr <- !is.na(beta_vec) & !is.na(expr_vec)
        if (sum(ok_corr) >= 5 && length(xaf1_row_idx) > 0) {
            try(
                {
                    ct <- suppressWarnings(cor.test(beta_vec[ok_corr], expr_vec[ok_corr], method = "spearman"))
                    corr_results[[length(corr_results) + 1]] <- data.frame(
                        cancer = c, probe = p, rho = as.numeric(ct$estimate), p_value = as.numeric(ct$p.value), n = sum(ok_corr),
                        stringsAsFactors = FALSE
                    )
                },
                silent = TRUE
            )
        } else {
            corr_results[[length(corr_results) + 1]] <- data.frame(
                cancer = c, probe = p, rho = NA_real_, p_value = NA_real_, n = sum(ok_corr), stringsAsFactors = FALSE
            )
        }
    } # end probe loop
} # end cancer loop

# 合并结果表
surv_df <- bind_rows(surv_results)
corr_df <- bind_rows(corr_results)

surv_df_clean <- surv_df %>%
    filter(is.finite(HR), is.finite(p_value))

corr_df_clean <- corr_df %>%
    filter(is.finite(rho), is.finite(p_value))

surv_df_best <- surv_df_clean %>%
    group_by(cancer, probe) %>% # 按癌种 + 探针分组
    slice_min(order_by = p_value, n = 1) %>% # 选取最小 p 值行
    ungroup()

anno_use_clean <- anno_use %>%
    rowwise() %>% # 按行处理
    mutate(
        UCSC_RefGene_Name = paste(unique(unlist(strsplit(UCSC_RefGene_Name, ";"))), collapse = ";"),
        UCSC_RefGene_Group = paste(unique(unlist(strsplit(UCSC_RefGene_Group, ";"))), collapse = ";")
    ) %>%
    ungroup()

surv_df_best <- merge(surv_df_best, anno_use_clean, by.x = "probe", by.y = "Name")
corr_df_clean <- merge(corr_df_clean, anno_use_clean, by.x = "probe", by.y = "Name")

# 保存输出
write_tsv(surv_df_best, file.path(result.path, "XAF1_probe_survival_results.tsv"))
write_tsv(corr_df_clean, file.path(result.path, "XAF1_probe_expression_correlation.tsv"))

message("Finished. Results saved to: ", base_dir)

head(surv_df_best)

plot_HR_volcano <- function(df) {
    library(ggplot2)
    library(ggrepel)
    library(dplyr)

    df <- df %>%
        mutate(adj.p = p.adjust(p_value, method = "BH")) |>
        mutate(
            log10p = -log10(adj.p),
            color_group = case_when(
                adj.p < 0.05 & HR > 1 ~ "pos",
                adj.p < 0.05 & HR < 1 ~ "neg",
                TRUE ~ "ns"
            )
        )
    max_finite <- max(df$log10p[is.finite(df$log10p)], na.rm = TRUE)
    df$log10p[!is.finite(df$log10p)] <- max_finite + 1

    ggplot(df, aes(x = HR, y = log10p)) +
        # --- 点 ---
        geom_point(
            aes(fill = color_group),
            shape = 21,
            color = "black",
            stroke = 0.25,
            size = 6
        ) +

        # --- 正相关显著点标签（右侧） ---
        geom_text_repel(
            data = subset(df, color_group == "pos"),
            aes(label = cancer, x = HR, y = log10p),
            box.padding = 0.5,
            force = 20,
            max.time = 10,
            fontface = "bold",
            nudge_x = 0.08,
            nudge_y = 0.02,
            max.overlaps = Inf,
            point.padding = 0.1,
            segment.curvature = -0.1,
            segment.ncp = 3,
            segment.angle = 20,
            direction = "y",
            hjust = "left",
            size = 6,
            segment.linetype = "dashed",
            segment.alpha = 0.5
        ) +

        # --- 负相关显著点标签（左侧） ---
        geom_text_repel(
            data = subset(df, color_group == "neg"),
            aes(label = cancer),
            box.padding = 0.3,
            fontface = "bold",
            nudge_x = -0.08,
            nudge_y = -0.02,
            max.overlaps = 100,
            point.padding = 0.1,
            segment.curvature = 0.1,
            segment.ncp = 3,
            segment.angle = 20,
            direction = "y",
            hjust = "right",
            size = 6,
            segment.linetype = "dashed",
            segment.alpha = 0.5
        ) +

        # --- 自定义颜色 ---
        scale_fill_manual(
            values = c("pos" = "#1F77B4", "neg" = "#009E73", "ns" = "white"),
            guide = "none"
        ) +

        # --- 坐标轴 ---
        scale_x_continuous(
            limits = c(-0.2, ceiling(max(df$HR)) + 0.2),
            breaks = seq(0, ceiling(max(df$HR)), by = 0.5),
            expand = c(0, 0)
        ) +
        geom_vline(xintercept = 1, color = "black", linewidth = 0.3) +
        geom_text(
            data = data.frame(y = seq(0, ceiling(max(df$log10p)), by = 2)),
            aes(x = 0.97, y = y, label = y), hjust = "right",
            inherit.aes = FALSE, size = 5
        ) +
        geom_segment(
            data = data.frame(y = seq(0, ceiling(max(df$log10p)), by = 2)),
            aes(x = 0.98, xend = 1, y = y, yend = y),
            inherit.aes = FALSE,
            color = "black",
            size = 0.3
        ) +

        # --- 样式 ---
        theme_classic(base_size = 14) +
        labs(
            x = expression(HR),
            y = expression(-log[10](adj.p))
        ) +
        theme(
            axis.line = element_line(size = 0.4),
            axis.ticks = element_line(size = 0.3),
            panel.grid = element_blank(),
            axis.line.y = element_blank(), # 隐藏左侧 y 轴线
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 13, color = "black")
        )
}

### 绘制甲基化与生存率的结果图 #####
for (i in unique(surv_df_best$probe)) {
    d1 <- surv_df_best |>
        filter(probe == i) |>
        as.data.frame()

    pdf(file.path(result.path, paste("TCGA mythelation450 XAF1", i, "survival valcano.pdf", sep = " ")), height = 7, width = 7)
    p1 <- plot_HR_volcano(d1) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title = paste("XAF1 ", i, " (", unique(d1$UCSC_RefGene_Group), ")", " HR", sep = ""),
            theme = theme(
                plot.title = element_text(
                    size = 16, face = "bold", hjust = 0.5,
                    margin = margin(b = 10)
                ),
                plot.margin = margin(20, 20, 20, 20)
            )
        )
    print(p1)
    dev.off()
}


corr_df_best <- corr_df_clean

plot_Rho_volcano <- function(df) {
    library(ggplot2)
    library(ggrepel)
    library(dplyr)

    df <- df %>%
        mutate(adj.p = p.adjust(p_value, method = "BH")) |>
        mutate(
            log10p = -log10(adj.p),
            color_group = case_when(
                adj.p < 0.05 & rho > 0 ~ "pos",
                adj.p < 0.05 & rho < 0 ~ "neg",
                TRUE ~ "ns"
            )
        )
    max_finite <- max(df$log10p[is.finite(df$log10p)], na.rm = TRUE)
    df$log10p[!is.finite(df$log10p)] <- max_finite + 1

    ggplot(df, aes(x = rho, y = log10p)) +
        # --- 点 ---
        geom_point(
            aes(fill = color_group),
            shape = 21,
            color = "black",
            stroke = 0.25,
            size = 6
        ) +

        # --- 正相关显著点标签（右侧） ---
        geom_text_repel(
            data = subset(df, color_group == "pos"),
            aes(label = cancer, x = rho, y = log10p),
            box.padding = 0.5,
            force = 20,
            max.time = 10,
            fontface = "bold",
            nudge_x = 0.08,
            nudge_y = 0.02,
            max.overlaps = Inf,
            point.padding = 0.1,
            segment.curvature = -0.1,
            segment.ncp = 3,
            segment.angle = 20,
            direction = "y",
            hjust = "left",
            size = 6,
            segment.linetype = "dashed",
            segment.alpha = 0.5
        ) +

        # --- 负相关显著点标签（左侧） ---
        geom_text_repel(
            data = subset(df, color_group == "neg"),
            aes(label = cancer, x = rho, y = log10p),
            box.padding = 0.3,
            fontface = "bold",
            nudge_x = -0.08,
            nudge_y = -0.02,
            max.overlaps = 100,
            point.padding = 0.1,
            segment.curvature = 0.1,
            segment.ncp = 3,
            segment.angle = 20,
            direction = "y",
            hjust = "right",
            size = 6,
            segment.linetype = "dashed",
            segment.alpha = 0.5
        ) +

        # --- 自定义颜色 ---
        scale_fill_manual(
            values = c("pos" = "#E63946", "neg" = "#457B9D", "ns" = "white"),
            guide = "none"
        ) +

        # --- 坐标轴 ---
        scale_x_continuous(
            limits = c(-1, 1),
            breaks = seq(-1, 1, by = 0.2),
            expand = c(0, 0)
        ) +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
        geom_text(
            data = data.frame(y = seq(0, ceiling(max(df$log10p)), by = 5)),
            aes(x = -0.03, y = y, label = y), hjust = "right",
            inherit.aes = FALSE, size = 5
        ) +
        geom_segment(
            data = data.frame(y = seq(0, ceiling(max(df$log10p)), by = 5)),
            aes(x = -0.02, xend = 0, y = y, yend = y),
            inherit.aes = FALSE,
            color = "black",
            size = 0.3
        ) +

        # --- 样式 ---
        theme_classic(base_size = 14) +
        labs(
            x = expression(rho),
            y = expression(-log[10](adj.p))
        ) +
        theme(
            axis.line = element_line(size = 0.4),
            axis.ticks = element_line(size = 0.3),
            panel.grid = element_blank(),
            axis.line.y = element_blank(), # 隐藏左侧 y 轴线
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 13, color = "black")
        )
}

### 绘制甲基化与XAF1表达的相关性图 #####
for (i in unique(corr_df_best$probe)) {
    d1 <- corr_df_best |>
        filter(probe == i) |>
        as.data.frame()

    pdf(file.path(result.path, paste("TCGA mythelation450 correlation of XAF1 and", i, ".pdf", sep = " ")), height = 7, width = 7)
    p1 <- plot_Rho_volcano(d1) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title = paste(expression(rho), " of XAF1 and ", i, " (", unique(d1$UCSC_RefGene_Group), ")", sep = ""),
            theme = theme(
                plot.title = element_text(
                    size = 16, face = "bold", hjust = 0.5,
                    margin = margin(b = 10)
                ),
                plot.margin = margin(20, 20, 20, 20)
            )
        )
    print(p1)
    dev.off()
}

### XAF1 expression and OS #######

surv_results_XAF1 <- list()
for (c in common_cancers) {
    message("Processing cancer: ", c)
    mpath <- methy_df$path[methy_df$cancer == c][1]
    tpath <- tpm_df$path[tpm_df$cancer == c][1]
    clpath <- clin_df$path[clin_df$cancer == c][1]

    # 读取数据（兼容 gz）
    methy <- read_tsv(mpath, col_names = TRUE, show_col_types = FALSE)
    tpm <- read_tsv(tpath, col_names = TRUE, show_col_types = FALSE)
    clin <- read_tsv(clpath, col_names = TRUE, show_col_types = FALSE)

    clin[1:5, 1:4]
    # 确认列名
    # methy 第一列为探针名（Composite Element REF）
    probe_col <- colnames(methy)[1]
    # tpm 第一列为基因名（hgnc_symbol）
    gene_col <- colnames(tpm)[1]
    # clinical 的样本列名可能为 "sample"（你示例是 sample）
    # 确保使用 OS.time 和 OS 字段
    if (!("OS.time" %in% colnames(clin) && "OS" %in% colnames(clin) && "sample" %in% colnames(clin))) {
        stop("临床文件需要包含列：sample, OS.time, OS。癌种: ", c)
    }

    # 找到 XAF1 在 tpm 行（exact match），若没有则跳过相关性分析
    xaf1_row_idx <- which(tpm[[gene_col]] %in% c("XAF1", "xaf1"))
    if (length(xaf1_row_idx) == 0) {
        message("Warning: XAF1 not found in tpm for ", c, " — correlation results will be NA.")
    }

    # 只保留在三类数据中同时存在的样本
    methy_samples <- colnames(methy)[-1]
    tpm_samples <- colnames(tpm)[-1]
    clinic_samples <- clin$sample

    common_samples <- intersect(intersect(methy_samples, tpm_samples), clinic_samples)
    if (length(common_samples) < 10) {
        message("Warning: common samples < 10 for ", c, " (", length(common_samples), "). Skipping if too few.")
    }

    # 提取 beta 向量（按 common_samples 顺序）
    XAF1_vec <- as.numeric(tpm[xaf1_row_idx, common_samples])
    names(XAF1_vec) <- common_samples

    # 若大部分为 NA，跳过
    if (sum(!is.na(XAF1_vec)) < 5) {
        for (split_method in c("mean", "Q75", "Q25")) {
            surv_results_XAF1[[length(surv_results_XAF1) + 1]] <- data.frame(
                cancer = c, probe = p, split = split_method,
                HR = NA_real_, p_value = NA_real_, n_high = sum(!is.na(XAF1_vec) & XAF1_vec > -Inf), n_low = NA_integer_,
                stringsAsFactors = FALSE
            )
        }
        next
    }

    # 提取 XAF1 表达（按 common_samples 顺序）
    expr_vec <- rep(NA_real_, length(common_samples))
    names(expr_vec) <- common_samples
    if (length(xaf1_row_idx) > 0) {
        expr_row <- tpm[xaf1_row_idx[1], common_samples, drop = TRUE]
        expr_vec <- as.numeric(expr_row)
        names(expr_vec) <- common_samples
    }

    # 临床数据按 sample
    clin_sub <- clin %>%
        filter(sample %in% common_samples) %>%
        arrange(match(sample, common_samples))

    if (nrow(clin_sub) != length(common_samples)) {
        # 如果 clinical 有重复或缺失，重新对齐
        clin_map <- clin %>% filter(sample %in% common_samples)
        clin_sub <- clin_map[match(common_samples, clin_map$sample), ]
    }

    # 确保 OS.time 和 OS 非空
    if (any(is.na(clin_sub$OS.time)) || any(is.na(clin_sub$OS))) {
        # 仍可继续但会跳过 Cox 若事件/时间不足
    }

    # 三种分割阈值
    thresholds <- c(
        median = median(expr_vec, na.rm = TRUE),
        Q75 = as.numeric(quantile(expr_vec, 0.75, na.rm = TRUE)),
        Q25 = as.numeric(quantile(expr_vec, 0.25, na.rm = TRUE))
    )

    for (split_method in names(thresholds)) {
        th <- thresholds[[split_method]]
        # 分组： high = 1 if beta > th, low = 0 if beta <= th
        group <- as.integer(expr_vec > th)
        # 只保留有观察值的样本（非 NA）
        ok <- !is.na(group) & !is.na(clin_sub$OS.time) & !is.na(clin_sub$OS)
        if (sum(ok) < 10) {
            # 样本太少，记录 NA
            surv_results_XAF1[[length(surv_results_XAF1) + 1]] <- data.frame(
                cancer = c, probe = p, split = split_method,
                HR = NA_real_, p_value = NA_real_, n_high = sum(group[ok] == 1, na.rm = TRUE),
                n_low = sum(group[ok] == 0, na.rm = TRUE), stringsAsFactors = FALSE
            )
            next
        }

        # 组内至少要有事件或样本量，否则 Cox 会失败
        n_high <- sum(group[ok] == 1, na.rm = TRUE)
        n_low <- sum(group[ok] == 0, na.rm = TRUE)
        if (n_high < 3 || n_low < 3) {
            surv_results_XAF1[[length(surv_results_XAF1) + 1]] <- data.frame(
                cancer = c, probe = p, split = split_method,
                HR = NA_real_, p_value = NA_real_, n_high = n_high, n_low = n_low, stringsAsFactors = FALSE
            )
            next
        }

        # 构建生存对象并拟合 Cox（group as factor）
        try(
            {
                df_cox <- data.frame(
                    time = clin_sub$OS.time,
                    status = clin_sub$OS,
                    group = group
                )[ok, , drop = FALSE]

                cox_model <- coxph(Surv(time, status) ~ group, data = df_cox)
                s <- summary(cox_model)
                beta_coef <- s$coefficients[, "coef"]
                hr <- exp(beta_coef)
                pval <- s$coefficients[, "Pr(>|z|)"] # Wald p-value

                surv_results_XAF1[[length(surv_results_XAF1) + 1]] <- data.frame(
                    cancer = c, probe = "XAF1", split = split_method,
                    HR = as.numeric(hr), p_value = as.numeric(pval),
                    n_high = n_high, n_low = n_low, stringsAsFactors = FALSE
                )
            },
            silent = TRUE
        )
    } # end split methods
}

surv_df_XAF1 <- bind_rows(surv_results_XAF1)

surv_df_XAF1_clean <- surv_df_XAF1 %>%
    filter(is.finite(HR), is.finite(p_value))

surv_df_XAF1_best <- surv_df_XAF1_clean %>%
    group_by(cancer) %>% # 按癌种 + 探针分组
    slice_min(order_by = p_value, n = 1) %>% # 选取最小 p 值行
    ungroup()

write_tsv(surv_df_XAF1_best, file.path(result.path, "XAF1_gene_survival_results.tsv"))


pdf(file.path(result.path, "TCGA XAF1 survival valcano.pdf"), height = 7, width = 7)
p1 <- plot_HR_volcano(surv_df_XAF1_best) +
    plot_layout(guides = "collect") +
    plot_annotation(
        title = paste("XAF1 HR", sep = ""),
        theme = theme(
            plot.title = element_text(
                size = 16, face = "bold", hjust = 0.5,
                margin = margin(b = 10)
            ),
            plot.margin = margin(20, 20, 20, 20)
        )
    )
print(p1)
dev.off()

### easier 预测 ICB response ########
library(easier)
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

cancer_use <- pan_cancer_set # Fig. 7f methylation panel: predefined 19-cancer collection.

easier_res <- list()

for (c in cancer_use) {
    message("Processing cancer: ", c)
    mpath <- methy_df$path[methy_df$cancer == c][1]
    tpath <- tpm_df$path[tpm_df$cancer == c][1]
    clpath <- clin_df$path[clin_df$cancer == c][1]
    countpath <- count_df$path[count_df$cancer == c][1]

    # 读取数据（兼容 gz）
    methy <- read_tsv(mpath, col_names = TRUE, show_col_types = FALSE)
    tpm <- read_tsv(tpath, col_names = TRUE, show_col_types = FALSE)
    clin <- read_tsv(clpath, col_names = TRUE, show_col_types = FALSE)
    count <- read_tsv(countpath, col_names = TRUE, show_col_types = FALSE)

    clin[1:5, 1:4]
    # 确认列名
    # methy 第一列为探针名（Composite Element REF）
    probe_col <- colnames(methy)[1]
    # tpm 第一列为基因名（hgnc_symbol）
    gene_col <- colnames(tpm)[1]
    # clinical 的样本列名可能为 "sample"（你示例是 sample）
    # 确保使用 OS.time 和 OS 字段
    if (!("OS.time" %in% colnames(clin) && "OS" %in% colnames(clin) && "sample" %in% colnames(clin))) {
        stop("临床文件需要包含列：sample, OS.time, OS。癌种: ", c)
    }

    # 找到 XAF1 在 tpm 行（exact match），若没有则跳过相关性分析
    xaf1_row_idx <- which(tpm[[gene_col]] %in% c("XAF1", "xaf1"))
    if (length(xaf1_row_idx) == 0) {
        message("Warning: XAF1 not found in tpm for ", c, " — correlation results will be NA.")
    }

    # 只保留在三类数据中同时存在的样本
    methy_samples <- colnames(methy)[-1]
    tpm_samples <- colnames(tpm)[-1]
    clinic_samples <- clin$sample

    common_samples <- intersect(intersect(methy_samples, tpm_samples), clinic_samples)
    if (length(common_samples) < 10) {
        message("Warning: common samples < 10 for ", c, " (", length(common_samples), "). Skipping if too few.")
    }

    if (c %in% c("COAD", "READ")) {
        cancer_type <- "CRC"
    } else if (c %in% c("ESCA")) {
        cancer_type <- "STAD"
    } else {
        cancer_type <- c
    }

    RNA_tpm <- tpm[, common_samples] |>
        as.data.frame()
    rownames(RNA_tpm) <- tpm$hgnc_symbol

    RNA_counts <- count[, common_samples] |>
        as.data.frame()
    rownames(RNA_counts) <- count$hgnc_symbol

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


    for (p in probes) {
        # 在 methy 中找到该 probe 行
        idx <- which(methy[[probe_col]] == p)
        if (length(idx) == 0) {
            # probe 不存在于该癌种甲基化数据，记录 NA（可选）
            # 记录在 surv_results/corr_results 中
            next
        }

        # 提取 beta 向量（按 common_samples 顺序）
        beta_vec <- as.numeric(methy[idx, common_samples])
        names(beta_vec) <- common_samples

        # 若大部分为 NA，跳过
        if (sum(!is.na(beta_vec)) < 5) {
            next
        }

        # 提取 XAF1 表达（按 common_samples 顺序）
        expr_vec <- rep(NA_real_, length(common_samples))
        names(expr_vec) <- common_samples
        if (length(xaf1_row_idx) > 0) {
            expr_row <- tpm[xaf1_row_idx[1], common_samples, drop = TRUE]
            expr_vec <- as.numeric(expr_row)
            names(expr_vec) <- common_samples
        }

        # 临床数据按 sample
        clin_sub <- clin %>%
            filter(sample %in% common_samples) %>%
            arrange(match(sample, common_samples))

        if (nrow(clin_sub) != length(common_samples)) {
            # 如果 clinical 有重复或缺失，重新对齐
            clin_map <- clin %>% filter(sample %in% common_samples)
            clin_sub <- clin_map[match(common_samples, clin_map$sample), ]
        }

        # 确保 OS.time 和 OS 非空
        if (any(is.na(clin_sub$OS.time)) || any(is.na(clin_sub$OS))) {
            # 仍可继续但会跳过 Cox 若事件/时间不足
        }

        # 三种分割阈值
        thresholds <- c(
            mean = mean(beta_vec, na.rm = TRUE),
            Q75  = as.numeric(quantile(beta_vec, 0.75, na.rm = TRUE)),
            Q25  = as.numeric(quantile(beta_vec, 0.25, na.rm = TRUE))
        )

        for (split_method in names(thresholds)) {
            th <- thresholds[[split_method]]
            # 分组： high = 1 if beta > th, low = 0 if beta <= th
            group <- as.integer(beta_vec > th)
            # 只保留有观察值的样本（非 NA）
            ok <- !is.na(group) & !is.na(clin_sub$OS.time) & !is.na(clin_sub$OS)
            if (sum(ok) < 10) {
                next
            }

            # 组内至少要有事件或样本量，否则 Cox 会失败
            n_high <- sum(group[ok] == 1, na.rm = TRUE)
            n_low <- sum(group[ok] == 0, na.rm = TRUE)
            if (n_high < 3 || n_low < 3) {
                next
            }

            # 构建生存对象并拟合 Cox（group as factor）
            try(
                {
                    XAF1_factor <- ifelse(beta_vec > th, "R", "NR")

                    XAF1_factor2 <- ifelse(beta_vec > th, "XAF1_mythelation_high", "XAF1_mythelation_low")

                    easier_df <- data.frame(easier_scores = easier_scores$easier_score, group = as.character(XAF1_factor2))
                    easier_df$thresholds <- split_method
                    easier_df$cancer <- c
                    easier_df$probe <- p
                    easier_res[[length(easier_res) + 1]] <- easier_df
                },
                silent = TRUE
            )
        } # end split methods
    } # end probe loop
}

easier_res_df <- bind_rows(easier_res)
easier_res_df <- easier_res_df |>
    na.omit()

head(easier_res_df)
unique(easier_res_df$group)
unique(easier_res_df$thresholds)
unique(easier_res_df$cancer)
unique(easier_res_df$probe)

library(effsize)
library(rstatix)
library(purrr)
library(tidyr)
# 假设 easier_res_df 已经存在
# 确保 group 是 factor，并指定顺序（high在前、low在后）
easier_res_df$group <- factor(easier_res_df$group,
    levels = c("XAF1_mythelation_high", "XAF1_mythelation_low")
)

# 诊断：查看每个组在每个子集的样本数量（便于定位为什么 NA）
count_diag <- easier_res_df %>%
    group_by(probe, cancer, thresholds, group) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = n, values_fill = 0)

# 查看有问题的组合（比如某一组为0）
problematic <- count_diag %>%
    filter(`XAF1_mythelation_high` == 0 | `XAF1_mythelation_low` == 0)
# 如果你想先看问题样本，可以 print(problematic) 或 write.csv(problematic, "problematic.csv")

# 稳健计算 cliff.delta 和 wilcox p-value：
all.score.df <- easier_res_df %>%
    group_by(probe, cancer, thresholds) %>%
    group_modify(~ {
        df <- .x
        # 计算每组样本数
        n_high <- sum(df$group == "XAF1_mythelation_high", na.rm = TRUE)
        n_low <- sum(df$group == "XAF1_mythelation_low", na.rm = TRUE)

        # 初始化返回行
        res <- tibble(
            cliff.delta = NA_real_,
            p.value = NA_real_,
            n_high = n_high,
            n_low = n_low
        )

        # 只有当两组都有样本时才计算
        if (n_high > 0 && n_low > 0) {
            # 确保 eager_scores 为数值
            df$easier_scores <- as.numeric(df$easier_scores)

            # 计算 cliff.delta（捕获可能的错误）
            cd <- tryCatch(
                {
                    effsize::cliff.delta(easier_scores ~ group, data = df)
                },
                error = function(e) NULL
            )

            if (!is.null(cd) && !is.na(cd$estimate)) {
                # cd$estimate 的类型可能是 named vector，取第一个数值
                res$cliff.delta <- as.numeric(cd$estimate)[1]
            }

            # 计算 wilcox p-value（捕获错误）
            w <- tryCatch(
                {
                    wilcox.test(easier_scores ~ group, data = df)
                },
                error = function(e) NULL
            )

            if (!is.null(w)) {
                res$p.value <- w$p.value
            }
        }

        res
    }, .keep = TRUE) %>%
    ungroup()

all.score.df <- all.score.df |>
    group_by(probe, cancer) |>
    slice_min(p.value) |>
    ungroup()

# 多重检验校正与显著性标记
all.score.df <- all.score.df %>%
    mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        signif = case_when(
            !is.na(p.adj) & p.adj < 0.05 ~ "adj.p < 0.05",
            TRUE ~ ""
        ),
        # 构建纵轴标签：probe + thresholds（根据偏好可改顺序）
        subtype = paste0(probe, "_", thresholds)
    )

library(ggplot2)
library(readr)

probe_order <- sort(unique(all.score.df$probe)) # 或自定义排序

all.score.df$probe <- factor(all.score.df$probe,
    levels = probe_order
)

filter_probe <- c("cg13823366", "cg10381036", "cg17753212")

all.score.df2 <- all.score.df %>%
    filter(!probe %in% filter_probe)


rho_tb <- read_delim(file.path(result.path, "XAF1_probe_expression_correlation.tsv"), delim = "\t", col_names = TRUE, show_col_types = FALSE)
rho_tb[1:5, 1:18]

all.score.df2 <- all.score.df2 %>%
    left_join(
        rho_tb %>%
            select(probe, cancer, rho, p_value) %>%
            dplyr::rename(rho.p.value = p_value),
        by = c("probe", "cancer")
    )
all.score.df2[1:5, 1:12]

library(ggstar)
p <- ggplot() +
    geom_tile(
        data = all.score.df2,
        aes(x = cancer, y = probe),
        fill = "#d9d9d9", color = "white", width = 0.65, height = 0.65
    ) +
    geom_tile(
        data = all.score.df2 %>% filter(!is.na(cliff.delta) & p.adj < 0.05),
        aes(x = cancer, y = probe, fill = cliff.delta),
        color = NA, width = 0.65, height = 0.65
    ) +
    geom_star(
        data = all.score.df2 %>% na.omit() %>% dplyr::filter(p.adj < 0.05),
        aes(x = cancer, y = probe, starshape = signif),
        size = 2, fill = "#d7301f", color = "#d7301f"
    ) +
    scale_fill_gradient2(
        low = "#984EA3", mid = "white", high = "#FF7F00", midpoint = 0,
        name = bquote(~Cliff ~ "'" ~ s ~ " Delta" ~ " (" ~ mythe^high / ~ mythe^low ~ ")")
    ) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "") +
    coord_equal(clip = "off") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face = "bold", color = "black", size = 9),
        axis.text.x = element_text(face = "bold", color = "black", size = 10, angle = 45, hjust = 0),
        legend.background = element_rect(color = "#737373")
    )

print(p)

pdf(file.path(result.path, "Fig7f_XAF1_methylation_EaSIeR_score.pdf"), height = 6, width = 10)
print(p)
dev.off()


library(ggplot2)
library(dplyr)
library(ggnewscale)

# 模拟数据
set.seed(1)
df <- all.score.df2

# 构造两个半格：上右半（cliff.delta）和下左半（rho）
half_tile <- function(data, side = "upper") {
    if (side == "upper") {
        data %>%
            mutate(
                x1 = as.numeric(factor(cancer)) - 0.5,
                x2 = as.numeric(factor(cancer)) + 0.5,
                y1 = as.numeric(factor(probe)) - 0.5,
                y2 = as.numeric(factor(probe)) + 0.5
            ) %>%
            rowwise() %>%
            mutate(coords = list(data.frame(
                x = c(x1, x2, x2),
                y = c(y2, y2, y1)
            ))) %>%
            unnest(coords)
    } else {
        data %>%
            mutate(
                x1 = as.numeric(factor(cancer)) - 0.5,
                x2 = as.numeric(factor(cancer)) + 0.5,
                y1 = as.numeric(factor(probe)) - 0.5,
                y2 = as.numeric(factor(probe)) + 0.5
            ) %>%
            rowwise() %>%
            mutate(coords = list(data.frame(
                x = c(x1, x1, x2),
                y = c(y1, y2, y1)
            ))) %>%
            unnest(coords)
    }
}

# 上右半：cliff.delta
upper <- half_tile(df, "upper") %>% mutate(value = cliff.delta, type = "Cliff Delta")

# 下左半：rho
lower <- half_tile(df, "lower") %>% mutate(value = rho, type = "rho")

# 合并
df_plot <- bind_rows(upper, lower) |>
    as.data.frame()

head(df_plot)

df_plot$value[which(df_plot$p.adj >= 0.05 & df_plot$type == "Cliff Delta")] <- NA
df_plot$value[which(df_plot$rho.p.value >= 0.05 & df_plot$type == "rho")] <- NA

center_df <- df_plot %>%
    distinct(cancer, probe, cliff.delta, rho, p.adj, rho.p.value, type) %>%
    mutate(
        x_center = as.numeric(factor(cancer)),
        y_center = as.numeric(factor(probe))
    )
head(center_df)

center_df <- center_df %>%
    mutate(
        signif = case_when(
            type == "Cliff Delta" & p.adj < 0.05 ~ "adj.p < 0.05",
            type == "rho" & rho.p.value < 0.05 ~ "adj.p < 0.05",
            TRUE ~ ""
        )
    )

# 绘图

library(ggplot2)
library(ggnewscale)
library(dplyr)
library(ggstar)
library(cowplot)

p <- ggplot() +
    # --- 第一层：Cliff Delta ---
    geom_polygon(
        data = df_plot %>% filter(type == "Cliff Delta"),
        aes(x = x, y = y, group = interaction(cancer, probe, type), fill = value),
        color = "white", size = 0.2
    ) +
    scale_fill_gradient2(
        low = "#984EA3", mid = "white", high = "#FF7F00", midpoint = 0,
        na.value = "#F0F0F0",
        name = bquote(~Cliff ~ "'" ~ s ~ " Delta" ~ " (" ~ mythe^high / ~ mythe^low ~ ")"),
        guide = guide_colorbar( # 🎯 独立设置 Cliff Delta 图例样式
            barwidth = 8,
            barheight = 0.6,
            title.position = "top",
            title.hjust = 0.5,
            title.theme = element_text(size = 9),
            label.theme = element_text(size = 8)
        )
    ) +
    new_scale_fill() + # 🔄 重置 fill，为 rho 启动新的颜色映射

    # --- 第二层：rho ---
    geom_polygon(
        data = df_plot %>% filter(type == "rho"),
        aes(x = x, y = y, group = interaction(cancer, probe, type), fill = value),
        color = "white", size = 0.2
    ) +
    scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
        na.value = "#F0F0F0",
        name = bquote(~rho ~ " (" ~ mythe / ~RNA ~ ")"),
        guide = guide_colorbar( # 🎯 独立设置 rho 图例样式
            barwidth = 4,
            barheight = 0.6,
            title.position = "top",
            title.hjust = 0.5,
            title.theme = element_text(size = 9),
            label.theme = element_text(size = 8)
        )
    ) +

    # --- 星号标记 ---
    geom_star(
        data = center_df %>% filter(p.adj < 0.05),
        aes(x = x_center + 0.2, y = y_center + 0.2),
        size = 2, color = "#D6604D", fill = "#D6604D"
    ) +
    geom_star(
        data = center_df %>% filter(rho.p.value < 0.05),
        aes(x = x_center - 0.2, y = y_center - 0.2),
        size = 2, color = "#2166AC", fill = "#2166AC"
    ) +

    # --- 坐标和主题 ---
    scale_x_continuous(
        breaks = 1:length(unique(df_plot$cancer)),
        labels = unique(df_plot$cancer),
        position = "top",
        expand = c(0, 0)
    ) +
    scale_y_continuous(
        breaks = 1:length(unique(df_plot$probe)),
        labels = unique(df_plot$probe),
        expand = c(0, 0)
    ) +
    coord_equal(clip = "off") +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face = "bold", color = "black", size = 15),
        axis.text.x = element_text(face = "bold", color = "black", size = 15, angle = 45, hjust = 0),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.x = unit(1, "cm") # 两个 legend 之间加点间距
    )

print(p)

pdf(file.path(result.path, "Fig7f_XAF1_methylation_EaSIeR_expression.pdf"), height = 7, width = 15)
print(p)
dev.off()
