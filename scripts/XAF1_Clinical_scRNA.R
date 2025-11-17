##############
##### SCP pipeline ###############
options(reticulate.conda_binary = "/home/bio/biosoft/miniconda3/bin/conda")
options(RETICULATE_PYTHON = "/home/bio/biosoft/miniconda3/envs/SCP_env/bin/python")

# ---- Packages: centralized at top of script ----
# 把所有 library() 放在脚本顶部，确保在脚本任何位置调用函数之前包已加载。
library(Seurat)
library(SCP)
library(reticulate)
library(readr)
library(harmony)
library(cowplot)
library(dplyr)
library(ggstar)
library(patchwork)
library(ggplot2)
library(GSVA)
library(tidyr)
library(pheatmap)
library(BiocParallel)
library(pbapply)
library(ggforce)
library(viridis)
library(scales)
# optional - keep commented if not needed
# library(msigdbr)
# -----------------------------------------------

### load data ##############
path_full <- "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/"

sce1 <- readRDS(paste(path_full, "data/IAV_and_COVID-19_data/GSE271284/230217_recast_after_qc.rds", sep = ""))
head(sce1@meta.data)
DefaultAssay(sce1) <- "RNA"

p1 <- CellDimPlot(
    srt = sce1, group.by = c(""), # stat.by = "Phase",
    reduction = "umap", theme_use = "theme_blank", label = T, label.size = 5, label_insitu = F, legend.position = "none", bg_color = NA
)
p1

p2 <- FeatureDimPlot(
    srt = sce1, features = c("CD19", "MS4A1", "CD3D", "CD3E", "CD8A", "CD4", "CD14", "LYZ", "NCAM1", "KLRD1"),
    reduction = "umap", theme_use = "theme_blank"
)
p2

# 1️⃣ 标准预处理
sce1.1 <- NormalizeData(sce1)
sce1.1 <- FindVariableFeatures(sce1.1, selection.method = "vst", nfeatures = 2000)
sce1.1 <- ScaleData(sce1.1, features = VariableFeatures(sce1.1))
sce1.1 <- RunPCA(sce1.1, features = VariableFeatures(sce1.1))

# 2️⃣ 使用 Harmony 整合
# 'group.by.vars' 是你想要去除批次效应的分组变量，这里是 "Patient.Id"
sce1.1 <- RunHarmony(
    object = sce1.1,
    "Patient.Id", sigma = 0.4, theta = 2, lambda = 1, plot_convergence = F
)

# 3️⃣ 下游分析使用 Harmony 校正后的结果
sce1.1 <- RunUMAP(sce1.1, reduction = "harmony", dims = 1:30)
sce1.1 <- FindNeighbors(sce1.1, reduction = "harmony", dims = 1:30)
sce1.1 <- FindClusters(sce1.1, resolution = 0.6)
DefaultAssay(sce)
View(sce1.1@meta.data)
p1 <- CellDimPlot(
    srt = sce1.1, group.by = c("seurat_clusters"), # stat.by = "Phase",
    reduction = "umap", theme_use = "theme_blank", label = T, label.size = 5, label_insitu = T, bg_color = NA
)
p1

p2 <- FeatureDimPlot(
    srt = sce1.1, features = c(
        "CD19", "MS4A1", "CD3D", "CD3E", "CD8A", "CD4",
        "CD14", "LYZ", "NCAM1", "KLRD1", "S100A8", "FCGR3A"
    ),
    reduction = "umap", theme_use = "theme_blank"
)
p2


sce1.1$condition <- factor(sce1.1$condition, levels = c(
    "healthy child",
    "no-pneumonia child",
    "asymptomatic child",
    "pneumonia child",
    "healthy adult",
    "mild adult",
    "severe adult",
    "IFNa severe adult"
))

sce1.1$condition_unified <- factor(sce1.1$condition_unified, levels = c(
    "healthy child",
    "asymptomatic child",
    "symptomatic child",
    "healthy adult",
    "mild adult",
    "severe adult"
))

indent2 <- as.character(sce1.1$seurat_clusters)
unique(indent2)

sce1.1@meta.data <- sce1@meta.data |>
    dplyr::mutate(cluster_anno = case_when(
        indent2 %in% c("0", "3", "10", "18", "21", "30", "31", "32", "36") ~ "B cell",
        indent2 %in% c("2", "6", "13", "17", "24", "28") ~ "CD4+ T",
        indent2 %in% c("4", "5", "7", "12", "16", "22", "25", "27") ~ "CD8+ T",
        indent2 %in% c("1", "8", "9", "15", "19", "20", "26", "34") ~ "Mono/Macro",
        TRUE ~ "Others"
    ))

sce1.1$condition_unified <- factor(sce1.1$condition_unified, levels = c(
    "healthy child",
    "asymptomatic child",
    "symptomatic child",
    "healthy adult",
    "mild adult",
    "severe adult"
))


p1 <- FeatureStatPlot(
    srt = sce1.1, group.by = "cluster_anno", , split.by = "condition_unified", bg.by = NULL, ncol = 2, plot_type = "box",
    stat.by = c("XAF1"), add_box = TRUE
)
p1

Mono_Macro <- subset(sce1.1, cells = rownames(sce1.1@meta.data)[which(sce1.1$cluster_anno %in% c("Mono/Macro"))])
CD8_T <- subset(sce1.1, cells = rownames(sce1.1@meta.data)[which(sce1.1$cluster_anno %in% c("CD8+ T"))])
CD4_T <- subset(sce1.1, cells = rownames(sce1.1@meta.data)[which(sce1.1$cluster_anno %in% c("CD4+ T"))])
B_cell <- subset(sce1.1, cells = rownames(sce1.1@meta.data)[which(sce1.1$cluster_anno %in% c("B cell"))])

unique(Mono_Macro$condition_unified)

p1 <- FeatureStatPlot(
    srt = Mono_Macro, group.by = "condition_unified", bg.by = NULL, ncol = 2, plot_type = "box",
    stat.by = c("XAF1"), add_box = TRUE, r, legend.position = "right", title = "XAF1 in Monocyte/Macropahge",
    comparisons = list(
        c("healthy child", "asymptomatic child"),
        c("asymptomatic child", "symptomatic child"),
        c("healthy adult", "mild adult"),
        c("mild adult", "severe adult")
    )
)
p1

p2 <- FeatureStatPlot(
    srt = CD8_T, group.by = "condition_unified", bg.by = NULL, ncol = 2, plot_type = "box",
    stat.by = c("XAF1"), add_box = TRUE, legend.position = "right", title = "XAF1 in CD8+ T",
    comparisons = list(
        c("healthy child", "asymptomatic child"),
        c("asymptomatic child", "symptomatic child"),
        c("healthy adult", "mild adult"),
        c("mild adult", "severe adult")
    )
)
p2

p3 <- FeatureStatPlot(
    srt = CD4_T, group.by = "condition_unified", bg.by = NULL, ncol = 2, plot_type = "box",
    stat.by = c("XAF1"), add_box = TRUE, legend.position = "right", title = "XAF1 in CD4+ T",
    comparisons = list(
        c("healthy child", "asymptomatic child"),
        c("asymptomatic child", "symptomatic child"),
        c("healthy adult", "mild adult"),
        c("mild adult", "severe adult")
    )
)
p3

p4 <- FeatureStatPlot(
    srt = B_cell, group.by = "condition_unified", bg.by = NULL, ncol = 2, plot_type = "box",
    stat.by = c("XAF1"), add_box = TRUE, legend.position = "right", title = "XAF1 in B cell",
    comparisons = list(
        c("healthy child", "asymptomatic child"),
        c("asymptomatic child", "symptomatic child"),
        c("healthy adult", "mild adult"),
        c("mild adult", "severe adult")
    )
)
p4

plot_grid(p1, p2, p3, p4, ncol = 2)

pdf(paste(path_full, "results/COVID-19_scRNA/XAF1 in cell types expression.pdf", sep = ""), height = 8, width = 12)
plot_grid(p1, p2, p3, p4, ncol = 2)
dev.off()

p1 <- CellDimPlot(
    srt = sce1.1, group.by = c("cluster_anno"), split.by = "condition_unified",
    reduction = "umap", theme_use = "theme_blank", label = T, label.size = 5, label_insitu = T, bg_color = NA
)
p1

sce1.2 <- subset(sce1.1, cells = rownames(sce1.1@meta.data)[which(sce1.1$cluster_anno %in% c("Mono/Macro", "CD8+ T", "CD4+ T", "B cell"))])


group.by <- "cluster_anno" # 细胞按什么分类

Reduction <- "umap" # 细胞用什么降维方法展示
population <- "condition_unified"


df <- data.frame(
    dim1 = sce1.2@reductions[[Reduction]]@cell.embeddings[, 1],
    dim2 = sce1.2@reductions[[Reduction]]@cell.embeddings[, 2],
    cluster = sce1.2@meta.data[[group.by]],
    population = sce1.2@meta.data[[population]],
    XAF1 = sce1.2@assays$RNA@data["XAF1", ]
)

setwd("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/COVID-19_scRNA/")
saveRDS(df, file = "df_for_cluster_umap.rds")

# 载入必要包
if (!requireNamespace("GSVA", quietly = TRUE)) BiocManager::install("GSVA")
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
# NOTE: package imports are centralized later in the script to avoid duplication

# ========== 参数（按需修改） ==========
seurat_obj <- sce1.2 # 替换为你脚本中的 Seurat 对象名
celltype_col <- "cluster_anno" # 细胞类型列名（替换为你实际列）
xaf1_gene <- "XAF1" # XAF1 的基因名（如果是 Ensembl，做映射）
n_top_celltypes <- 4 # 取样本数最多的 4 个细胞类型做比较
min_geneset_size <- 10 # 过滤太小的通路
kcdf_choice <- "Gaussian" # 对 log/normalized 表达用 Gaussian；若 counts 严重非正态可试 "Poisson"
# =====================================

# 1) 获取 KEGG 基因集（使用 msigdbr）或使用你已有的 prepare_DB 产物
# 说明：msigdbr 的 subcategory 为 "CP:KEGG"
kegg_df <- PrepareDB(species = "Homo_sapiens", db = c("KEGG"))
kegg_df <- kegg_df$Homo_sapiens$KEGG

head(kegg_df$TERM2GENE)
# kegg_df 包含 columns: gs_name, gene_symbol 等
kegg_list <- split(kegg_df$TERM2GENE$symbol, kegg_df$TERM2GENE$Term)
head(kegg_list)

# 2) 从 Seurat 对象取表达矩阵（genes x cells），使用 normalized data（slot="data"）
expr_mat <- as.matrix(GetAssayData(seurat_obj, slot = "data")) # genes (rows) x cells (cols)
# 确认行名为 SYMBOL
head(rownames(expr_mat))

# 3) 清洗/映射基因名，保证 gene sets 与 expr rownames 匹配（大写/去版本等）
# 例如去掉 ENSG 后缀: 如果你的行名是 ENSG... .1 需先处理
# 如果是 SYMBOL，通常无需更改。这里做一个安全检查并报告匹配率。
genes_in_expr <- rownames(expr_mat)
# 统计每个基因集的匹配数
match_counts <- sapply(kegg_list, function(gs) sum(gs %in% genes_in_expr))
summary(match_counts)
# 过滤掉匹配基因过少的集合
kegg_list_filt <- kegg_list[match_counts >= min_geneset_size]
length(kegg_list_filt) # 有多少通路被保留

# 4) 运行 ssGSEA（GSVA）得到 pathway scores（pathways x cells）
# 注意：这步可能较慢 / 占内存，建议在多线程或大内存环境下运行。
# 可选：设置 BiocParallel 以加速（例如 MulticoreParam）
# library(BiocParallel)
# detect cores or set manually
bpp <- BiocParallel::MulticoreParam(workers = 90) # 根据机器核数调整

ssgseapar <- ssgseaParam(expr_mat, kegg_list_filt)

ssgsea_results <- gsva(
    ssgseapar,
    BPPARAM = bpp
)
ssgsea_results[1:5, 1:5]
ssgsea_scores <- ssgsea_results # pathways (rows) x cells (cols)
# ssgsea_scores: matrix (npathways x ncells)

# 5) 把 XAF1 表达取出（同样按 cells）
# 确保 XAF1 在表达矩阵里
if (!(xaf1_gene %in% rownames(expr_mat))) stop("找不到 XAF1 基因，请确认基因名或先做 ID 映射。")
xaf1_vec <- expr_mat[xaf1_gene, ] # numeric vector named by cell

# 6) 选择四个主要细胞类型（样本数最多的前4）
ct_tab <- table(seurat_obj[[celltype_col]][, 1]) # 如果是 dataframe-like metadata
ct_tab_sorted <- sort(ct_tab, decreasing = TRUE)
major_cts <- names(ct_tab_sorted)[1:n_top_celltypes]
major_cts

# 7) 在每个细胞类型内对每条通路和 XAF1 做 Spearman 相关（逐通路）
## library(pbapply)
cells_by_ct <- lapply(major_cts, function(ct) {
    which(seurat_obj[[celltype_col]][, 1] == ct) # 索引（如果 metadata 是单列）
})
names(cells_by_ct) <- major_cts

# 结果容器
res_list <- list()

for (ct in major_cts) {
    idx_cells <- cells_by_ct[[ct]]
    cell_names <- colnames(expr_mat)[idx_cells]

    # 对每条通路做 cor.test
    res_ct <- pbapply::pbsapply(rownames(ssgsea_scores), function(pathway) {
        scores <- ssgsea_scores[pathway, cell_names]
        xaf <- as.numeric(xaf1_vec[cell_names])
        # remove NA
        ok <- !is.na(scores) & !is.na(xaf)
        if (sum(ok) < 4) {
            return(c(rho = NA, p = NA, n = sum(ok)))
        } else {
            ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
            c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
        }
    })
    res_ct <- t(res_ct)
    res_ct <- as.data.frame(res_ct, stringsAsFactors = FALSE)
    res_ct$adj.p <- p.adjust(res_ct$p, method = "BH")
    res_ct$pathway <- rownames(res_ct)
    res_ct$celltype <- ct
    res_list[[ct]] <- res_ct
}

res_all <- do.call(rbind, res_list)
colnames(res_all)[1] <- "rho"

head(res_all)

# 整理列类型
res_all <- res_all %>%
    mutate(
        rho = as.numeric(rho),
        p = as.numeric(p),
        adj.p = as.numeric(adj.p),
        n = as.integer(n)
    )

head(kegg_df$TERM2NAME)

res_all <- merge(res_all, kegg_df$TERM2NAME, by.x = "pathway", by.y = "Term", all.x = TRUE)
# 9) 导出结果表
write.csv(res_all, file = "XAF1_pathway_spearman_by_celltype.csv", row.names = FALSE)

res_all <- read.csv("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/COVID-19_scRNA/XAF1_pathway_spearman_by_celltype.csv",header = TRUE)

# 10) 判断细胞类型特异性 —— 两类方法（A）基于显著性交叉比较（显著/非显著分布）
# 简单筛选：在某细胞类型 q < 0.05 且在其它类型 q >= 0.05 -> 细胞类型特异性
res_sig <- res_all %>%
    filter(!is.na(adj.p)) %>%
    group_by(Name) %>%
    summarise(
        n_sig = sum(adj.p < 0.05),
        sig_ct = paste(celltype[adj.p < 0.05], collapse = ";"),
        .groups = "drop"
    ) %>%
    arrange(desc(n_sig))
# pathways where n_sig ==1 are potentially cell-type-specific
pathways_ct_specific <- res_sig %>% filter(n_sig <= 3)

# 11) 更严格的比较：两两比较相关系数差异（Fisher z）
# 函数：比较两组 rho 的差异
fisher_z_test <- function(r1, n1, r2, n2) {
    if (is.na(r1) | is.na(r2) | n1 < 4 | n2 < 4) {
        return(NA)
    }
    z1 <- atanh(r1)
    z2 <- atanh(r2)
    se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z <- (z1 - z2) / se
    p <- 2 * pnorm(-abs(z))
    return(p)
}

# 对每个 pathway 做 major_cts 两两比较
pairwise_results <- list()
for (pathway in rownames(ssgsea_scores)) {
    tabp <- res_all %>% filter(pathway == !!pathway & celltype %in% major_cts)
    # if not all present, next
    if (nrow(tabp) < 2) next
    combos <- combn(major_cts, 2)
    for (i in seq_len(ncol(combos))) {
        a <- combos[1, i]
        b <- combos[2, i]
        r1 <- tabp$rho[tabp$celltype == a]
        n1 <- tabp$n[tabp$celltype == a]
        r2 <- tabp$rho[tabp$celltype == b]
        n2 <- tabp$n[tabp$celltype == b]
        pval <- fisher_z_test(r1, n1, r2, n2)
        pairwise_results[[length(pairwise_results) + 1]] <- data.frame(pathway = pathway, Name = unique(tabp$Name), ct1 = a, ct2 = b, p_diff = pval, r1 = r1, r2 = r2, n1 = n1, n2 = n2)
    }
}
pairwise_df <- do.call(rbind, pairwise_results)
# 多重检验校正（按比较或按 pathway）
pairwise_df <- pairwise_df %>% mutate(q_diff = p.adjust(p_diff, method = "BH"))
head(pairwise_df)
write.csv(pairwise_df, file = "XAF1_pathway_pairwise_corr_diff.csv", row.names = FALSE)

# 12) 可视化：取 top N 显著相关的 pathway（按任一 celltype q）
topN <- 40
# 选基于任一 celltype 最小 q 排序
path_order <- res_all %>%
    group_by(Name) %>%
    summarise(min_q = min(adj.p, na.rm = TRUE)) %>%
    arrange(min_q) %>%
    pull(Name)
top_paths <- head(path_order, topN)

res_plot <- res_all %>% filter(Name %in% top_paths)

## Consolidated imports removed here (moved to top of script)


p <- ggplot() +
    geom_tile(
        data = res_plot,
        aes(x = celltype, y = Name),
        fill = "#d9d9d9", color = "white", width = 0.65, height = 0.65
    ) +
    geom_tile(
        data = res_plot %>% filter(adj.p < 0.05),
        aes(x = celltype, y = Name, fill = rho),
        color = NA, width = 0.65, height = 0.65
    ) +
    geom_star(
        data = res_plot %>% na.omit() %>% dplyr::filter(adj.p < 0.05),
        aes(x = celltype, y = Name),
        size = 1.5, fill = "#d7301f", color = "#d7301f"
    ) +
    scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
        name = bquote(~rho ~ "(" ~ "XAF1 & pathway activity" ~ ")")
    ) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "") +
    coord_equal(clip = "off") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face = "bold", color = "black", size = 11),
        axis.text.x = element_text(face = "bold", color = "black", size = 11, angle = 45, hjust = 0),
        legend.background = element_rect(color = "#737373")
    ) +
    plot_annotation(
        title = paste("TOP40 significant pathway", sep = ""),
        theme = theme(
            plot.title = element_text(
                size = 16, face = "bold", hjust = 0.5,
                margin = margin(b = 10)
            ),
            plot.margin = margin(20, 20, 20, 20)
        )
    )

pdf("XAF1_pathway_corr_heatmap_top40.pdf", width = 8, height = 12)
print(p)
dev.off()




### 切换 plot_star_cloud 环境 ####

df_all <- readRDS("df_for_cluster_umap.rds")
summary(df$XAF1)

remove_outlier <- function(x, y, proportion = 0.9) {
    stopifnot(length(x) == length(y))

    distance_to_centroid <- (x - mean(x))^2 + (y - mean(y))^2

    kept_site <- order(distance_to_centroid)[1:ceiling(length(x) * proportion)]

    return(kept_site)
}

unique(df$population)

for (i in levels(df_all$population)) {
    df <- df_all[df_all$population == i, , drop = F]
    p <- ggplot(
        df_all,
        aes(x = dim1, y = dim2)
    ) +
        stat_density_2d(
            geom = "raster",
            aes(fill = after_stat(density)),
            contour = FALSE, n = 200
        ) +
        viridis::scale_fill_viridis(option = "magma", direction = 1, alpha = 1) +
        coord_cartesian(expand = F) +
        theme_void() +
        geom_point(
            data = df |> dplyr::filter(XAF1 == 0),
            aes(x = dim1, y = dim2), color = "#FFFFFFAA",
            size = 0.005, alpha = 0.05
        ) +
        geom_point(
            data = df |> dplyr::filter(XAF1 > 0),
            aes(x = dim1, y = dim2, color = XAF1),
            size = 0.2, alpha =
            ) +
        scale_color_gradientn(
            colors = c("#FFFFFFAA", "green"), values = scales::rescale(c(0, 5)), # 对应表达量区间 limits = c(0, 5), # 限制范围（超出会截断） oob = scales::squish # 超出范围的值挤到边界 )
        ) +
        theme(legend.position = "none")

    df2 <- lapply(unique(df_all$cluster), function(x) {
        new_df <- df_all[df_all$cluster == x, , drop = F]

        return(new_df[remove_outlier(new_df$dim1, new_df$dim2), , drop = F])
    }) |> dplyr::bind_rows()

    p1 <- p + ggforce::geom_mark_hull(
        data = df2,
        aes(x = dim1, y = dim2, group = cluster, label = cluster),
        linewidth = 0.8,
        concavity = 5,
        color = "white",
        con.colour = "white",
        con.type = "straight",
        label.fontsize = 24,
        label.colour = "white",
        label.fill = "black",
        expand = unit(3, "mm"),
        linetype = "f8"
    ) +
        plot_annotation(
            title = i,
            theme = theme(
                plot.title = element_text(
                    size = 16, face = "bold", hjust = 0.5,
                    margin = margin(b = 10)
                ),
                plot.margin = margin(20, 20, 20, 20)
            )
        )
    pdf(paste(i, "COVID19_Major_Immune_Cell_Types_UMAP.pdf"), width = 10, height = 10)
    print(p)
    dev.off()

    # ----------------- 新增：基于 res_all 的汇总与可视化 -----------------
    # 1) 每个细胞类型中，与 XAF1 显著相关的 pathway 数量（adj.p < 0.05）
    sig_res <- res_all %>% filter(!is.na(adj.p) & adj.p < 0.05)

    counts_df <- sig_res %>% group_by(celltype) %>% summarise(n_sig = n()) %>% arrange(desc(n_sig))

    if (nrow(counts_df) > 0) {
        p_counts <- ggplot(counts_df, aes(x = reorder(celltype, -n_sig), y = n_sig)) +
            geom_col(fill = "#2b83ba") +
            geom_text(aes(label = n_sig), vjust = -0.3, size = 3.5) +
            labs(x = "Cell type", y = "Number of significant pathways (adj.p < 0.05)", title = "XAF1-significant KEGG pathways by cell type") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

        pdf("XAF1_sig_pathway_counts_by_celltype.pdf", width = 6, height = 4)
        print(p_counts)
        dev.off()
    } else {
        message("No significant pathways (adj.p < 0.05) found in res_all -- skipping counts plot.")
    }

    # 2) 四个细胞类型共同显著相关的 pathway（在所有 major_cts 中均 adj.p < 0.05）
    # 确保 major_cts 变量可用（脚本前面已定义），否则使用独立计算的前四类
    if (!exists("major_cts")) {
        major_cts <- unique(res_all$celltype)[1:4]
    }

    shared_names <- sig_res %>% filter(celltype %in% major_cts) %>% group_by(Name) %>% summarise(n_ct = n_distinct(celltype)) %>% filter(n_ct == length(major_cts)) %>% pull(Name)

    if (length(shared_names) > 0) {
        mat_shared <- res_all %>% filter(Name %in% shared_names & celltype %in% major_cts) %>% select(Name, celltype, rho) %>% tidyr::pivot_wider(names_from = celltype, values_from = rho)
        matm <- as.data.frame(mat_shared)
        rownames(matm) <- matm$Name
        matm$Name <- NULL
        matm_mat <- as.matrix(matm)

        # order pathways by mean rho
        ord <- order(rowMeans(matm_mat, na.rm = TRUE), decreasing = TRUE)
        matm_mat <- matm_mat[ord, , drop = FALSE]

        if (nrow(matm_mat) > 0) {
            pdf("XAF1_shared_pathways_heatmap_all4.pdf", width = 6 + 0.4 * nrow(matm_mat), height = 6 + 0.15 * nrow(matm_mat))
            pheatmap::pheatmap(matm_mat,
                color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(50),
                cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE,
                main = paste0("Pathways significant in all ", length(major_cts), " major cell types"),
                fontsize_row = 8
            )
            dev.off()
        }
    } else {
        message("No pathways are significant across all major cell types.")
    }

    # 3) 非全部共同显著（partial / cell-type-specific）pathways：dotplot（color = rho, size = -log10(adj.p) ）
    noncommon_tbl <- sig_res %>% filter(celltype %in% major_cts) %>% group_by(Name) %>% summarise(n_ct = n_distinct(celltype), min_q = min(adj.p, na.rm = TRUE)) %>% filter(n_ct < length(major_cts)) %>% arrange(n_ct, min_q)
    noncommon_names <- noncommon_tbl$Name

    if (length(noncommon_names) > 0) {
        topN_noncommon <- 40
        top_noncommon <- head(noncommon_names, topN_noncommon)

        df_dot <- res_all %>% filter(Name %in% top_noncommon & celltype %in% major_cts) %>%
            mutate(adj.p = as.numeric(adj.p)) %>%
            mutate(size_val = ifelse(!is.na(adj.p) & adj.p > 0, pmin(-log10(adj.p), 50), 0),
                   sig_flag = ifelse(!is.na(adj.p) & adj.p < 0.05, TRUE, FALSE))

        # keep order of pathways by number of cell types (more specific first)
        pathway_order <- noncommon_tbl %>% filter(Name %in% top_noncommon) %>% arrange(n_ct, min_q) %>% pull(Name)

        p_dot <- ggplot(df_dot, aes(x = celltype, y = factor(Name, levels = rev(pathway_order)))) +
            geom_point(aes(size = size_val, color = rho, alpha = sig_flag)) +
            scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0, na.value = "grey80") +
            scale_size_continuous(range = c(0.5, 6), name = expression(-log[10](adj.p))) +
            scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = FALSE) +
            labs(x = "Cell type", y = "Pathway (non-shared significant)", title = paste0("Top ", length(pathway_order), " pathways significant in subset of cell types")) +
            theme_bw() +
            theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1))

        pdf("XAF1_noncommon_pathways_dotplot_top40.pdf", width = 8, height = max(4, 0.2 * length(pathway_order)))
        print(p_dot)
        dev.off()
    } else {
        message("No partially-shared significant pathways to plot (all significant pathways are either unique or shared by all).")
    }

    # ----------------- 新增结束 -----------------

    ### 切换 plot_star_cloud 环境 ####
p <- ggplot(
    df_all,
    aes(x = dim1, y = dim2)
) +
    stat_density_2d(
        geom = "raster",
        aes(fill = after_stat(density)),
        contour = FALSE, n = 200
    ) +
    viridis::scale_fill_viridis(option = "magma", direction = 1, alpha = 1) +
    coord_cartesian(expand = F) +
    theme_void() +
    geom_point(
        data = df_all |> dplyr::filter(XAF1 == 0),
        aes(x = dim1, y = dim2), color = "#FFFFFFAA",
        size = 0.005, alpha = 0.05
    ) +
    geom_point(
        data = df_all |> dplyr::filter(XAF1 > 0),
        aes(x = dim1, y = dim2, color = XAF1),
        size = 0.005, alpha = 0.4
    ) +
    scale_color_gradientn(
        colors = c("#FFFFFFAA", "green"), values = scales::rescale(c(0, 5)), # 对应表达量区间 limits = c(0, 5), # 限制范围（超出会截断） oob = scales::squish # 超出范围的值挤到边界 )
    )

df2 <- lapply(unique(df_all$cluster), function(x) {
    new_df <- df_all[df_all$cluster == x, , drop = F]

    return(new_df[remove_outlier(new_df$dim1, new_df$dim2), , drop = F])
}) |> dplyr::bind_rows()

p1 <- p + ggforce::geom_mark_hull(
    data = df2,
    aes(x = dim1, y = dim2, group = cluster, label = cluster),
    linewidth = 0.8,
    concavity = 5,
    color = "white",
    con.colour = "white",
    con.type = "straight",
    label.fontsize = 24,
    label.colour = "white",
    label.fill = "black",
    expand = unit(3, "mm"),
    linetype = "f8"
) +
    plot_annotation(
        title = "UMAP of Major Immune Cell Types",
        theme = theme(
            plot.title = element_text(
                size = 16, face = "bold", hjust = 0.5,
                margin = margin(b = 10)
            ),
            plot.margin = margin(20, 20, 20, 20)
        )
    )

pdf("COVID19_Major_Immune_Cell_Types_UMAP.pdf", width = 10, height = 10)
print(p1)
dev.off()

png("COVID19_Major_Immune_Cell_Types_UMAP.png", width = 10, height = 10, units = "in", res = 300)
print(p1)
dev.off()
