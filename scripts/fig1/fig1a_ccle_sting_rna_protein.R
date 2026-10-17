#!/usr/bin/env Rscript

# Fig. 1a: CCLE STING/TMEM173 RNA-protein discordance.

library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)

project_root <- Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd())
source(file.path(project_root, "scripts", "common", "config.R"))

da_path <- xaf1_env_path("XAF1_CCLE_DATA_DIR", xaf1_path("data", "CCLE"))
res_path <- xaf1_output_dir("fig1")
list.files(da_path)

protein_da <- read_csv(file.path(da_path, "protein_quant_current_normalized.csv.gz"), col_names = TRUE)
head(protein_da)
colnames(protein_da)
protein_da <- protein_da[, -c(7:48)]
protein_sample <- apply(as.data.frame(colnames(protein_da)[-c(1:6)]), 1, function(x) strsplit(x, "_")[[1]][1])
protein_sample_df <- data.frame(sample = colnames(protein_da)[-c(1:6)], CellLine = protein_sample)
protein_sample_df <- protein_sample_df[!duplicated(protein_sample_df$CellLine), ]
head(protein_sample_df)
protein_da <- protein_da[, c("Gene_Symbol", protein_sample_df$sample)]
dim(protein_da)
protein_da <- protein_da[!duplicated(protein_da$Gene_Symbol), ]

mRNA_da <- read_csv(file.path(da_path, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"), col_names = TRUE)
head(mRNA_da)
mRNA_da <- mRNA_da %>%
    filter(IsDefaultEntryForModel == "Yes")

mRNA_model_info <- read_csv(file.path(da_path, "Model.csv"), col_names = TRUE)
head(mRNA_model_info)
mRNA_model_info <- as.data.frame(mRNA_model_info)
mRNA_model_info <- mRNA_model_info |>
    filter(ModelID %in% mRNA_da$ModelID)

union_cell <- intersect(protein_sample_df$CellLine, mRNA_model_info$StrippedCellLineName)
# union_cell <- union_cell[union_cell != "A549"] # 去掉 A549 这个细胞系，因为它的蛋白质数据缺失较多

length(union_cell)
length(unique(union_cell))
sum(is.na(protein_da$Gene_Symbol) | protein_da$Gene_Symbol == "")
sum(duplicated(protein_da$Gene_Symbol))
protein_da_use <- protein_da |>
    filter(!is.na(Gene_Symbol) & Gene_Symbol != "") |>
    column_to_rownames("Gene_Symbol") |>
    select(all_of(protein_sample_df$sample[protein_sample_df$CellLine %in% union_cell]))

head(protein_da_use)

protein_sample_df <- protein_sample_df |> filter(CellLine %in% union_cell)

protein_da_use <- protein_da_use[, protein_sample_df$sample]
colnames(protein_da_use) <- protein_sample_df$CellLine

mRNA_da_use <- mRNA_da |>
    left_join(mRNA_model_info[, c("ModelID", "StrippedCellLineName")], by = "ModelID") |>
    as.data.frame() |>
    filter(StrippedCellLineName %in% union_cell) |>
    column_to_rownames("StrippedCellLineName") |>
    t()

mRNA_da_use <- mRNA_da_use[-c(1:6), ] # 去掉前6行非数据行
mRNA_da_use_num <- as.data.frame(
    apply(mRNA_da_use, 2, as.numeric),
    row.names = rownames(mRNA_da_use)
)
mRNA_da_use_num[1:5, 1:5]
dim(mRNA_da_use_num)
new.row <- apply(as.data.frame(rownames(mRNA_da_use_num)), 1, function(x) strsplit(x, " ")[[1]][2])
new.row <- gsub("\\(|\\)", "", new.row)

head(new.row)
entrez_ids <- new.row

symbol_map <- mapIds(
    org.Hs.eg.db,
    keys = entrez_ids,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
)

sum(duplicated(rownames(mRNA_da_use_num)))
symbol_map[symbol_map == "STING1"] <- "TMEM173"
symbol_map[is.na(symbol_map)] <- entrez_ids[is.na(symbol_map)]

rownames(mRNA_da_use_num) <- symbol_map

union_gene <- intersect(rownames(protein_da_use), rownames(mRNA_da_use_num))

length(union_gene)

protein_da_final <- protein_da_use[union_gene, union_cell]
mRNA_da_final <- mRNA_da_use_num[union_gene, union_cell]

dim(protein_da_final)
dim(mRNA_da_final)

protein_da_final[1:5, 1:5]
mRNA_da_final[1:5, 1:5]

which(union_gene == "TMEM173")

rna_vec <- as.numeric(mRNA_da_final["TMEM173", ])
prot_vec <- as.numeric(protein_da_final["TMEM173", ])

rna_z <- scale(rna_vec)
prot_z <- scale(prot_vec)

decouple_score <- rna_z - prot_z

result_df <- data.frame(
    CellLine = colnames(mRNA_da_final),
    RNA = rna_vec,
    Protein = prot_vec,
    RNA_z = as.numeric(rna_z),
    Protein_z = as.numeric(prot_z),
    Score = as.numeric(decouple_score)
)

result_df <- result_df[order(-result_df$Score), ]
result_df <- na.omit(result_df)
head(result_df, 10)


library(dplyr)
library(ggplot2)
library(ggrepel)

# 1) 选出两类极端点
top_hi_mrna_lo_prot <- result_df %>%
    arrange(desc(Score)) %>%
    slice_head(n = 20) %>%
    mutate(group = "High mRNA / Low protein (Top 20)")

top_lo_mrna_hi_prot <- result_df %>%
    arrange(Score) %>%
    slice_head(n = 20) %>%
    mutate(group = "Low mRNA / High protein (Top 20)")

# 2) 给全体点打标签分组（默认 Others）
plot_df <- result_df %>%
    mutate(group = "Others") %>%
    # 把两类极端点覆盖上去
    left_join(
        bind_rows(
            top_hi_mrna_lo_prot %>% dplyr::select(CellLine, group),
            top_lo_mrna_hi_prot %>% dplyr::select(CellLine, group)
        ),
        by = "CellLine",
        suffix = c("", ".new")
    ) %>%
    mutate(group = coalesce(group.new, group)) %>%
    dplyr::select(-group.new)

# 3) 绘图
p <- ggplot(plot_df, aes(x = RNA_z, y = Protein_z)) +
    geom_point(aes(color = group), alpha = 0.75, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_text_repel(
        data = plot_df %>% filter(group != "Others"),
        aes(label = CellLine, color = group),
        size = 3.2,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.25,
        segment.alpha = 0.5,
        show.legend = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(
        values = c(
            "Others" = "grey70",
            "High mRNA / Low protein (Top 20)" = "#d73027",
            "Low mRNA / High protein (Top 20)" = "#4575b4"
        )
    ) +
    theme_classic() +
    labs(x = "TMEM173 mRNA (z-score)", y = "TMEM173 protein (z-score)")

pdf(file.path(res_path, "TMEM173_mRNA_protein_scatter.pdf"),
    width = 9, height = 6
)
print(p)
dev.off()

head(plot_df)
write.table(plot_df, file.path(res_path, "TMEM173_mRNA_protein_scatter_data.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


library(dplyr)
library(ComplexHeatmap)
library(circlize)

# 1) 取要展示的细胞系（按你希望的顺序：先High mRNA/Low protein，再Low mRNA/High protein）
cells_hi <- top_hi_mrna_lo_prot$CellLine
cells_lo <- top_lo_mrna_hi_prot$CellLine
cells_all <- unique(c(cells_hi, cells_lo)) # 去重并保持顺序

# 2) 从 result_df 提取 z-score，构建 2 x N 矩阵
mat_df <- result_df %>%
    filter(CellLine %in% cells_all) %>%
    mutate(CellLine = factor(CellLine, levels = cells_all)) %>%
    arrange(CellLine) %>%
    dplyr::select(CellLine, RNA_z, Protein_z)

mat <- rbind(
    RNA = mat_df$RNA_z,
    Protein = mat_df$Protein_z
)
colnames(mat) <- as.character(mat_df$CellLine)

# 3) 给列做一个分组注释（属于哪一类极端点）
grp <- ifelse(colnames(mat) %in% cells_hi, "High mRNA / Low protein",
    ifelse(colnames(mat) %in% cells_lo, "Low mRNA / High protein", "Other")
)
grp <- factor(grp, levels = c("High mRNA / Low protein", "Low mRNA / High protein"))

ha <- rowAnnotation(
    Group = grp,
    col = list(
        Group = c(
            "High mRNA / Low protein" = "#d73027", # 红色
            "Low mRNA / High protein" = "#4575b4" # 蓝色
        )
    )
)

# 4) 颜色映射：z-score 常用 -2~2 截断更稳
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

cell_lines1 <- c(
    "LS513",
    "PANC0203",
    "NCIH196",
    "A704",
    "NCIH747",
    "KP2",
    "HSC4",
    "DETROIT562",
    "OV56",
    "HCC827",
    "JHOS2",
    "HCC2218",
    "KPL1",
    "KYM1",
    "HUH1",
    "NCIH446",
    "SUDHL4",
    "DMS114",
    "JHH7",
    "RKO",
    "THP1",
    "TF1",
    "ASPC1",
    "IALM",
    "WM115",
    "JHUEM2",
    "WM793",
    "NCIH2286",
    "NCIH647",
    "SKES1",
    "CAL27",
    "OVCAR8",
    "RH30",
    "HCC44",
    "SNU685",
    "MKN7",
    "NUGC3",
    "KMS12BM",
    "RH41",
    "JEKO1"
)

# 5) 画热图
ht1 <- Heatmap(
    t(mat[, cell_lines1]),
    name = "z-score",
    col = col_fun,
    left_annotation = ha,
    # top_annotation = ha,
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_row_dend = FALSE,
    column_names_rot = 45,
    row_names_side = "left",
    heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2))
)

pdf(file.path(res_path, "TMEM173_mRNA_protein_heatmap.pdf"), width = 6, height = 8)
draw(ht1)
dev.off()

############# top10 #################

# 1) 选出两类极端点
top_hi_mrna_lo_prot <- result_df %>%
    arrange(desc(Score)) %>%
    slice_head(n = 10) %>%
    mutate(group = "High mRNA / Low protein (Top 10)")

top_lo_mrna_hi_prot <- result_df %>%
    arrange(Score) %>%
    slice_head(n = 10) %>%
    mutate(group = "Low mRNA / High protein (Top 10)")

# 2) 给全体点打标签分组（默认 Others）
plot_df <- result_df %>%
    mutate(group = "Others") %>%
    # 把两类极端点覆盖上去
    left_join(
        bind_rows(
            top_hi_mrna_lo_prot %>% dplyr::select(CellLine, group),
            top_lo_mrna_hi_prot %>% dplyr::select(CellLine, group)
        ),
        by = "CellLine",
        suffix = c("", ".new")
    ) %>%
    mutate(group = coalesce(group.new, group)) %>%
    dplyr::select(-group.new)

# 3) 绘图
p <- ggplot(plot_df, aes(x = RNA_z, y = Protein_z)) +
    geom_point(aes(color = group), alpha = 0.75, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_text_repel(
        data = plot_df %>% filter(group != "Others"),
        aes(label = CellLine, color = group),
        size = 3.2,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.25,
        segment.alpha = 0.5,
        show.legend = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(
        values = c(
            "Others" = "grey70",
            "High mRNA / Low protein (Top 10)" = "#d73027",
            "Low mRNA / High protein (Top 10)" = "#4575b4"
        )
    ) +
    theme_classic() +
    labs(x = "TMEM173 mRNA (z-score)", y = "TMEM173 protein (z-score)")

pdf(file.path(res_path, "TMEM173_mRNA_protein_scatter(top10).pdf"),
    width = 9, height = 6
)
print(p)
dev.off()

head(plot_df)
write.table(plot_df, file.path(res_path, "TMEM173_mRNA_protein_scatter_data(top10).txt"), sep = "\t", row.names = FALSE, quote = FALSE)


library(dplyr)
library(ComplexHeatmap)
library(circlize)

# 1) 取要展示的细胞系（按你希望的顺序：先High mRNA/Low protein，再Low mRNA/High protein）
cells_hi <- top_hi_mrna_lo_prot$CellLine
cells_lo <- top_lo_mrna_hi_prot$CellLine
cells_all <- unique(c(cells_hi, cells_lo)) # 去重并保持顺序

# 2) 从 result_df 提取 z-score，构建 2 x N 矩阵
mat_df <- result_df %>%
    filter(CellLine %in% cells_all) %>%
    mutate(CellLine = factor(CellLine, levels = cells_all)) %>%
    arrange(CellLine) %>%
    dplyr::select(CellLine, RNA_z, Protein_z)

mat <- rbind(
    RNA = mat_df$RNA_z,
    Protein = mat_df$Protein_z
)
colnames(mat) <- as.character(mat_df$CellLine)

# 3) 给列做一个分组注释（属于哪一类极端点）
grp <- ifelse(colnames(mat) %in% cells_hi, "High mRNA / Low protein",
    ifelse(colnames(mat) %in% cells_lo, "Low mRNA / High protein", "Other")
)
grp <- factor(grp, levels = c("High mRNA / Low protein", "Low mRNA / High protein"))

ha <- rowAnnotation(
    Group = grp,
    col = list(
        Group = c(
            "High mRNA / Low protein" = "#d73027", # 红色
            "Low mRNA / High protein" = "#4575b4" # 蓝色
        )
    )
)

# 4) 颜色映射：z-score 常用 -2~2 截断更稳
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

cell_lines2 <- c(
    "LS513",
    "PANC0203",
    "NCIH747",
    "KP2",
    "DETROIT562",
    "HCC827",
    "KPL1",
    "KYM1",
    "NCIH446",
    "DMS114",
    "TF1",
    "ASPC1",
    "IALM",
    "WM115",
    "NCIH2286",
    "NCIH647",
    "SKES1",
    "SNU685",
    "MKN7",
    "RH41"
)

# 5) 画热图
ht1 <- Heatmap(
    t(mat[, cell_lines2]),
    name = "z-score",
    col = col_fun,
    left_annotation = ha,
    # top_annotation = ha,
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_row_dend = FALSE,
    column_names_rot = 45,
    row_names_side = "left",
    heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2))
)

pdf(file.path(res_path, "TMEM173_mRNA_protein_heatmap(top10).pdf"), width = 6, height = 8)
draw(ht1)
dev.off()





########

protein_mRNA_cor <- read_excel(file.path(da_path, "Table_S4_Protein_RNA_Correlation_and_Enrichments.xlsx"), sheet = 2)
head(protein_mRNA_cor)
dim(protein_mRNA_cor)

# =========================
# 1) 输入：你的相关性数据 df
# =========================
# df <- read.csv("rna_protein_spearman.csv", stringsAsFactors = FALSE)

# 如果你的表里还有一列标注方法（如 Method / Type），只保留 Spearman：
# df <- df %>% filter(Method == "Spearman")

# 统一列名（如不需要可删）
df2 <- protein_mRNA_cor %>%
    transmute(
        Gene  = as.character(`Gene Symbol`),
        Value = as.numeric(Spearman)
    ) %>%
    filter(!is.na(Value)) %>%
    # 限制范围到 [-1, 1]（如果你希望保留超界值可改为 pmin/pmax）
    filter(Value >= -1, Value <= 1)
head(df2)
# =========================
# 2) 构建 0.05 宽度的 bins
# =========================
bin_width <- 0.05
breaks <- seq(-1, 1, by = bin_width)

# 让每个值落入区间 [a, b)；include.lowest=TRUE 保证 -1 包含
df_bins <- df2 %>%
    mutate(
        bin = cut(Value, breaks = breaks, include.lowest = TRUE, right = FALSE)
    )

hist_df <- df_bins %>%
    count(bin, name = "n") %>%
    mutate(
        # 提取区间左右端点
        bin_left  = as.numeric(sub("^\\[([^,]+),.*$", "\\1", gsub("\\(|\\]", "[", as.character(bin)))),
        bin_right = bin_left + bin_width,
        bin_mid   = (bin_left + bin_right) / 2
    ) %>%
    arrange(bin_left)

# =========================
# 3) TMEM173 的定位：找到它的 spearman 值、所属 bin、该柱高度
# =========================
target_gene <- "TMEM173" # 你也可以换成 "STING1" 等
tmem_row <- df2 %>%
    filter(Gene == target_gene) %>%
    slice(1)

if (nrow(tmem_row) == 0) {
    warning(paste0("在 df 中没有找到基因：", target_gene, "。将不会标注。"))
} else {
    tmem_cor <- tmem_row$Value

    # 找 TMEM173 落在哪个 bin
    tmem_bin_left <- breaks[max(which(breaks <= tmem_cor))]
    # 若 tmem_cor == 1，max(which(breaks<=1)) 会得到最后一个break，需往前一个 bin
    if (isTRUE(all.equal(tmem_cor, 1))) tmem_bin_left <- 1 - bin_width

    tmem_mid <- tmem_bin_left + bin_width / 2

    tmem_count <- hist_df %>%
        filter(abs(bin_left - tmem_bin_left) < 1e-12) %>%
        pull(n)

    if (length(tmem_count) == 0) tmem_count <- 0
}

# =========================
# 4) 画图
# =========================
blue_nice <- "#2F80ED" # 好看的蓝色（你可改成 "#1F77B4" 等）

p <- ggplot(hist_df, aes(x = bin_mid, y = n)) +
    geom_col(width = bin_width * 0.95, fill = blue_nice, color = NA) +
    # 0.5 的虚线
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.6) +
    scale_x_continuous(
        limits = c(-1, 1),
        breaks = c(-1, -0.5, 0, 0.5, 1),
        labels = c("-1", "-0.5", "0", "0.5", "1")
    ) +
    labs(
        x = "Spearman correlation (mRNA vs protein)",
        y = "Number of genes"
    ) +
    theme_classic(base_size = 13) +
    theme(
        axis.line = element_line(linewidth = 0.7),
        axis.ticks = element_line(linewidth = 0.7),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 11)
    ) +
    xlim(-0.6, 1)

# =========================
# 5) 标注 TMEM173：牵引线起点 (x=实际cor, y=该柱高度)，指向右上方标签
# =========================
if (exists("tmem_cor") && !is.na(tmem_cor)) {
    # 标签放置位置：稍微往右、往上
    label_x <- min(1, tmem_cor + 0.03)
    label_y <- tmem_count + max(hist_df$n, na.rm = TRUE) * 0.03 + 1

    p <- p +
        # 牵引线：从 (tmem_cor, tmem_count) 指到 (label_x, label_y)
        geom_segment(
            aes(x = tmem_cor, y = tmem_count, xend = label_x, yend = label_y),
            inherit.aes = FALSE,
            linewidth = 0.6
        ) +
        # 标签文字
        annotate(
            "text",
            x = label_x, y = label_y,
            label = target_gene,
            hjust = 0, vjust = 0.5,
            size = 4
        ) +
        # 可选：在柱顶点加个小点，强调“这个柱”
        geom_point(
            aes(x = tmem_mid, y = tmem_count),
            inherit.aes = FALSE,
            size = 2
        )
}

print(p)

pdf(file.path(res_path, "protein_mRNA_correlation_histogram.pdf"), width = 6, height = 4)
print(p)
dev.off()
