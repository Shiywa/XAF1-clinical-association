library(BiocParallel)
library(tictoc)
library(GEOquery)
library(readr)
library(dplyr)
library(tibble)
library(tidyverse)
library(survival)
library(survminer)

res_path <- "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/ICB_real"
## load ICB data ###########

ex1 <- read_tsv("/home/Toshiba4/kidney/data/ICB/GSE91061/GSE91061_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")

ex1 <- as.data.frame(ex1)
ex1[1:5, 1:5]

anno1 <- read_tsv("/home/Toshiba4/kidney/data/ICB/GSE78220/Human.GRCh38.p13.annot.tsv.gz")
anno1 <- as.data.frame(anno1)
head(anno1)

clinical_da1 <- read.csv("/home/Toshiba4/kidney/data/ICB/GSE91061/bms038_clinical_data.csv", header = T)
head(clinical_da1)
table(clinical_da1$SampleType)
length(intersect(meta1$patient, clinical_da1$PatientID))
clinical_da1 <- clinical_da1[match(meta1$patient, clinical_da1$PatientID), ]
dim(clinical_da1)

exp1 <- ex1 |>
    inner_join(anno1[, c("GeneID", "Symbol")], by = "GeneID") %>%
    group_by(Symbol) %>%
    summarise(across(-GeneID, max, na.rm = TRUE)) %>%
    column_to_rownames(var = "Symbol")
exp1[1:5, 1:5]

saveRDS(exp1, file = "/home/Toshiba4/kidney/data/ICB/GSE91061/melanoma_bulk_expression_tpm.rds")
meta1 <- readRDS("/home/Toshiba4/kidney/data/ICB/GSE91061/melanoma_bulk_metadata.rds")

###

ex2 <- read_tsv("/home/Toshiba4/kidney/data/ICB/GSE78220/GSE78220_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
ex2 <- as.data.frame(ex2)
ex2[1:5, 1:5]

exp2 <- ex2 |>
    inner_join(anno1[, c("GeneID", "Symbol")], by = "GeneID") %>%
    group_by(Symbol) %>%
    summarise(across(-GeneID, \(x) max(x, na.rm = TRUE))) %>%
    column_to_rownames(var = "Symbol")

exp2[1:5, 1:5]
meta2 <- readRDS("/home/Toshiba4/kidney/data/ICB/GSE78220/melanoma_bulk_metadata.rds")

common_sample2 <- intersect(colnames(exp2), meta2$geo_accession)

exp2 <- exp2[, common_sample2]

saveRDS(exp2, file = "/home/Toshiba4/kidney/data/ICB/GSE78220/melanoma_bulk_expression_tpm.rds")
###

exp3_count <- readRDS("/home/Toshiba4/kidney/data/ICB/bladder/bladder_bulk_expression.rds")

exp3_count[1:5, 1:5]
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
    countData = exp3_count,
    colData = data.frame(
        condition = rep("A", ncol(exp3_count)),
        row.names = colnames(exp3_count)
    ),
    design = ~1 # 不做差异分析，只做归一化
)
vsd <- varianceStabilizingTransformation(dds)
exp3 <- assay(vsd)

saveRDS(exp3, file = "/home/Toshiba4/kidney/data/ICB/bladder/bladder_bulk_expression_vst.rds")

meta3 <- readRDS("/home/Toshiba4/kidney/data/ICB/bladder/bladder_bulk_metadata.rds")

### exp1 GSE91061 XAF1 生存分析 ##########

XAF1_exp_df <- data.frame(geo_accession = colnames(exp1), XAF1 = as.numeric(exp1["XAF1", ]))
head(XAF1_exp_df)

meta_clinial1 <- meta1 |>
    select("title", "geo_accession", patient) |>
    left_join(clinical_da1, by = c("patient" = "PatientID"))

surv_data <- meta_clinial1 %>%
    inner_join(XAF1_exp_df, by = "geo_accession")

head(surv_data)
# head(surv_data)
# 2. 构建生存对象
# OS 生存分析
surv_data$OS_event <- as.numeric(surv_data$OS_SOR) # 0 = alive, 1 = dead
surv_data$OS_time <- as.numeric(surv_data$OS) # 天数
surv_data$OS_group <- ifelse(surv_data$XAF1 >= median(surv_data$XAF1, na.rm = TRUE),
    "XAF1 high", "XAF1 low"
)

# 3. Cox 回归模型
# OS
cox_OS <- coxph(Surv(OS_time, OS_event) ~ OS_group, data = surv_data)
# summary(cox_OS)


# OS Kaplan-Meier曲线
fit_OS <- survfit(Surv(OS_time, OS_event) ~ OS_group, data = surv_data)

p1 <- ggsurvplot(fit_OS,
    pval = TRUE, conf.int = FALSE,
    risk.table = TRUE, # Add risk table
    risk.table.col = "strata", # Change risk table color by groups
    linetype = "strata", # Change line type by groups
    surv.median.line = "hv", # Specify median survival
    ggtheme = theme_bw(), # Change ggplot2 theme
    palette = c("#d94936", "#4fb0ca"),
    title = "Riaz et al.\n(Melanoma, n = 51)"
)

plotp <- p1$plot + theme_classic() + theme(plot.margin = margin(l = 90), legend.position = "top")
tablep <- p1$table + theme_classic() + theme(legend.position = "none") + scale_y_discrete(labels = c("Low", "High"))

p2 <- aplot::gglist(gglist = list(plotp, tablep), ncol = 1, heights = c(1.5, 0.3))
pdf(file.path(res_path, paste("1.XAF1_OS_KM_curve (GSE91061).pdf", sep = "")), height = 6, width = 6)
print(p2)
dev.off()

###  exp2 GSE78220 XAF1 生存分析 ##########

XAF1_exp_df2 <- data.frame(geo_accession = colnames(exp2), XAF1 = as.numeric(exp2["XAF1", ]))
head(XAF1_exp_df2)

surv_data <- meta2 %>%
    inner_join(XAF1_exp_df2, by = "geo_accession")
# 2. 构建生存对象
# OS 生存分析
surv_data$OS_event <- as.numeric(surv_data$OS_SOR) # 0 = alive, 1 = dead
surv_data$OS_time <- as.numeric(surv_data$OS) # 天数
surv_data$OS_group <- ifelse(surv_data$XAF1 >= median(surv_data$XAF1, na.rm = TRUE),
    "XAF1 high", "XAF1 low"
)

# OS Kaplan-Meier曲线
fit_OS <- survfit(Surv(OS_time, OS_event) ~ OS_group, data = surv_data)

p1 <- ggsurvplot(fit_OS,
    pval = TRUE, conf.int = FALSE,
    risk.table = TRUE, # Add risk table
    risk.table.col = "strata", # Change risk table color by groups
    linetype = "strata", # Change line type by groups
    surv.median.line = "hv", # Specify median survival
    ggtheme = theme_bw(), # Change ggplot2 theme
    palette = c("#d94936", "#4fb0ca"),
    title = "Hugo et al.\n(Melanoma, n = 26)"
)

plotp <- p1$plot + theme_classic() + theme(plot.margin = margin(l = 90), legend.position = "top")
tablep <- p1$table + theme_classic() + theme(legend.position = "none") + scale_y_discrete(labels = c("Low", "High"))
p2 <- aplot::gglist(gglist = list(plotp, tablep), ncol = 1, heights = c(1.5, 0.3))
pdf(file.path(res_path, paste("1.XAF_OS_KM_curve (GSE78220).pdf", sep = "")), height = 6, width = 6)
print(p2)
dev.off()

###  exp3 bladder XAF1 生存分析 ##########

XAF1_exp_df3 <- data.frame(geo_accession = colnames(exp3), XAF1 = as.numeric(exp3["XAF1", ]))
head(XAF1_exp_df3)

合并临床信息
surv_data <- meta3 %>%
    rownames_to_column(var = "geo_accession") %>%
    inner_join(XAF1_exp_df3, by = "geo_accession")
# 2. 构建生存对象
# OS 生存分析
surv_data$OS_event <- as.numeric(surv_data$censOS) # 0 = alive, 1 = dead
surv_data$OS_time <- as.numeric(surv_data$os) # 天数
surv_data$OS_group <- ifelse(surv_data$XAF1 >= median(surv_data$XAF1, na.rm = TRUE),
    "XAF1 high", "XAF1 low"
)

# OS Kaplan-Meier曲线
fit_OS <- survfit(Surv(OS_time, OS_event) ~ OS_group, data = surv_data)

p1 <- ggsurvplot(fit_OS,
    pval = TRUE, conf.int = FALSE,
    risk.table = TRUE, # Add risk table
    risk.table.col = "strata", # Change risk table color by groups
    linetype = "strata", # Change line type by groups
    surv.median.line = "hv", # Specify median survival
    ggtheme = theme_bw(), # Change ggplot2 theme
    palette = c("#d94936", "#4fb0ca"),
    title = "Mariathasan et al.\n(Bladder cancer, n = 348)"
)

plotp <- p1$plot + theme_classic() + theme(plot.margin = margin(l = 90), legend.position = "top")
tablep <- p1$table + theme_classic() + theme(legend.position = "none") + scale_y_discrete(labels = c("Low", "High"))

p2 <- aplot::gglist(gglist = list(plotp, tablep), ncol = 1, heights = c(1.5, 0.3))
pdf(file.path(res_path, paste("1.XAF1_OS_KM_curve (bladder).pdf", sep = "")), height = 6, width = 6)
print(p2)
dev.off()
