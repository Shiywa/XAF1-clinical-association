library(readr)
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)

### 450K probes ########

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno[1:5, 1:33]

anno_use <- anno[, c("chr", "pos", "strand", "Name", "Probe_rs", "Probe_maf", "CpG_rs", "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Enhancer", "Regulatory_Feature_Group", "DHS")]
anno_use[1:5, ]

anno_use <- anno_use |>
    as.data.frame() |>
    filter(grepl("XAF1", UCSC_RefGene_Name))

probes <- unique(as.character(anno_use$Name))

### GSE 168739  #######

beta1 <- read_csv("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE168739/GSE168739_betas.csv.gz", col_names = T)
beta1[1:5, 1:5]

gset1 <- getGEO("GSE168739", GSEMatrix = TRUE, getGPL = FALSE, destdir = "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE168739/", AnnotGPL = F)
meta1 <- pData(gset1[[1]])

beta1 <- beta1[, c(1, grep("EPICOVID", colnames(beta1)))]
which(beta1$ID_REF %in% probes)
dim(beta1)

summary(as.numeric(meta1$`age:ch1`))

meta1 <- meta1 %>%
    mutate(
        age = as.numeric(`age:ch1`),
        age_group = cut_number(age, n = 4),
        age_group2 = cut_number(age, n = 3)
    )

# 提取样本名交集
common_samples <- intersect(colnames(beta1), meta1$description)

beta_long <- beta1 %>%
    filter(ID_REF %in% probes) %>%
    pivot_longer(
        cols = -ID_REF,
        names_to = "sample",
        values_to = "beta"
    ) %>%
    left_join(
        meta1 %>%
            mutate(sample = meta1$description) %>%
            select(sample, age, age_group),
        by = "sample"
    )

library(ggplot2)

b.list <- list()

for (i in probes) {
    # 举例：画第一个探针
    probe_to_plot <- i

    unique(meta1$age_group)

    my_comparisons <- list(
        c("[19,35]", "(35,42]"),
        c("[19,35]", "(42,50]"),
        c("[19,35]", "(50,61]"),
        c("(35,42]", "(50,61]"),
        c("(35,42]", "(42,50]"),
        c("(42,50]", "(50,61]")
    )

    p <- ggboxplot(beta_long %>% filter(ID_REF == probe_to_plot),
        x = "age_group", y = "beta",
        color = "age_group", palette = "lancet",
        bxp.errorbar = T, bxp.errorbar.width = 0.1,
        add = "jitter", short.panel.labs = FALSE, nrow = 1
    ) +
        labs(
            title = probe_to_plot,
            x = "Age (years)",
            y = expression(beta)
        ) +
        # scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
        expand_limits(y = 0) +
        stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
        theme(
            axis.text.y = element_text(size = 15),
            axis.title.y = element_text(size = 15)
        ) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    b.list[[i]] <- p
}

setwd("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/methylation/age")

pdf("1. GSE168739 age.pdf", height = 12, width = 24)
plot_grid(plotlist = b.list, ncol = 5)
dev.off()

beta_long2 <- beta1 %>%
    filter(ID_REF %in% probes) %>%
    pivot_longer(
        cols = -ID_REF,
        names_to = "sample",
        values_to = "beta"
    ) %>%
    left_join(
        meta1 %>%
            mutate(sample = meta1$description) %>%
            select(sample, age, age_group2),
        by = "sample"
    )

library(ggplot2)

b.list2 <- list()

for (i in probes) {
    # 举例：画第一个探针
    probe_to_plot <- i

    unique(meta1$age_group2)

    my_comparisons <- list(
        c("[19,37]", "(37,47]"),
        c("[19,37]", "(47,61]"),
        c("(37,47]", "(47,61]")
    )

    p <- ggboxplot(beta_long2 %>% filter(ID_REF == probe_to_plot),
        x = "age_group2", y = "beta",
        color = "age_group2", palette = "lancet",
        bxp.errorbar = T, bxp.errorbar.width = 0.1,
        add = "jitter", short.panel.labs = FALSE, nrow = 1
    ) +
        labs(
            title = probe_to_plot,
            x = "Age (years)",
            y = expression(beta)
        ) +
        expand_limits(y = 0) +
        stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
        theme(
            axis.text.y = element_text(size = 15),
            axis.title.y = element_text(size = 15)
        ) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    b.list2[[i]] <- p
}

pdf("1. GSE168739 age 2.pdf", height = 12, width = 24)
plot_grid(plotlist = b.list2, ncol = 5)
dev.off()

age.list <- list()

for (i in probes) {
    probe_to_plot <- i
    df_probe <- beta_long %>%
        filter(ID_REF == probe_to_plot) %>%
        mutate(age = as.numeric(age))

    cor_res <- cor.test(df_probe$beta, df_probe$age, method = "spearman")

    cor_res$estimate # rho
    cor_res$p.value # p值

    p <- ggplot(df_probe, aes(x = age, y = beta)) +
        geom_point(color = "#1F77B4", size = 2, alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, color = "#E6550D", linetype = "solid") +
        labs(
            title = probe_to_plot,
            x = "Age (years)",
            y = expression(beta)
        ) +
        annotate(
            "text",
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 5, color = "black",
            label = paste0(
                "Spearman ρ = ",
                sprintf("%.2f", cor_res$estimate),
                "\nP = ", formatC(cor_res$p.value, format = "e", digits = 2)
            )
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_blank()
        )

    age.list[[i]] <- p
}

pdf("1. GSE168739 age cor.pdf", height = 12, width = 24)
plot_grid(plotlist = age.list, ncol = 5)
dev.off()

### GSE167202 ######

beta2 <- read_delim("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE167202/GSE167202_ProcessedBetaValues.txt.gz", col_names = T, delim = "\t")
beta2 <- as.data.frame(beta2)
dim(beta2)
beta2[1:10, 1:5]
beta2 <- beta2[, -c(2:22)]

anno2 <- read_delim("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE167202/GPL29753_platform.txt.gz", col_names = T, skip = 22, delim = "\t")
anno2[1:5, 1:15]
anno_use2 <- anno2 |>
    as.data.frame() |>
    filter(grepl("XAF1", Gene_NCBI)) |>
    filter(Probe_Source == "EPIC")

probes2 <- unique(anno_use2$ID)

gset2 <- getGEO("GSE167202", GSEMatrix = TRUE, getGPL = FALSE, destdir = "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE167202/", AnnotGPL = F)
meta2 <- pData(gset2[[1]])
table(meta2$`covid_status:ch1`)
table(meta2$`ed_disposition:ch1`)
dim(meta2)
table(meta2$`characteristics_ch1.29`)

meta2 <- meta2 %>%
    mutate(
        covid_status = as.character(`covid_status:ch1`),
        state = as.character(`ed_disposition:ch1`),
        severity = recode(`characteristics_ch1.29`,
            "severity score: 4" = "severity score: 3"
        )
    ) |>
    filter(covid_status != "other infection") |>
    mutate(severity = recode(severity,
        "severity score: NA" = "Healthy",
        "severity score: 1" = "S1",
        "severity score: 2" = "S2",
        "severity score: 3" = "S3"
    )) |>
    mutate(severity = factor(severity, levels = c("Healthy", "S1", "S2", "S3")))

table(meta2$state)

common_samples <- intersect(colnames(beta2), meta2$description.1)

beta_long2 <- beta2 %>%
    filter(ID_REF %in% probes2) %>%
    pivot_longer(
        cols = -ID_REF,
        names_to = "sample",
        values_to = "beta"
    ) %>%
    left_join(
        meta2 %>%
            mutate(sample = meta2$description.1) %>%
            select(sample, severity),
        by = "sample"
    )
table(beta_long2$ID_REF)
beta_long2 <- na.omit(beta_long2)
b.list3 <- list()

for (i in probes2) {
    # 举例：画第一个探针
    probe_to_plot <- i

    my_comparisons <- list(
        c("Healthy", "S1"),
        c("Healthy", "S2"),
        c("Healthy", "S3"),
        c("S1", "S2"),
        c("S1", "S3"),
        c("S2", "S3")
    )

    pd1 <- beta_long2 %>% filter(ID_REF == probe_to_plot)
    pd1 <- merge(pd1, anno_use2, by.x = "ID_REF", by.y = "ID")

    p <- ggboxplot(pd1,
        x = "severity", y = "beta",
        color = "severity", palette = "lancet",
        bxp.errorbar = T, bxp.errorbar.width = 0.1,
        add = "jitter", short.panel.labs = FALSE, nrow = 1
    ) +
        labs(
            title = paste(probe_to_plot, unique(pd1$Feat_NCBI)),
            x = "Severity",
            y = expression(beta)
        ) +
        # scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
        expand_limits(y = 0) +
        stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
        theme(
            axis.text.y = element_text(size = 15),
            axis.title.y = element_text(size = 15)
        ) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    b.list3[[i]] <- p
}

setwd("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/methylation/severity")

pdf("1. GSE167202 age.pdf", height = 16, width = 24)
plot_grid(plotlist = b.list3, ncol = 5)
dev.off()

### GSE246337 ######

beta3 <- read_csv("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE246337/GSE246337_betas.csv.gz", col_names = T)
beta3 <- as.data.frame(beta3)
dim(beta3)
colnames(beta3)[1] <- "ID_REF"
beta3[1:5, 1:5]

gset3 <- getGEO("GSE246337", GSEMatrix = TRUE, getGPL = FALSE, destdir = "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE246337/", AnnotGPL = F)
meta3 <- pData(gset3[[1]])
unique(meta3$description)

meta3 <- meta3 |>
    mutate(age = as.numeric(recode(`age:ch1`, ">89" = "89"))) |>
    mutate(age_group = cut_number(age, n = 4))

summary(as.numeric(meta3$age))
table(meta3$age_group)

anno3 <- read_delim("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/methylation/GSE246337/GPL33022-233683.txt", col_names = T, skip = 50)
anno3 <- as.data.frame(anno3)
anno3[1:5, 1:35]

anno_use3 <- anno3 |>
    as.data.frame() |>
    filter(grepl("XAF1", UCSC_RefGene_Name)) |>
    select("ID", "SPOT_ID", "Name", "UCSC_RefGene_Group", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island") |>
    rowwise() %>% # 按行处理
    mutate(UCSC_RefGene_Group = gsub("_1|_2|_3|_4|_5|_6|_7|_8", "", UCSC_RefGene_Group)) |>
    mutate(
        UCSC_RefGene_Name = paste(unique(unlist(strsplit(UCSC_RefGene_Name, ";"))), collapse = ";"),
        UCSC_RefGene_Group = paste(unique(unlist(strsplit(UCSC_RefGene_Group, ";"))), collapse = ";")
    ) %>%
    ungroup()

probes3 <- anno_use3$ID

common_samples <- intersect(colnames(beta3), meta3$description)

beta_long3 <- beta3 %>%
    filter(ID_REF %in% probes3) %>%
    pivot_longer(
        cols = -ID_REF,
        names_to = "sample",
        values_to = "beta"
    ) %>%
    left_join(
        meta3 %>%
            mutate(sample = meta3$description) %>%
            select(sample, age, age_group),
        by = "sample"
    )

b.list4 <- list()

for (i in probes3) {
    # 举例：画第一个探针
    probe_to_plot <- i

    unique(meta3$age_group)

    my_comparisons <- list(
        c("[18.1,34.2]", "(34.2,54]"),
        c("[18.1,34.2]", "(54,72]"),
        c("[18.1,34.2]", "(72,89]"),
        c("(34.2,54]", "(54,72]"),
        c("(34.2,54]", "(72,89]"),
        c("(54,72]", "(72,89]")
    )

    pd1 <- beta_long3 %>% filter(ID_REF == probe_to_plot)
    pd1 <- merge(pd1, anno_use3, by.x = "ID_REF", by.y = "ID")

    p <- ggboxplot(pd1,
        x = "age_group", y = "beta",
        color = "age_group", palette = "lancet",
        bxp.errorbar = T, bxp.errorbar.width = 0.1,
        add = "jitter", short.panel.labs = FALSE, nrow = 1
    ) +
        labs(
            title = paste(probe_to_plot, unique(pd1$UCSC_RefGene_Group)),
            x = "Age (years)",
            y = expression(beta)
        ) +
        # scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
        expand_limits(y = 0) +
        stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
        theme(
            axis.text.y = element_text(size = 15),
            axis.title.y = element_text(size = 15)
        ) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
        theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
    b.list4[[i]] <- p
}

setwd("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/methylation/age")

pdf("2. GSE246337 age.pdf", height = 12, width = 30)
plot_grid(plotlist = b.list4, ncol = 6)
dev.off()

age.list2 <- list()

for (i in probes3) {
    probe_to_plot <- i
    df_probe <- beta_long3 %>%
        filter(ID_REF == probe_to_plot) %>%
        mutate(age = as.numeric(age))

    cor_res <- cor.test(df_probe$beta, df_probe$age, method = "spearman")

    cor_res$estimate # rho
    cor_res$p.value # p值

    p <- ggplot(df_probe, aes(x = age, y = beta)) +
        geom_point(color = "#1F77B4", size = 2, alpha = 0.6) +
        geom_smooth(method = "lm", se = TRUE, color = "#E6550D", linetype = "solid") +
        labs(
            title = probe_to_plot,
            x = "Age (years)",
            y = expression(beta)
        ) +
        annotate(
            "text",
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 5, color = "black",
            label = paste0(
                "Spearman ρ = ",
                sprintf("%.2f", cor_res$estimate),
                "\nP = ", formatC(cor_res$p.value, format = "e", digits = 2)
            )
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_blank()
        )

    age.list2[[i]] <- p
}

pdf("2. GSE246337 age cor.pdf", height = 12, width = 24)
plot_grid(plotlist = age.list, ncol = 5)
dev.off()
