setwd("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/IAV_and_COVID-19_data")
library(GEOquery)

library(GEOquery)
library(limma)
library(idmap2)
library(SCP)
# Annotation packages for Ensembl -> Symbol mapping
suppressMessages({
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    message("Package 'AnnotationDbi' not installed. Install if mapping via org.Hs.eg.db is desired.")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("Package 'org.Hs.eg.db' not installed. 'biomaRt' will be used as a fallback if available.")
  }
})

# Helper: convert Ensembl IDs to gene symbols. Removes version suffixes like .1 by default.
ensembl_to_symbol <- function(ids, remove_version = TRUE) {
  if (remove_version) {
    ids_clean <- sub("\\..*$", "", as.character(ids))
  } else {
    ids_clean <- as.character(ids)
  }

  # Try AnnotationDbi + org.Hs.eg.db first
  if (requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    suppressMessages(library(AnnotationDbi))
    suppressMessages(library(org.Hs.eg.db))
    # select may return multiple rows for one key; we'll match back to the input order
    res <- tryCatch(
      AnnotationDbi::select(org.Hs.eg.db, keys = unique(ids_clean), keytype = "ENSEMBL", columns = c("SYMBOL", "ENSEMBL")),
      error = function(e) NULL
    )
    if (!is.null(res) && nrow(res) > 0) {
      # ensure one-to-one mapping by taking the first SYMBOL for each ENSEMBL
      res <- res[!duplicated(res$ENSEMBL), , drop = FALSE]
      symbol_vec <- res$SYMBOL[match(ids_clean, res$ENSEMBL)]
      return(symbol_vec)
    }
  }

  # Fallback: biomaRt
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    suppressMessages(library(biomaRt))
    ensembl <- tryCatch(
      biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl"),
      error = function(e) NULL
    )
    if (!is.null(ensembl)) {
      res2 <- tryCatch(
        biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ids_clean, mart = ensembl),
        error = function(e) NULL
      )
      if (!is.null(res2) && nrow(res2) > 0) {
        res2 <- res2[!duplicated(res2$ensembl_gene_id), , drop = FALSE]
        symbol_vec <- res2$hgnc_symbol[match(ids_clean, res2$ensembl_gene_id)]
        return(symbol_vec)
      }
    }
  }

  # If all fails, return NA vector
  warning("Could not map Ensembl IDs to symbols: required packages or internet access for biomaRt may be missing. Returning NA for unmapped IDs.")
  return(rep(NA_character_, length(ids)))
}

# Additional packages used later in the script
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(readr)
library(ggpmisc)
library(patchwork)
library(idmap1)
library(idmap2)
library(tibble)
### 1. IAV GSE61821 #####

gset1 <- getGEO("GSE61821", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE61821", AnnotGPL = F)
meta1 <- pData(gset1[[1]])

ex1 <- exprs(gset1[[1]])
ex1 <- normalizeBetweenArrays(ex1)
boxplot(ex1, outline = FALSE, notch = T, col = factor(meta1$geo_accession), las = 2) ## 查看表达水平

index <- gset1[[1]]@annotation

ids <- getIDs("gpl10558")
head(ids)

ex1 <- as.data.frame(ex1)

exp1 <- ex1 |>
  dplyr::mutate(probe_id = rownames(ex1)) %>%
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

exp1[1:5, 1:5]
exp1_max <- exp1 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp1_max) <- exp1_max$symbol
exp1_max <- exp1_max[, -c(1:2)]

#### 1.1 all #####

used_feature <- c("XAF1", "IFNB1", "TMEM173", "CD8A")
plist <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(exp1_max[i, meta1$geo_accession, drop = T]),
    severity = meta1$`severity:ch1`,
    timepoint = meta1$`timepoint:ch1`
  ) |>
    mutate(timepoint = recode(timepoint, "day_0" = "Acute", "day_28" = "Follow Up")) |>
    mutate(Group = paste(severity, "\n(", timepoint, ")", sep = "")) |>
    mutate(
      timepoint = factor(timepoint, levels = c("Acute", "Follow Up")),
      severity = factor(severity, levels = c("OFI", "Mild", "Moderate", "Severe")),
      Group = factor(Group, levels = c("OFI\n(Acute)", "OFI\n(Follow Up)", "Mild\n(Acute)", "Mild\n(Follow Up)", "Moderate\n(Acute)", "Moderate\n(Follow Up)", "Severe\n(Acute)", "Severe\n(Follow Up)"))
    )

  group_n <- data3 %>%
    group_by(Group) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(Group, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$Group)

  my_comparisons <- list(
    c("OFI\n(Acute)", "Mild\n(Acute)"),
    c("Mild\n(Acute)", "Mild\n(Follow Up)"),
    c("Mild\n(Acute)", "Moderate\n(Acute)"),
    c("Moderate\n(Acute)", "Moderate\n(Follow Up)"),
    c("Moderate\n(Acute)", "Severe\n(Acute)")
  )
  p <- ggboxplot(data3,
    x = "Group", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    )
  p
  plist[[i]] <- p
}


pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/1. GSE61821 XAF1 and other gene expresion boxplot.pdf", height = 14, width = 16)
plot_grid(plotlist = plist, ncol = 2)
dev.off()

used_feature2 <- c("XAF1", "TMEM173", "CD8A")
plist2 <- list()
for (i in used_feature2) {
  data3 <- data.frame(
    exp = unlist(exp1_max[i, meta1$geo_accession, drop = T]),
    severity = meta1$`severity:ch1`,
    timepoint = meta1$`timepoint:ch1`
  ) |>
    mutate(timepoint = recode(timepoint, "day_0" = "Acute", "day_28" = "Follow Up")) |>
    mutate(Group = paste(severity, "\n(", timepoint, ")", sep = "")) |>
    mutate(
      timepoint = factor(timepoint, levels = c("Acute", "Follow Up")),
      severity = factor(severity, levels = c("OFI", "Mild", "Moderate", "Severe")),
      Group = factor(Group, levels = c("OFI\n(Acute)", "OFI\n(Follow Up)", "Mild\n(Acute)", "Mild\n(Follow Up)", "Moderate\n(Acute)", "Moderate\n(Follow Up)", "Severe\n(Acute)", "Severe\n(Follow Up)"))
    )

  data3 <- data3 |>
    filter(timepoint %in% c("Acute"))

  group_n <- data3 %>%
    group_by(Group) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(Group, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$Group)

  my_comparisons <- list(
    c("OFI\n(Acute)", "Mild\n(Acute)"),
    c("Mild\n(Acute)", "Moderate\n(Acute)"),
    c("Moderate\n(Acute)", "Severe\n(Acute)")
  )
  p <- ggboxplot(data3,
    x = "Group", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist2[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/1. GSE61821 XAF1 and other gene expresion boxplot2.pdf", height = 12, width = 10)
plot_grid(plotlist = plist2, ncol = 2)
dev.off()


used_feature2 <- c("IFNB1", "TMEM173", "CD8A")
plist2 <- list()
for (i in used_feature2) {
  data3 <- data.frame(
    exp = unlist(exp1_max["XAF1", meta1$geo_accession, drop = T]),
    exp2 = unlist(exp1_max[i, meta1$geo_accession, drop = T]),
    severity = meta1$`severity:ch1`,
    timepoint = meta1$`timepoint:ch1`
  ) |>
    mutate(timepoint = recode(timepoint, "day_0" = "Acute", "day_28" = "Follow Up")) |>
    mutate(Group = paste(severity, "\n(", timepoint, ")", sep = "")) |>
    mutate(
      timepoint = factor(timepoint, levels = c("Acute", "Follow Up")),
      severity = factor(severity, levels = c("OFI", "Mild", "Moderate", "Severe")),
      Group = factor(Group, levels = c("OFI\n(Acute)", "OFI\n(Follow Up)", "Mild\n(Acute)", "Mild\n(Follow Up)", "Moderate\n(Acute)", "Moderate\n(Follow Up)", "Severe\n(Acute)", "Severe\n(Follow Up)"))
    )

  # data3 <- data3 |>
  #   filter(severity %in% c("Mild","Moderate") & timepoint %in% c("Acute"))

  p5 <- ggplot(data3, aes(x = exp, y = exp2)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    stat_poly_eq(
      aes(
        label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~") # 加入 p.value
      ),
      formula = y ~ x,
      parse = TRUE,
      output.type = "expression", # 以可解析的数学公式输出
      label.x.npc = "left",
      label.y.npc = 0.95
    ) +
    xlab("XAF1 expression") +
    ylab(paste(i, "expression", sep = " ")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = "black")
    )
  plist2[[i]] <- p5
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/1. GSE61821 correlation of XAF1 and other genes.pdf", height = 12, width = 12)
plot_grid(plotlist = plist2, ncol = 2)
dev.off()

### 2. IAV GSE90732 ####

gset2 <- getGEO("GSE90732", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE90732/", AnnotGPL = F)

ex2 <- exprs(gset2[[1]])
ex2 <- normalizeBetweenArrays(ex2)
boxplot(ex2, outline = FALSE, notch = T, col = factor(meta1$geo_accession), las = 2) ## 查看表达水平

index <- gset2[[1]]@annotation

ids <- getIDs("gpl10558")
head(ids)

ex2 <- as.data.frame(ex2)

exp2 <- ex2 |>
  dplyr::mutate(probe_id = rownames(ex2)) %>%
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

exp2[1:5, 1:5]
exp2_max <- exp2 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp2_max) <- exp2_max$symbol
exp2_max <- exp2_max[, -c(1:2)]

####

meta2 <- pData(gset2[[1]])
length(unique(meta2$`subject id:ch1`))

meta2 <- meta2 |>
  arrange(`subject id:ch1`, `time point:ch1`)

meta2 <- meta2[-grep("Replicate 2", meta2$title), ]

exp2_max <- exp2_max[, meta2$geo_accession]

used_feature <- c("XAF1", "IFNB1", "TMEM173", "CD8A")
plist3 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(exp2_max[i, meta2$geo_accession, drop = T]),
    sample = meta2$`subject id:ch1`,
    timepoint = meta2$`time point:ch1`
  ) |>
    mutate(timepoint = factor(timepoint, levels = c("baseline", "24 hrs", "48 hrs", "72 hrs", "96 hrs")))

  group_n <- data3 %>%
    group_by(timepoint) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(timepoint, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$timepoint)

  my_comparisons <- list(
    c("baseline", "24 hrs"),
    c("24 hrs", "48 hrs"),
    c("48 hrs", "72 hrs"),
    c("72 hrs", "96 hrs")
  )
  p <- ggboxplot(data3,
    x = "timepoint", y = "exp",
    color = "timepoint", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", paired = T) +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

  p
  plist3[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/2. GSE90732 XAF1 and other gene expresion boxplot.pdf", height = 14, width = 12)
plot_grid(plotlist = plist3, ncol = 2)
dev.off()

###

used_feature2 <- c("IFNB1", "TMEM173", "CD8A")
plist4 <- list()
for (i in used_feature2) {
  data3 <- data.frame(
    exp = unlist(exp2_max["XAF1", meta2$geo_accession, drop = T]),
    exp2 = unlist(exp2_max[i, meta2$geo_accession, drop = T]),
    sample = meta2$`subject id:ch1`,
    timepoint = meta2$`time point:ch1`
  ) |>
    mutate(timepoint = factor(timepoint, levels = c("baseline", "24 hrs", "48 hrs", "72 hrs", "96 hrs")))

  p5 <- ggplot(data3, aes(x = exp, y = exp2)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    stat_poly_eq(
      aes(
        label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~") # 加入 p.value
      ),
      formula = y ~ x,
      parse = TRUE,
      output.type = "expression", # 以可解析的数学公式输出
      label.x.npc = "left",
      label.y.npc = 0.95
    ) +
    xlab("XAF1 expression") +
    ylab(paste(i, "expression", sep = " ")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = "black")
    )
  plist4[[i]] <- p5
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/2. GSE90732 correlation of XAF1 and other genes.pdf", height = 12, width = 12)
plot_grid(plotlist = plist4, ncol = 2)
dev.off()

### 3. IAV GSE101702 #####

gset3 <- getGEO("GSE101702", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE101702/", AnnotGPL = F)

ex3 <- exprs(gset3[[1]])
ex3 <- normalizeBetweenArrays(ex3)
boxplot(ex3, outline = FALSE, notch = T, col = factor(meta1$geo_accession), las = 2) ## 查看表达水平

index <- gset3[[1]]@annotation

gpl <- read_delim("./GSE101702/GPL21185-21174.txt", delim = "\t", col_names = T, skip = 18)

ids <- gpl[, c(1, 6)]
colnames(ids) <- c("probe_id", "symbol")

ids <- as.data.frame(ids) |>
  na.omit()

ex3 <- as.data.frame(ex3)

exp3 <- ex3 |>
  dplyr::mutate(probe_id = rownames(ex3)) %>%
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

exp3[1:5, 1:5]
exp3_max <- exp3 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp3_max) <- exp3_max$symbol
exp3_max <- exp3_max[, -c(1:2)]

meta3 <- pData(gset3[[1]])
length(unique(meta2$`subject id:ch1`))

meta3 <- meta3 |>
  mutate(severity = recode(`severity:ch1`, "hlty_ctrl" = "Ctrl", "flu_mod" = "Moderate", "flu_svre" = "Severe")) |>
  mutate(severity = factor(severity, levels = c("Ctrl", "Moderate", "Severe")))

exp3_max <- exp3_max[, meta3$geo_accession]

plist5 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(exp3_max[i, meta3$geo_accession, drop = T]),
    severity = meta3$severity
  ) |>
    mutate(severity = factor(severity, levels = c("Ctrl", "Moderate", "Severe")))

  group_n <- data3 %>%
    group_by(severity) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(severity, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$severity)

  my_comparisons <- list(
    c("Ctrl", "Moderate"),
    c("Moderate", "Severe")
  )
  p <- ggboxplot(data3,
    x = "severity", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

  p
  plist5[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/3. GSE101702 XAF1 and other gene expresion boxplot.pdf", height = 10, width = 9)
plot_grid(plotlist = plist5, ncol = 2)
dev.off()

####

plist6 <- list()
for (i in used_feature2) {
  data3 <- data.frame(
    exp = unlist(exp3_max["XAF1", meta3$geo_accession, drop = T]),
    exp2 = unlist(exp3_max[i, meta3$geo_accession, drop = T])
  )

  p5 <- ggplot(data3, aes(x = exp, y = exp2)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    stat_poly_eq(
      aes(
        label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~") # 加入 p.value
      ),
      formula = y ~ x,
      parse = TRUE,
      output.type = "expression", # 以可解析的数学公式输出
      label.x.npc = "left",
      label.y.npc = 0.95
    ) +
    xlab("XAF1 expression") +
    ylab(paste(i, "expression", sep = " ")) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = "black")
    )
  plist6[[i]] <- p5
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/3. GSE101702 correlation of XAF1 and other genes.pdf", height = 12, width = 12)
plot_grid(plotlist = plist6, ncol = 2)
dev.off()

### 4. IAV GSE111368 #####

gset4 <- getGEO("GSE111368", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE111368/", AnnotGPL = F)

ex4 <- exprs(gset4[[1]])
ex4 <- normalizeBetweenArrays(ex4)
boxplot(ex4, outline = FALSE, notch = T, col = factor(meta1$geo_accession), las = 2) ## 查看表达水平

index <- gset4[[1]]@annotation

ids <- getIDs("gpl10558")
head(ids)

ex4 <- as.data.frame(ex4)

exp4 <- ex4 |>
  dplyr::mutate(probe_id = rownames(ex4)) %>%
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

exp4[1:5, 1:5]
exp4_max <- exp4 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp4_max) <- exp4_max$symbol
exp4_max <- exp4_max[, -c(1:2)]

####

meta4 <- pData(gset4[[1]])
length(unique(meta2$`subject id:ch1`))
unique(meta4$`timepoint:ch1`)

exp4_max <- exp4_max[, meta4$geo_accession]

used_feature <- c("XAF1", "TMEM173", "CD8A")
plist7 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(exp4_max[i, meta4$geo_accession, drop = T]),
    severity = meta4$`t1severity:ch1`,
    timepoint = meta4$`timepoint:ch1`
  ) |>
    mutate(Group = paste(severity, "\n(", timepoint, ")", sep = "")) |>
    mutate(
      timepoint = factor(timepoint, levels = c("HC", "T1", "T2", "T3")),
      severity = factor(severity, levels = c("HC", "1", "2", "3")),
      Group = factor(Group, levels = c("HC\n(HC)", "1\n(T1)", "1\n(T2)", "1\n(T3)", "2\n(T1)", "2\n(T2)", "2\n(T3)", "3\n(T1)", "3\n(T2)", "3\n(T3)"))
    )

  group_n <- data3 %>%
    group_by(Group) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(Group, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$Group)

  my_comparisons <- list(
    c("HC\n(HC)", "1\n(T1)"),
    c("1\n(T1)", "2\n(T1)"),
    c("2\n(T1)", "3\n(T1)")
  )
  p <- ggboxplot(data3,
    x = "Group", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist7[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/4. GSE111368 XAF1 and other gene expresion boxplot.pdf", height = 14, width = 13)
plot_grid(plotlist = plist7, ncol = 2)
dev.off()

plist8 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(exp4_max[i, meta4$geo_accession, drop = T]),
    severity = meta4$`t1severity:ch1`,
    timepoint = meta4$`timepoint:ch1`
  ) |>
    mutate(Group = paste(severity, "\n(", timepoint, ")", sep = "")) |>
    mutate(
      timepoint = factor(timepoint, levels = c("HC", "T1", "T2", "T3")),
      severity = factor(severity, levels = c("HC", "1", "2", "3")),
      Group = factor(Group, levels = c("HC\n(HC)", "1\n(T1)", "1\n(T2)", "1\n(T3)", "2\n(T1)", "2\n(T2)", "2\n(T3)", "3\n(T1)", "3\n(T2)", "3\n(T3)"))
    )

  data3 <- data3 |>
    filter(timepoint %in% c("HC", "T1"))

  group_n <- data3 %>%
    group_by(Group) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(Group, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$Group)

  my_comparisons <- list(
    c("HC\n(HC)", "1\n(T1)"),
    c("1\n(T1)", "2\n(T1)"),
    c("2\n(T1)", "3\n(T1)")
  )
  p <- ggboxplot(data3,
    x = "Group", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist8[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/4. GSE111368 XAF1 and other gene expresion boxplot2.pdf", height = 10, width = 9)
plot_grid(plotlist = plist8, ncol = 2)
dev.off()

### 5. COVID-19 GSE202805 #####

gset5 <- getGEO("GSE202805", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE202805/", AnnotGPL = F)

ex5 <- read_delim("./GSE202805/GSE202805_genecounts_20220406.txt", delim = "\t", col_names = T)
ex5 <- as.data.frame(ex5)
head(ex5)
ex5[1:5, 1:8]
dim(ex5)

ex5$new_id <- apply(ex5, 1, function(x) strsplit(x["Geneid"], split = "\\.")[[1]][1])
head(ex5$new_id)

# Map Ensembl IDs (without version) to gene symbols
ex5$symbol <- ensembl_to_symbol(ex5$new_id)

ex5$symbol[which(is.na(ex5$symbol))] <- ex5$new_id[which(is.na(ex5$symbol))]

# Build expression matrix: keep Length, select symbol and sample columns (those starting with 'R')
# and aggregate counts by symbol while summarising a representative gene length (median).
exp5 <- ex5 |>
  dplyr::select(symbol, Length, starts_with("R")) |>
  # ensure numeric for Length and sample columns
  dplyr::mutate(Length = as.numeric(Length), dplyr::across(starts_with("R"), ~ as.numeric(.x))) |>
  # aggregate identical symbols: median Length and sum counts
  dplyr::group_by(symbol) %>%
  dplyr::summarise(Length = median(Length, na.rm = TRUE), dplyr::across(starts_with("R"), ~ sum(.x, na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  as.data.frame()

# Remove genes with all-zero expression across all samples (exclude Length column)
if (nrow(exp5) > 0) {
  sample_cols <- grep("^R", colnames(exp5), value = TRUE)
  expr_vals <- exp5[, sample_cols, drop = FALSE]
  nonzero_mask <- rowSums(expr_vals, na.rm = TRUE) > 0
  exp5 <- exp5[nonzero_mask, , drop = FALSE]
}

# Set rownames to symbol and keep Length column for TPM calculation
if (nrow(exp5) > 0) {
  rownames(exp5) <- exp5$symbol
  # keep Length column alongside sample columns for now
  # drop the symbol column
  exp5 <- exp5[, !(colnames(exp5) %in% c("symbol")), drop = FALSE]
}

# Clean sample column names (remove .bam suffix) but leave Length intact
if (nrow(exp5) > 0) {
  sample_cols <- setdiff(colnames(exp5), "Length")
  colnames(exp5)[colnames(exp5) %in% sample_cols] <- gsub("\\.bam$", "", sample_cols)
}

meta5 <- pData(gset5[[1]])
meta5 <- meta5[, c("title", "disease severity:ch1", "group:ch1")]
colnames(meta5) <- c("Sample", "severity", "group")
meta5$severity[which(meta5$group == "convalescent")] <- "convalescent"
meta5$severity[which(meta5$severity == "severe_and_icu")] <- "severe"
meta5$severity <- factor(meta5$severity, levels = c("healthy", "mild", "moderate", "severe", "convalescent"))

# Reorder sample columns to match meta5; ensure Length column is excluded when subsetting
if (nrow(exp5) > 0) {
  samples_in_meta <- meta5$Sample
  # keep only samples present both in exp5 and meta5
  samples_present <- intersect(samples_in_meta, colnames(exp5))
  if (length(samples_present) == 0) {
    stop("No matching samples between exp5 columns and meta5$Sample")
  }
  counts_mat <- as.matrix(exp5[, samples_present, drop = FALSE])

  # compute TPM
  # lengths in kilobases
  lengths_kb <- exp5$Length / 1000
  # avoid division by zero lengths
  lengths_kb[lengths_kb == 0 | is.na(lengths_kb)] <- NA
  # RPK: counts / length_kb (result keeps NA where length unknown)
  rpk <- sweep(counts_mat, 1, lengths_kb, FUN = "/")
  # column sums of RPK
  rpk_sum <- colSums(rpk, na.rm = TRUE)
  # TPM: rpk / (rpk_sum / 1e6)
  tpm <- sweep(rpk, 2, rpk_sum / 1e6, FUN = "/")
  # replace NaN/Inf produced by zero division with 0
  tpm[!is.finite(tpm)] <- 0

  exp5_tpm <- as.data.frame(tpm)
  # ensure rownames preserved
  rownames(exp5_tpm) <- rownames(exp5)
}

used_feature <- c("XAF1", "STING1", "CD8A")
plist9 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(exp5_tpm[i, meta5$Sample, drop = T]),
    severity = meta5$severity
  ) |>
    mutate(severity = factor(severity, levels = c("healthy", "mild", "moderate", "severe", "convalescent")))

  group_n <- data3 %>%
    group_by(severity) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(severity, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$severity)

  my_comparisons <- list(
    c("healthy", "mild"),
    c("mild", "moderate"),
    c("moderate", "severe"),
    c("severe", "convalescent")
  )
  p <- ggboxplot(data3,
    x = "severity", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist9[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/5. GSE202805 XAF1 and other gene expresion boxplot.pdf", height = 10, width = 10)
plot_grid(plotlist = plist9, ncol = 2)
dev.off()

### 6. COVID-19 GSE227080 #####
gset6 <- getGEO("GSE227080", GSEMatrix = TRUE, getGPL = FALSE, destdir = "/home/Toshiba1/Lab_Ma/yuanyi/XAF1_clinical/GSE227080/", AnnotGPL = F)

meta6 <- pData(gset6[[1]])
meta6 <- meta6[, c("geo_accession", "severity:ch1")]
colnames(meta6) <- c("Sample", "severity")
meta6 <- meta6 |>
  mutate(severity = recode(severity, "NEG" = "healthy", "MILD" = "mild", "MOD_SEV" = "moderate_severe")) |>
  mutate(severity = factor(severity, levels = c("healthy", "mild", "moderate_severe")))

ex6 <- exprs(gset6[[1]])

used_feature <- c("TMEM173", "CD8A")
plist10 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex6[i, meta6$Sample, drop = T]),
    severity = meta6$severity
  )

  group_n <- data3 %>%
    group_by(severity) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(severity, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$severity)

  my_comparisons <- list(
    c("healthy", "mild"),
    c("mild", "moderate_severe")
  )
  p <- ggboxplot(data3,
    x = "severity", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist10[[i]] <- p
}
pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/6. GSE227080 TMEM173 and CD8A expresion boxplot.pdf", height = 4, width = 8)
plot_grid(plotlist = plist10, ncol = 2)
dev.off()

### 7. COVID-19 GSE215262 #####

gset7 <- getGEO("GSE215262", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE215262/", AnnotGPL = F)

meta7 <- pData(gset7[[1]])

ex7 <- read.table("./GSE215262/GSE215262_All.counts.txt", header = T, row.names = 1, sep = "\t", check.names = F)

meta7 <- data.frame(
  Sample = colnames(ex7)
) |>
  mutate(severity = ifelse(grepl("Normal", Sample), "healthy", "infected")) |>
  mutate(severity = factor(severity, levels = c("healthy", "infected")))

used_feature <- c("XAF1", "TMEM173", "CD8A")

plist11 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex7[i, meta7$Sample, drop = T]),
    severity = meta7$severity
  )
  group_n <- data3 %>%
    group_by(severity) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(severity, "\n(n=", n, ")"))

  labels_vec <- setNames(group_n$label, group_n$severity)
  my_comparisons <- list(
    c("healthy", "infected")
  )
  p <- ggboxplot(data3,
    x = "severity", y = "exp",
    color = "severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist11[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/7. GSE215262 XAF1 and other gene expresion boxplot.pdf", height = 8, width = 8)
plot_grid(plotlist = plist11, ncol = 2)
dev.off()

### 8. COVID-19 GSE152075 #####

gset8 <- getGEO("GSE152075", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE152075/", AnnotGPL = F)
meta8 <- pData(gset8[[1]])
meta8 <- meta8[, c("title", "sars-cov-2 positivity:ch1", "gender:ch1", "age:ch1")]
colnames(meta8) <- c("Sample", "status", "gender", "age")
meta8 <- meta8 |>
  mutate(status = recode(status, "neg" = "healthy", "pos" = "infected")) |>
  mutate(status = factor(status, levels = c("healthy", "infected"))) |>
  mutate(gender = recode(gender, "M" = "Male", "F" = "Female", "not collected" = "Unknown")) |>
  mutate(gender = factor(gender, levels = c("Male", "Female", "Unknown")))

ex8 <- read.table("./GSE152075/GSE152075_raw_counts_GEO.txt", header = T, row.names = 1, sep = " ", check.names = F)

# Preserve original counts
ex8_counts <- as.matrix(ex8)

# Match samples between counts and meta8
samples_present <- intersect(colnames(ex8_counts), meta8$Sample)
if (length(samples_present) == 0) {
  stop("No matching samples between ex8 columns and meta8$Sample")
}

# Subset and reorder counts to match metadata
ex8_counts <- ex8_counts[, samples_present, drop = FALSE]

# Remove genes with all-zero counts
keep_mask <- rowSums(ex8_counts, na.rm = TRUE) > 0
ex8_counts <- ex8_counts[keep_mask, , drop = FALSE]

# Coerce to integer storage mode expected by DESeq2
storage.mode(ex8_counts) <- "integer"

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("Package 'DESeq2' is required for VST. Install with BiocManager::install('DESeq2')")
}
suppressMessages(library(DESeq2))

# Prepare colData for DESeq2 (rownames must match column names of count matrix)
coldata <- meta8[match(colnames(ex8_counts), meta8$Sample), , drop = FALSE]
rownames(coldata) <- coldata$Sample

# Build DESeqDataSet and compute VST. Use a simple design (~1) because we only want transformation.
dds <- DESeqDataSetFromMatrix(countData = ex8_counts, colData = coldata, design = ~1)

# Try vst directly; if it fails (rare), run DESeq first then vst
vst_res <- tryCatch(
  {
    vst(dds, blind = TRUE)
  },
  error = function(e) {
    message("vst() failed directly (", e$message, "); running DESeq() then vst() as fallback")
    dds2 <- DESeq(dds, quiet = TRUE)
    vst(dds2, blind = TRUE)
  }
)

ex8_vst <- as.data.frame(assay(vst_res))
# Ensure column order follows meta8
ex8_vst <- ex8_vst[, coldata$Sample, drop = FALSE]

used_feature <- c("XAF1", "STING1", "CD8A")
plist12 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex8_vst[i, meta8$Sample, drop = T]),
    status = meta8$status
  )
  group_n <- data3 %>%
    group_by(status) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(status, "\n(n=", n, ")"))
  labels_vec <- setNames(group_n$label, group_n$status)
  my_comparisons <- list(
    c("healthy", "infected")
  )
  p <- ggboxplot(data3,
    x = "status", y = "exp",
    color = "status", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist12[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/8. GSE152075 XAF1 and other gene expresion boxplot.pdf", height = 8, width = 8)
plot_grid(plotlist = plist12, ncol = 2)
dev.off()

meta8_2 <- meta8 |>
  filter(gender != "Unknown")

plist13 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex8_vst[i, meta8_2$Sample, drop = T]),
    status = meta8_2$status,
    gender = meta8_2$gender
  ) |>
    mutate(group = paste(status, gender, sep = "_")) |>
    mutate(group = factor(group, levels = c("healthy_Male", "healthy_Female", "infected_Male", "infected_Female")))
  group_n <- data3 %>%
    group_by(group) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(group, "\n(n=", n, ")"))
  labels_vec <- setNames(group_n$label, group_n$group)
  my_comparisons <- list(
    c("healthy_Male", "infected_Male"),
    c("healthy_Female", "infected_Female"),
    c("infected_Male", "infected_Female")
  )
  p <- ggboxplot(data3,
    x = "group", y = "exp",
    color = "group", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist13[[i]] <- p
}

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/8. GSE152075 XAF1 and other gene expresion boxplot2.pdf", height = 8, width = 12)
plot_grid(plotlist = plist13, ncol = 2)
dev.off()

# Create meta8_3: parse numeric age where possible and create age_stage
# Strategy: parse a numeric age (first integer found), keep only rows with numeric age for binning,
# for infected samples create ~balanced bins using quantiles (nbins determined by number of infected samples),
# label healthy samples as 'healthy'. This keeps groups reasonably balanced without requiring exact equality.
meta8_3 <- meta8_2 |>
  mutate(
    age_raw = as.character(age),
    age_num = suppressWarnings(as.numeric(sub("^.*?(\\d{1,3}).*$", "\\1", age_raw)))
  )

# Determine infected ages and choose number of bins based on sample size
infected_idx <- which(meta8_3$status == "infected" & !is.na(meta8_3$age_num))
n_inf <- length(infected_idx)
nbins <- if (n_inf >= 60) 4 else if (n_inf >= 30) 3 else if (n_inf >= 15) 2 else if (n_inf >= 5) 2 else 1

if (n_inf > 0 && nbins > 1) {
  # compute quantile breaks; if breaks collapse (due to identical ages) reduce nbins
  q <- unique(quantile(meta8_3$age_num[infected_idx], probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
  while (length(q) < (nbins + 1) && nbins > 1) {
    nbins <- nbins - 1
    q <- unique(quantile(meta8_3$age_num[infected_idx], probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
  }
  if (nbins > 1) {
    # create human-readable labels from the quantile breaks, e.g. infected_20-29
    # lower and upper bounds for each bin
    lower_bounds <- as.integer(floor(q[-length(q)]))
    upper_bounds <- as.integer(ceiling(q[-1]))
    inf_labels <- paste0("infected_", lower_bounds, "-", upper_bounds)
    # cut infected ages using these breaks and labels
    cuts <- rep(NA_character_, nrow(meta8_3))
    cuts[infected_idx] <- as.character(cut(meta8_3$age_num[infected_idx], breaks = q, include.lowest = TRUE, labels = inf_labels))
    meta8_3$age_stage <- ifelse(meta8_3$status == "healthy", "healthy", cuts)
    # ensure factor with stable levels: healthy first then age-interval labels
    meta8_3$age_stage <- factor(meta8_3$age_stage, levels = c("healthy", inf_labels))
  } else {
    # single infected bin
    meta8_3$age_stage <- ifelse(meta8_3$status == "healthy", "healthy", "infected_all_ages")
    meta8_3$age_stage <- factor(meta8_3$age_stage, levels = c("healthy", "infected_all_ages"))
  }
} else {
  # not enough infected age info: mark all infected as single group
  meta8_3$age_stage <- ifelse(meta8_3$status == "healthy", "healthy", "infected_all_ages")
  meta8_3$age_stage <- factor(meta8_3$age_stage, levels = c("healthy", "infected_all_ages"))
}
table(meta8_3$age_stage)
plist14 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex8_vst[i, meta8_3$Sample, drop = T]),
    age_stage = meta8_3$age_stage
  )
  group_n <- data3 %>%
    group_by(age_stage) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(age_stage, "\n(n=", n, ")"))
  labels_vec <- setNames(group_n$label, group_n$age_stage)
  my_comparisons <- list(
    c("healthy", "infected_2-40"),
    c("healthy", "infected_40-56"),
    c("healthy", "infected_56-71"),
    c("healthy", "infected_71-90"),
    c("infected_2-40", "infected_40-56"),
    c("infected_2-40", "infected_56-71"),
    c("infected_2-40", "infected_71-90")
  )

  p <- ggboxplot(data3,
    x = "age_stage", y = "exp",
    color = "age_stage", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist14[[i]] <- p
}
pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/8. GSE152075 XAF1 and other gene expresion boxplot3.pdf", height = 10, width = 14)
plot_grid(plotlist = plist14, ncol = 2)
dev.off()

### 9. COVID-19 GSE157103 #####
ex9 <- read_delim("./GSE157103/GSE157103_genes.tpm.tsv", delim = "\t", col_names = T)
ex9 <- as.data.frame(ex9)
rownames(ex9) <- ex9[, 1]
ex9 <- ex9[, -1]

gset9 <- getGEO("GSE157103", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE157103/", AnnotGPL = F)
meta9 <- pData(gset9[[1]])
meta9 <- meta9[, c("title", "Sex:ch1", "disease state:ch1", "icu:ch1", "age (years):ch1")]
colnames(meta9) <- c("Sample", "Sex", "disease", "icu", "age")
meta9 <- meta9 |>
  mutate(disease = recode(disease, "COVID-19" = "infected", "non-COVID-19" = "healthy")) |>
  mutate(disease = factor(disease, levels = c("healthy", "infected"))) |>
  mutate(Severity = ifelse(disease == "healthy", "healthy", ifelse(icu == "yes", "ICU", "non ICU"))) |>
  mutate(Severity = factor(Severity, levels = c("healthy", "non ICU", "ICU")))
meta9$Sample2 <- colnames(ex9)

used_feature <- c("XAF1", "STING1", "CD8A")
plist15 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex9[i, meta9$Sample2, drop = T]),
    Severity = meta9$Severity
  )
  group_n <- data3 %>%
    group_by(Severity) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(Severity, "\n(n=", n, ")"))
  labels_vec <- setNames(group_n$label, group_n$Severity)
  my_comparisons <- list(
    c("healthy", "non ICU"),
    c("non ICU", "ICU"),
    c("healthy", "ICU")
  )
  p <- ggboxplot(data3,
    x = "Severity", y = "exp",
    color = "Severity", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist15[[i]] <- p
}
pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/9. GSE157103 XAF1 and other gene expresion boxplot.pdf", height = 8, width = 8)
plot_grid(plotlist = plist15, ncol = 2)
dev.off()

meta9_2 <- meta9 |>
  filter(Sex != "unknown") |>
  mutate(Sex = factor(Sex, levels = c("male", "female")))

plist16 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex9[i, meta9_2$Sample2, drop = T]),
    disease = meta9_2$disease,
    Sex = meta9_2$Sex
  ) |>
    mutate(group = paste(disease, Sex, sep = "_")) |>
    mutate(group = factor(group, levels = c("healthy_male", "healthy_female", "infected_male", "infected_female")))

  group_n <- data3 %>%
    group_by(group) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(group, "\n(n=", n, ")"))
  labels_vec <- setNames(group_n$label, group_n$group)
  my_comparisons <- list(
    c("healthy_male", "infected_male"),
    c("healthy_female", "infected_female"),
    c("infected_male", "infected_female")
  )
  p <- ggboxplot(data3,
    x = "group", y = "exp",
    color = "group", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist16[[i]] <- p
}
pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/9. GSE157103 XAF1 and other gene expresion boxplot2.pdf", height = 8, width = 11)
plot_grid(plotlist = plist16, ncol = 2)
dev.off()

# Create meta9_3: parse numeric age where possible and create age_stage
# Strategy: parse a numeric age (first integer found), keep only rows with numeric age for binning,
# for infected samples create ~balanced bins using quantiles (nbins determined by number of infected samples),
# label healthy samples as 'healthy'. This keeps groups reasonably balanced without requiring exact equality.
meta9_3 <- meta9_2 |>
  filter(age != ":") |>
  mutate(
    age_raw = as.character(age),
    age_num = suppressWarnings(as.numeric(sub("^.*?(\\d{1,3}).*$", "\\1", age_raw)))
  )
# Determine infected ages and choose number of bins based on sample size
infected_idx <- which(meta9_3$disease == "infected")
n_inf <- length(infected_idx)
nbins <- if (n_inf >= 60) 4 else if (n_inf >= 30) 3 else if (n_inf >= 15) 2 else if (n_inf >= 5) 2 else 1
if (n_inf > 0 && nbins > 1) {
  # compute quantile breaks; if breaks collapse (due to identical ages) reduce nbins
  q <- unique(quantile(meta9_3$age_num[infected_idx], probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
  while (length(q) < (nbins + 1) && nbins > 1) {
    nbins <- nbins - 1
    q <- unique(quantile(meta9_3$age_num[infected_idx], probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
  }
  if (nbins > 1) {
    # create human-readable labels from the quantile breaks, e.g. infected_20-29
    # lower and upper bounds for each bin
    lower_bounds <- as.integer(floor(q[-length(q)]))
    upper_bounds <- as.integer(ceiling(q[-1]))
    inf_labels <- paste0("infected_", lower_bounds, "-", upper_bounds)
    # cut infected ages using these breaks and labels
    cuts <- rep(NA_character_, nrow(meta9_3))
    cuts[infected_idx] <- as.character(cut(meta9_3$age_num[infected_idx], breaks = q, include.lowest = TRUE, labels = inf_labels))
    meta9_3$age_stage <- ifelse(meta9_3$disease == "healthy", "healthy", cuts)
    # ensure factor with stable levels: healthy first then age-interval labels
    meta9_3$age_stage <- factor(meta9_3$age_stage, levels = c("healthy", inf_labels))
  } else {
    # single infected bin
    meta9_3$age_stage <- ifelse(meta9_3$disease == "healthy", "healthy", "infected_all_ages")
    meta9_3$age_stage <- factor(meta9_3$age_stage, levels = c("healthy", "infected_all_ages"))
  }
} else {
  # not enough infected age info: mark all infected as single group
  meta9_3$age_stage <- ifelse(meta9_3$disease == "healthy", "healthy", "infected_all_ages")
  meta9_3$age_stage <- factor(meta9_3$age_stage, levels = c("healthy", "infected_all_ages"))
}
table(meta9_3$age_stage)

plist17 <- list()
for (i in used_feature) {
  data3 <- data.frame(
    exp = unlist(ex9[i, meta9_3$Sample2, drop = T]),
    age_stage = meta9_3$age_stage
  )
  group_n <- data3 %>%
    group_by(age_stage) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(age_stage, "\n(n=", n, ")"))
  labels_vec <- setNames(group_n$label, group_n$age_stage)
  my_comparisons <- list(
    c("healthy", "infected_21-51"),
    c("healthy", "infected_50-62"),
    c("healthy", "infected_62-74"),
    c("healthy", "infected_74-89"),
    c("infected_21-51", "infected_50-62"),
    c("infected_21-51", "infected_62-74"),
    c("infected_21-51", "infected_74-89")
  )
  p <- ggboxplot(data3,
    x = "age_stage", y = "exp",
    color = "age_stage", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab(paste(i, "expression", sep = " ")) +
    scale_x_discrete(labels = labels_vec) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  p
  plist17[[i]] <- p
}
pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/9. GSE157103 XAF1 and other gene expresion boxplot3.pdf", height = 10, width = 14)
plot_grid(plotlist = plist17, ncol = 2)
dev.off()

### integrate all plots into one PDF #####

#### IAV GSE101702 ######

gset3 <- getGEO("GSE101702", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE101702/", AnnotGPL = F)

ex3 <- exprs(gset3[[1]])

qx <- quantile(ex3, c(0, 0.25, 0.5, 0.75, 0.99, 1))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex3[ex3 <= 0] <- NA
  ex3 <- log2(ex3 + 1)
}

ex3 <- normalizeBetweenArrays(ex3)
boxplot(ex3, outline = FALSE, notch = T, col = factor(meta1$geo_accession), las = 2) ## 查看表达水平

index <- gset3[[1]]@annotation

gpl <- read_delim("./GSE101702/GPL21185-21174.txt", delim = "\t", col_names = T, skip = 18)

ids <- gpl[, c(1, 6)]
colnames(ids) <- c("probe_id", "symbol")

ids <- as.data.frame(ids) |>
  na.omit()

ex3 <- as.data.frame(ex3)

exp3 <- ex3 |>
  dplyr::mutate(probe_id = rownames(ex3)) %>%
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

exp3[1:5, 1:5]
exp3_max <- exp3 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp3_max) <- exp3_max$symbol
exp3_max <- exp3_max[, -c(1:2)]

meta3 <- pData(gset3[[1]])
length(unique(meta2$`subject id:ch1`))

meta3 <- meta3 |>
  mutate(severity = recode(`severity:ch1`, "hlty_ctrl" = "Healthy", "flu_mod" = "Moderate", "flu_svre" = "Severe")) |>
  mutate(severity = factor(severity, levels = c("Healthy", "Moderate", "Severe")))

exp3_max <- exp3_max[, meta3$geo_accession]

data3 <- data.frame(
  exp = unlist(exp3_max["XAF1", meta3$geo_accession, drop = T]),
  severity = meta3$severity
) |>
  mutate(Group = factor(severity, levels = c("Healthy", "Moderate", "Severe")))

group_n <- data3 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))

labels_vec <- setNames(group_n$label, group_n$Group)

data3$label <- factor(labels_vec[data3$Group], levels = labels_vec)

head(data3)

data3$Accession <- "GSE101702"
data3$Virus <- "IAV"
data3 <- data3 |>
  mutate(stage = recode(Group, "Moderate" = "Non-severe"))

head(data3)

#### GSE111368 IAV ####

gset4 <- getGEO("GSE111368", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE111368/", AnnotGPL = F)

ex4 <- exprs(gset4[[1]])

qx <- quantile(ex4, c(0, 0.25, 0.5, 0.75, 0.99, 1))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex4[ex4 <= 0] <- NA
  ex4 <- log2(ex4 + 1)
}

ex4 <- normalizeBetweenArrays(ex4)


index <- gset4[[1]]@annotation

ids <- getIDs("gpl10558")
head(ids)

ex4 <- as.data.frame(ex4)
ex4[1:5, 1:5]

# intersect(rownames(ex4), ids$probe_id)
exp4 <- ex4 |>
  dplyr::mutate(probe_id = rownames(ex4)) %>%
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

dim(exp4)
exp4[1:5, 355:360]

exp4_max <- exp4 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp4_max) <- exp4_max$symbol
exp4_max <- exp4_max[, -c(1:2)]

####

meta4 <- pData(gset4[[1]])
unique(meta4$`timepoint:ch1`)

exp4_max <- exp4_max[, meta4$geo_accession]
exp4_max[1:5, 1:5]

data4 <- data.frame(
  exp = unlist(exp4_max["XAF1", meta4$geo_accession, drop = T]),
  severity = meta4$`t1severity:ch1`,
  timepoint = meta4$`timepoint:ch1`
) |>
  filter(timepoint %in% c("HC", "T1")) |>
  mutate(Group = recode(severity, `HC` = "Healthy", `1` = "Mild", `2` = "Mild", `3` = "Severe")) |>
  mutate(
    Group = factor(Group, levels = c("Healthy", "Mild", "Severe"))
  ) |>
  select(exp, Group)

group_n <- data4 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))

labels_vec <- setNames(group_n$label, group_n$Group)

data4$label <- factor(labels_vec[data4$Group], levels = labels_vec)

head(data4)

data4$Accession <- "GSE111368"
data4$Virus <- "IAV"
data4 <- data4 |>
  mutate(stage = recode(Group, "Mild" = "Non-severe", "Severe" = "Severe"))


#### GSE161731 SARS-CoV-2 #######

gset6 <- getGEO("GSE161731", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE161731/", AnnotGPL = F)

meta6 <- pData(gset6[[1]])
table(meta6$`cohort:ch1`)
meta6 <- meta6 |>
  filter(`cohort:ch1` %in% c("COVID-19", "healthy"))

ex6 <- read_tsv("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/IAV_and_COVID-19_data/GSE161731/GSE161731_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")

ex6 <- as.data.frame(ex6)
rownames(ex6) <- as.character(ex6$GeneID)
ex6 <- ex6[, -1]
common_sample2 <- intersect(colnames(ex6), meta6$geo_accession)

exp6 <- ex6[, common_sample2]
meta6 <- meta6[common_sample2, ]
table(meta6$`cohort:ch1`)
meta6$`time_since_onset:ch1`[which(meta6$`time_since_onset:ch1` == "NA")] <- "Healthy"
table(meta6$`time_since_onset:ch1`)

anno2 <- read_tsv("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/IAV_and_COVID-19_data/GSE161731/Human.GRCh38.p13.annot.tsv.gz")
anno2 <- as.data.frame(anno2)
entrze_id2 <- as.character(anno2$GeneID[which(anno2$Symbol == "XAF1")])

data6 <- data.frame(
  exp = unlist(exp6[entrze_id2, meta6$geo_accession, drop = T]),
  severity = meta6$`time_since_onset:ch1`
) |>
  mutate(Group = recode(severity, `early` = "Mild", `middle` = "Severe", `late` = "Severe")) |>
  mutate(
    Group = factor(Group, levels = c("Healthy", "Mild", "Severe"))
  ) |>
  select(exp, Group)

group_n <- data6 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))

labels_vec <- setNames(group_n$label, group_n$Group)

data6$label <- factor(labels_vec[data6$Group], levels = labels_vec)

data6$Accession <- "GSE161731"
data6$Virus <- "SARS-CoV-2"

data6 <- data6 |>
  mutate(stage = recode(Group, "Mild" = "Non-severe", "Severe" = "Severe"))

head(data6)

#### GSE51808 DENV #######

gset7 <- getGEO("GSE51808", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE51808/", AnnotGPL = F)

ex7 <- exprs(gset7[[1]])

qx <- quantile(ex7, c(0, 0.25, 0.5, 0.75, 0.99, 1))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex7[ex7 <= 0] <- NA
  ex7 <- log2(ex7 + 1)
}

ex7 <- normalizeBetweenArrays(ex7)

ex7[1:5, 1:5]

index <- gset7[[1]]@annotation

gpl <- read_delim("./GSE51808/GPL13158.annot.gz", delim = "\t", col_names = T, skip = 27)

ids <- gpl[, c(1, 3)]
colnames(ids) <- c("probe_id", "symbol")

ids <- as.data.frame(ids) |>
  na.omit()

ex7 <- as.data.frame(ex7)

exp7 <- ex7 |>
  dplyr::mutate(probe_id = rownames(ex7)) %>%
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

exp7[1:5, 1:5]
exp7_max <- exp7 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp7_max) <- exp7_max$symbol
exp7_max <- exp7_max[, -c(1:2)]

meta7 <- pData(gset7[[1]])
table(meta7$`status:ch1`)

exp7_max <- exp7_max[, meta7$geo_accession]

exp7_max[1:5, 1:5]

data7 <- data.frame(
  exp = unlist(exp7_max["XAF1", meta7$geo_accession, drop = T]),
  severity = meta7$`status:ch1`
) |>
  mutate(Group = recode(severity, `control` = "Healthy", `convalescent` = "Convalescent", `DF` = "DF", `DHF` = "DHF")) |>
  filter(Group != "Convalescent") |>
  mutate(
    Group = factor(Group, levels = c("Healthy", "DF", "DHF"))
  ) |>
  select(exp, Group)

group_n <- data7 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))

labels_vec <- setNames(group_n$label, group_n$Group)

data7$label <- factor(labels_vec[data7$Group], levels = labels_vec)

head(data7)

data7$Accession <- "GSE51808"
data7$Virus <- "DENV"

data7 <- data7 |>
  mutate(stage = recode(Group, "DF" = "Non-severe", "DHF" = "Severe"))


#### GSE105450 RSV #######

gset10 <- getGEO("GSE105450", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE105450/", AnnotGPL = F)

ex10 <- exprs(gset10[[1]])

qx <- quantile(ex10, c(0, 0.25, 0.5, 0.75, 0.99, 1))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex10[ex10 <= 0] <- NA
  ex10 <- log2(ex10 + 1)
}

ex10[1:5, 1:5]

index <- gset10[[1]]@annotation
gpl <- getIDs(index)
head(gpl)

gpl <- gpl[, c("probe_id", "symbol")]

ex10 <- as.data.frame(ex10)
exp10 <- ex10 |>
  dplyr::mutate(probe_id = rownames(ex10)) %>%
  inner_join(gpl, by = "probe_id")

exp10 <- exp10 %>%
  dplyr::select(probe_id, symbol, everything())

exp10[1:5, 1:5]

exp10_max <- exp10 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp10_max) <- exp10_max$symbol
exp10_max <- exp10_max[, -c(1:2)]

meta10 <- pData(gset10[[1]])
table(meta10$`condition:ch1`)
summary(as.numeric(meta10$`age (mos):ch1`))
meta10$severity <- meta10$`condition:ch1`
table(meta10$severity)

exp10_max <- exp10_max[, meta10$geo_accession]
exp10_max[1:5, 1:5]

exp10_max[is.na(exp10_max)] <- 0

data10 <- data.frame(
  exp = unlist(exp10_max["XAF1", meta10$geo_accession, drop = T]),
  Severity = meta10$severity
) |>
  mutate(Group = recode(Severity, `RSV Outpatient` = "OP", `RSV Inpatient` = "IP")) |>
  mutate(
    Group = factor(Group, levels = c("Healthy", "OP", "IP"))
  ) |>
  select(exp, Group)

group_n <- data10 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))

labels_vec <- setNames(group_n$label, group_n$Group)

data10$label <- factor(labels_vec[data10$Group], levels = labels_vec)

head(data10)

data10$Accession <- "GSE105450"
data10$Virus <- "RSV"
data10 <- data10 |>
  mutate(stage = recode(Group, "OP" = "Non-severe", "IP" = "Severe"))

#### GSE157103 SARS-CoV-2 #####
ex9 <- read_delim("./GSE157103/GSE157103_genes.tpm.tsv", delim = "\t", col_names = T)
ex9 <- as.data.frame(ex9)
rownames(ex9) <- ex9[, 1]
ex9 <- ex9[, -1]

gset9 <- getGEO("GSE157103", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE157103/", AnnotGPL = F)
meta9 <- pData(gset9[[1]])
meta9 <- meta9[, c("title", "Sex:ch1", "disease state:ch1", "icu:ch1", "age (years):ch1")]
colnames(meta9) <- c("Sample", "Sex", "disease", "icu", "age")
meta9 <- meta9 |>
  mutate(disease = recode(disease, "COVID-19" = "infected", "non-COVID-19" = "Healthy")) |>
  mutate(disease = factor(disease, levels = c("Healthy", "infected"))) |>
  mutate(Severity = ifelse(disease == "Healthy", "Healthy", ifelse(icu == "yes", "ICU", "non-ICU"))) |>
  mutate(Severity = factor(Severity, levels = c("Healthy", "non-ICU", "ICU")))

meta9$Sample2 <- colnames(ex9)

data9 <- data.frame(
  exp = unlist(ex9["XAF1", meta9$Sample2, drop = T]),
  Severity = meta9$Severity
) |>
  mutate(Group = Severity) |>
  select(exp, Group)

group_n <- data9 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))
labels_vec <- setNames(group_n$label, group_n$Group)

data9$label <- factor(labels_vec[data9$Group], levels = labels_vec)

head(data9)

data9$Accession <- "GSE157103"
data9$Virus <- "SARS-CoV-2"
data9 <- data9 |>
  mutate(stage = recode(Group, "non-ICU" = "Non-severe", "ICU" = "Severe"))

#### GSE94892 DENV #######

gset8 <- getGEO("GSE94892", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE94892/", AnnotGPL = F)

meta8 <- pData(gset8[[1]])
table(meta8$`cohort:ch1`)


ex8 <- read_tsv("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/IAV_and_COVID-19_data/GSE94892/GSE94892_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")

ex8 <- as.data.frame(ex8)
ex8[1:5, 1:5]
rownames(ex8) <- as.character(ex8$GeneID)
ex8 <- ex8[, -1]
common_sample2 <- intersect(colnames(ex8), meta8$geo_accession)

exp8 <- ex8[, common_sample2]
meta8 <- meta8[common_sample2, ]
table(meta8$`disease:ch1`)

anno2 <- read_tsv("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/IAV_and_COVID-19_data/GSE161731/Human.GRCh38.p13.annot.tsv.gz")
anno2 <- as.data.frame(anno2)
entrze_id2 <- as.character(anno2$GeneID[which(anno2$Symbol == "XAF1")])

data8 <- data.frame(
  exp = unlist(exp8[entrze_id2, meta8$geo_accession, drop = T]),
  severity = meta8$`disease:ch1`
) |>
  mutate(Group = recode(severity, `Other Febrile Illness` = "OFI", `Dengue fever` = "DF", `Dengue hemorrhagic fever` = "DHF", `Dengue shock syndrome` = "DSS")) |>
  filter(Group %in% c("OFI", "DF", "DHF")) |>
  mutate(
    Group = factor(Group, levels = c("OFI", "DF", "DHF"))
  ) |>
  select(exp, Group)

group_n <- data8 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")"))

labels_vec <- setNames(group_n$label, group_n$Group)

data8$label <- factor(labels_vec[data8$Group], levels = labels_vec)

data8$Accession <- "GSE94892"
data8$Virus <- "DENV"

data8 <- data8 |>
  mutate(stage = recode(Group, "OFI" = "Healthy", "DF" = "Non-severe", "DHF" = "Severe"))

head(data8)

#### GSE69606 RSV ########

gset11 <- getGEO("GSE69606", GSEMatrix = TRUE, getGPL = FALSE, destdir = "./GSE69606/", AnnotGPL = F)

ex11 <- exprs(gset11[[1]])

qx <- quantile(ex11, c(0, 0.25, 0.5, 0.75, 0.99, 1))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex11[ex11 <= 0] <- NA
  ex11 <- log2(ex11 + 1)
}

ex11[1:5, 1:5]

index <- gset11[[1]]@annotation
gpl <- getIDs(index)
head(gpl)

ex11 <- as.data.frame(ex11)
exp11 <- ex11 |>
  dplyr::mutate(probe_id = rownames(ex11)) %>%
  inner_join(gpl, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())

exp11[1:5, 1:5]

exp11_max <- exp11 %>%
  group_by(symbol) %>%
  slice_max(
    order_by = rowMeans(across(where(is.numeric))), # 按所有数值列平均排序
    with_ties = FALSE # 只保留最大的一行
  ) %>%
  ungroup() |>
  as.data.frame()

rownames(exp11_max) <- exp11_max$symbol
exp11_max <- exp11_max[, -c(1:2)]

meta11 <- pData(gset11[[1]])
table(meta11$`disease:ch1`)
meta11$`disease:ch1`[which(meta11$`acute or recovery:ch1` == "recovery")] <- "Recovery"

exp11_max <- exp11_max[, meta11$geo_accession]
exp11_max[1:5, 1:5]

exp11_max[is.na(exp11_max)] <- 0

data11 <- data.frame(
  exp = unlist(exp11_max["XAF1", meta11$geo_accession, drop = T]),
  severity = meta11$`disease:ch1`
) |>
  mutate(Group = recode(severity, `Recovery` = "Healthy", `mild` = "Mild", `moderate` = "Moderate", `severe` = "Severe")) |>
  filter(Group != "Moderate") |>
  mutate(
    Group = factor(Group, levels = c("Healthy", "Mild", "Severe"))
  ) |>
  select(exp, Group)

data11 <- data11[-c(2, 12), ]

group_n <- data11 %>%
  group_by(Group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(Group, "\n(n=", n, ")")) %>%
  select(label, Group)

labels_vec <- setNames(group_n$label, group_n$Group)

data11$label <- factor(labels_vec[data11$Group], levels = labels_vec)
head(data11)

data11$Accession <- "GSE69606"
data11$Virus <- "RSV"
data11 <- data11 |>
  mutate(stage = recode(Group, "Mild" = "Non-severe", "Severe" = "Severe"))


#### plot all ###########

head(data3)
table(data10$Group)

plot_list <- list(data3, data4, data6, data7, data10, data9, data8, data11)

for (i in 1:8) {
  plot_df <- plot_list[[i]]

  if (i == 1) {
    my_comparisons <- list(
      c("Healthy", "Moderate"),
      c("Moderate", "Severe")
    )
  } else if (i == 4) {
    my_comparisons <- list(
      c("Healthy", "DF"),
      c("DF", "DHF")
    )
  } else if (i == 5) {
    my_comparisons <- list(
      c("Healthy", "OP"),
      c("OP", "IP")
    )
  } else if (i == 6) {
    my_comparisons <- list(
      c("Healthy", "non-ICU"),
      c("non-ICU", "ICU")
    )
  } else if (i == 7) {
    my_comparisons <- list(
      c("OFI", "DF"),
      c("DF", "DHF")
    )
  } else {
    my_comparisons <- list(
      c("Healthy", "Mild"),
      c("Mild", "Severe")
    )
  }

  p <- ggboxplot(plot_df,
    x = "Group", y = "exp",
    color = "Group", palette = "lancet",
    bxp.errorbar = T, bxp.errorbar.width = 0.1,
    add = "jitter", short.panel.labs = FALSE, nrow = 1
  ) +
    ylab("XAF1 normalized expression") +
    scale_x_discrete(labels = levels(plot_df$label)) + # 替换 x 轴显示
    expand_limits(y = 0) +
    stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
    theme(
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    plot_annotation(
      title = paste(unique(plot_df$Accession), unique(plot_df$Virus), sep = ": "),
      theme = theme(
        plot.title = element_text(
          size = 16, face = "bold", hjust = 0.5,
          margin = margin(b = 10)
        ),
        plot.margin = margin(20, 20, 20, 20)
      )
    )
  pdf(paste("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/0. ", unique(plot_df$Accession), "-", unique(plot_df$Virus), ".pdf", sep = ""), height = 5, width = 5)
  print(p)
  dev.off()
}


# 依赖包
library(dplyr)
library(ggplot2)
library(tidyr)

# --- 可选：设置随机种子（如果后续用了随机分布的点） ---
set.seed(123)

# --- 1) 缩放函数（默认 z-score 标准化） ---
# 也支持 method = "minmax"
scale_exp_df <- function(df, exp_col = "exp", method = c("z", "minmax"), dataset_name = NULL) {
  method <- match.arg(method)
  if (!(exp_col %in% names(df))) stop("找不到 exp 列: ", exp_col)
  if (method == "z") {
    df <- df %>% mutate(scaled = as.numeric(scale(.data[[exp_col]])))
  } else {
    rng_min <- min(df[[exp_col]], na.rm = TRUE)
    rng_max <- max(df[[exp_col]], na.rm = TRUE)
    df <- df %>% mutate(scaled = (.data[[exp_col]] - rng_min) / (rng_max - rng_min))
  }
  if (!is.null(dataset_name)) df$dataset <- dataset_name
  df
}

# --- 2) 应用到三个数据集 ---
data3s <- scale_exp_df(data3, exp_col = "exp", method = "z", dataset_name = unique(data3$Accession))
data4s <- scale_exp_df(data4, exp_col = "exp", method = "z", dataset_name = unique(data4$Accession))
data6s <- scale_exp_df(data6, exp_col = "exp", method = "z", dataset_name = unique(data6$Accession))
data7s <- scale_exp_df(data7, exp_col = "exp", method = "z", dataset_name = unique(data7$Accession))
data10s <- scale_exp_df(data10, exp_col = "exp", method = "z", dataset_name = unique(data10$Accession))
data9s <- scale_exp_df(data9, exp_col = "exp", method = "z", dataset_name = unique(data9$Accession))
data11s <- scale_exp_df(data11, exp_col = "exp", method = "z", dataset_name = unique(data11$Accession))
data8s <- scale_exp_df(data8, exp_col = "exp", method = "z", dataset_name = unique(data8$Accession))

head(data8s)

# 检查一下缩放结果（可选）
summary(data3s$scaled)
summary(data4s$scaled)
summary(data9s$scaled)
summary(data6s$scaled)
summary(data7s$scaled)
summary(data10s$scaled)

# --- 3) 合并 ---
combined <- bind_rows(data3s, data4s, data6s, data7s, data10s, data9s, data11s, data8s)
unique(data10s$stage)
# 确保 Group 有一个稳定的顺序（若需要按 appearance 排序）
# 如果 Group 已是 factor 且顺序正确，可以跳过这步
combined <- combined %>% mutate(stage = as.character(stage))
group_levels <- c("Healthy", "Non-severe", "Severe") # appearance order
combined <- combined %>% mutate(stage = factor(stage, levels = group_levels))

head(combined)
# --- 4) 按 stage 和 dataset 汇总（计算 mean、sd、n、se、95% CI） ---
summary_df <- combined %>%
  group_by(Virus, dataset, stage) %>%
  summarise(
    mean_scaled = mean(scaled, na.rm = TRUE),
    sd_scaled = sd(scaled, na.rm = TRUE),
    n = sum(!is.na(scaled)),
    .groups = "drop"
  ) %>%
  mutate(
    se = sd_scaled / sqrt(pmax(n, 1)),
    tval = qt(0.975, df = pmax(n - 1, 1)),
    ci = tval * se,
    lower = mean_scaled - ci,
    upper = mean_scaled + ci
  )

summary_df <- summary_df |>
  as.data.frame() |>
  mutate(x = recode(stage, "Healthy" = 0, "Non-severe" = 1, "Severe" = 2)) |>
  mutate(x = as.numeric(x))

summary_df <- summary_df %>%
  mutate(stage = factor(stage, levels = c("Healthy", "Non-severe", "Severe")))
my_colors <- c(
  "#E64B35FF", "#F08065FF", # 红色系
  "#4DBBD5FF", "#76D1E0FF", # 蓝色系
  "#00A087FF", "#33B59FFF", # 绿色系
  "#9B59B6FF", "#B07AC9FF" # 紫色系
)

summary_df$label <- paste(summary_df$Virus, summary_df$dataset, sep = ":")

# --- 5) 绘图：每个 dataset 一条线 + ribbon（95% CI） ---
p <- ggplot(summary_df, aes(x = x, y = mean_scaled, group = dataset, color = label, fill = label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.18, color = NA,
    position = position_dodge(width = 0.15)
  ) +
  geom_line(position = position_dodge(width = 0.15), size = 0.8) +
  geom_point(position = position_dodge(width = 0.15), size = 1.6) +
  theme_classic(base_size = 12) +
  labs(x = "Group", y = "XAF1 mRNAS Scaled exp (z-score)", title = "Trend of scaled expression by disease severity") +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = 0:2, labels = levels(summary_df$stage)) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(axis.text = element_text(size = 18, color = "black"))

# 显示并保存
print(p)
ggsave("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/0. trend_scaled_exp_by_severity_ribbon.pdf", p, width = 6, height = 4)


### KEGG pathway activity in infected samples #########

kegg_df <- PrepareDB(species = "Homo_sapiens", db = c("KEGG"))
kegg_df <- kegg_df$Homo_sapiens$KEGG

head(kegg_df$TERM2GENE)
# kegg_df 包含 columns: gs_name, gene_symbol 等
kegg_list <- split(kegg_df$TERM2GENE$symbol, kegg_df$TERM2GENE$Term)
head(kegg_list)

exp3_max[1:5, 1:5]

bpp <- BiocParallel::MulticoreParam(workers = 90) # 根据机器核数调整

ssgseapar3 <- ssgseaParam(as.matrix(exp3_max), kegg_list)

ssgsea_results3 <- gsva(
  ssgseapar3,
  BPPARAM = bpp
)

ssgsea_results3[1:5, 1:5]

xaf1_vec3 <- exp3_max["XAF1", ]

res_ct3 <- pbapply::pbsapply(rownames(ssgsea_results3), function(pathway) {
  scores <- ssgsea_results3[pathway, ]
  xaf <- as.numeric(xaf1_vec3)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct3 <- t(res_ct3)
res_ct3 <- as.data.frame(res_ct3, stringsAsFactors = FALSE)
res_ct3$adj.p <- p.adjust(res_ct3$p, method = "BH")
res_ct3$pathway_ID <- rownames(res_ct3)

colnames(res_ct3)[1] <- "rho"

head(res_ct3)

# 整理列类型
res_ct3 <- res_ct3 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct3)
res_ct3 <- merge(res_ct3, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct3$Accession <- "GSE101702"
res_ct3$Virus <- "IAV"

####

ssgseapar4 <- ssgseaParam(as.matrix(exp4_max), kegg_list)

ssgsea_results4 <- gsva(
  ssgseapar4,
  BPPARAM = bpp
)

ssgsea_results4[1:5, 1:5]
xaf1_vec4 <- exp4_max["XAF1", ]

res_ct4 <- pbapply::pbsapply(rownames(ssgsea_results4), function(pathway) {
  scores <- ssgsea_results4[pathway, ]
  xaf <- as.numeric(xaf1_vec4)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct4 <- t(res_ct4)
res_ct4 <- as.data.frame(res_ct4, stringsAsFactors = FALSE)
res_ct4$adj.p <- p.adjust(res_ct4$p, method = "BH")
res_ct4$pathway_ID <- rownames(res_ct4)

colnames(res_ct4)[1] <- "rho"

head(res_ct4)

# 整理列类型
res_ct4 <- res_ct4 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct4)
res_ct4 <- merge(res_ct4, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct4$Accession <- "GSE111368"
res_ct4$Virus <- "IAV"

####
anno2$GeneID <- as.character(anno2$GeneID)

exp6_symbol <- exp6 |>
  rownames_to_column(var = "GeneID") %>%
  inner_join(anno2[, c("GeneID", "Symbol")], by = "GeneID") %>%
  dplyr::select(GeneID, Symbol, everything())

head(exp6_symbol)
exp6_symbol <- exp6_symbol[!duplicated(exp6_symbol$Symbol), ]
rownames(exp6_symbol) <- exp6_symbol$Symbol
exp6_symbol <- exp6_symbol[, -c(1:2)]

ssgseapar6 <- ssgseaParam(as.matrix(exp6_symbol), kegg_list)

ssgsea_results6 <- gsva(
  ssgseapar6,
  BPPARAM = bpp
)

ssgsea_results6[1:5, 1:5]
xaf1_vec6 <- exp6_symbol["XAF1", ]

res_ct6 <- pbapply::pbsapply(rownames(ssgsea_results6), function(pathway) {
  scores <- ssgsea_results6[pathway, ]
  xaf <- as.numeric(xaf1_vec6)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct6 <- t(res_ct6)
res_ct6 <- as.data.frame(res_ct6, stringsAsFactors = FALSE)
res_ct6$adj.p <- p.adjust(res_ct6$p, method = "BH")
res_ct6$pathway_ID <- rownames(res_ct6)

colnames(res_ct6)[1] <- "rho"

head(res_ct6)

# 整理列类型
res_ct6 <- res_ct6 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct6)
res_ct6 <- merge(res_ct6, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct6$Accession <- "GSE161731"
res_ct6$Virus <- "SARS-CoV-2"

####

ssgseapar7 <- ssgseaParam(as.matrix(exp7_max), kegg_list)

ssgsea_results7 <- gsva(
  ssgseapar7,
  BPPARAM = bpp
)

ssgsea_results7[1:5, 1:5]
xaf1_vec7 <- exp7_max["XAF1", ]

res_ct7 <- pbapply::pbsapply(rownames(ssgsea_results7), function(pathway) {
  scores <- ssgsea_results7[pathway, ]
  xaf <- as.numeric(xaf1_vec7)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct7 <- t(res_ct7)
res_ct7 <- as.data.frame(res_ct7, stringsAsFactors = FALSE)
res_ct7$adj.p <- p.adjust(res_ct7$p, method = "BH")
res_ct7$pathway_ID <- rownames(res_ct7)

colnames(res_ct7)[1] <- "rho"

head(res_ct7)

# 整理列类型
res_ct7 <- res_ct7 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct7)
res_ct7 <- merge(res_ct7, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct7$Accession <- "GSE51808"
res_ct7$Virus <- "DENV"

####
exp8_symbol <- exp8 |>
  rownames_to_column(var = "GeneID") %>%
  inner_join(anno2[, c("GeneID", "Symbol")], by = "GeneID") %>%
  dplyr::select(GeneID, Symbol, everything())

head(exp8_symbol)
exp8_symbol <- exp8_symbol[!duplicated(exp8_symbol$Symbol), ]
rownames(exp8_symbol) <- exp8_symbol$Symbol
exp8_symbol <- exp8_symbol[, -c(1:2)]

ssgseapar8 <- ssgseaParam(as.matrix(exp8_symbol), kegg_list)

ssgsea_results8 <- gsva(
  ssgseapar8,
  BPPARAM = bpp
)

ssgsea_results8[1:5, 1:5]
xaf1_vec8 <- exp8_symbol["XAF1", ]

res_ct8 <- pbapply::pbsapply(rownames(ssgsea_results8), function(pathway) {
  scores <- ssgsea_results8[pathway, ]
  xaf <- as.numeric(xaf1_vec8)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct8 <- t(res_ct8)
res_ct8 <- as.data.frame(res_ct8, stringsAsFactors = FALSE)
res_ct8$adj.p <- p.adjust(res_ct8$p, method = "BH")
res_ct8$pathway_ID <- rownames(res_ct8)

colnames(res_ct8)[1] <- "rho"

head(res_ct8)

# 整理列类型
res_ct8 <- res_ct8 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct8)
res_ct8 <- merge(res_ct8, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct8$Accession <- "GSE94892"
res_ct8$Virus <- "DENV"

####

ssgseapar9 <- ssgseaParam(as.matrix(ex9), kegg_list)

ssgsea_results9 <- gsva(
  ssgseapar9,
  BPPARAM = bpp
)

ssgsea_results9[1:5, 1:5]
xaf1_vec9 <- ex9["XAF1", ]

res_ct9 <- pbapply::pbsapply(rownames(ssgsea_results9), function(pathway) {
  scores <- ssgsea_results9[pathway, ]
  xaf <- as.numeric(xaf1_vec9)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct9 <- t(res_ct9)
res_ct9 <- as.data.frame(res_ct9, stringsAsFactors = FALSE)
res_ct9$adj.p <- p.adjust(res_ct9$p, method = "BH")
res_ct9$pathway_ID <- rownames(res_ct9)

colnames(res_ct9)[1] <- "rho"
head(res_ct9)

# 整理列类型
res_ct9 <- res_ct9 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct9)
res_ct9 <- merge(res_ct9, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct9$Accession <- "GSE157103"
res_ct9$Virus <- "SARS-CoV-2"

####

ssgseapar10 <- ssgseaParam(as.matrix(exp10_max), kegg_list)

ssgsea_results10 <- gsva(
  ssgseapar10,
  BPPARAM = bpp
)

ssgsea_results10[1:5, 1:5]
xaf1_vec10 <- exp10_max["XAF1", ]

res_ct10 <- pbapply::pbsapply(rownames(ssgsea_results10), function(pathway) {
  scores <- ssgsea_results10[pathway, ]
  xaf <- as.numeric(xaf1_vec10)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct10 <- t(res_ct10)
res_ct10 <- as.data.frame(res_ct10, stringsAsFactors = FALSE)
res_ct10$adj.p <- p.adjust(res_ct10$p, method = "BH")
res_ct10$pathway_ID <- rownames(res_ct10)

colnames(res_ct10)[1] <- "rho"

head(res_ct10)

# 整理列类型
res_ct10 <- res_ct10 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct10)
res_ct10 <- merge(res_ct10, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct10$Accession <- "GSE105450"
res_ct10$Virus <- "RSV"

####

exp11_max <- exp11_max[, meta11$geo_accession]
exp11_max <- exp11_max[, -c(2, 12)]

ssgseapar11 <- ssgseaParam(as.matrix(exp11_max), kegg_list)

head(exp11_max)

ssgsea_results11 <- gsva(
  ssgseapar11,
  BPPARAM = bpp
)

ssgsea_results11[1:5, 1:5]
xaf1_vec11 <- exp11_max["XAF1", ]

res_ct11 <- pbapply::pbsapply(rownames(ssgsea_results11), function(pathway) {
  scores <- ssgsea_results11[pathway, ]
  xaf <- as.numeric(xaf1_vec11)
  # remove NA
  ok <- !is.na(scores) & !is.na(xaf)
  if (sum(ok) < 4) {
    return(c(rho = NA, p = NA, n = sum(ok)))
  } else {
    ctst <- suppressWarnings(cor.test(xaf[ok], scores[ok], method = "spearman", exact = FALSE))
    c(rho = ctst$estimate, p = ctst$p.value, n = sum(ok))
  }
})
res_ct11 <- t(res_ct11)
res_ct11 <- as.data.frame(res_ct11, stringsAsFactors = FALSE)
res_ct11$adj.p <- p.adjust(res_ct11$p, method = "BH")
res_ct11$pathway_ID <- rownames(res_ct11)

colnames(res_ct11)[1] <- "rho"

head(res_ct11)

# 整理列类型
res_ct11 <- res_ct11 %>%
  mutate(
    rho = as.numeric(rho),
    p = as.numeric(p),
    adj.p = as.numeric(adj.p),
    n = as.integer(n)
  )

head(res_ct11)
res_ct11 <- merge(res_ct11, kegg_df$TERM2NAME, by.x = "pathway_ID", by.y = "Term", all.x = TRUE)
res_ct11$Accession <- "GSE69606"
res_ct11$Virus <- "RSV"

####
head(res_ct7)

rho_all_res <- rbind(res_ct3, res_ct4, res_ct6, res_ct7, res_ct10, res_ct9, res_ct8, res_ct11)

write.csv(rho_all_res, "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/0. XAF1_vs_KEGG_pathway_correlation_all_datasets.csv", row.names = F)

rho_all_res <- rho_all_res |>
  arrange(Virus, Accession) |>
  as.data.frame()

length(rho_all_res$Name)
unique(rho_all_res$Accession)
# 统计每个细胞类型中显著相关的pathway数量
res_count <- rho_all_res %>%
  filter(adj.p < 0.05) %>%
  group_by(Virus, Accession, n) %>%
  summarise(sig_pathway_count = n_distinct(Name)) %>%
  arrange(desc(sig_pathway_count))

res_count$label <- paste0(res_count$Accession, "\n(n=", res_count$n, ")")
# 绘制条形图
p1 <- ggplot(res_count, aes(
  x = reorder(Accession, -sig_pathway_count),
  y = sig_pathway_count, fill = Virus
)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sig_pathway_count), vjust = -0.3, size = 4) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = res_count$label, breaks = res_count$Accession) +
  labs(
    x = "", y = "Number of significant pathways",
    title = "Significant pathways correlated with XAF1"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 10)
  )
p1

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/0.XAF1_pathway_sig_count_by_celltype.pdf", width = 9, height = 6)
print(p1)
dev.off()

res_sig <- rho_all_res %>%
  filter(!is.na(adj.p)) %>%
  group_by(Name) %>%
  summarise(
    n_sig = sum(adj.p < 0.05),
    sig_ct = paste(Accession[adj.p < 0.05], collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig))

# 提取4个细胞类型都显著的pathway
common_sig <- res_sig %>% filter(n_sig >= 6)

View(res_sig)
res_plot <- rho_all_res %>%
  filter(Name %in% common_sig$Name)

# 12) 可视化：取 top N 显著相关的 pathway（按任一 celltype q）
topN2 <- 30
# 选基于任一 celltype 最小 q 排序
path_order2 <- res_plot %>%
  group_by(Name) %>%
  summarise(min_q = min(adj.p, na.rm = TRUE)) %>%
  arrange(min_q) %>%
  pull(Name)
top_paths <- head(path_order2, topN2)

res_plot <- rho_all_res %>% filter(Name %in% top_paths)

res_plot$label <- paste0(res_plot$Virus, ":", res_plot$Accession, "\n(n=", res_plot$n, ")")

res_plot$Accession <- factor(res_plot$Accession, levels = unique(res_plot$Accession))

res_plot$Name <- factor(res_plot$Name, levels = rev(sort(unique(res_plot$Name))))


KEGG_Class <- read.table("/home/bio/biosoft/GSEA_ref/KEGG_pathway_Class.txt", header = T, sep = "\t")
head(KEGG_Class)

res_plot2 <- res_plot %>%
  left_join(KEGG_Class, by = c("Name" = "pathway_name")) |>
  mutate(
    Name = factor(Name, levels = rev(intersect(KEGG_Class$pathway_name, Name))),
    Main = factor(Main, levels = intersect(KEGG_Class$Main, Main)),
    Class = factor(Class, levels = intersect(KEGG_Class$Class, Class))
  )

head(res_plot2)


library(RColorBrewer)

main_cols <- setNames(brewer.pal(length(unique(res_plot2$Main)), "Set2"), unique(res_plot2$Main))
class_cols <- setNames(brewer.pal(length(unique(res_plot2$Class)), "Set3"), unique(res_plot2$Class))

annot_df <- res_plot2 |>
  select("Name", "Main", "Class") |>
  mutate(Name = as.character(Name))
head(annot_df)
annot_df <- annot_df[!duplicated(annot_df$Name), ]
annot_df$Name <- factor(annot_df$Name, levels = levels(res_plot2$Name))

p1 <- annot_df %>%
  ggplot() +
  geom_tile(aes(x = 1, y = Name, fill = Main), color = NA, width = 0.5, height = 0.5) +
  scale_fill_manual(values = main_cols) +
  coord_equal() +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "bold", color = "#000000", size = 13),
    legend.background = element_rect(color = "#737373")
  )
p1

p2 <- annot_df %>%
  ggplot() +
  geom_tile(aes(x = 1, y = Name, fill = Class), color = NA, width = 0.5, height = 0.5) +
  scale_fill_manual(values = class_cols) +
  coord_equal() +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.background = element_rect(color = "#737373")
  )
p2

# 绘制热图
p3 <- ggplot() +
  geom_tile(
    data = res_plot2,
    aes(x = Accession, y = Name),
    fill = "#d9d9d9", color = "white", width = 0.65, height = 0.65
  ) +
  geom_tile(
    data = res_plot2 %>% filter(adj.p < 0.05),
    aes(x = Accession, y = Name, fill = rho),
    color = NA, width = 0.65, height = 0.65
  ) +
  geom_star(
    data = res_plot2 %>% filter(adj.p < 0.05),
    aes(x = Accession, y = Name),
    size = 1.5, fill = "#d7301f", color = "#d7301f"
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D", midpoint = 0,
    name = expression(rho ~ "(XAF1 & pathway activity)")
  ) +
  scale_x_discrete(position = "top", labels = res_plot$label, breaks = res_plot$Accession) +
  labs(x = "", y = "") +
  coord_equal(clip = "off") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      face = "bold", color = "black", size = 9,
      angle = 45, hjust = 0
    ),
    legend.background = element_rect(color = "#737373")
  )
p3


# 添加标题到 patchwork 图中
p_combine <- p1 + p2 + p3 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Top30 Common significant pathways (in all 4 Virus)",
    theme = theme(
      plot.title = element_text(
        size = 16, face = "bold", hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.margin = margin(20, 20, 20, 20)
    )
  )

p_combine

pdf("/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/IAV_and_COVID-19/0.XAF1_common_pathway_corr_heatmap_top30.pdf", width = 13, height = 15)
print(p_combine)
dev.off()
