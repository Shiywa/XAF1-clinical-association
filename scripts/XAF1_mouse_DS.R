library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(lme4) # lmer
library(lmerTest) # 给 lmer 提供 p 值
library(emmeans) # 辅助对比
library(ggplot2)
library(survival)
library(survminer)
library(ordinal)
library(gridExtra)
library(grid)

ds_file <- "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/mouse/disease_Score/disease_score.xlsx"
res_path <- "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/mouse/disease_score"


ds_EV71 <- read_excel(ds_file, sheet = 1)
ds_HSV_brain <- read_excel(ds_file, sheet = 2)
df_HSV_tail <- read_excel(ds_file, sheet = 3)

ds_EV71 <- as.data.frame(ds_EV71)
ds_HSV_brain <- as.data.frame(ds_HSV_brain)
df_HSV_tail <- as.data.frame(df_HSV_tail)

### HSV brain infection disease score analysis ######

df <- ds_HSV_brain

# 1. 创建唯一老鼠 ID
df$MouseID <- paste0("M", seq_len(nrow(df)))

# 2. 解析 group（从第一列）
df$Group_label <- df[[1]]

# 3. 去掉原始 group 列
df2 <- df %>% select(-1)

# 4. wide → long
df_long <- df2 %>%
    pivot_longer(
        cols = -c(MouseID, Group_label),
        names_to = "Day",
        values_to = "Score"
    )

# 5. Day 转成数值
df_long$Day <- as.numeric(df_long$Day)

# 6. Alive 状态（Score = NA → 死亡）
df_long <- df_long %>%
    mutate(Alive = ifelse(is.na(Score), 0, 1)) |>
    as.data.frame()

df_long$Group <- ifelse(df_long$Group_label == "WT (n=13)", "WT", "XAF1_KO")

df_long$MouseID <- as.factor(df_long$MouseID)
df_long$Group <- as.factor(df_long$Group)

df_long$Day <- as.numeric(df_long$Day) # 若 day 原本就是数值，这一步无害

df_mod <- df_long %>% filter(!is.na(Score))
df_mod$Score <- factor(df_mod$Score,
    levels = sort(unique(df_mod$Score)),
    ordered = TRUE
)
str(df_mod$Score)


HSV_brain_clmm <- clmm(
    Score ~ Group * Day + (1 | MouseID),
    data = df_mod, link = "logit"
)
summary(HSV_brain_clmm) # 固定效应估计值与 p 值（lmerTest）

# ---------------------------------------------------
# 2. 绘制折线图
# ---------------------------------------------------
summary_df <- df_mod %>%
    mutate(Score = as.numeric(as.character(Score))) %>%
    group_by(Group_label, Day) %>%
    summarise(
        Mean = mean(Score, na.rm = TRUE),
        SD = sd(Score, na.rm = TRUE),
        N = n(),
        SE = SD / sqrt(N),
        Lower = Mean - 1.96 * SE,
        Upper = Mean + 1.96 * SE,
        .groups = "drop"
    )
unique(summary_df$Group_label)

model_table <- data.frame(
    Term = c("Group × Day"),
    Estimate = c("beta = 0.46"),
    Pvalue = c("p = 8.12e-05")
)

tb <- tableGrob(model_table,
    rows = NULL,
    theme = ttheme_minimal(
        core = list(bg_params = list(fill = "white", alpha = 0.7)),
        colhead = list(bg_params = list(fill = "white", alpha = 0.7))
    )
)

p <- ggplot(summary_df, aes(x = Day, y = Mean, color = Group_label, group = Group_label)) +

    # CI 区域（推荐，有 publication 质感）
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group_label),
        alpha = 0.15, color = NA
    ) +

    # 均值折线
    geom_line(size = 1.2) +

    # 每天的误差棒
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
        width = 0.2, size = 0.5
    ) +
    # 均值点
    scale_color_manual(values = c("WT (n=13)" = "#000000", "Xaf1-/-(n=10)" = "#E64B35")) +
    scale_fill_manual(values = c("WT (n=13)" = "#000000", "Xaf1-/-(n=10)" = "#E64B35")) +
    labs(
        x = "Day post infection",
        y = "Mean Clinical Score",
        title = "Daily Mean Disease Scores with 95% CI",
        subtitle = "WT (n=13) vs Xaf1-/-(n=10)"
    ) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 23, color = "black"),
        axis.title = element_text(size = 25)
    ) +
    labs(color = "Group", fill = "Group") +
    annotation_custom(
        grob = tb,
        xmin = -0.4, xmax = 5,
        ymin = 4.5, ymax = Inf
    ) +
    coord_cartesian(clip = "off")

pdf(file.path(res_path, "1.XAF1_KO_HSV_brain_Disease_Score_Lineplot.pdf"), width = 10, height = 8)
print(p)
dev.off()

### EV71 infection disease score analysis ######

df2 <- ds_EV71

# 1. 创建唯一老鼠 ID
df2$MouseID <- paste0("M", seq_len(nrow(df2)))
# 2. 解析 group（从第一列）
df2$Group_label <- df2[[1]]
# 3. 去掉原始 group 列
df2b <- df2 %>% select(-1)
# 4. wide → long
df2_long <- df2b %>%
    pivot_longer(
        cols = -c(MouseID, Group_label),
        names_to = "Day",
        values_to = "Score"
    )
# 5. Day 转成数值
df2_long$Day <- as.numeric(df2_long$Day)
# 6. Alive 状态（Score = NA → 死亡）
df2_long <- df2_long %>%
    mutate(Alive = ifelse(is.na(Score), 0, 1)) |>
    as.data.frame()

df2_long$Group <- ifelse(df2_long$Group_label == "WT (n=10)", "WT", "XAF1_KO")

df2_long$MouseID <- as.factor(df2_long$MouseID)
df2_long$Group <- as.factor(df2_long$Group)

df2_long$Day <- as.numeric(df2_long$Day) # 若 day 原本就是数值，这一步无害

df2_mod <- df2_long %>% filter(!is.na(Score))
df2_mod$Score <- factor(df2_mod$Score,
    levels = sort(unique(df2_mod$Score)),
    ordered = TRUE
)
str(df2_mod$Score)
table(df2_mod$Score, df2_mod$Group, df2_mod$Day)

EV71_clmm <- clmm(
    Score ~ Group * Day + (1 | MouseID),
    data = df2_mod |>
        filter(Day >= 3), link = "logit"
)
summary(EV71_clmm) # 固定效应估计值与 p 值

# ---------------------------------------------------
# 2. 绘制折线图
# ---------------------------------------------------
summary_df2 <- df2_mod %>%
    mutate(Score = as.numeric(as.character(Score))) %>%
    group_by(Group_label, Day) %>%
    summarise(
        Mean = mean(Score, na.rm = TRUE),
        SD = sd(Score, na.rm = TRUE),
        N = n(),
        SE = SD / sqrt(N),
        Lower = Mean - 1.96 * SE,
        Upper = Mean + 1.96 * SE,
        .groups = "drop"
    )
unique(summary_df2$Group_label)
model_table2 <- data.frame(
    Term = c("Group × Day"),
    Estimate = c("beta = 3.5469"),
    Pvalue = c("p = 6.96e-06")
)
tb2 <- tableGrob(model_table2,
    rows = NULL,
    theme = ttheme_minimal(
        core = list(bg_params = list(fill = "white", alpha = 0.7)),
        colhead = list(bg_params = list(fill = "white", alpha = 0.7))
    )
)
p2 <- ggplot(summary_df2, aes(x = Day, y = Mean, color = Group_label, group = Group_label)) +

    # CI 区域（推荐，有 publication 质感）
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group_label),
        alpha = 0.15, color = NA
    ) +

    # 均值折线
    geom_line(size = 1.2) +

    # 每天的误差棒
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
        width = 0.2, size = 0.5
    ) +
    scale_x_continuous(labels = c(0, 5, 10), breaks = c(0, 5, 10)) +
    # 均值点
    scale_color_manual(values = c("WT (n=10)" = "#000000", "Xaf1-/-(n=11)" = "#E64B35")) +
    scale_fill_manual(values = c("WT (n=10)" = "#000000", "Xaf1-/-(n=11)" = "#E64B35")) +
    labs(
        x = "Day post infection",
        y = "Mean Clinical Score",
        title = "Daily Mean Disease Scores with 95% CI",
        subtitle = "WT (n=10) vs Xaf1-/-(n=11)"
    ) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 23, color = "black"),
        axis.title = element_text(size = 25)
    ) +
    labs(color = "Group", fill = "Group") +
    annotation_custom(
        grob = tb2,
        xmin = -0.4, xmax = 4.5,
        ymin = 4.5, ymax = Inf
    ) +
    coord_cartesian(clip = "off")

pdf(file.path(res_path, "2.XAF1_KO_EV71_Disease_Score_Lineplot.pdf"), width = 10, height = 8)
print(p2)
dev.off()

### HSV tail infection disease score analysis ######

df3 <- df_HSV_tail
# 1. 创建唯一老鼠 ID
df3$MouseID <- paste0("M", seq_len(nrow(df3)))
# 2. 解析 group（从第一列）
df3$Group_label <- df3[[1]]
# 3. 去掉原始 group 列
df3b <- df3 %>% select(-1)
# 4. wide → long
df3_long <- df3b %>%
    pivot_longer(
        cols = -c(MouseID, Group_label),
        names_to = "Day",
        values_to = "Score"
    )
# 5. Day 转成数值
df3_long$Day <- as.numeric(df3_long$Day)
# 6. Alive 状态（Score = NA → 死亡）
df3_long <- df3_long %>%
    mutate(Alive = ifelse(is.na(Score), 0, 1)) |>
    as.data.frame()
df3_long$Group <- ifelse(df3_long$Group_label == "WT (n=12)", "WT", "XAF1_KO")
df3_long$MouseID <- as.factor(df3_long$MouseID)
df3_long$Group <- as.factor(df3_long$Group)
df3_long$Day <- as.numeric(df3_long$Day) # 若 day 原本就是数值，这一步无害
df3_mod <- df3_long %>% filter(!is.na(Score))
df3_mod$Score <- factor(df3_mod$Score,
    levels = sort(unique(df3_mod$Score)),
    ordered = TRUE
)
str(df3_mod$Score)
table(df3_mod$Score, df3_mod$Group, df3_mod$Day)
HSV_tail_clmm <- clmm(
    Score ~ Group * Day + (1 | MouseID),
    data = df3_mod, link = "logit"
)
summary(HSV_tail_clmm) # 固定效应估计值与 p 值
# ---------------------------------------------------
# 2. 绘制折线图
# ---------------------------------------------------
summary_df3 <- df3_mod %>%
    mutate(Score = as.numeric(as.character(Score))) %>%
    group_by(Group_label, Day) %>%
    summarise(
        Mean = mean(Score, na.rm = TRUE),
        SD = sd(Score, na.rm = TRUE),
        N = n(),
        SE = SD / sqrt(N),
        Lower = Mean - 1.96 * SE,
        Upper = Mean + 1.96 * SE,
        .groups = "drop"
    )
unique(summary_df3$Group_label)
model_table3 <- data.frame(
    Term = c("Group × Day"),
    Estimate = c("beta = 1.59856"),
    Pvalue = c("p = 1.10e-07")
)
tb3 <- tableGrob(model_table3,
    rows = NULL,
    theme = ttheme_minimal(
        core = list(bg_params = list(fill = "white", alpha = 0.7)),
        colhead = list(bg_params = list(fill = "white", alpha = 0.7))
    )
)
p3 <- ggplot(summary_df3, aes(x = Day, y = Mean, color = Group_label, group = Group_label)) +
    # CI 区域（推荐，有 publication 质感）
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Group_label),
        alpha = 0.15, color = NA
    ) +
    # 均值折线
    geom_line(size = 1.2) +
    # 每天的误差棒
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
        width = 0.2, size = 0.5
    ) +
    scale_x_continuous(labels = c(0, 5, 10), breaks = c(0, 5, 10)) +
    # 均值点
    scale_color_manual(values = c("WT (n=12)" = "#000000", "Xaf1-/- (n=12)" = "#E64B35")) +
    scale_fill_manual(values = c("WT (n=12)" = "#000000", "Xaf1-/- (n=12)" = "#E64B35")) +
    labs(
        x = "Day post infection",
        y = "Mean Clinical Score",
        title = "Daily Mean Disease Scores with 95% CI",
        subtitle = "WT (n=12) vs Xaf1-/-(n=12)"
    ) +
    ylim(-0.5, 6.5) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 23, color = "black"),
        axis.title = element_text(size = 25)
    ) +
    labs(color = "Group", fill = "Group") +
    annotation_custom(
        grob = tb3,
        xmin = -0.1, xmax = 6,
        ymin = 5.2, ymax = Inf
    ) +
    coord_cartesian(clip = "off")

pdf(file.path(res_path, "3.XAF1_KO_HSV_tail_Disease_Score_Lineplot.pdf"), width = 10, height = 8)
print(p3)
dev.off()
