suppressPackageStartupMessages({
    library(Gviz)
    library(GenomicRanges)
    library(rtracklayer)
    library(AnnotationDbi)
    library(dplyr)
    library(stringr)
})

# =========================
# 0) 你需要填的路径
# =========================
genome <- "hg38"

res_path <- "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/results/CLIP"
bw_path <- "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/Clip/CS2601026_人_clip/2_Peaks/2.7_bw"
narrowPeak_path <- "/home/Toshiba1/Lab_Ma/yuanyi/XAF1-clinical-association/data/Clip/CS2601026_人_clip/2_Peaks/2.1_peaks"
bw_ip <- c(
    HSV_IP_rep1   = file.path(bw_path, "IP-3-HSV-1.bw"),
    HSV_IP_rep2   = file.path(bw_path, "IP-4-HSV-2.bw")
)
bw_in <- c(
    HSV_input_rep1 = file.path(bw_path, "Input-3-HSV-1.bw"),
    HSV_input_rep2 = file.path(bw_path, "Input-4-HSV-1.bw")
)

peaks_narrowPeak1 <- file.path(narrowPeak_path, "IP-3-HSV-1_macs_peaks.narrowPeak")

peaks_narrowPeak2 <- file.path(narrowPeak_path, "IP-4-HSV-2_macs_peaks.narrowPeak")

summits_bed1 <- file.path(narrowPeak_path, "IP-3-HSV-1_macs_summits.bed")

summits_bed2 <- file.path(narrowPeak_path, "IP-4-HSV-2_macs_summits.bed")

# 你要展示的基因（人：TMEM173 是 STING1）
gene_symbol <- "STING1" # 或 "STING1"（但更常用 TMEM173）

# =========================
# 1) 获取基因坐标 + 转录本结构
#    方案A：TxDb UCSC knownGene
# =========================
### 构建 KCNMA1 gene model
gtf <- read.delim("/home/bio/biosoft/reference/GRCh38/gtf/STING1_exons.gtf",
    header = FALSE, comment.char = "#", stringsAsFactors = FALSE
)

colnames(gtf) <- c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
head(gtf)
# -------------------------
# 2. 筛选 KCNMA1 的 exon
# -------------------------
gtf_exon <- gtf %>%
    filter(feature == "exon" & grepl("ENSG00000184584", attribute))

get_attr <- function(x, attr) {
    pattern <- paste0(attr, ' "?([^";]+)"?')
    val <- str_match(x, pattern)[, 2]
    trimws(val)
}

gtf_exon$gene_id <- get_attr(gtf_exon$attribute, "gene_id")
gtf_exon$transcript_id <- get_attr(gtf_exon$attribute, "transcript_id")
gtf_exon$exon_id <- get_attr(gtf_exon$attribute, "exon_id")
gtf_exon$gene_name <- get_attr(gtf_exon$attribute, "gene_name")
gtf_exon$width <- abs(gtf_exon$start - gtf_exon$end)

head(gtf_exon)
gtf_exon$chromosome <- paste0("chr", gtf_exon$chromosome)
STING1_geneModels <- gtf_exon %>%
    select(
        chromosome,
        start,
        end,
        width,
        strand,
        feature,
        gene = gene_id,
        exon = exon_id,
        transcript = transcript_id,
        symbol = gene_name
    ) %>%
    distinct() # 去重
head(STING1_geneModels)
table(STING1_geneModels$exon)

grep(STING1_geneModels$transcript, pattern = "ENST00000330794", value = T)

STING1_geneModels <- STING1_geneModels %>%
    filter(transcript %in% c(
        "ENST00000330794"
    ))

STING1_grtrack <- GeneRegionTrack(STING1_geneModels,
    genome = "hg38",
    chromosome = "chr5", name = "STING1 transcripts"
)

############

# =========================
# 2) bigWig 信号 DataTrack
# =========================
make_bw_track <- function(bw_path, track_name) {
    DataTrack(
        range = bw_path,
        genome = genome,
        chromosome = chr,
        name = track_name,
        type = "histogram", # 或 "l"
        col.histogram = "#2C7FB8",
        fill.histogram = "#2C7FB8",
        baseline = 0,
        ylim = c(0, 60)
    )
}

ip_tracks <- lapply(names(bw_ip), function(nm) make_bw_track(bw_ip[[nm]], nm))
in_tracks <- lapply(names(bw_in), function(nm) make_bw_track(bw_in[[nm]], nm))

# 统一一下 input 的颜色（可选）
for (i in seq_along(in_tracks)) {
    displayPars(in_tracks[[i]]) <- list(col.histogram = "#A6A6A6", fill.histogram = "#A6A6A6")
}

# =========================
# 3) peaks / summits AnnotationTrack
# =========================
gr_peaks1 <- rtracklayer::import(peaks_narrowPeak1) # narrowPeak 也能 import 成 GRanges

gr_peaks2 <- rtracklayer::import(peaks_narrowPeak2)

chr <- "chr5"
peakTrack1 <- AnnotationTrack(gr_peaks1,
    genome = genome,
    chromosome = chr,
    name = "MACS peaks",
    fill = "#FDB863",
    col = NA
)
peaksTrack2 <- AnnotationTrack(gr_peaks2,
    genome = genome,
    chromosome = chr,
    name = "MACS peaks",
    fill = "#FDB863",
    col = NA
)

# =========================
# 4) Axis + Ideogram（ideogram可选）
# =========================
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = genome, chromosome = chr)

# =========================
# 5) 输出
# =========================
tracks <- c(
    list(itrack, gtrack),
    ip_tracks[[1]], in_tracks[[1]],
    ip_tracks[[2]], in_tracks[[2]],
    STING1_grtrack
)

pdf(file.path(res_path, "HSV_seCLIP_TMEM173_Gviz.pdf"), width = 7, height = 5)
plotTracks(tracks,
    chromosome = "chr5", from = 139475533 - 200, to = 139482935 + 200,
    transcriptAnnotation = "transcript",
    background.panel = "white",
    col.axis = "black",
    col.title = "black",
    main = paste0("HSV seCLIP at ", gene_symbol)
)
dev.off()


library(rtracklayer)
library(GenomicRanges)

# 定义你关心的区域
region <- GRanges(
    "chr5",
    IRanges(
        139475533 - 200,
        139482935 + 200
    )
)

# 读取该区间的bigWig信号
bw_data <- import(file.path(bw_path, "IP-4-HSV-2.bw"), which = region)

bw_data
head(bw_data)

bw_data_strong <- bw_data[score(bw_data) > 0]
bw_data_strong

bw_data_strong_merged <- reduce(bw_data_strong)
bw_data_strong_merged
df <- as.data.frame(bw_data_strong_merged)
write.csv(df,
    file.path(res_path, "HSV_rep2_strong_binding_regions.csv"),
    row.names = FALSE, quote = FALSE
)
