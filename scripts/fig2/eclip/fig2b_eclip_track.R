#!/usr/bin/env Rscript

# Fig. 2b: XAF1 eCLIP signal across the human STING1/TMEM173 locus.
# Inputs are processed bigWig tracks plus a small GTF containing the displayed
# STING1 transcript. Raw-read processing is in run_eclip_preprocessing.sh.

suppressPackageStartupMessages({
    library(Gviz)
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
    library(stringr)
})

project_root <- Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd())
source(file.path(project_root, "scripts", "common", "config.R"))

genome <- "hg38"
chromosome <- "chr5"
plot_from <- 139475533L - 200L
plot_to <- 139482935L + 200L

bigwig_dir <- xaf1_env_path(
    "XAF1_ECLIP_BIGWIG_DIR",
    xaf1_path("data", "eclip", "bigwig")
)
peak_dir <- xaf1_env_path(
    "XAF1_ECLIP_PEAK_DIR",
    xaf1_path("data", "eclip", "peaks")
)
sting_gtf <- xaf1_env_path(
    "XAF1_STING1_GTF",
    xaf1_path("data", "eclip", "STING1_exons.gtf"),
    must_exist = TRUE
)
output_dir <- xaf1_output_dir("fig2")

track_files <- list(
    rep1 = c(IP = "IP-3-HSV-1.bw", input = "Input-3-HSV-1.bw"),
    rep2 = c(IP = "IP-4-HSV-2.bw", input = "Input-4-HSV-1.bw")
)
replicate_to_plot <- Sys.getenv("XAF1_ECLIP_REPLICATE", unset = "rep2")
if (!replicate_to_plot %in% names(track_files)) {
    stop("XAF1_ECLIP_REPLICATE must be one of: ", paste(names(track_files), collapse = ", "))
}

selected <- file.path(bigwig_dir, track_files[[replicate_to_plot]])
names(selected) <- names(track_files[[replicate_to_plot]])
invisible(lapply(selected, xaf1_require_file, label = "eCLIP bigWig"))

gtf <- read.delim(
    sting_gtf,
    header = FALSE,
    comment.char = "#",
    stringsAsFactors = FALSE
)
colnames(gtf) <- c(
    "chromosome", "source", "feature", "start", "end", "score",
    "strand", "frame", "attribute"
)

get_attr <- function(x, attribute) {
    pattern <- paste0(attribute, ' "?([^";]+)"?')
    trimws(str_match(x, pattern)[, 2])
}

sting_models <- gtf |>
    filter(feature == "exon", grepl("ENSG00000184584", attribute)) |>
    mutate(
        gene = get_attr(attribute, "gene_id"),
        transcript = get_attr(attribute, "transcript_id"),
        exon = get_attr(attribute, "exon_id"),
        symbol = get_attr(attribute, "gene_name"),
        chromosome = ifelse(startsWith(chromosome, "chr"), chromosome, paste0("chr", chromosome)),
        width = abs(start - end)
    ) |>
    filter(transcript == "ENST00000330794") |>
    select(chromosome, start, end, width, strand, feature, gene, exon, transcript, symbol) |>
    distinct()

if (nrow(sting_models) == 0L) {
    stop("The requested STING1 transcript ENST00000330794 was not found in ", sting_gtf)
}

gene_track <- GeneRegionTrack(
    sting_models,
    genome = genome,
    chromosome = chromosome,
    name = "STING1 transcript"
)

make_signal_track <- function(path, label, colour) {
    DataTrack(
        range = path,
        genome = genome,
        chromosome = chromosome,
        name = label,
        type = "histogram",
        col.histogram = colour,
        fill.histogram = colour,
        baseline = 0,
        ylim = c(0, 60)
    )
}

ip_track <- make_signal_track(selected[["IP"]], paste0("XAF1 IP (", replicate_to_plot, ")"), "#2C7FB8")
input_track <- make_signal_track(selected[["input"]], paste0("SMInput (", replicate_to_plot, ")"), "#A6A6A6")

# Peak files are the source of the three binding regions annotated during
# figure assembly. The panel itself uses labels instead of a separate peak row.
peak_file <- file.path(
    peak_dir,
    if (replicate_to_plot == "rep1") "IP-3-HSV-1_macs_peaks.narrowPeak" else "IP-4-HSV-2_macs_peaks.narrowPeak"
)
xaf1_require_file(peak_file, "MACS2 narrowPeak file")

tracks <- list(
    IdeogramTrack(genome = genome, chromosome = chromosome),
    GenomeAxisTrack(),
    ip_track,
    input_track,
    gene_track
)

output_file <- file.path(output_dir, "Fig2b_STING1_eCLIP_track.pdf")
pdf(output_file, width = 5, height = 4)
plotTracks(
    tracks,
    chromosome = chromosome,
    from = plot_from,
    to = plot_to,
    transcriptAnnotation = "transcript",
    background.panel = "white",
    col.axis = "black",
    col.title = "black"
)
dev.off()

message("Created ", output_file)
