# Software dependencies

Package versions should be recorded when the final manuscript outputs are
frozen. The list below describes the software interfaces used by the curated
scripts.

## R

The analyses use R 4.4-compatible releases of these packages:

- Data and plotting: `readr`, `readxl`, `dplyr`, `tidyr`, `tibble`, `ggplot2`,
  `ggrepel`, `patchwork`, `cowplot`, `viridis`, `ggpubr`, `ggstar`,
  `ggnewscale`, `ComplexHeatmap`, and `circlize`.
- Statistical analysis: `ordinal`, `survival`, `rstatix`, and `effsize`.
- Genomics: `AnnotationDbi`, `org.Hs.eg.db`, `GEOquery`, `limma`, `GSVA`,
  `BiocParallel`, `SCP`, `idmap2`, `Seurat`, `harmony`, `Gviz`,
  `GenomicRanges`, `rtracklayer`, `biomaRt`, `AnnotationHub`, `easier`, and
  `matrixStats`.
- Methylation: `minfi` and
  `IlluminaHumanMethylation450kanno.ilmn12.hg19`.

The EaSIeR scripts also require the matching `opt_models.rds` and
`opt_xtrain_stats.rds` objects identified by `XAF1_EASIER_ROOT`.

## Python and cross-species IDR tools

- Python 3.10 or newer.
- `numpy>=2.0,<3` from `scripts/fig2/idr/requirements.txt`.
- NCBI Datasets CLI (analysis snapshot used 18.39.0).
- MAFFT (analysis snapshot used 7.526); `run_fig2h.sh` invokes L-INS-i through
  `--localpair --maxiterate 1000` for the final retained set.
- The source distribution of IUPred2A, used in `long` mode.

## eCLIP command-line tools

The preprocessing workflow checks for `FastQC`, `MultiQC`, `fastp`, `STAR`,
`samtools`, `MACS2`, Java/Picard, and `awk`. Sample bigWig generation also
requires deepTools `bamCoverage`; bedGraph conversion requires UCSC
`bedGraphToBigWig` when enabled.

## Protected data

CCLE/DepMap files, TCGA/GDC derivatives, eCLIP reads and tracks, the disease
score workbook, downloaded GEO matrices, reference indexes, and EaSIeR model
objects are not stored in Git. Their locations are supplied through
`config/environment` using `config/environment.example` as the template.
