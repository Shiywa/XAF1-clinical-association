# XAF1 figure-analysis code

This repository contains the analysis code that directly supports the
computational panels in the XAF1 manuscript. The branch was reorganized on
2026-10-17 to remove unrelated exploratory and side-project scripts, replace
workstation-specific paths with configuration variables, and make the link
between each manuscript panel and its code explicit.

本仓库仅整理 XAF1 稿件中计算分析图所对应的代码。2026-10-17 的整理版删除了与
定稿图无关的探索性脚本和其他项目代码，将个人电脑绝对路径改为环境变量，并为每个
图面板标明入口脚本、输入数据和主要输出。

## Repository layout

```text
config/environment.example       input and software path template
docs/FIGURE_SCRIPT_MAP.tsv       machine-readable panel-to-script map
docs/DEPENDENCIES.md             software and R/Python dependencies
scripts/common/                  shared path helpers
scripts/fig1/                    CCLE RNA/protein analysis
scripts/fig2/eclip/              eCLIP processing and locus track
scripts/fig2/idr/                mammalian XAF1 IDR conservation analysis
scripts/fig4_fig5_ed5/           longitudinal mouse disease scores
scripts/fig6_ed7/                viral bulk and single-cell cohorts
scripts/fig7_ed8/                TCGA methylation and EaSIeR analyses
```

Raw data, reference indexes, intermediate files, and generated results are
excluded from Git. They remain under `data/`, `work/`, or `results/` in a local
checkout.

## Configuration

Run commands from the repository root. Copy the template, replace the paths
needed by the selected panel, and export them into the shell:

```bash
cp config/environment.example config/environment
set -a
source config/environment
set +a
```

Only variables used by the chosen script need to be set. Most R scripts also
provide a repository-relative default under `data/` when practical.

## Figure entry points

| Figure panel | Analysis | Entry point |
|---|---|---|
| Fig. 1a | CCLE STING1/TMEM173 RNA and protein | `scripts/fig1/fig1a_ccle_sting_rna_protein.R` |
| Fig. 2b | XAF1 eCLIP processing and STING1 locus track | `scripts/fig2/eclip/run_eclip_preprocessing.sh`; `scripts/fig2/eclip/fig2b_eclip_track.R` |
| Fig. 2h | XAF1 IDR conservation across non-human mammals | `scripts/fig2/idr/download_ncbi_orthologues.sh`; `scripts/fig2/idr/run_fig2h.sh` |
| Fig. 4f, Fig. 5e, ED Fig. 5e | Longitudinal mouse disease scores | `scripts/fig4_fig5_ed5/fig4f_fig5e_ed5e_disease_score_clmm.R` |
| Fig. 6k-l, ED Fig. 7i | Viral-cohort XAF1 expression and KEGG associations | `scripts/fig6_ed7/fig6kl_ed7i_bulk_viral.R` |
| ED Fig. 7j-k | GSE271284 single-cell analysis | `scripts/fig6_ed7/ed7jk_scRNA_xaf1.R` |
| Fig. 7e, ED Fig. 8c | Pan-cancer EaSIeR analysis | `scripts/fig7_ed8/fig7e_ed8c_easier_pancancer.R` |
| Fig. 7f | Pan-cancer methylation, expression and EaSIeR | `scripts/fig7_ed8/fig7f_methylation_easier.R` |

The complete input/output map is in
[`docs/FIGURE_SCRIPT_MAP.tsv`](docs/FIGURE_SCRIPT_MAP.tsv).

## Important scope decisions

- The reorganized Fig. 2h pipeline queries the NCBI Datasets mammalian
  orthologue set, chooses one protein per non-human taxon, applies alignment
  coverage and identity quality control, calculates IUPred2 long-mode scores,
  and projects them to human XAF1 coordinates. The current analysis retained
  220 non-human mammalian orthologues; the count is written from the run data
  and may change when NCBI records are updated.
- Fig. 7f uses 19 cancer types after excluding OV for the methylation display.
  Fig. 7e and ED Fig. 8c use the 20-cancer expression/EaSIeR collection,
  including OV. The two scopes are intentionally recorded separately.
- ED Fig. 7k displays cell-level distributions but exports donor-level summary
  values and Wilcoxon tests so cells from one donor are not counted as
  independent biological replicates.

## Typical commands

```bash
Rscript scripts/fig1/fig1a_ccle_sting_rna_protein.R
bash scripts/fig2/eclip/run_eclip_preprocessing.sh
Rscript scripts/fig2/eclip/fig2b_eclip_track.R
bash scripts/fig2/idr/download_ncbi_orthologues.sh
bash scripts/fig2/idr/run_fig2h.sh
Rscript scripts/fig4_fig5_ed5/fig4f_fig5e_ed5e_disease_score_clmm.R
Rscript scripts/fig6_ed7/fig6kl_ed7i_bulk_viral.R
Rscript scripts/fig6_ed7/ed7jk_scRNA_xaf1.R
Rscript scripts/fig7_ed8/fig7e_ed8c_easier_pancancer.R
Rscript scripts/fig7_ed8/fig7f_methylation_easier.R
```

Use `bash scripts/check_repository.sh` for syntax and portability checks that
do not require the protected input datasets.


