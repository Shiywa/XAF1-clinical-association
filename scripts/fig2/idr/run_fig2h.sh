#!/usr/bin/env bash
set -euo pipefail

# Reproduce Fig. 2h after download_ncbi_orthologues.sh has completed.

PROJECT_ROOT="${XAF1_PROJECT_ROOT:-$(pwd)}"
PYTHON_BIN="${XAF1_PYTHON_BIN:-python3}"
RSCRIPT_BIN="${XAF1_RSCRIPT_BIN:-Rscript}"
MAFFT_BIN="${XAF1_MAFFT_BIN:-mafft}"
IUPRED_SOURCE="${XAF1_IUPRED_SOURCE:?Set XAF1_IUPRED_SOURCE to the saezlab/iupred source directory}"

INPUT_DIR="${PROJECT_ROOT}/data/fig2/idr/ncbi"
RESULT_DIR="${PROJECT_ROOT}/results/fig2/idr"
FIGURE_DIR="${PROJECT_ROOT}/results/fig2"

mkdir -p "${RESULT_DIR}" "${FIGURE_DIR}"

"${PYTHON_BIN}" "${PROJECT_ROOT}/scripts/fig2/idr/build_xaf1_idr_data.py" \
    --package-dir "${INPUT_DIR}/packages" \
    --gene-report "${INPUT_DIR}/xaf1_mammals_gene_report.jsonl" \
    --iupred-source "${IUPRED_SOURCE}" \
    --mafft "${MAFFT_BIN}" \
    --outdir "${RESULT_DIR}" \
    --min-human-coverage 0.65 \
    --min-idr-coverage 0.60 \
    --min-identity 0.20

"${RSCRIPT_BIN}" "${PROJECT_ROOT}/scripts/fig2/idr/plot_xaf1_idr.R" \
    "${RESULT_DIR}/position_summary_human_coordinates.tsv" \
    "${FIGURE_DIR}"

echo "Created Fig. 2h files in ${FIGURE_DIR}"
