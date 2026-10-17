#!/usr/bin/env bash
set -euo pipefail

# Download the current NCBI Datasets mammalian ortholog set for human XAF1.
# The resulting directory is intentionally ignored by Git.

PROJECT_ROOT="${XAF1_PROJECT_ROOT:-$(pwd)}"
DATASETS_BIN="${XAF1_DATASETS_BIN:-datasets}"
OUT_DIR="${PROJECT_ROOT}/data/fig2/idr/ncbi"

mkdir -p "${OUT_DIR}/packages"

"${DATASETS_BIN}" summary gene gene-id 54739 \
    --ortholog mammals \
    --as-json-lines \
    > "${OUT_DIR}/xaf1_mammals_gene_report.jsonl"

"${DATASETS_BIN}" download gene gene-id 54739 \
    --ortholog mammals \
    --include protein,rna \
    --filename "${OUT_DIR}/packages/batch_01.zip"

echo "Downloaded XAF1 mammalian orthologues to ${OUT_DIR}"
