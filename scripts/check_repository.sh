#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="${XAF1_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RSCRIPT_BIN="${XAF1_RSCRIPT_BIN:-Rscript}"
PYTHON_BIN="${XAF1_PYTHON_BIN:-python3}"

cd "${PROJECT_ROOT}"

while IFS= read -r file; do
    "${RSCRIPT_BIN}" -e 'parse(file=commandArgs(TRUE)[1])' "${file}" >/dev/null
done < <(find scripts -type f -name '*.R' -print | sort)

while IFS= read -r file; do
    "${PYTHON_BIN}" -m py_compile "${file}"
done < <(find scripts -type f -name '*.py' -print | sort)

while IFS= read -r file; do
    bash -n "${file}"
done < <(find scripts -type f -name '*.sh' -print | sort)

if rg -n '/home/|/Users/|Toshiba[0-9]*' scripts config README.md docs --glob '!check_repository.sh'; then
    echo "Absolute workstation path found." >&2
    exit 1
fi

git diff --check
echo "Repository syntax and portability checks passed."
