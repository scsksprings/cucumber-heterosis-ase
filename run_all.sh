#!/usr/bin/env bash
# Regenerate phenotype statistics with the original R script, then figures and Tables S1–S8.
set -euo pipefail
cd "$(dirname "$0")"

export REPO_ROOT="$PWD"
export FIG_OUT="${FIG_OUT:-$REPO_ROOT/figures/output}"
export PHENO_OUT="${PHENO_OUT:-$FIG_OUT/phenotype}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$FIG_OUT/.matplotlib-cache}"
PYTHON_BIN="${PYTHON:-python3}"
RSCRIPT_BIN="${RSCRIPT:-Rscript}"

if ! command -v "$RSCRIPT_BIN" >/dev/null 2>&1; then
  echo "Rscript is required. Install base R >= 4.3 or set RSCRIPT to its executable path." >&2
  exit 1
fi
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Python is required. Set PYTHON to a Python >= 3.10 executable." >&2
  exit 1
fi

export MPLBACKEND=Agg
mkdir -p "$PHENO_OUT" "$MPLCONFIGDIR"
"$PYTHON_BIN" scripts/validate_release.py
# These c22 CSV files are generated intermediates, not missing source inputs.
# Preserve the original R analysis; do not substitute extracts from frozen result tables.
"$RSCRIPT_BIN" --vanilla analysis/phenotype/heterosis_rcbd_lm.R \
  data/phenotype/block_observations.csv "$PHENO_OUT/c22"
for s in fig1 fig2 fig3 fig4 fig5 fig6 figS6 figS_others figS4_atlas; do
  echo "Rendering $s"
  "$PYTHON_BIN" "figures/scripts/$s.py"
done
"$PYTHON_BIN" tables/build_tables.py
"$PYTHON_BIN" scripts/validate_release.py --built
echo "Done: $FIG_OUT and $REPO_ROOT/tables/output/"
