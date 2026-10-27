#!/usr/bin/env bash
# Run from bundled fragment counts; all output is separate from the paper tables.
set -euo pipefail
cd "$(dirname "$0")/.."
profile="${1:-common3}"
model="${2:-genewise}"
case "$profile" in all|common3|common4) ;; *) echo 'Profile must be all, common3 or common4' >&2; exit 2;; esac
case "$model" in genewise|shared) ;; *) echo 'Model must be genewise or shared' >&2; exit 2;; esac
py="${PYTHON:-python3}"
if [[ "$model" == genewise ]]; then
  "$py" pipeline/scripts/fit_ase_beta_binomial.py \
    --manifest data/rnaseq/sample_metadata.tsv --counts-dir "data/ase/$profile" \
    --output-dir "work/ase/$profile/genewise" \
    --group-min-depth 30 --replicate-min-depth 5 --min-covered-replicates 2 \
    --fdr 0.05 --min-allelic-ratio 2 --diff-min-log2-ratio 1
else
  args=()
  if [[ "$profile" != all ]]; then args+=(--same-sites); fi
  "$py" pipeline/scripts/calibrate_ase_shared_dispersion.py \
    --manifest data/rnaseq/sample_metadata.tsv --counts-dir "data/ase/$profile" \
    --output-dir "work/ase/$profile/shared" --workers "${WORKERS:-4}" \
    --draws 199999 --null-checks 500 --seed 20261013 "${args[@]}"
fi
