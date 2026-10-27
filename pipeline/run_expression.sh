#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p work/expression/inheritance
Rscript --vanilla pipeline/scripts/run_deseq2_expression.R \
  --counts data/rnaseq/star_gene_counts.tsv --metadata data/rnaseq/sample_metadata.tsv \
  --output-dir work/expression/deseq2
Rscript --vanilla pipeline/scripts/run_fruit_stage_interaction.R \
  --counts data/rnaseq/star_gene_counts.tsv --metadata data/rnaseq/sample_metadata.tsv \
  --output work/expression/fruit_stage_interactions.tsv \
  --summary work/expression/fruit_stage_summary.tsv
"${PYTHON:-python3}" pipeline/scripts/classify_expression_inheritance.py \
  --contrasts work/expression/deseq2/deseq2_all_contrasts.tsv \
  --output work/expression/inheritance/gene_results.tsv \
  --summary work/expression/inheritance/summary.tsv
"${PYTHON:-python3}" pipeline/scripts/reproduce_pca.py
"${PYTHON:-python3}" scripts/validate_expression.py
