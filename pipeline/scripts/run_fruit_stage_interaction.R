#!/usr/bin/env Rscript

# Formally test whether F1 departure from the parental log-expression mean
# changes among fruit 3, 6 and 9 DAP.  This avoids inferring an interaction
# merely because one stage is significant and another is not.

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop("Missing argument: ", flag)
  args[[idx + 1]]
}

count_path <- get_arg("--counts")
metadata_path <- get_arg("--metadata")
output_path <- get_arg("--output")
summary_path <- get_arg("--summary")

counts <- fread(count_path)
gene_ids <- counts[[1]]
counts <- as.matrix(counts[, -1])
rownames(counts) <- gene_ids
metadata <- as.data.frame(fread(metadata_path))
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(counts), , drop = FALSE]

fruit_stages <- c("fruit_3d", "fruit_6d", "fruit_9d")
crosses <- list(MD = "D", MS = "S")
all_results <- list()

tidy_result <- function(result, hybrid, paternal, comparison) {
  frame <- as.data.frame(result)
  frame$gene_id <- rownames(frame)
  frame$hybrid <- hybrid
  frame$paternal_parent <- paternal
  frame$comparison <- comparison
  frame[, c(
    "gene_id", "hybrid", "paternal_parent", "comparison", "baseMean",
    "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  )]
}

for (hybrid in names(crosses)) {
  paternal <- crosses[[hybrid]]
  keep_samples <- metadata$stage %in% fruit_stages &
    metadata$genotype %in% c("C", paternal, hybrid)
  meta <- metadata[keep_samples, , drop = FALSE]
  meta$stage <- relevel(factor(meta$stage, levels = fruit_stages), ref = "fruit_3d")
  meta$genotype <- relevel(
    factor(meta$genotype, levels = c("C", paternal, hybrid)), ref = "C"
  )
  matrix <- counts[, rownames(meta), drop = FALSE]
  keep_genes <- rowSums(matrix >= 10) >= 3
  dds <- DESeqDataSetFromMatrix(
    countData = matrix[keep_genes, , drop = FALSE],
    colData = meta,
    design = ~ stage + genotype + stage:genotype
  )
  dds <- DESeq(dds, quiet = TRUE)
  coefficient_names <- resultsNames(dds)

  interaction_name <- function(genotype, stage) {
    candidates <- c(
      paste0("stage", stage, ".genotype", genotype),
      paste0("genotype", genotype, ".stage", stage)
    )
    found <- candidates[candidates %in% coefficient_names]
    if (!length(found)) {
      stop(
        "Cannot locate interaction coefficient for ", genotype, " at ", stage,
        ". Available: ", paste(coefficient_names, collapse = ",")
      )
    }
    found[[1]]
  }

  departure_change_vector <- function(stage) {
    vector <- setNames(rep(0, length(coefficient_names)), coefficient_names)
    vector[interaction_name(hybrid, stage)] <- 1
    vector[interaction_name(paternal, stage)] <- -0.5
    vector
  }

  change_6_vs_3 <- departure_change_vector("fruit_6d")
  change_9_vs_3 <- departure_change_vector("fruit_9d")
  change_9_vs_6 <- change_9_vs_3 - change_6_vs_3
  contrasts <- list(
    fruit_6d_vs_fruit_3d = change_6_vs_3,
    fruit_9d_vs_fruit_3d = change_9_vs_3,
    fruit_9d_vs_fruit_6d = change_9_vs_6
  )
  for (comparison in names(contrasts)) {
    result <- results(
      dds,
      contrast = contrasts[[comparison]],
      independentFiltering = TRUE
    )
    all_results[[paste(hybrid, comparison, sep = "__")]] <-
      tidy_result(result, hybrid, paternal, comparison)
  }
}

result_table <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(result_table, output_path, sep = "\t", na = "NA")

summary <- result_table[, .(
  tested = sum(!is.na(padj)),
  significant_interactions = sum(
    padj < 0.05 & abs(log2FoldChange) >= 1, na.rm = TRUE
  ),
  positive = sum(padj < 0.05 & log2FoldChange >= 1, na.rm = TRUE),
  negative = sum(padj < 0.05 & log2FoldChange <= -1, na.rm = TRUE)
), by = .(hybrid, paternal_parent, comparison)]
fwrite(summary, summary_path, sep = "\t", na = "NA")
