#!/usr/bin/env Rscript

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
output_dir <- get_arg("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- fread(count_path)
gene_ids <- counts[[1]]
counts <- as.matrix(counts[, -1])
rownames(counts) <- gene_ids
metadata <- as.data.frame(fread(metadata_path))
rownames(metadata) <- metadata$sample_id
metadata <- metadata[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(metadata), colnames(counts)))

stages <- c("seedling_leaf", "fruiting_leaf", "fruit_3d", "fruit_6d", "fruit_9d")
all_results <- list()
normalized <- matrix(NA_real_, nrow = nrow(counts), ncol = ncol(counts),
                     dimnames = dimnames(counts))

# A single transformation across all 75 libraries is used only for global PCA
# and sample-level QC. Stage-specific DE tests below retain stage-specific size
# factors so that developmental shifts do not dominate within-stage contrasts.
global_keep <- rowSums(counts >= 10) >= 3
global_meta <- metadata
global_meta$stage <- factor(global_meta$stage)
global_meta$genotype <- factor(global_meta$genotype)
global_dds <- DESeqDataSetFromMatrix(counts[global_keep, , drop = FALSE], global_meta,
                                    ~ stage + genotype)
global_vst <- assay(vst(global_dds, blind = TRUE))

tidy_result <- function(result, stage, comparison) {
  frame <- as.data.frame(result)
  frame$gene_id <- rownames(frame)
  frame$stage <- stage
  frame$comparison <- comparison
  frame[, c("gene_id", "stage", "comparison", "baseMean", "log2FoldChange",
            "lfcSE", "stat", "pvalue", "padj")]
}

for (stage_name in stages) {
  samples <- metadata$sample_id[metadata$stage == stage_name]
  stage_counts <- counts[, samples, drop = FALSE]
  stage_meta <- metadata[samples, , drop = FALSE]
  keep <- rowSums(stage_counts >= 10) >= 3
  dds <- DESeqDataSetFromMatrix(stage_counts[keep, , drop = FALSE], stage_meta, ~ genotype)
  dds$genotype <- relevel(factor(dds$genotype), ref = "C")
  design(dds) <- ~ genotype
  dds <- DESeq(dds, quiet = TRUE)
  normalized[rownames(dds), samples] <- counts(dds, normalized = TRUE)

  simple <- list(
    D_vs_C = c("genotype", "D", "C"),
    S_vs_C = c("genotype", "S", "C"),
    MD_vs_C = c("genotype", "MD", "C"),
    MD_vs_D = c("genotype", "MD", "D"),
    MS_vs_C = c("genotype", "MS", "C"),
    MS_vs_S = c("genotype", "MS", "S")
  )
  for (comparison in names(simple)) {
    all_results[[paste(stage_name, comparison, sep = "__")]] <-
      tidy_result(results(dds, contrast = simple[[comparison]], independentFiltering = TRUE),
                  stage_name, comparison)
  }

  coefficient_names <- resultsNames(dds)
  contrast_md <- setNames(rep(0, length(coefficient_names)), coefficient_names)
  contrast_md["genotype_MD_vs_C"] <- 1
  contrast_md["genotype_D_vs_C"] <- -0.5
  all_results[[paste(stage_name, "MD_vs_midparent", sep = "__")]] <-
    tidy_result(results(dds, contrast = contrast_md, independentFiltering = TRUE),
                stage_name, "MD_vs_midparent")

  contrast_ms <- setNames(rep(0, length(coefficient_names)), coefficient_names)
  contrast_ms["genotype_MS_vs_C"] <- 1
  contrast_ms["genotype_S_vs_C"] <- -0.5
  all_results[[paste(stage_name, "MS_vs_midparent", sep = "__")]] <-
    tidy_result(results(dds, contrast = contrast_ms, independentFiltering = TRUE),
                stage_name, "MS_vs_midparent")
}

result_table <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
fwrite(result_table, file.path(output_dir, "deseq2_all_contrasts.tsv"), sep = "\t", na = "NA")
fwrite(data.table(gene_id = rownames(normalized), normalized),
       file.path(output_dir, "deseq2_normalized_counts.tsv"), sep = "\t", na = "NA")
fwrite(data.table(gene_id = rownames(global_vst), global_vst),
       file.path(output_dir, "vst_blind_all_samples.tsv"), sep = "\t", na = "NA")

summary <- result_table[, .(
  tested = sum(!is.na(padj)),
  significant = sum(padj < 0.05 & abs(log2FoldChange) >= 1, na.rm = TRUE),
  up = sum(padj < 0.05 & log2FoldChange >= 1, na.rm = TRUE),
  down = sum(padj < 0.05 & log2FoldChange <= -1, na.rm = TRUE)
), by = .(stage, comparison)]
fwrite(summary, file.path(output_dir, "deseq2_contrast_summary.tsv"), sep = "\t")
