# Data semantics

`data_dictionary.tsv` lists each bundled tabular file, its row/column counts and
column names. Full Excel workbooks preserve the original supplementary schema;
`tables/output/` contains the condensed bilingual submission versions.

- `sample_id`: biological RNA library ID; `library_id`: sequencing-provider ID.
  `genotype` is C, D, S, MD or MS. Each genotype/context has three libraries.
- `stage`: `seedling_leaf`, `fruiting_leaf`, `fruit_3d`, `fruit_6d`, `fruit_9d`.
  Fruit stages are days after anthesis in the updated manuscript. Historical source
  comments saying DAP do not assert that pollination was performed.
- Phenotype `plot` / `block`: randomized block. Weight is recorded in kg; model
  yield and fruit mass are g, fruit number is per plant. The analysis unit is the
  block × genotype summary, with three blocks, not individual harvested fruits.
- Expression `log2FoldChange_midparent` / `D`: deviation from the parental log mean,
  not the arithmetic expression mean. `padj` / `q_BH` is BH adjusted P in the stated
  family. Failure to cross a threshold is not evidence of equal expression.
- Enrichment `background_n`, `selected_n`, `term_background_n`, `overlap_n` are
  respectively M, n, K, k of the hypergeometric test. `q_global_sensitivity` spans
  all 20 families and 11,608 tests, including zero-hit terms.
- ASE `maternal_fragments` and `paternal_fragments`: informative paired-read fragments.
  A fragment is counted once within a gene; ambiguous fragments are separate.
  Maternal fraction is C/(C+paternal). Shared and gene-wise model estimates can
  differ from a simple pooled count ratio.
- `delta_log2_ratio_MS_minus_MD`: log2(C/S) in MS minus log2(C/D) in MD. Positive
  values indicate relatively greater C bias in MS.
- `differential_ase`: statistical/effect-size call for a cross-hybrid event.
  `testable=False`, a missing row, and `differential_ase=False` are distinct states.
  Full gene-wise/shared files have different eligibility denominators and should
  be joined on `gene_id,stage`, never by row order.
- `all`, `common3`, `common4`: cross-specific SNPs, exact common SNPs with ≥3 parent
  library support, and exact common SNPs with ≥4 support. S6's historical
  `TableS20_complete_MLE_ASE` and `TableS21_complete_shared_ASE` are within-hybrid
  **all-profile** tables; use the explicitly named differential files for shared-site
  cross-hybrid inference.
- qPCR `target_ct`, `reference_ct`, `delta_ct`: Ct values and target minus Actin2.
  `relative_expression=2^(-delta_delta_ct)` assumes 100% amplification efficiency;
  C's geometric mean is one within gene/context. `replicate` is biological.
  Pairwise inference uses ΔCt, not untransformed relative expression.
- qPCR effect counts are 540 possible parental effects / 270 log-midparent effects,
  with 524 / 262 finite RNA-seq matches retained. The 24 excluded effects concern
  four gene-contexts across three genes. `S8_qPCR_primer_sequences.tsv` includes
  27 targets plus Actin2 (28 distinct primer target labels).

Many empty numeric cells mean not testable/not available rather than zero. Preserve
these states when importing into spreadsheets. `results/` is the supplied scientific
snapshot; scripts write new calculations to separate output directories.
