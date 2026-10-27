# Validation of the prepared 2026-10-26 repository

The author source and manuscript were refreshed at **2026-10-26b (`投稿文件_20261026b.zip` and `GitHub仓库_20261026b.zip`; code received at 2026-10-26 16:50 UTC / 2026-10-27 00:50 Asia/Shanghai)**.
The figure/table presentation code was updated; all statistical inputs and original
R/statistical model code remain unchanged. Records below distinguish the prior
statistical validation (still applicable by hashes) from the refreshed presentation
checks. `validation/` contains author-provided reports; `docs/validation/` contains
checks performed on this integrated repository. No remote GitHub Actions run is claimed.

Current update/build status is recorded in `docs/validation/update_20261026b.json`
(completed successfully). Full read processing and simulations are
not rerun for this presentation-only update.

| Check | Observed result |
|---|---|
| Isolated ZIP build | Prepared archive extracted into a new directory outside the project; `bash run_all.sh` completed including checksum/data/output checks |
| Refreshed figure/table reconstruction | Six main figures, six supplementary figures, three S4 page SVGs and PNGs, 27 gene-panel PDFs and eight condensed workbooks; the changed figure scripts and independent ZIP full build both passed |
| Refreshed submission workbooks | Compare every S1–S8 cell, effective style and rich-text run against the new submission; detailed totals are recorded in `validation/presentation_update.json`  |
| STAR input | 24,317 genes × 75 libraries; unique IDs, integer nonnegative counts, 25 genotype/context groups with n=3 |
| DESeq2/inheritance | **177,202 gene × genotype × context rows** compared with the frozen snapshot; 5,179 expression departure calls retained |
| Fruit interactions | **110,994 gene × contrast records** compared with the frozen snapshot |
| PCA | Both subsets and their top-500 gene sets reproduced; axis sign flips allowed |
| qPCR-matched RNA-seq effects | All **810 potential effects** checked against newly fitted DESeq2 results, including unavailable values |
| qPCR from Ct | 27 genes, 135 panels, 2,025 biological observations; Welch/Games–Howell results and letters matched |
| qPCR permutations | 99,999 gene-block permutations; 524/262 matched effects; ρ=0.5218998254/0.4506512454; both Holm P=0.00004 |
| Enrichment | All **11,608 tests / 20 families**, memberships, P, local BH and global BH reproduced; 83 local-q significant terms across the complete family set |
| Complete differential ASE exports | All four tables matched TableS13 q values and calls by gene/context; common3 jointly supports 8 events and common4 supports 7 |
| Alignment/QC dry run | Full alignment/QC DAG resolved using actual external read/reference paths and a new output directory; no jobs executed |
| Workflow parsing | All four Snakemake entry points parsed and listed their rules successfully (7.32.3); Python, shell and YAML syntax checks passed |
| ASE scientific regression suite | **9 passed**: direction, balanced alleles, differential effects, dispersion recovery, conditional simulation, allele/group symmetry, ambiguous exons and common-site orientation/support |

The statistical checks listed above were completed for the prior prepared release and their input hashes were confirmed unchanged for this update. New table postprocessing changes 33 text cells and applies 649 rich-text rules, without changing statistical values.

The executable checks are `scripts/validate_release.py`, `scripts/validate_expression.py`,
`pipeline/scripts/reproduce_qpcr.py`, `pipeline/scripts/reproduce_enrichment.py` and `tests/`.
Machine-readable records are in [validation/](validation/).

Statistical environment: Python 3.10.12; R 4.3.1; DESeq2 1.40.2; data.table 1.16.0.
Rendering environment: Python 3.12.12 / base R 4.3.1, using the exact direct packages in
`requirements-figures.txt`; the author used Python 3.12.14 / R 4.3.3. Any PDF pixel differences are reported explicitly.
Direct Python package versions are recorded in `requirements-tested.txt`.

Scientific data in the supplied result workbooks and TSVs were preserved. The phenotype CSV line endings were normalized to LF without changing values. Only
server-specific path prefixes in VCF history headers were normalized; all VCF
variant/genotype records were checked unchanged. The included main/figure scripts
retain the supplied manuscript presentation. PDF bytes are not asserted identical
because timestamps and font rendering may differ.

Scope limits: the raw-read alignment/GATK/WASP pipeline and full 199,999-draw shared
ASE simulation were not rerun during packaging. The frozen statistical code,
fragment inputs and full outputs are supplied, with focused scientific tests.
The full S2 sensitivity/Fieller/decomposition output remains a publication snapshot;
the base-R scripts reproduce the principal phenotype analysis, not every S2 sheet.

Current presentation verification matched **114,155 cells**, including **72,094 numbers**, **859 Boolean values**, effective styles and **649 rich-text cells**. The unchanged table postprocessor was verified byte-idempotent in the preceding revision. The source ZIP also completed an independent full build and reproduced all 14 PNGs exactly.

All **14 PNGs** (including the three 600-dpi S4 pages) are pixel-identical to the new submission at their original resolution. All **12 main/supplementary PDFs** are pixel-identical when rendered at 72 dpi. Figure 6 has six panel PDFs; S4 contains three pages and 27 individual gene-panel PDFs.

The preceding revision’s PAL/WRKY33 log-midparent values, Games–Howell comparisons, four-PAL ranges and 9-day-fruit expression directions were checked against S8. Details are in `validation/update_20261026a.json`. Statistical source/input hashes are unchanged, so the previously completed statistical reconstruction remains applicable.

The b revision corrects color-bar placement in Fig. 3/Fig. S6 PDF and SVG exports, restores full gene identifiers in Figs. 3/S2, and describes Fig. 3C as a −4.5 to 4.5 color range. Numerical plotting tables are compared with the previous revision.

Current b-revision vector checks: all three revised SVG exports match the submission after normalizing metadata and internal IDs; all five affected panel PDFs match at 72 dpi. The Fig. 3B, Fig. 3C and Fig. S6 plotting tables are exactly unchanged. Fig. 3C contains 138 effects across 23 genes; maximum |effect| = 4.333042647, so no displayed value exceeds the −4.5 to 4.5 range. Details are in `validation/update_20261026b.json`.
