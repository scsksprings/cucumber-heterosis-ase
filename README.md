# Cucumber yield heterosis and allele-specific expression

Code, analysis inputs and complete processed results accompanying the **2026-10-26b manuscript revision (`投稿文件_20261026b.zip`) and the matching updated code archive**:

> 杨宗辉，李利斌，孟昭娟，祝利霞，曹齐卫. 同一母本两黄瓜杂交组合产量杂种优势与等位基因特异表达分析.
> Yang Z, Li L, Meng Z, Zhu L, Cao Q. *Yield heterosis and allele-specific expression in two cucumber hybrids sharing a maternal parent.*

Parents: MC7 (**C**, maternal), DRTJY2 (**D**) and SJ11 (**S**).
Hybrids: MC7×DRTJY2 (**MD**) and MC7×SJ11 (**MS**).
The study includes 75 RNA-seq libraries: five genotypes × five contexts × three biological replicates.
BioProject: **PRJCA077972**. See [data availability](docs/DATA_AVAILABILITY.md).

[GitHub repository](https://github.com/scsksprings/cucumber-heterosis-ase) · [GitHub Actions](https://github.com/scsksprings/cucumber-heterosis-ase/actions)

[中文使用与上传说明](docs/UPLOAD_ZH.md) · [Analysis workflow](pipeline/README.md) · [Validation](docs/VALIDATION.md) · [File/column inventory](docs/data_dictionary.tsv)

## Rebuild the manuscript figures and tables

Use Linux, macOS, or WSL with **Python 3.12** and **R ≥4.3** (base R is sufficient for this command). The exact author rendering dependencies are separate from the historical analysis environment.

```bash
python3.12 -m venv .venv-figures
source .venv-figures/bin/activate
python -m pip install -r requirements-figures.txt
bash run_all.sh
```

This checks the supplied data, recomputes the randomized-block phenotype statistics,
rebuilds **Figs. 1–6 and S1–S6**, and writes **eight condensed supplementary workbooks**.
Outputs are `figures/output/` and `tables/output/`; the complete published tables remain in `results/`.
Fig. 6 contains four selected-gene qRT-PCR profiles (A–D) and two RNA-seq/qPCR correlations (E–F). Fig. S4 is a three-page atlas of 27 genes (nine genes × five contexts per page), with three 600-dpi PNGs, three editable SVGs and 27 individual gene-panel PDFs. Figs. 3 and S2 retain full `CsaV3_` gene identifiers. The color bars in Figs. 3 and S6 use the revised PDF/SVG positioning; the Fig. 3C color scale spans −4.5 to 4.5 and contains all displayed effects. Main and supplementary figures also export per-panel PDFs. SVGs declare CJK font fallbacks. The full build was exercised locally; PDF binary hashes
can vary with timestamps, font libraries and rendering versions.

For the **Python 3.10 statistical-analysis environment**, install
`requirements-tested.txt` in a separate environment. It contains dependencies for the statistical modules and tests. The author figure pins use NumPy 2; the frozen statistical environment uses NumPy 1.26. Do not install both pinned environments together.
It records the tested environment, rather than a platform-independent lock of every transitive dependency.

## Recompute the statistics

```bash
# In a separate Python 3.10 analysis environment:
python -m pip install -r requirements-tested.txt
python pipeline/scripts/reproduce_qpcr.py
python pipeline/scripts/reproduce_enrichment.py
# Requires R packages DESeq2 1.40.2 and data.table:
bash pipeline/run_expression.sh
# Gene-wise ASE from the supplied common-SNP fragment counts:
bash pipeline/run_ase_models.sh common3 genewise
# Shared-dispersion simulation; substantially more CPU time:
WORKERS=4 bash pipeline/run_ase_models.sh common3 shared
```

Results are written under `work/`. The [pipeline guide](pipeline/README.md) describes all-SNP
and stricter-support profiles, external sequencing inputs and the staged Snakemake workflows.
The packaged qPCR module uses the **27-gene STAR-based manuscript analysis**: 524 parental-effect
pairs (Spearman ρ=0.521900) and 262 log-midparent pairs (ρ=0.450651), with 99,999
stratified gene-block permutations. The separate 29-gene BMK-count analysis is a different version.

## Repository contents

| Path | Contents |
|---|---|
| `analysis/phenotype/` | Base-R randomized-block ANOVA, LSD letters and two-sided heterosis contrasts |
| `data/phenotype/` | Fifteen block × genotype records |
| `data/rnaseq/` | STAR unstranded counts for 24,317 genes × 75 libraries; sample metadata |
| `data/ase/` | Fragment counts for 30 F1 libraries in each of three SNP profiles |
| `pipeline/scripts/` | DESeq2, inheritance classification, PCA, SNP selection, fragment counting, both ASE models, qPCR and enrichment |
| `pipeline/workflows/` | Staged STAR/QC, expression, GATK and WASP Snakemake workflows |
| `pipeline/vendor/WASP/` | Five mapping modules from the frozen analysis, with upstream license |
| `results/S1.xlsx`–`S8.xlsx` | Complete manuscript result workbooks; these are the inputs to the condensed table builder |
| `results/S1/`–`S8/` | Machine-readable complete tests and supporting data |
| `figures/scripts/` | Main and supplementary figure scripts |
| `figures/inputs/`, `figures/fonts/` | Supplied fruit-image canvas and licensed Chinese fonts |
| `tables/` | Condensed submission-table builder, bilingual labels and lossless rich-text/terminology postprocessing |
| `validation/`, `REPRODUCTION.md` | Author-supplied presentation validation and the original R-to-figure input/output chain |
| `scripts/`, `tests/` | Release/data checks, expression checks and scientific regression tests |
| `docs/` | Data dictionary, provenance, validation and upload guide |

## Results map

| Table | Content |
|---|---|
| S1 | Sample manifest, alignment QC, PCA scores and selected genes |
| S2 | Phenotypes, ANOVA, heterosis, sensitivity analyses, Fieller intervals, harvests and product decomposition |
| S3 | Complete expression departures, leaf overlap and enrichment tests, backgrounds and annotations |
| S4 | Complete fruit-stage interaction tests and enrichment |
| S5 | Parent-informative/shared SNP selection and common SNP VCFs at parent support ≥3 and ≥4 |
| S6 | Within-hybrid ASE, complete cross-hybrid tests for both models at both support thresholds, dispersion and 40 null-diagnostic families |
| S7 | Eight jointly supported ASE candidates, site consistency, total expression and representative fragment counts |
| S8 | 27-gene Ct values, primers, Welch/Games–Howell tests, matched effects and exclusions |

The four `results/S6/differential_ASE_*_all.tsv.gz` files contain **all returned cross-hybrid
records**, including nonsignificant results. The gene-wise and shared models have different
eligibility rules; missing rows are not equivalent to nonsignificant tests. See [data semantics](docs/DATA_DICTIONARY.md).

Sample types: SL, seedling leaf; FL, fruiting-stage leaf; F3/F6/F9, fruits at 3/6/9 days after anthesis
(parthenocarpic fruits, without pollination). Directory names such as `fruit_3d` are analysis labels.
Expression departure is `D = log2(μF1) − [log2(μC) + log2(μP)]/2`.
ASE fractions refer to maternal C; cross-hybrid effects are MS minus MD.

## Validation and scope

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests -q
python scripts/validate_release.py
# After bash run_all.sh:
python scripts/validate_release.py --built
```

The release records local reconstruction and comparison in [VALIDATION.md](docs/VALIDATION.md).
The latest table builder applies 33 text-cell updates and 649 rich-text rules while preserving numerical values. The original R script generates all four `c22_*.csv` intermediates before Figure 1; they are not missing source data.
`PYTHON`, `RSCRIPT`, `FIG_OUT` and `PHENO_OUT` are supported overrides (use absolute output paths).
The full 227 GB read alignment/variant-calling analysis and full shared-dispersion simulations
were not rerun during repository preparation. Their source and processed inputs/results are supplied.
The full S2 Fieller/sensitivity/decomposition results remain a validated publication snapshot;
the base-R entry point recomputes the principal phenotype statistics used by Fig. 1.

## Citation and licenses

Use [CITATION.cff](CITATION.cff) for the software authors and manuscript title.
No publication DOI or final journal citation has been assigned in this repository.
Original project code is MIT licensed ([LICENSE](LICENSE)). Vendored WASP code uses Apache 2.0;
fonts use SIL OFL 1.1. See [third-party notices](docs/THIRD_PARTY.md).
