# Reproduction notes / 复现与来源说明

## Integration note / 整合说明

本文件对应 `投稿文件_20261026b.zip` 及配套代码：图6包含A–F六个面板，图S4为3页、每页9个基因。
本版另同步图3、图S6的PDF/SVG色标位置及图3、图S2的完整基因编号；图3C色标范围为−4.5～4.5，未截断显示数据。
本文件记录作者新提交的展示层修订；其中“已验证项目”是作者提供的验证记录。
本仓库另已补齐上游代码、STAR计数、ASE片段计数和完整组合间检验。
本地本次验证以 [docs/VALIDATION.md](docs/VALIDATION.md) 为准。
作者图表依赖完整保存在 `requirements-figures.txt`（Python 3.12）；历史统计环境
单独记录在 `requirements-tested.txt`（Python 3.10）。原CSV换行和两个VCF的服务器
路径已做发布规范化，其余97个作者输入文件保持原始字节校验值。

## Scope / 范围

This source tree integrates the reviewed figure typography and table presentation with the user-uploaded repository. The paper title remains:

**同一母本两黄瓜杂交组合产量杂种优势与等位基因特异表达分析**

本次同步的是展示层、运行入口和来源说明。统计输入数值、结果、果实照片/位置及字体文件与上传仓库一致；三项发布规范化差异见上方整合说明。`analysis/phenotype/heterosis_rcbd_lm.R` 的统计计算代码原样保留。RNA-seq 比对、ASE 模型拟合、富集、qPCR 检验等仍使用仓库已有结果；本次未重新执行这些分析。

## The c22 files are generated outputs / c22 文件是生成物

`c22_A.csv`、`c22_M.csv`、`c22_H.csv`、`c22_Mult.csv` 均由原始脚本 **`analysis/phenotype/heterosis_rcbd_lm.R`** 生成。原上传仓库的 `run_all.sh` 已包含这一调用，因此新检出源码目录中不预置这些 CSV 是正常情况，不应判断为缺少研究数据。

此前单独的作图源包曾用冻结结果表的兼容视图替代这些 CSV。本源目录恢复原始生成流程，不包含或依赖 `figures/frozen_phenotype/`。

仓库中的另一个脚本 `analysis/phenotype/heterosis_plot_stats.R` 输出后缀为 `_plot_level.csv`、`_anova.csv`、`_means_letters.csv` 和 `_heterosis_tests.csv`；它不是当前 `run_all.sh` 的 c22 生成入口。两份原始 R 脚本均保留，不混用其输出。

```bash
mkdir -p figures/output/phenotype
Rscript --vanilla analysis/phenotype/heterosis_rcbd_lm.R \
  data/phenotype/block_observations.csv figures/output/phenotype/c22
```

| Input | Original producer | Generated output | Consumer |
|---|---|---|---|
| `data/phenotype/block_observations.csv` (15 block × genotype rows) | `analysis/phenotype/heterosis_rcbd_lm.R` | `c22_A.csv` (9 ANOVA rows) | Figure 1 B–D genotype P labels |
| Same input and producer | Same | `c22_M.csv` (15 adjusted means, SE, n and LSD labels) | Figure 1 B–D means, errors and letter groups |
| Same input and producer | Same | `c22_H.csv` (6 heterosis-test rows) | Reproducible phenotype statistics; Figure 1 does not read this file |
| Same input and producer | Same | `c22_Mult.csv` (2 product-reference rows) | Figure 1 E product-reference decomposition |
| `figures/inputs/fruit_canvas.png` and `fruit_pos.json` | Unmodified reviewed photo/layout | Figure 1 A | `figures/scripts/fig1.py` |
| Frozen `results/S1`–`results/S8` and workbooks | Uploaded processed research outputs | Figures 2–6 and S1–S6 | `figures/scripts/` |
| Same frozen result workbooks/TSVs | `tables/build_tables.py` then `tables/postprocess_tables.py` | Reviewed Tables S1–S8 | `tables/output/` |

`yield_expected_from_component_MPs` 是原始脚本的既有字段名，保留以保证接口兼容；在论文和图中它表示“乘积参照值”，不作为新的统计期望或机制贡献率。

Table S2 and Figure S5 continue to read the uploaded complete decomposition, Fieller intervals and sensitivity results. They are not regenerated from the abbreviated four-file c22 output. Upstream RNA-seq/ASE code and model inputs are included in this integrated release; see `pipeline/README.md` for the independently runnable statistical modules.

## Running / 运行

Author rendering environment: base R >= 4.3, Python 3.12 and `requirements-figures.txt`.
The historical analysis environment uses Python 3.10 with `requirements-tested.txt`.

```bash
python3 -m pip install -r requirements-figures.txt
bash run_all.sh
```

`run_all.sh` runs the unchanged R analysis before Figure 1. It then renders all main and supplementary figures and builds the reviewed supplementary workbooks. Table typography is part of the build entry point rather than a separate unrecorded manual step.

Optional overrides:

- `RSCRIPT`: absolute path to the Rscript executable.
- `PYTHON`: Python executable.
- `FIG_OUT`: absolute figure-output directory, default `figures/output/`.
- `PHENO_OUT`: absolute R-output directory, default `FIG_OUT/phenotype/`; the same path is used by R and Figure 1.
- `MPLCONFIGDIR`: writable Matplotlib cache location, default under `FIG_OUT/`.

`tables/output/` remains the supplementary-workbook destination. The generated figure/table directories are ignored by git. A fresh checkout should run the producer; do not manually reconstruct the c22 files from spreadsheet values.

## Author-provided verification / 作者提供的验证记录

1. The original R script was executed with **base R 4.3.3**, using the unchanged uploaded input records. It produced all four c22 CSV files.
2. Its 9 ANOVA rows, 15 adjusted-mean/SE rows and 2 product-reference rows agree with the formerly substituted frozen inputs to numerical precision. The maximum adjusted-mean difference was 4.55 × 10⁻¹² g; all LSD letters and sample counts are identical. The original ANOVA output also includes `n_obs`.
3. Figure 1 was rendered from these actual R outputs. Its 3504 × 4334 PNG is **pixel-identical** to the reviewed delivered Figure 1; all displayed text is identical. Photographs, colors, coordinates and uncertainty values were not changed.
4. A complete `bash run_all.sh` finished with exit code 0. All 14 final PNGs (11 main/supplementary figures and three S4 pages) are pixel-identical to the delivered figures; all three S4 atlas pages were checked at 600 dpi. The eight rebuilt workbook packages match the independently verified table builds apart from creation/modification timestamps.
5. The table build reproduces all 33 reviewed text-cell edits and 649 rich-text formatted cells. All 72,953 numeric/Boolean cells (72,094 numbers and 859 Boolean values), original identifiers, effective styles and worksheet structures agree with the delivered workbooks. A second table postprocessing pass is byte-idempotent.
6. `validation/input_manifest.json` lists all 100 original data/result/photo/font files and their SHA256 values. All still match the uploaded repository after the complete build.

The R script emits its existing messages about discarding automatically inferred row names. CSVs are written with `row.names = FALSE`; these messages do not change their cell values or column definitions. No statistical code was altered to suppress them.

Detailed checks are in `validation/phenotype_R_check.json`, `validation/figure1_R_check.json`, `validation/full_build_check.json` and `validation/table_build_check.json`. Software versions and the unchanged R-script hash are in `validation/build_environment.json`.

## Preserved conventions / 保留的展示规则

The figure scripts match the final reviewed figure source: 苗期叶片/结瓜期叶片; 单株瓜数/单瓜重; F₁ and mathematical subscripts; italic gene symbols and IDs when used as gene labels; upright descriptive suffixes, proteins and pathway annotations; and the Figure 4 term “Parent-informative SNPs”. Historical input names such as `context`, `diagnostic`, `concordance` or `yield_expected_from_component_MPs` remain unchanged where they are database fields/file interfaces.

## Repository status / 仓库状态

The original author archive contained no Git history. The integrated repository uses the `main` branch at [https://github.com/scsksprings/cucumber-heterosis-ase](https://github.com/scsksprings/cucumber-heterosis-ase). The reports in `docs/validation/` document local checks completed before publication; live hosted checks are available under [GitHub Actions](https://github.com/scsksprings/cucumber-heterosis-ase/actions).
