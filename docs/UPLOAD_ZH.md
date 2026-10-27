# 使用与上传 GitHub

本仓库对应作者提供的 **20261026b修改稿及配套的新版代码包**。qPCR采用论文中的27个基因、
STAR计数、524/262对效应和ρ=0.522/0.451，不混入另一套29基因的百迈客counts结果。

## 本次整理内容

- 同步b版图3/图S6的PDF与SVG色标定位修正、图3/图S2完整基因编号，以及图3C色标范围说明。

- 补齐原占位 `pipeline/`：STAR/QC、GATK/WASP、DESeq2、共同SNP计数、两种ASE模型，
  以及可独立运行的qPCR与富集复算模块。
- 增加24,317基因×75文库的STAR计数矩阵和3种SNP范围各30个F1文库的片段计数。
- 补齐两个ASE模型、两档亲本支持门槛的全量组合间检验结果，包含不显著事件。
- 同步图6：四个基因的qRT-PCR面板A–D和两幅相关性图E–F；图S4改为3页，每页9个基因，提供600 dpi PNG、SVG及27个基因面板PDF。
- 合并新版图表术语、数学符号/基因名排版、SVG字体回退、图S4逐页SVG和27个基因分面板、33项表格文字修订及649个富文本单元格；
- 保留作者修改后的图表布局和科学结果；重生成S1–S8与投稿版本逐单元格一致。
- 增加引用元数据、依赖版本、自动检查、数据字典、来源追踪和SHA256校验清单。
- 新运算输出集中在 `work/`、`figures/output/`、`tables/output/`，由Git忽略。

详细验证结果见 [VALIDATION.md](VALIDATION.md)。原ZIP和投稿Word/PDF保留在项目
`artical/` 下，论文正文不作为代码仓库的上传内容。

## 在本地复现

重绘最新版图表使用Python 3.12与R ≥4.3，在仓库根目录运行：

```bash
python3.12 -m venv .venv-figures
source .venv-figures/bin/activate
python -m pip install -r requirements-figures.txt
bash run_all.sh
```

上述命令重算图1表型统计，并从完整结果表重绘全部主图、附图和8份精简附表。
qPCR、富集和表达模型的独立复算命令见 [pipeline/README.md](../pipeline/README.md)。
新版图表的固定依赖见`requirements-figures.txt`；统计分析单独使用Python 3.10和
`requirements-tested.txt`，避免NumPy 1与2的依赖冲突。
全量reads重比对、变异检测和共享离散度大规模模拟不属于默认图表构建。

## 上传前本地检查

```bash
python scripts/validate_release.py
```

该命令检查样本数、输入计数、qPCR口径、完整ASE表和文件校验和。
完整报告见 `docs/validation/`。日后修改受校验文件后，先审查修改，再运行
`python scripts/build_release.py` 更新校验清单与发布压缩包。

## GitHub仓库与后续更新

正式仓库：[scsksprings/cucumber-heterosis-ase](https://github.com/scsksprings/cucumber-heterosis-ase)。
默认分支为 `main`；SSH推送地址为 `git@github.com:scsksprings/cucumber-heterosis-ase.git`。

下载源码：

```bash
git clone https://github.com/scsksprings/cucumber-heterosis-ase.git
cd cucumber-heterosis-ase
```

维护者修改代码或数据后，应完成相应验证，再更新校验清单、提交和推送：

```bash
python scripts/build_release.py
python scripts/validate_release.py
git add .
git commit -m "Update manuscript code and documentation"
git push origin main
```

仅在本仓库目录执行上述命令。图表输出和原始测序数据不纳入Git提交。
若修改作图代码，应先使用 `requirements-figures.txt` 对应环境运行 `bash run_all.sh`；
如已有旧版校验清单，先用 `python scripts/build_release.py` 更新清单后再构建。

推送后的运行状态见 [GitHub Actions](https://github.com/scsksprings/cucumber-heterosis-ase/actions)。
`docs/validation/` 中的报告记录发布前的本地验证，云端执行状态以Actions页面为准。
自动检查包含ASE测试、qPCR、富集和图表构建，不启动原始reads分析。
正式发表后可补充最终文章DOI和GSA的CRA编号，并为对应版本建立Release。
