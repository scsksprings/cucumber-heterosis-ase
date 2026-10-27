# Analysis workflow

Run commands from the repository root. `results/` is the manuscript snapshot;
new analysis outputs belong in `work/`. The default figure build needs only the
bundled processed results. Read [validation scope](../docs/VALIDATION.md) before
interpreting a reanalysis as identical to the paper.

## Dependencies

Use a separate **Python 3.10** environment for this historical pipeline. The author’s
new figure pins in `requirements-figures.txt` require Python 3.12 / NumPy 2 and belong
to a separate rendering environment; do not overlay them onto this environment.

Python dependencies: `python -m pip install -r pipeline/requirements.txt`.
The tested direct versions are in `requirements-tested.txt` at the repository root.
Expression analysis used R 4.3.1, DESeq2 1.40.2 and data.table 1.16.0.
For the historical R environment, use Bioconductor 3.17 with R 4.3 and confirm
`packageVersion("DESeq2")` before comparing numerical results. The package version
used in the paper matters; installation commands are not run automatically.

Read-level processing additionally requires STAR 2.7.11b, GATK 4.6.2.0,
FastQC 0.12.1, MultiQC 1.15, samtools, bcftools (including `+setGT`), tabix,
gffread, Java and zcat. The vendored WASP modules require pysam and PyTables.
The upstream WASP commit and deterministic duplicate-selection modification are
documented in [vendor/WASP/README.md](vendor/WASP/README.md).

## Expression and PCA from bundled counts

```bash
bash pipeline/run_expression.sh
```

The input is STAR GeneCounts **column 2, unstranded**, stored as
`data/rnaseq/star_gene_counts.tsv`. Models use counts ≥10 in at least three samples.
The script fits within-context `~ genotype` contrasts, a joint blind VST for PCA,
fruit `~ stage + genotype + stage:genotype` interactions, and expression inheritance
classes. PCA selects the 500 most variable genes separately within each displayed
subset of the same joint VST. Automatic checks compare effects, q values, calls,
PCA axes (allowing SVD sign flips), and all qPCR-matched RNA-seq effects with the paper.
The raw class names retain historical terminology; `publication_class` in S3
maps `additive` to `no_detected_logmean_departure`, `high_parent_dominance` to
`high_parent_like`, and corresponding low/above/below categories.

## qPCR from bundled Ct observations

```bash
python pipeline/scripts/reproduce_qpcr.py
```

This reconstructs ΔCt, C-calibrated relative expression, Welch ANOVA, Games–Howell
comparisons and letter assignments for **27 genes / 135 panels / 2,025 biological
observations**. It checks matched qPCR effects against independent Ct-derived effects,
then runs 99,999 whole-gene profile permutations within batch and identical coverage
strata. Serialized effects are retained for the permutation after the Ct check to
preserve exact rank ties. RNA-seq matched effects come from S8; their DESeq2
reconstruction is independently checked by `run_expression.sh`.

The two tests use a two-sided Monte Carlo P value with a +1 correction and a minimum
of 2×10⁻⁵. Their Holm-adjusted P values are 4×10⁻⁵. Effect pairs are not independent
biological replicates. A reduced `--permutations` value is for smoke testing only.

## GO/KEGG from bundled annotation and gene sets

```bash
python pipeline/scripts/reproduce_enrichment.py
```

This independently reconstructs all **11,608 tests in 20 families** from term-to-gene
mapping, selected sets and their testable backgrounds. It includes zero-hit terms,
requires at least five annotated background members per term, and compares memberships,
one-sided hypergeometric P, within-family BH and global BH with the supplied results.
The four shared-direction families use 16,938 jointly testable leaf genes before
database-specific annotation filtering. The input annotation is the supplied
Chinese Long v3 analysis snapshot, not a live annotation database query.

## ASE models from bundled fragment counts

```bash
bash pipeline/run_ase_models.sh common3 genewise
WORKERS=4 bash pipeline/run_ase_models.sh common3 shared
```

Profiles are `all` (cross-specific diagnostic SNPs), `common3` (2,348 shared,
uniquely assigned exonic SNPs) and `common4` (2,273 SNPs). Run either model for
any profile, for example `bash pipeline/run_ase_models.sh common4 genewise`.
`all` differential output from the gene-wise script compares different SNP sets;
**the manuscript's cross-hybrid evidence uses common3/common4 only**.

Each profile contains 30 F1 gene-level count files. The gene-wise model uses group
depth ≥30, at least two replicates with depth ≥5, FDR 0.05, within-hybrid allelic
ratio ≥2 and between-hybrid absolute Δlog2 ratio ≥1. The shared model uses the
frozen script's cross-fitted dispersions and conditional simulation, including
0.5×/1×/2× dispersion sensitivity, 199,999 draws, 500 null checks and seed 20261013.
It can take hours; it is not part of the default build or CI. Results and checkpoints
are isolated by profile/model under `work/ase/`. Reusing a checkpoint after changing
inputs or model parameters requires a new output directory.

All cross-hybrid returned records are supplied as:

| Profile | Gene-wise | Shared |
|---|---:|---:|
| common3 | 3,640 records; 8 calls | 3,632 records; 35 calls |
| common4 | 3,600 records; 7 calls | 3,594 records; 31 calls |

Different eligibility rules explain the row-count differences. The models jointly
support 8 common3 events and 7 common4 events; the eighth is not testable at the
stricter threshold. These are gene × context events.

## Processing sequencing reads

Place the matching clean FASTQs at the paths in `data/rnaseq/sample_metadata.tsv`,
or edit those paths in a separate local manifest. Supply the Chinese Long v3 FASTA
and GFF3. Copy `pipeline/config.yaml` to `pipeline/config.local.yaml` and set
`manifest`, `reference_fasta`, `reference_gff3` and `output_root`; keep `project_root`
pointing at this repository so the code paths resolve. Paths are relative to the
repository root unless absolute. Use paths without whitespace for these legacy
Snakemake shell commands.

The workflows are **staged**, not a single cross-file DAG. Complete each stage in order:

```bash
snakemake -s pipeline/workflows/Snakefile --configfile pipeline/config.local.yaml --cores 24 --resources mem_mb=120000
snakemake -s pipeline/workflows/expression.smk --configfile pipeline/config.local.yaml --cores 4
snakemake -s pipeline/workflows/variants.smk --configfile pipeline/config.local.yaml --cores 24 --resources mem_mb=120000
snakemake -s pipeline/workflows/wasp_ase.smk --configfile pipeline/config.local.yaml --cores 24 --resources mem_mb=120000
```

Add `--dry-run` before execution to check the stage's inputs and planned jobs.
Alignment/QC creates the reference GTF and STAR index; variant calling marks parent
duplicates before SplitNCigarReads and HaplotypeCaller. WASP filters allele-swap
remapping and uses deterministic fragment duplicate removal. Actual hard filters
also include MQRankSum <−12.5 and ReadPosRankSum <−8 when annotations are present;
the full thresholds are in `variants.smk`.

To prepare and recount shared SNPs after those stages, with the default `output_root`:

```bash
mkdir -p work/upstream/inputs
cp data/rnaseq/sample_metadata.tsv work/upstream/inputs/sample_manifest.tsv
python pipeline/scripts/final_ase_postprocess.py prepare --run-root work/upstream
python pipeline/scripts/final_ase_postprocess.py recount --run-root work/upstream --workers 4
```

If a local manifest was used, copy that same manifest instead. Shared-site selection
requires identical coordinates, REF/ALT and maternal/paternal orientation, parental
support and unique exonic gene assignment. Count each read pair once per gene.
Then call the model scripts directly with `--counts-dir work/upstream/common_support3/counts`
(or support4), using the options shown in `run_ase_models.sh`. The model runner by
default always reads the **bundled** counts, so it cannot silently switch to newly
recounted inputs.
