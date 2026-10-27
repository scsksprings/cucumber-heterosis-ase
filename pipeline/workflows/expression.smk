from pathlib import Path
import pandas as pd

configfile: "pipeline/config.yaml"

ROOT = Path(config["project_root"])
OUT = Path(config["output_root"])
MANIFEST = pd.read_csv(config["manifest"], sep="\t")
SAMPLES = MANIFEST.sample_id.tolist()


rule all:
    input:
        OUT / "expression/count_matrix.tsv",
        OUT / "expression/deseq2/deseq2_all_contrasts.tsv",
        OUT / "expression/deseq2/deseq2_normalized_counts.tsv",
        OUT / "expression/deseq2/vst_blind_all_samples.tsv",
        OUT / "expression/deseq2/deseq2_contrast_summary.tsv"


rule count_matrix:
    input:
        script=ROOT / "pipeline/scripts/build_star_count_matrix.py",
        counts=expand(str(OUT / "alignment/{sample}/{sample}.ReadsPerGene.out.tab"), sample=SAMPLES),
        manifest=config["manifest"]
    output:
        matrix=OUT / "expression/count_matrix.tsv",
        metadata=OUT / "expression/sample_metadata.tsv"
    params:
        alignment=OUT / "alignment"
    log:
        OUT / "logs/expression/count_matrix.log"
    shell:
        "mkdir -p $(dirname {output.matrix}) $(dirname {log}) && "
        "python {input.script} --manifest {input.manifest} --alignment-dir {params.alignment} "
        "--output-matrix {output.matrix} --output-metadata {output.metadata} > {log} 2>&1"


rule deseq2_expression:
    input:
        script=ROOT / "pipeline/scripts/run_deseq2_expression.R",
        matrix=OUT / "expression/count_matrix.tsv",
        metadata=OUT / "expression/sample_metadata.tsv"
    output:
        contrasts=OUT / "expression/deseq2/deseq2_all_contrasts.tsv",
        normalized=OUT / "expression/deseq2/deseq2_normalized_counts.tsv",
        vst=OUT / "expression/deseq2/vst_blind_all_samples.tsv",
        summary=OUT / "expression/deseq2/deseq2_contrast_summary.tsv"
    params:
        outdir=OUT / "expression/deseq2"
    log:
        OUT / "logs/expression/deseq2.log"
    shell:
        "mkdir -p {params.outdir} $(dirname {log}) && "
        "Rscript {input.script} --counts {input.matrix} --metadata {input.metadata} "
        "--output-dir {params.outdir} > {log} 2>&1"
