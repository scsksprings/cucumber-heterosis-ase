from pathlib import Path
import pandas as pd

configfile: "pipeline/config.yaml"

ROOT = Path(config["project_root"])
OUT = Path(config["output_root"])
MANIFEST = pd.read_csv(config["manifest"], sep="\t").set_index("sample_id")
F1 = MANIFEST[MANIFEST["genotype"].isin(["MD", "MS"])].index.tolist()
REF = OUT / "reference/ChineseLong_genome_v3.fa"
DICT = OUT / "reference/ChineseLong_genome_v3.dict"
WASP = ROOT / "pipeline/vendor/WASP/mapping"


def cross_for(sample):
    return "CD" if MANIFEST.loc[sample, "genotype"] == "MD" else "CS"


rule all:
    input:
        expand(str(OUT / "ase/site_counts/{sample}.ase.table"), sample=F1),
        expand(str(OUT / "ase/gene_counts/{sample}.gene_allele_counts.tsv.gz"), sample=F1),
        expand(str(OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam.bai"), sample=F1),
        OUT / "ase/statistics/ase_gene_results.tsv",
        OUT / "ase/statistics/differential_ase_MD_vs_MS.tsv",
        OUT / "ase/statistics/ase_stage_summary.tsv",
        OUT / "ase/statistics/ase_stage_patterns.tsv",
        OUT / "ase/statistics/ase_stage_intersection_summary.tsv",
        OUT / "ase/statistics/ase_direction_tests.tsv"


rule wasp_snp_text:
    input:
        script=ROOT / "pipeline/scripts/vcf_to_wasp_snp_dir.py",
        vcf=lambda w: str(OUT / "variants/informative" / f"{w.cross}.informative_snps.vcf.gz"),
        tbi=lambda w: str(OUT / "variants/informative" / f"{w.cross}.informative_snps.vcf.gz.tbi")
    output:
        done=OUT / "ase/wasp_snps/{cross}/.done",
        summary=OUT / "ase/wasp_snps/{cross}/conversion_summary.tsv"
    params:
        outdir=lambda w: str(OUT / "ase/wasp_snps" / w.cross),
    log:
        OUT / "logs/ase/wasp_snps/{cross}.log"
    shell:
        "python {input.script} --vcf {input.vcf} --output-dir {params.outdir} > {log} 2>&1"


rule filter_unique_proper_pairs:
    input:
        bam=OUT / "alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai=OUT / "alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        bam=OUT / "ase/wasp/{sample}/{sample}.q30.bam",
        bai=OUT / "ase/wasp/{sample}/{sample}.q30.bam.bai"
    threads: 3
    log:
        OUT / "logs/ase/wasp/{sample}.q30.log"
    shell:
        "mkdir -p $(dirname {output.bam}) $(dirname {log}) && "
        "samtools view -@ {threads} -b -q 30 -f 2 -F 2308 {input.bam} -o {output.bam} 2> {log} && "
        "samtools quickcheck -v {output.bam} && samtools index -@ {threads} {output.bam}"


rule find_intersecting_snps:
    input:
        script=WASP / "find_intersecting_snps.py",
        bam=OUT / "ase/wasp/{sample}/{sample}.q30.bam",
        bai=OUT / "ase/wasp/{sample}/{sample}.q30.bam.bai",
        snps=lambda w: str(OUT / "ase/wasp_snps" / cross_for(w.sample) / ".done")
    output:
        keep=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.keep.bam",
        remap=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.to.remap.bam",
        fq1=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.remap.fq1.gz",
        fq2=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.remap.fq2.gz"
    params:
        outdir=lambda w: str(OUT / "ase/wasp" / w.sample / "intersect"),
        snpdir=lambda w: str(OUT / "ase/wasp_snps" / cross_for(w.sample))
    log:
        OUT / "logs/ase/wasp/{sample}.find_intersecting.log"
    shell:
        "mkdir -p {params.outdir} $(dirname {log}) && "
        "python {input.script} --is_paired_end --is_sorted --snp_dir {params.snpdir} "
        "--output_dir {params.outdir} {input.bam} > {log} 2>&1"


rule remap_swapped_reads:
    input:
        fq1=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.remap.fq1.gz",
        fq2=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.remap.fq2.gz",
        index=OUT / "reference/star_index/SA"
    output:
        bam=OUT / "ase/wasp/{sample}/remap/{sample}.Aligned.out.bam"
    threads: 8
    params:
        index=OUT / "reference/star_index",
        prefix=lambda w: str(OUT / "ase/wasp" / w.sample / "remap" / f"{w.sample}.")
    log:
        OUT / "logs/ase/wasp/{sample}.remap.log"
    shell:
        "mkdir -p $(dirname {output.bam}) $(dirname {log}) && "
        "STAR --runThreadN {threads} --genomeDir {params.index} "
        "--readFilesIn {input.fq1} {input.fq2} --readFilesCommand zcat "
        "--outFileNamePrefix {params.prefix} --outSAMtype BAM Unsorted "
        "--outSAMattributes NH HI AS nM MD XS --outFilterMultimapNmax 10 "
        "--outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 "
        "--alignIntronMax 20000 --alignMatesGapMax 20000 > {log} 2>&1"


rule filter_remapped_reads:
    input:
        script=WASP / "filter_remapped_reads.py",
        original=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.to.remap.bam",
        remapped=OUT / "ase/wasp/{sample}/remap/{sample}.Aligned.out.bam"
    output:
        bam=OUT / "ase/wasp/{sample}/{sample}.remap.keep.bam"
    log:
        OUT / "logs/ase/wasp/{sample}.filter_remapped.log"
    shell:
        "python {input.script} {input.original} {input.remapped} {output.bam} > {log} 2>&1"


rule merge_wasp_reads:
    input:
        untouched=OUT / "ase/wasp/{sample}/intersect/{sample}.q30.keep.bam",
        remap_keep=OUT / "ase/wasp/{sample}/{sample}.remap.keep.bam"
    output:
        bam=OUT / "ase/wasp/{sample}/{sample}.wasp.sorted.bam"
    threads: 4
    params:
        tmp=lambda w: str(OUT / "ase/wasp" / w.sample / f"{w.sample}.wasp.unsorted.bam")
    shell:
        "samtools merge -@ {threads} -f {params.tmp} {input.untouched} {input.remap_keep} && "
        "samtools sort -@ {threads} -o {output.bam} {params.tmp} && "
        "samtools quickcheck -v {output.bam}"


rule unbiased_duplicate_filter:
    input:
        script=WASP / "rmdup_pe.py",
        bam=OUT / "ase/wasp/{sample}/{sample}.wasp.sorted.bam"
    output:
        bam=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.bam"
    log:
        OUT / "logs/ase/wasp/{sample}.rmdup.log"
    shell:
        "python {input.script} {input.bam} {output.bam} > {log} 2>&1"


rule sort_index_wasp_bam:
    input:
        bam=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.bam"
    output:
        bam=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam",
        bai=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam.bai"
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input.bam} && "
        "samtools quickcheck -v {output.bam} && samtools index -@ {threads} {output.bam}"


rule f1_heterozygous_vcf:
    """Represent every parental-informative SNP as heterozygous in the F1 sample.

    ASEReadCounter requires a genotype column whose sample name matches the BAM
    read group.  The parental VCF has no F1 columns, so retain one placeholder
    column, set every informative locus to 0/1, and rename it to the F1 sample.
    Maternal/paternal orientation remains available in MAT_ALLELE/PAT_ALLELE.
    """
    input:
        vcf=lambda w: str(OUT / "variants/informative" / f"{cross_for(w.sample)}.informative_snps.vcf.gz"),
        tbi=lambda w: str(OUT / "variants/informative" / f"{cross_for(w.sample)}.informative_snps.vcf.gz.tbi")
    output:
        vcf=OUT / "ase/f1_vcf/{sample}.informative_het.vcf.gz",
        tbi=OUT / "ase/f1_vcf/{sample}.informative_het.vcf.gz.tbi"
    params:
        tmp=lambda w: str(OUT / "ase/f1_vcf" / f"{w.sample}.setgt.vcf.gz"),
        rename=lambda w: str(OUT / "ase/f1_vcf" / f"{w.sample}.samples.txt")
    shell:
        "mkdir -p $(dirname {output.vcf}) && "
        "bcftools view -s C1 {input.vcf} -Ou | "
        "bcftools +setGT -Oz -o {params.tmp} -- -t a -n c:'0/1' && "
        "printf 'C1\\t{wildcards.sample}\\n' > {params.rename} && "
        "bcftools reheader -s {params.rename} -o {output.vcf} {params.tmp} && "
        "bcftools index -f -t {output.vcf} && "
        "rm -f {params.tmp} {params.rename}"


rule ase_read_counter:
    input:
        bam=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam",
        bai=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam.bai",
        vcf=OUT / "ase/f1_vcf/{sample}.informative_het.vcf.gz",
        tbi=OUT / "ase/f1_vcf/{sample}.informative_het.vcf.gz.tbi",
        ref=REF,
        fai=str(REF) + ".fai",
        dictionary=DICT
    output:
        table=OUT / "ase/site_counts/{sample}.ase.table"
    params:
        mq=config["ase_filters"]["mapping_quality"],
        bq=config["ase_filters"]["base_quality"]
    log:
        OUT / "logs/ase/counts/{sample}.log"
    shell:
        "mkdir -p $(dirname {output.table}) $(dirname {log}) && "
        "gatk ASEReadCounter -R {input.ref} -I {input.bam} -V {input.vcf} "
        "-O {output.table} --read-filter MappingQualityReadFilter "
        "--minimum-mapping-quality {params.mq} "
        "--min-base-quality {params.bq} --min-depth-of-non-filtered-base 1 "
        "--count-overlap-reads-handling COUNT_FRAGMENTS_REQUIRE_SAME_BASE "
        "--output-format TABLE > {log} 2>&1"


rule gene_fragment_allele_counter:
    input:
        script=ROOT / "pipeline/scripts/count_gene_alleles.py",
        bam=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam",
        bai=OUT / "ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam.bai",
        vcf=lambda w: str(OUT / "variants/informative" / f"{cross_for(w.sample)}.informative_snps.vcf.gz"),
        tbi=lambda w: str(OUT / "variants/informative" / f"{cross_for(w.sample)}.informative_snps.vcf.gz.tbi"),
        gtf=OUT / "reference/ChineseLong_gene_v3.gtf"
    output:
        genes=OUT / "ase/gene_counts/{sample}.gene_allele_counts.tsv.gz",
        sites=OUT / "ase/gene_counts/{sample}.site_fragment_counts.tsv.gz"
    params:
        mq=config["ase_filters"]["mapping_quality"],
        bq=config["ase_filters"]["base_quality"],
        unambiguous=("--exclude-multigene-snps" if config["ase_filters"].get("exclude_multigene_snps", False) else "")
    log:
        OUT / "logs/ase/counts/{sample}.gene_fragments.log"
    shell:
        "python {input.script} --bam {input.bam} --vcf {input.vcf} --gtf {input.gtf} "
        "--output {output.genes} --site-output {output.sites} "
        "--min-mapping-quality {params.mq} --min-base-quality {params.bq} {params.unambiguous} > {log} 2>&1"


rule beta_binomial_ase:
    input:
        script=ROOT / "pipeline/scripts/fit_ase_beta_binomial.py",
        counts=expand(str(OUT / "ase/gene_counts/{sample}.gene_allele_counts.tsv.gz"), sample=F1),
        manifest=config["manifest"]
    output:
        genes=OUT / "ase/statistics/ase_gene_results.tsv",
        differential=OUT / "ase/statistics/differential_ase_MD_vs_MS.tsv",
        summary=OUT / "ase/statistics/ase_stage_summary.tsv"
    params:
        counts_dir=OUT / "ase/gene_counts",
        outdir=OUT / "ase/statistics",
        group_depth=config["ase_filters"]["gene_group_min_depth"],
        replicate_depth=config["ase_filters"]["gene_replicate_min_depth"],
        min_replicates=config["ase_filters"]["gene_min_covered_replicates"],
        fdr=config["ase_filters"]["fdr"],
        ratio=config["ase_filters"]["min_allelic_ratio"],
        diff=config["ase_filters"]["differential_min_log2_ratio"]
    log:
        OUT / "logs/ase/beta_binomial.log"
    shell:
        "python {input.script} --manifest {input.manifest} --counts-dir {params.counts_dir} "
        "--output-dir {params.outdir} --group-min-depth {params.group_depth} "
        "--replicate-min-depth {params.replicate_depth} "
        "--min-covered-replicates {params.min_replicates} --fdr {params.fdr} "
        "--min-allelic-ratio {params.ratio} --diff-min-log2-ratio {params.diff} "
        "> {log} 2>&1"


rule summarize_ase_dynamics:
    input:
        script=ROOT / "pipeline/scripts/summarize_ase_dynamics.py",
        ase=rules.beta_binomial_ase.output.genes
    output:
        patterns=OUT / "ase/statistics/ase_stage_patterns.tsv",
        intersections=OUT / "ase/statistics/ase_stage_intersection_summary.tsv",
        direction=OUT / "ase/statistics/ase_direction_tests.tsv"
    log:
        OUT / "logs/ase/ase_dynamics.log"
    shell:
        "python {input.script} --ase {input.ase} --patterns {output.patterns} "
        "--intersections {output.intersections} --direction-tests {output.direction} "
        "> {log} 2>&1"
