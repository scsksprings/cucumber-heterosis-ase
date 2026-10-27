from pathlib import Path
import pandas as pd

configfile: "pipeline/config.yaml"

ROOT = Path(config["project_root"])
OUT = Path(config["output_root"])
SOURCE_REF = Path(config["reference_fasta"])
MANIFEST = pd.read_csv(config["manifest"], sep="\t").set_index("sample_id")
PARENT_SAMPLES = MANIFEST[MANIFEST["genotype"].isin(["C", "D", "S"])].index.tolist()
REF = OUT / "reference/ChineseLong_genome_v3.fa"
DICT = OUT / "reference/ChineseLong_genome_v3.dict"


rule all:
    input:
        OUT / "variants/filtered/parents.snp.pass.vcf.gz.tbi",
        OUT / "variants/informative/CD.informative_snps.vcf.gz.tbi",
        OUT / "variants/informative/CS.informative_snps.vcf.gz.tbi",
        OUT / "variants/informative/informative_snp_summary.tsv"


rule local_reference:
    input:
        SOURCE_REF
    output:
        REF
    shell:
        "cp -p {input} {output}"


rule reference_fai:
    input:
        REF
    output:
        str(REF) + ".fai"
    shell:
        "samtools faidx {input}"


rule reference_dict:
    input:
        REF
    output:
        DICT
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"


rule mark_parent_duplicates:
    """Mark optical/PCR duplicate candidates before RNA-seq variant calling."""
    input:
        bam=OUT / "alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai=OUT / "alignment/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        bam=OUT / "variants/dedup/{sample}.dedup.bam",
        bai=OUT / "variants/dedup/{sample}.dedup.bam.bai",
        metrics=OUT / "variants/dedup/{sample}.markduplicates.metrics.txt"
    threads: 2
    resources:
        mem_mb=8000
    log:
        OUT / "logs/variants/markduplicates/{sample}.log"
    shell:
        "mkdir -p $(dirname {output.bam}) $(dirname {log}) && "
        "gatk --java-options '-Xmx7g' MarkDuplicates "
        "-I {input.bam} -O {output.bam} -M {output.metrics} "
        "--CREATE_INDEX false --VALIDATION_STRINGENCY SILENT "
        "> {log} 2>&1 && "
        "samtools index -@ {threads} -o {output.bai} {output.bam} && "
        "samtools quickcheck -v {output.bam}"


rule split_n_cigar:
    input:
        bam=rules.mark_parent_duplicates.output.bam,
        bai=rules.mark_parent_duplicates.output.bai,
        ref=REF,
        fai=str(REF) + ".fai",
        dictionary=DICT
    output:
        bam=OUT / "variants/split/{sample}.split.bam",
        bai=OUT / "variants/split/{sample}.split.bai"
    threads: 2
    resources:
        mem_mb=6000
    log:
        OUT / "logs/variants/split/{sample}.log"
    shell:
        "mkdir -p $(dirname {output.bam}) $(dirname {log}) && "
        "gatk --java-options '-Xmx5g' SplitNCigarReads -R {input.ref} -I {input.bam} "
        "-O {output.bam} > {log} 2>&1 && "
        "samtools index -@ {threads} -o {output.bai} {output.bam} && "
        "samtools quickcheck -v {output.bam}"


rule haplotype_caller:
    input:
        bam=OUT / "variants/split/{sample}.split.bam",
        bai=OUT / "variants/split/{sample}.split.bai",
        ref=REF,
        fai=str(REF) + ".fai",
        dictionary=DICT
    output:
        vcf=OUT / "variants/gvcf/{sample}.g.vcf.gz",
        tbi=OUT / "variants/gvcf/{sample}.g.vcf.gz.tbi"
    threads: 2
    resources:
        mem_mb=6000
    log:
        OUT / "logs/variants/haplotypecaller/{sample}.log"
    shell:
        "mkdir -p $(dirname {output.vcf}) $(dirname {log}) && "
        "gatk --java-options '-Xmx5g' HaplotypeCaller -R {input.ref} -I {input.bam} "
        "-O {output.vcf} -ERC GVCF --dont-use-soft-clipped-bases true "
        "--native-pair-hmm-threads {threads} --standard-min-confidence-threshold-for-calling 20 > {log} 2>&1"


rule combine_parent_gvcfs:
    input:
        vcfs=expand(str(OUT / "variants/gvcf/{sample}.g.vcf.gz"), sample=PARENT_SAMPLES),
        tbis=expand(str(OUT / "variants/gvcf/{sample}.g.vcf.gz.tbi"), sample=PARENT_SAMPLES),
        ref=REF,
        fai=str(REF) + ".fai",
        dictionary=DICT
    output:
        vcf=OUT / "variants/joint/parents.combined.g.vcf.gz",
        tbi=OUT / "variants/joint/parents.combined.g.vcf.gz.tbi"
    resources:
        mem_mb=24000
    params:
        variants=lambda w, input: " ".join(f"-V {v}" for v in input.vcfs)
    log:
        OUT / "logs/variants/combine_gvcfs.log"
    shell:
        "mkdir -p $(dirname {output.vcf}) $(dirname {log}) && "
        "gatk --java-options '-Xmx22g' CombineGVCFs -R {input.ref} {params.variants} "
        "-O {output.vcf} > {log} 2>&1"


rule genotype_parent_gvcfs:
    input:
        vcf=OUT / "variants/joint/parents.combined.g.vcf.gz",
        tbi=OUT / "variants/joint/parents.combined.g.vcf.gz.tbi",
        ref=REF,
        fai=str(REF) + ".fai",
        dictionary=DICT
    output:
        vcf=OUT / "variants/joint/parents.raw.vcf.gz",
        tbi=OUT / "variants/joint/parents.raw.vcf.gz.tbi"
    resources:
        mem_mb=16000
    log:
        OUT / "logs/variants/genotype_gvcfs.log"
    shell:
        "gatk --java-options '-Xmx14g' GenotypeGVCFs -R {input.ref} -V {input.vcf} "
        "-O {output.vcf} > {log} 2>&1"


rule select_snps:
    input:
        vcf=OUT / "variants/joint/parents.raw.vcf.gz",
        tbi=OUT / "variants/joint/parents.raw.vcf.gz.tbi",
        ref=REF,
        fai=str(REF) + ".fai",
        dictionary=DICT
    output:
        vcf=OUT / "variants/filtered/parents.snp.raw.vcf.gz",
        tbi=OUT / "variants/filtered/parents.snp.raw.vcf.gz.tbi"
    log:
        OUT / "logs/variants/select_snps.log"
    shell:
        "mkdir -p $(dirname {output.vcf}) && "
        "gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include SNP "
        "--restrict-alleles-to BIALLELIC -O {output.vcf} > {log} 2>&1"


rule hard_filter_snps:
    input:
        vcf=OUT / "variants/filtered/parents.snp.raw.vcf.gz",
        tbi=OUT / "variants/filtered/parents.snp.raw.vcf.gz.tbi",
        ref=REF,
        fai=str(REF) + ".fai",
        dictionary=DICT
    output:
        vcf=OUT / "variants/filtered/parents.snp.filtered.vcf.gz",
        tbi=OUT / "variants/filtered/parents.snp.filtered.vcf.gz.tbi"
    log:
        OUT / "logs/variants/filter_snps.log"
    shell:
        "gatk VariantFiltration -R {input.ref} -V {input.vcf} -O {output.vcf} "
        "--filter-name QD2 --filter-expression 'QD < 2.0' "
        "--filter-name QUAL30 --filter-expression 'QUAL < 30.0' "
        "--filter-name FS60 --filter-expression 'FS > 60.0' "
        "--filter-name SOR3 --filter-expression 'SOR > 3.0' "
        "--filter-name MQ40 --filter-expression 'MQ < 40.0' "
        "--filter-name MQRankSum12.5 --filter-expression 'MQRankSum < -12.5' "
        "--filter-name ReadPosRankSum8 --filter-expression 'ReadPosRankSum < -8.0' "
        "> {log} 2>&1"


rule retain_pass_snps:
    input:
        vcf=OUT / "variants/filtered/parents.snp.filtered.vcf.gz",
        tbi=OUT / "variants/filtered/parents.snp.filtered.vcf.gz.tbi"
    output:
        vcf=OUT / "variants/filtered/parents.snp.pass.vcf.gz",
        tbi=OUT / "variants/filtered/parents.snp.pass.vcf.gz.tbi"
    shell:
        "bcftools view -f PASS,. -Oz -o {output.vcf} {input.vcf} && "
        "tabix -f -p vcf {output.vcf}"


rule informative_parent_snps:
    input:
        script=ROOT / "pipeline/scripts/select_informative_parent_snps.py",
        vcf=OUT / "variants/filtered/parents.snp.pass.vcf.gz",
        tbi=OUT / "variants/filtered/parents.snp.pass.vcf.gz.tbi",
        manifest=config["manifest"]
    output:
        cd_vcf=OUT / "variants/informative/CD.informative_snps.vcf.gz",
        cd_tbi=OUT / "variants/informative/CD.informative_snps.vcf.gz.tbi",
        cs_vcf=OUT / "variants/informative/CS.informative_snps.vcf.gz",
        cs_tbi=OUT / "variants/informative/CS.informative_snps.vcf.gz.tbi",
        summary=OUT / "variants/informative/informative_snp_summary.tsv"
    params:
        min_dp=config["ase_filters"]["parent_min_depth"],
        min_gq=config["ase_filters"]["genotype_quality"],
        min_af=config["ase_filters"]["parent_min_allele_fraction"],
        min_support=config["ase_filters"]["parent_min_supporting_samples"]
    log:
        OUT / "logs/variants/informative_snps.log"
    shell:
        "mkdir -p $(dirname {output.cd_vcf}) && "
        "python {input.script} --vcf {input.vcf} --manifest {input.manifest} "
        "--output-dir $(dirname {output.cd_vcf}) --min-dp {params.min_dp} "
        "--min-gq {params.min_gq} --min-allele-fraction {params.min_af} "
        "--min-supporting-samples {params.min_support} > {log} 2>&1"
