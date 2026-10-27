#!/usr/bin/env python3
"""Select fixed, opposite-homozygous parent SNPs for C-D and C-S ASE."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd
import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-dp", type=int, default=10)
    parser.add_argument("--min-gq", type=float, default=30)
    parser.add_argument("--min-allele-fraction", type=float, default=0.95)
    parser.add_argument("--min-supporting-samples", type=int, default=2)
    return parser.parse_args()


def high_quality_homozygote(call, min_dp: int, min_gq: float, min_af: float):
    gt = call.get("GT")
    dp = call.get("DP")
    gq = call.get("GQ")
    ad = call.get("AD")
    if gt is None or dp is None or gq is None or ad is None:
        return None
    if dp < min_dp or gq < min_gq or len(gt) != 2 or gt[0] is None or gt[1] is None:
        return None
    if gt[0] != gt[1] or gt[0] not in (0, 1):
        return "conflict"
    allele = gt[0]
    total = sum(x for x in ad if x is not None)
    if total < min_dp or allele >= len(ad) or ad[allele] / total < min_af:
        return None
    return allele


def consensus(record, samples, args):
    calls = [
        high_quality_homozygote(
            record.samples[s], args.min_dp, args.min_gq, args.min_allele_fraction
        )
        for s in samples
    ]
    if "conflict" in calls:
        return None, 0, "heterozygous_high_quality_call"
    alleles = [x for x in calls if x in (0, 1)]
    if len(alleles) < args.min_supporting_samples:
        return None, len(alleles), "insufficient_parent_support"
    if len(set(alleles)) != 1:
        return None, len(alleles), "discordant_parent_calls"
    return alleles[0], len(alleles), "pass"


def main() -> None:
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    by_genotype = {
        genotype: manifest.loc[manifest.genotype == genotype, "sample_id"].tolist()
        for genotype in ("C", "D", "S")
    }

    source = pysam.VariantFile(args.vcf)
    missing = sorted(set(sum(by_genotype.values(), [])) - set(source.header.samples))
    if missing:
        raise SystemExit(f"Samples absent from VCF: {','.join(missing)}")

    source.header.info.add("MAT_ALLELE", 1, "Integer", "Maternal C allele index (0=REF,1=ALT)")
    source.header.info.add("PAT_ALLELE", 1, "Integer", "Paternal allele index (0=REF,1=ALT)")
    source.header.info.add("MAT_N", 1, "Integer", "Number of high-quality C samples supporting consensus")
    source.header.info.add("PAT_N", 1, "Integer", "Number of high-quality paternal samples supporting consensus")
    source.header.info.add("ASE_CROSS", 1, "String", "Cross for which SNP is informative")

    paths = {
        "CD": outdir / "CD.informative_snps.vcf.gz",
        "CS": outdir / "CS.informative_snps.vcf.gz",
    }
    writers = {cross: pysam.VariantFile(path, "wz", header=source.header) for cross, path in paths.items()}
    counters = {cross: Counter() for cross in paths}
    total = 0
    for record in source:
        total += 1
        c_allele, c_n, c_reason = consensus(record, by_genotype["C"], args)
        for cross, paternal in (("CD", "D"), ("CS", "S")):
            counters[cross]["pass_vcf_snps_examined"] += 1
            if c_reason != "pass":
                counters[cross][f"reject_C_{c_reason}"] += 1
                continue
            p_allele, p_n, p_reason = consensus(record, by_genotype[paternal], args)
            if p_reason != "pass":
                counters[cross][f"reject_{paternal}_{p_reason}"] += 1
                continue
            if c_allele == p_allele:
                counters[cross]["reject_same_parent_allele"] += 1
                continue
            copied = record.copy()
            copied.info["MAT_ALLELE"] = c_allele
            copied.info["PAT_ALLELE"] = p_allele
            copied.info["MAT_N"] = c_n
            copied.info["PAT_N"] = p_n
            copied.info["ASE_CROSS"] = cross
            writers[cross].write(copied)
            counters[cross]["informative_snps"] += 1

    source.close()
    for writer in writers.values():
        writer.close()
    for path in paths.values():
        pysam.tabix_index(str(path), preset="vcf", force=True)

    categories = sorted(set().union(*(counter.keys() for counter in counters.values())))
    rows = []
    for cross in ("CD", "CS"):
        row = {"cross": cross, **{category: counters[cross][category] for category in categories}}
        row.update(
            min_dp=args.min_dp,
            min_gq=args.min_gq,
            min_allele_fraction=args.min_allele_fraction,
            min_supporting_samples=args.min_supporting_samples,
        )
        rows.append(row)
    pd.DataFrame(rows).fillna(0).to_csv(
        outdir / "informative_snp_summary.tsv", sep="\t", index=False
    )
    print(f"Examined {total} PASS biallelic SNPs")
    for cross in ("CD", "CS"):
        print(f"{cross}: {counters[cross]['informative_snps']} informative SNPs")


if __name__ == "__main__":
    main()
