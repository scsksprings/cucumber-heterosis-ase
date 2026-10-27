#!/usr/bin/env python3
"""Convert a biallelic VCF to WASP's one-text-file-per-contig SNP format."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pysam


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    handles = {}
    counts = {}
    vcf = pysam.VariantFile(args.vcf)
    try:
        for record in vcf:
            if len(record.ref) != 1 or len(record.alts or ()) != 1 or len(record.alts[0]) != 1:
                continue
            if record.contig not in handles:
                handles[record.contig] = gzip.open(
                    output / f"{record.contig}.snps.txt.gz", "wt"
                )
                counts[record.contig] = 0
            handles[record.contig].write(f"{record.pos} {record.ref} {record.alts[0]}\n")
            counts[record.contig] += 1
    finally:
        vcf.close()
        for handle in handles.values():
            handle.close()

    with (output / "conversion_summary.tsv").open("w") as handle:
        handle.write("contig\tn_snps\n")
        for contig, count in counts.items():
            handle.write(f"{contig}\t{count}\n")
    (output / ".done").write_text("complete\n")
    print(f"Wrote {sum(counts.values())} SNPs on {len(counts)} contigs to {output}")


if __name__ == "__main__":
    main()
