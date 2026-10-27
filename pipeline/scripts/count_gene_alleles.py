#!/usr/bin/env python3
"""Count maternal/paternal RNA fragments per gene without double-counting SNPs.

Each read pair (query name) contributes at most one observation to a gene even
when it overlaps several informative SNPs. Fragments with discordant allele
observations within a gene are marked ambiguous and excluded.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path

import pysam
from intervaltree import IntervalTree


ATTRIBUTE = re.compile(r'([A-Za-z0-9_.-]+)\s+"([^"]+)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--site-output", required=True)
    parser.add_argument("--min-mapping-quality", type=int, default=30)
    parser.add_argument("--min-base-quality", type=int, default=20)
    parser.add_argument("--exclude-multigene-snps", action="store_true",
                        help="Exclude SNPs overlapping exons of more than one gene.")
    return parser.parse_args()


def open_text(path: str, mode: str):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def exon_trees(gtf: str):
    trees = defaultdict(IntervalTree)
    gene_names = {}
    with open_text(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(ATTRIBUTE.findall(fields[8]))
            gene = attrs.get("gene_id") or attrs.get("gene")
            if not gene:
                continue
            start0, end = int(fields[3]) - 1, int(fields[4])
            trees[fields[0]].addi(start0, end, gene)
            gene_names[gene] = attrs.get("gene_name", gene)
    return trees, gene_names


def usable_alignment(read, min_mq: int) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_qcfail
        or read.is_duplicate
        or not read.is_proper_pair
        or read.mapping_quality < min_mq
    )


def main() -> None:
    args = parse_args()
    trees, gene_names = exon_trees(args.gtf)
    bam = pysam.AlignmentFile(args.bam, "rb")
    vcf = pysam.VariantFile(args.vcf)
    output = Path(args.output)
    site_output = Path(args.site_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    site_output.parent.mkdir(parents=True, exist_ok=True)

    with open_text(str(output), "wt") as gene_handle, open_text(str(site_output), "wt") as site_handle:
        gene_writer = csv.writer(gene_handle, delimiter="\t", lineterminator="\n")
        site_writer = csv.writer(site_handle, delimiter="\t", lineterminator="\n")
        gene_writer.writerow([
            "gene_id", "gene_name", "maternal_fragments", "paternal_fragments",
            "ambiguous_fragments", "total_informative_fragments", "n_informative_snps"
        ])
        site_writer.writerow([
            "contig", "position", "ref", "alt", "gene_ids", "maternal_fragments",
            "paternal_fragments", "other_fragments", "total_pass_fragments"
        ])

        for contig in bam.references:
            if contig not in trees or contig not in vcf.header.contigs:
                continue
            gene_calls = defaultdict(dict)
            gene_snps = Counter()
            for record in vcf.fetch(contig):
                if len(record.ref) != 1 or len(record.alts or ()) != 1 or len(record.alts[0]) != 1:
                    continue
                genes = {interval.data for interval in trees[contig].at(record.pos - 1)}
                if not genes:
                    continue
                if args.exclude_multigene_snps and len(genes) != 1:
                    continue
                mat_index = record.info.get("MAT_ALLELE")
                pat_index = record.info.get("PAT_ALLELE")
                if mat_index not in (0, 1) or pat_index not in (0, 1):
                    continue
                alleles = (record.ref.upper(), record.alts[0].upper())
                maternal_base, paternal_base = alleles[mat_index], alleles[pat_index]
                site_calls = {}
                for column in bam.pileup(
                    contig, record.pos - 1, record.pos, truncate=True,
                    stepper="all", min_base_quality=0, max_depth=1_000_000
                ):
                    if column.reference_pos != record.pos - 1:
                        continue
                    for pileup_read in column.pileups:
                        read = pileup_read.alignment
                        qpos = pileup_read.query_position
                        if pileup_read.is_del or pileup_read.is_refskip or qpos is None:
                            continue
                        if not usable_alignment(read, args.min_mapping_quality):
                            continue
                        qualities = read.query_qualities
                        if qualities is None or qualities[qpos] < args.min_base_quality:
                            continue
                        base = read.query_sequence[qpos].upper()
                        call = 0 if base == maternal_base else 1 if base == paternal_base else -1
                        previous = site_calls.get(read.query_name)
                        if previous is None:
                            site_calls[read.query_name] = call
                        elif previous != call:
                            site_calls[read.query_name] = -1

                site_counts = Counter(site_calls.values())
                site_writer.writerow([
                    contig, record.pos, record.ref, record.alts[0], ",".join(sorted(genes)),
                    site_counts[0], site_counts[1], site_counts[-1], len(site_calls)
                ])
                for gene in genes:
                    gene_snps[gene] += 1
                    calls = gene_calls[gene]
                    for query_name, call in site_calls.items():
                        previous = calls.get(query_name)
                        if previous is None:
                            calls[query_name] = call
                        elif previous != call:
                            calls[query_name] = -1

            for gene in sorted(gene_calls):
                counts = Counter(gene_calls[gene].values())
                gene_writer.writerow([
                    gene, gene_names.get(gene, gene), counts[0], counts[1], counts[-1],
                    counts[0] + counts[1], gene_snps[gene]
                ])

    bam.close()
    vcf.close()
    print(output)


if __name__ == "__main__":
    main()
