#!/usr/bin/env python3
"""Shared-site preparation/recount extracted from the frozen final ASE run."""
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import argparse, hashlib, json, subprocess, sys
import pandas as pd
import pysam
from count_gene_alleles import exon_trees

def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False)+"\n")
    tmp.replace(path)

def snp_key(rec):
    return rec.contig, rec.pos, rec.ref, tuple(rec.alts or ())

def prepare_shared(root):
    trees, _ = exon_trees(str(root/"reference/ChineseLong_gene_v3.gtf"))
    vcfdir = root/"variants/informative"
    summaries = []
    for cross in ("CD", "CS"):
        c = Counter()
        with pysam.VariantFile(vcfdir/f"{cross}.informative_snps.vcf.gz") as vcf:
            for r in vcf:
                c["total"] += 1
                genes = {x.data for x in trees[r.contig].at(r.pos-1)}
                c["no_exonic_gene" if not genes else "single_gene" if len(genes) == 1 else "multiple_genes"] += 1
        summaries.append({"cross": cross, **c})
    cs = {}
    with pysam.VariantFile(vcfdir/"CS.informative_snps.vcf.gz") as source:
        for r in source:
            key = snp_key(r)
            if key in cs:
                raise ValueError(f"Duplicate SNP in CS: {key}")
            cs[key] = dict(r.info)
    records = []
    for support in (3, 4):
        out = root/f"common_support{support}"
        out.mkdir(exist_ok=True)
        target = out/"shared_informative_snps.vcf.gz"
        counters = Counter()
        with pysam.VariantFile(vcfdir/"CD.informative_snps.vcf.gz") as source:
            header = source.header.copy()
            header.info.add("CS_MAT_N", 1, "Integer", "C supporting libraries in C-S comparison")
            header.info.add("CS_PAT_N", 1, "Integer", "S supporting libraries in C-S comparison")
            temp = out/"shared_informative_snps.tmp.vcf.gz"
            with pysam.VariantFile(temp, "wz", header=header) as writer:
                seen = set()
                for r in source:
                    key = snp_key(r)
                    if key in seen:
                        raise ValueError(f"Duplicate SNP in CD: {key}")
                    seen.add(key)
                    if key not in cs:
                        continue
                    counters["exact_allele_intersection"] += 1
                    other = cs[key]
                    if any(r.info.get(k) != other.get(k) for k in ("MAT_ALLELE", "PAT_ALLELE")):
                        counters["orientation_mismatch"] += 1
                        continue
                    if min(r.info["MAT_N"], r.info["PAT_N"], other["MAT_N"], other["PAT_N"]) < support:
                        counters["insufficient_parent_support"] += 1
                        continue
                    genes = {x.data for x in trees[r.contig].at(r.pos-1)}
                    if len(genes) != 1:
                        counters["no_exonic_gene" if not genes else "multiple_genes"] += 1
                        continue
                    r = r.copy()
                    r.translate(header)
                    r.info["CS_MAT_N"], r.info["CS_PAT_N"] = other["MAT_N"], other["PAT_N"]
                    writer.write(r)
                    counters["retained_same_allele_single_gene_snps"] += 1
            temp.replace(target)
        pysam.tabix_index(str(target), preset="vcf", force=True)
        records.append({"minimum_parent_support": support, **counters})
        if counters["retained_same_allele_single_gene_snps"] == 0:
            raise ValueError(f"No shared SNPs pass support {support}")
    pd.DataFrame(summaries).fillna(0).to_csv(root/"variants/gene_assignment_audit.tsv", sep="\t", index=False)
    pd.DataFrame(records).fillna(0).to_csv(root/"variants/shared_snp_audit.tsv", sep="\t", index=False)
    print(pd.DataFrame(records).to_string(index=False), flush=True)

def recount_shared(root, workers):
    manifest = pd.read_csv(root/"inputs/sample_manifest.tsv", sep="\t")
    jobs = [(s, support) for s in manifest.loc[manifest.genotype.isin(["MD", "MS"]), "sample_id"] for support in (3,4)]
    script = Path(__file__).with_name("count_gene_alleles.py")

    def run(job):
        sample, support = job
        out = root/f"common_support{support}/counts"
        out.mkdir(parents=True, exist_ok=True)
        genes = out/f"{sample}.gene_allele_counts.tsv.gz"
        sites = out/f"{sample}.site_fragment_counts.tsv.gz"
        marker = out/f"{sample}.complete.json"
        bam = root/f"ase/wasp/{sample}/{sample}.wasp.rmdup.sorted.bam"
        vcf = root/f"common_support{support}/shared_informative_snps.vcf.gz"
        h = hashlib.sha256(script.read_bytes()+vcf.read_bytes())
        h.update(f"{bam}:{bam.stat().st_size}:{bam.stat().st_mtime_ns}".encode())
        digest = h.hexdigest()
        if marker.exists() and genes.exists() and sites.exists() and json.loads(marker.read_text()).get("digest") == digest:
            return f"cached {sample} support={support}"
        log = root/f"logs/recount_{sample}_support{support}.log"
        with log.open("a") as handle:
            subprocess.run([sys.executable, str(script), "--bam", str(bam), "--vcf", str(vcf),
                "--gtf", str(root/"reference/ChineseLong_gene_v3.gtf"), "--output", str(genes),
                "--site-output", str(sites), "--min-mapping-quality", "30", "--min-base-quality", "20",
                "--exclude-multigene-snps"], stdout=handle, stderr=subprocess.STDOUT, check=True)
        g, s = pd.read_csv(genes, sep="\t"), pd.read_csv(sites, sep="\t")
        if g.gene_id.duplicated().any() or s.gene_ids.astype(str).str.contains(",").any():
            raise ValueError(f"Ambiguous recount: {sample}")
        write_json(marker, dict(digest=digest, gene_rows=len(g), site_rows=len(s)))
        return f"done {sample} support={support}"
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for message in pool.map(run, jobs):
            print(message, flush=True)

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=["prepare","recount"])
    p.add_argument("--run-root",required=True,type=Path)
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    (a.run_root/"logs").mkdir(parents=True,exist_ok=True)
    if a.action=="prepare": prepare_shared(a.run_root)
    else: recount_shared(a.run_root,a.workers)
