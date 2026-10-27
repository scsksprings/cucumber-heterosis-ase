#!/usr/bin/env python3
"""Build a single unstranded gene-count matrix from STAR GeneCounts files."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def build_matrix(manifest: pd.DataFrame, alignment_dir: Path) -> pd.DataFrame:
    """Require complete, matching STAR gene sets before assembling counts."""
    if "sample_id" not in manifest or manifest.empty:
        raise ValueError("Manifest must contain at least one sample_id")
    samples = manifest.sample_id
    if samples.isna().any() or samples.astype(str).str.strip().eq("").any():
        raise ValueError("Manifest contains missing or empty sample_id values")
    if samples.duplicated().any():
        raise ValueError(f"Duplicate sample_id: {samples[samples.duplicated()].tolist()}")

    series = []
    reference_genes = None
    for sample in samples:
        path = alignment_dir / sample / f"{sample}.ReadsPerGene.out.tab"
        frame = pd.read_csv(path, sep="\t", header=None, dtype=str)
        if frame.shape[1] != 4:
            raise ValueError(f"{sample}: STAR GeneCounts must have exactly four columns")
        frame.columns = ["gene_id", "unstranded", "forward", "reverse"]
        expected_summary = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]
        if frame.gene_id.iloc[:4].tolist() != expected_summary:
            raise ValueError(f"{sample}: missing or malformed STAR summary rows")
        frame = frame.iloc[4:].copy()
        if frame.empty or frame.gene_id.isna().any() or frame.gene_id.duplicated().any():
            raise ValueError(f"{sample}: gene IDs must be nonempty and unique")
        counts = frame[["unstranded", "forward", "reverse"]].apply(pd.to_numeric, errors="raise")
        values = counts.to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any() or (values % 1 != 0).any():
            raise ValueError(f"{sample}: counts must be finite nonnegative integers")
        genes = pd.Index(frame.gene_id, name="gene_id")
        if reference_genes is None:
            reference_genes = genes
        else:
            missing = reference_genes.difference(genes)
            extra = genes.difference(reference_genes)
            if len(missing) or len(extra):
                raise ValueError(
                    f"{sample}: gene set differs from {samples.iloc[0]} "
                    f"(missing={len(missing)}, extra={len(extra)}); "
                    "check truncation and reference annotation; missing rows are not zero counts"
                )
        column = pd.Series(counts.unstranded.to_numpy(dtype="int64"), index=genes, name=sample)
        series.append(column.reindex(reference_genes))
    return pd.concat(series, axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--alignment-dir", required=True)
    parser.add_argument("--output-matrix", required=True)
    parser.add_argument("--output-metadata", required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    matrix = build_matrix(manifest, Path(args.alignment_dir))
    if list(matrix.columns) != manifest.sample_id.tolist():
        raise RuntimeError("Count-matrix columns no longer match manifest order")
    Path(args.output_matrix).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_metadata).parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output_matrix, sep="\t")
    manifest.to_csv(args.output_metadata, sep="\t", index=False)
    print(f"{matrix.shape[0]} genes x {matrix.shape[1]} samples")


if __name__ == "__main__":
    main()
