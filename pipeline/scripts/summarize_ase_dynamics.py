#!/usr/bin/env python3
"""Summarize stable/stage-specific ASE and test allelic-bias proportions."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact
from statsmodels.stats.multitest import multipletests


STAGES = ["seedling_leaf", "fruiting_leaf", "fruit_3d", "fruit_6d", "fruit_9d"]


def bh(values: pd.Series) -> np.ndarray:
    result = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        result[valid] = multipletests(values[valid], method="fdr_bh")[1]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ase", required=True)
    parser.add_argument("--patterns", required=True)
    parser.add_argument("--intersections", required=True)
    parser.add_argument("--direction-tests", required=True)
    args = parser.parse_args()

    ase = pd.read_csv(args.ase, sep="\t")
    ase["stage"] = pd.Categorical(ase.stage, STAGES, ordered=True)
    rows = []
    for (genotype, gene_id), frame in ase.groupby(["genotype", "gene_id"], observed=True):
        frame = frame.sort_values("stage")
        tested = frame.loc[frame.testable == True]  # noqa: E712
        significant = frame.loc[frame.ase_significant == True]  # noqa: E712
        sig_stages = significant.stage.astype(str).tolist()
        directions = significant.ase_direction.tolist()
        rows.append({
            "genotype": genotype,
            "gene_id": gene_id,
            "gene_name": frame.gene_name.dropna().iloc[0] if frame.gene_name.notna().any() else gene_id,
            "testable_stages": ",".join(tested.stage.astype(str)),
            "n_testable_stages": len(tested),
            "ase_stages": ",".join(sig_stages),
            "n_ase_stages": len(significant),
            "ase_stage_pattern": "|".join(sig_stages) if sig_stages else "none",
            "stable_ase_all_five_stages": len(significant) == 5,
            "single_stage_ase": len(significant) == 1,
            "ase_directions": ",".join(directions),
            "direction_consistent_across_ase_stages": len(set(directions)) <= 1,
            "max_abs_log2_ratio": (
                significant.log2_maternal_paternal_ratio.abs().max()
                if len(significant) else np.nan
            ),
        })
    patterns = pd.DataFrame(rows)

    intersections = (
        patterns[patterns.n_ase_stages > 0]
        .groupby(["genotype", "ase_stage_pattern", "n_ase_stages"], observed=True)
        .size().rename("n_genes").reset_index()
        .sort_values(["genotype", "n_genes"], ascending=[True, False])
    )

    direction_rows = []
    for (genotype, stage), frame in ase.groupby(["genotype", "stage"], observed=True):
        significant = frame[frame.ase_significant == True]  # noqa: E712
        maternal = int((significant.ase_direction == "maternal_C_bias").sum())
        paternal = int(significant.ase_direction.str.startswith("paternal_").sum())
        total = maternal + paternal
        pvalue = binomtest(maternal, total, 0.5).pvalue if total else np.nan
        direction_rows.append({
            "test": "within_genotype_maternal_fraction_vs_half",
            "genotype": genotype, "stage": str(stage),
            "maternal_bias": maternal, "paternal_bias": paternal,
            "total_significant_ase": total,
            "maternal_fraction": maternal / total if total else np.nan,
            "odds_ratio": np.nan, "pvalue": pvalue,
        })

    for stage in STAGES:
        stage_rows = [row for row in direction_rows if row["stage"] == stage]
        by_genotype = {row["genotype"]: row for row in stage_rows}
        if set(by_genotype) != {"MD", "MS"}:
            continue
        md, ms = by_genotype["MD"], by_genotype["MS"]
        table = [
            [md["maternal_bias"], md["paternal_bias"]],
            [ms["maternal_bias"], ms["paternal_bias"]],
        ]
        if sum(map(sum, table)):
            odds_ratio, pvalue = fisher_exact(table)
        else:
            odds_ratio, pvalue = np.nan, np.nan
        direction_rows.append({
            "test": "MD_vs_MS_maternal_fraction",
            "genotype": "MD_vs_MS", "stage": stage,
            "maternal_bias": np.nan, "paternal_bias": np.nan,
            "total_significant_ase": np.nan, "maternal_fraction": np.nan,
            "odds_ratio": odds_ratio, "pvalue": pvalue,
        })
    direction = pd.DataFrame(direction_rows)
    direction["padj"] = bh(direction.pvalue)

    for path in (args.patterns, args.intersections, args.direction_tests):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    patterns.to_csv(args.patterns, sep="\t", index=False)
    intersections.to_csv(args.intersections, sep="\t", index=False)
    direction.to_csv(args.direction_tests, sep="\t", index=False)
    print(
        f"ASE dynamics: genes={len(patterns)}; "
        f"genes_with_ASE={int((patterns.n_ase_stages > 0).sum())}; "
        f"stable_all_five={int(patterns.stable_ase_all_five_stages.sum())}"
    )


if __name__ == "__main__":
    main()
