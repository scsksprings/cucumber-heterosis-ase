#!/usr/bin/env python3
"""Classify F1 expression inheritance using uniform DESeq2 contrasts.

The five requested biological classes are assigned conservatively. A gene must
first differ from the mid-parent value at both FDR and effect-size thresholds
to be called nonadditive. Dominance calls additionally require the F1 to differ
from the lower/higher parent in the expected direction while not showing a
threshold-sized significant difference from the matching parent.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


STAGE_ORDER = ["seedling_leaf", "fruiting_leaf", "fruit_3d", "fruit_6d", "fruit_9d"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contrasts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-abs-log2fc", type=float, default=1.0)
    args = parser.parse_args()

    data = pd.read_csv(args.contrasts, sep="\t")
    frames = []
    for genotype, paternal in (("MD", "D"), ("MS", "S")):
        names = {
            "parent": f"{paternal}_vs_C",
            "hybrid_C": f"{genotype}_vs_C",
            "hybrid_paternal": f"{genotype}_vs_{paternal}",
            "midparent": f"{genotype}_vs_midparent",
        }
        subset = data[data.comparison.isin(names.values())].copy()
        values = subset.pivot(
            index=["gene_id", "stage"], columns="comparison",
            values=["log2FoldChange", "padj"],
        )
        values.columns = [
            f"{metric}_{next(key for key, value in names.items() if value == comparison)}"
            for metric, comparison in values.columns
        ]
        values = values.reset_index()
        values["genotype"] = genotype
        values["paternal_parent"] = paternal

        def significant(label: str) -> pd.Series:
            return (
                (values[f"padj_{label}"] < args.fdr)
                & (values[f"log2FoldChange_{label}"].abs() >= args.min_abs_log2fc)
            )

        mp_sig = significant("midparent")
        hc_sig = significant("hybrid_C")
        hp_sig = significant("hybrid_paternal")
        parent_lfc = values["log2FoldChange_parent"]
        paternal_high = parent_lfc > 0

        mp_tested = values["padj_midparent"].notna()
        values["inheritance_class"] = "not_tested"
        values.loc[mp_tested, "inheritance_class"] = "additive"
        values.loc[mp_sig, "inheritance_class"] = "ambiguous_nonadditive"

        over = (
            mp_sig & hc_sig & hp_sig
            & (values.log2FoldChange_hybrid_C > 0)
            & (values.log2FoldChange_hybrid_paternal > 0)
        )
        under = (
            mp_sig & hc_sig & hp_sig
            & (values.log2FoldChange_hybrid_C < 0)
            & (values.log2FoldChange_hybrid_paternal < 0)
        )
        values.loc[over, "inheritance_class"] = "overdominance"
        values.loc[under, "inheritance_class"] = "underdominance"

        high_diff_sig = np.where(paternal_high, hp_sig, hc_sig)
        high_diff_lfc = np.where(
            paternal_high,
            values.log2FoldChange_hybrid_paternal,
            values.log2FoldChange_hybrid_C,
        )
        low_diff_sig = np.where(paternal_high, hc_sig, hp_sig)
        low_diff_lfc = np.where(
            paternal_high,
            values.log2FoldChange_hybrid_C,
            values.log2FoldChange_hybrid_paternal,
        )
        high_dom = (
            mp_sig & ~over & ~under
            & ~high_diff_sig & low_diff_sig & (low_diff_lfc > 0)
        )
        low_dom = (
            mp_sig & ~over & ~under
            & ~low_diff_sig & high_diff_sig & (high_diff_lfc < 0)
        )
        values.loc[high_dom, "inheritance_class"] = "high_parent_dominance"
        values.loc[low_dom, "inheritance_class"] = "low_parent_dominance"
        values["midparent_tested"] = mp_tested
        values["midparent_significant"] = mp_sig
        values["C_parent_is_higher"] = ~paternal_high
        values["classification_fdr"] = args.fdr
        values["classification_min_abs_log2fc"] = args.min_abs_log2fc
        # Preserve legacy machine-readable classes for downstream
        # compatibility, but provide a claim-safe interpretation.  A
        # non-significant difference is not an equivalence test.
        interpretation = {
            "additive": "no_detected_parental_log_mean_departure",
            "high_parent_dominance": "high_parent_like_logmean_departure",
            "low_parent_dominance": "low_parent_like_logmean_departure",
            "overdominance": "above_both_parents_logmean_departure",
            "underdominance": "below_both_parents_logmean_departure",
            "ambiguous_nonadditive": "other_parental_logmean_departure",
            "not_tested": "not_tested",
        }
        values["claim_safe_interpretation"] = values.inheritance_class.map(
            interpretation
        )
        frames.append(values)

    result = pd.concat(frames, ignore_index=True)
    result["stage"] = pd.Categorical(result.stage, STAGE_ORDER, ordered=True)
    result = result.sort_values(["genotype", "stage", "gene_id"])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False)
    summary = (
        result.groupby(
            ["genotype", "stage", "inheritance_class"],
            observed=True,
        ).size().rename("n_genes").reset_index()
    )
    summary.to_csv(args.summary, sep="\t", index=False)
    print(
        f"Classified {len(result)} genotype-stage-gene rows; "
        f"nonadditive={int(result.midparent_significant.sum())}"
    )


if __name__ == "__main__":
    main()
