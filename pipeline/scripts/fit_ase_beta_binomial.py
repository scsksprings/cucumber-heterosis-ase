#!/usr/bin/env python3
"""Replicate-aware beta-binomial ASE and differential-ASE analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln, logit
from scipy.stats import binomtest, chi2
from statsmodels.stats.multitest import multipletests


STAGE_ORDER = ["seedling_leaf", "fruiting_leaf", "fruit_3d", "fruit_6d", "fruit_9d"]


def log_bb(y, n, mu, phi):
    a, b = mu * phi, (1 - mu) * phi
    return np.sum(
        gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1)
        + betaln(y + a, n - y + b) - betaln(a, b)
    )


def numerical_hessian(function, point, step=1e-4):
    point = np.asarray(point, dtype=float)
    size = len(point)
    hessian = np.zeros((size, size))
    f0 = function(point)
    for i in range(size):
        ei = np.zeros(size); ei[i] = step
        hessian[i, i] = (function(point + ei) - 2 * f0 + function(point - ei)) / step**2
        for j in range(i):
            ej = np.zeros(size); ej[j] = step
            value = (
                function(point + ei + ej) - function(point + ei - ej)
                - function(point - ei + ej) + function(point - ei - ej)
            ) / (4 * step**2)
            hessian[i, j] = hessian[j, i] = value
    return hessian


def robust_minimize(function, start, bounds):
    result = minimize(function, start, method="L-BFGS-B", bounds=bounds)
    if not result.success or not np.isfinite(result.fun):
        result = minimize(function, start, method="Powell", bounds=bounds,
                          options={"maxiter": 5000, "xtol": 1e-8, "ftol": 1e-10})
    return result


def fit_single_group(y, n):
    pooled = np.clip(y.sum() / n.sum(), 1e-5, 1 - 1e-5)

    def alt_obj(par):
        return -log_bb(y, n, expit(par[0]), np.exp(par[1]))

    def null_obj(par):
        return -log_bb(y, n, 0.5, np.exp(par[0]))

    alt = robust_minimize(alt_obj, [logit(pooled), np.log(100)], [(-10, 10), (-7, 14)])
    null = robust_minimize(null_obj, [np.log(100)], [(-7, 14)])
    if not alt.success or not null.success:
        return pooled, np.nan, np.nan, np.nan, False
    mu, phi = expit(alt.x[0]), np.exp(alt.x[1])
    lrt = max(0.0, 2 * (-alt.fun + null.fun))
    try:
        covariance = np.linalg.inv(numerical_hessian(alt_obj, alt.x))
        eta_se = np.sqrt(covariance[0, 0]) if covariance[0, 0] > 0 else np.nan
    except (np.linalg.LinAlgError, ValueError):
        eta_se = np.nan
    return mu, phi, chi2.sf(lrt, 1), eta_se, True


def fit_two_groups(y, n, group):
    pooled = np.clip(y.sum() / n.sum(), 1e-5, 1 - 1e-5)
    m0 = np.clip(y[group == 0].sum() / n[group == 0].sum(), 1e-5, 1 - 1e-5)
    m1 = np.clip(y[group == 1].sum() / n[group == 1].sum(), 1e-5, 1 - 1e-5)

    def alt_obj(par):
        mus = expit(par[0] + par[1] * group)
        phi = np.exp(par[2])
        a, b = mus * phi, (1 - mus) * phi
        return -np.sum(
            gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1)
            + betaln(y + a, n - y + b) - betaln(a, b)
        )

    def null_obj(par):
        return -log_bb(y, n, expit(par[0]), np.exp(par[1]))

    start = [logit(m0), logit(m1) - logit(m0), np.log(100)]
    alt = robust_minimize(alt_obj, start, [(-10, 10), (-20, 20), (-7, 14)])
    null = robust_minimize(null_obj, [logit(pooled), np.log(100)], [(-10, 10), (-7, 14)])
    if not alt.success or not null.success:
        return m0, m1, np.nan, np.nan, False
    mu0 = expit(alt.x[0])
    mu1 = expit(alt.x[0] + alt.x[1])
    lrt = max(0.0, 2 * (-alt.fun + null.fun))
    return mu0, mu1, np.exp(alt.x[2]), chi2.sf(lrt, 1), True


def bh(values):
    values = np.asarray(values, dtype=float)
    result = np.full(len(values), np.nan)
    ok = np.isfinite(values)
    if ok.any():
        result[ok] = multipletests(values[ok], method="fdr_bh")[1]
    return result


def log2_ratio(mu):
    mu = np.clip(mu, 1e-6, 1 - 1e-6)
    return np.log2(mu / (1 - mu))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--counts-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--group-min-depth", type=int, default=30)
    parser.add_argument("--replicate-min-depth", type=int, default=5)
    parser.add_argument("--min-covered-replicates", type=int, default=2)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-allelic-ratio", type=float, default=2.0)
    parser.add_argument("--diff-min-log2-ratio", type=float, default=1.0)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    manifest = manifest[manifest.genotype.isin(["MD", "MS"])].copy()
    count_frames = []
    for row in manifest.itertuples(index=False):
        path = Path(args.counts_dir) / f"{row.sample_id}.gene_allele_counts.tsv.gz"
        frame = pd.read_csv(path, sep="\t")
        frame["sample_id"] = row.sample_id
        frame["genotype"] = row.genotype
        frame["stage"] = row.stage
        frame["replicate"] = row.replicate
        count_frames.append(frame)
    counts = pd.concat(count_frames, ignore_index=True)
    counts["total"] = counts.maternal_fragments + counts.paternal_fragments

    rows = []
    grouped = counts.groupby(["genotype", "stage", "gene_id", "gene_name"], sort=False)
    for (genotype, stage, gene, name), frame in grouped:
        frame = frame.set_index("sample_id").reindex(
            manifest.loc[(manifest.genotype == genotype) & (manifest.stage == stage), "sample_id"]
        ).fillna(0)
        y = frame.maternal_fragments.to_numpy(dtype=float)
        n = (frame.maternal_fragments + frame.paternal_fragments).to_numpy(dtype=float)
        ambiguous = frame.ambiguous_fragments.to_numpy(dtype=float)
        snp_counts = frame.n_informative_snps.to_numpy(dtype=float)
        covered = int((n >= args.replicate_min_depth).sum())
        testable = n.sum() >= args.group_min_depth and covered >= args.min_covered_replicates
        covered_mask = n >= args.replicate_min_depth
        replicate_signs = np.sign(2 * y[covered_mask] - n[covered_mask])
        non_tied_signs = replicate_signs[replicate_signs != 0]
        direction_consistent = (
            len(non_tied_signs) > 0 and len(np.unique(non_tied_signs)) == 1
        )
        replicate_ratios = [
            log2_ratio(mat / depth) if depth > 0 else np.nan
            for mat, depth in zip(y, n)
        ]
        observed_snp_counts = snp_counts[snp_counts > 0]
        result = {
            "genotype": genotype, "stage": stage, "gene_id": gene, "gene_name": name,
            "maternal_total": int(y.sum()), "paternal_total": int(n.sum() - y.sum()),
            "ambiguous_total": int(ambiguous.sum()),
            "total_depth": int(n.sum()), "covered_replicates": covered,
            "replicate_depths": ",".join(str(int(x)) for x in n), "testable": testable,
            "replicate_log2_ratios": ",".join(
                "NA" if not np.isfinite(x) else f"{x:.6g}" for x in replicate_ratios
            ),
            "replicate_direction_consistent": direction_consistent,
            "informative_snps_min": (
                int(observed_snp_counts.min()) if len(observed_snp_counts) else 0
            ),
            "informative_snps_max": (
                int(observed_snp_counts.max()) if len(observed_snp_counts) else 0
            ),
        }
        if testable:
            mu, phi, pvalue, eta_se, converged = fit_single_group(y, n)
            result.update(
                maternal_fraction=mu, log2_maternal_paternal_ratio=log2_ratio(mu),
                log2_ratio_se=eta_se / np.log(2),
                beta_binomial_concentration=phi, pvalue=pvalue, model_converged=converged,
                pooled_binomial_pvalue=binomtest(int(y.sum()), int(n.sum()), 0.5).pvalue,
            )
        rows.append(result)
    ase = pd.DataFrame(rows)
    ase["padj"] = np.nan
    for _, idx in ase.groupby(["genotype", "stage"]).groups.items():
        ase.loc[idx, "padj"] = bh(ase.loc[idx, "pvalue"])
    cutoff = np.log2(args.min_allelic_ratio)
    ase["ase_significant"] = (
        ase.testable & (ase.padj < args.fdr)
        & (ase.log2_maternal_paternal_ratio.abs() >= cutoff)
    )
    ase["ase_direction"] = "balanced"
    ase.loc[ase.ase_significant & (ase.log2_maternal_paternal_ratio > 0), "ase_direction"] = "maternal_C_bias"
    ase.loc[ase.ase_significant & (ase.genotype == "MD") & (ase.log2_maternal_paternal_ratio < 0), "ase_direction"] = "paternal_D_bias"
    ase.loc[ase.ase_significant & (ase.genotype == "MS") & (ase.log2_maternal_paternal_ratio < 0), "ase_direction"] = "paternal_S_bias"

    diff_rows = []
    testable = ase[ase.testable].copy()
    for (stage, gene), pair in testable.groupby(["stage", "gene_id"]):
        if set(pair.genotype) != {"MD", "MS"}:
            continue
        frames = []
        for genotype in ("MD", "MS"):
            sub = counts[(counts.stage == stage) & (counts.genotype == genotype) & (counts.gene_id == gene)]
            sub = sub.set_index("sample_id").reindex(
                manifest.loc[(manifest.genotype == genotype) & (manifest.stage == stage), "sample_id"]
            ).fillna(0)
            frames.append(sub)
        y = np.concatenate([x.maternal_fragments.to_numpy(dtype=float) for x in frames])
        n = np.concatenate([(x.maternal_fragments + x.paternal_fragments).to_numpy(dtype=float) for x in frames])
        group = np.array([0] * 3 + [1] * 3)
        mu_md, mu_ms, phi, pvalue, converged = fit_two_groups(y, n, group)
        diff_rows.append({
            "stage": stage, "gene_id": gene, "gene_name": pair.gene_name.iloc[0],
            "MD_maternal_fraction": mu_md, "MS_maternal_fraction": mu_ms,
            "MD_log2_ratio": log2_ratio(mu_md), "MS_log2_ratio": log2_ratio(mu_ms),
            "delta_log2_ratio_MS_minus_MD": log2_ratio(mu_ms) - log2_ratio(mu_md),
            "beta_binomial_concentration": phi, "pvalue": pvalue, "model_converged": converged,
        })
    differential = pd.DataFrame(diff_rows)
    if not differential.empty:
        differential["padj"] = np.nan
        for _, idx in differential.groupby("stage").groups.items():
            differential.loc[idx, "padj"] = bh(differential.loc[idx, "pvalue"])
        differential["differential_ase"] = (
            (differential.padj < args.fdr)
            & (differential.delta_log2_ratio_MS_minus_MD.abs() >= args.diff_min_log2_ratio)
        )

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    ase.to_csv(outdir / "ase_gene_results.tsv", sep="\t", index=False)
    differential.to_csv(outdir / "differential_ase_MD_vs_MS.tsv", sep="\t", index=False)
    summary = (
        ase.groupby(["genotype", "stage"], observed=True)
        .agg(
            genes_with_counts=("gene_id", "size"), testable_genes=("testable", "sum"),
            significant_ase_genes=("ase_significant", "sum"),
            maternal_bias=("ase_direction", lambda x: (x == "maternal_C_bias").sum()),
            paternal_bias=("ase_direction", lambda x: x.isin(["paternal_D_bias", "paternal_S_bias"]).sum()),
        ).reset_index()
    )
    summary["ase_fraction_of_testable"] = summary.significant_ase_genes / summary.testable_genes
    summary["stage"] = pd.Categorical(summary.stage, STAGE_ORDER, ordered=True)
    summary.sort_values(["genotype", "stage"]).to_csv(outdir / "ase_stage_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
