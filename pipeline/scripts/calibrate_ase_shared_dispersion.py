#!/usr/bin/env python3
"""Cross-fitted shared-dispersion sensitivity analysis with parametric simulation.

This is a companion to, not a silent replacement for, the legacy gene-wise MLE.
Each gene's dispersion is estimated from OTHER genes in the same genotype/stage.
Beta-binomial null simulations calibrate a quasi-score statistic conditional on
that dispersion. Differential ASE uses a fitted common null mean (parametric
bootstrap, not an exact test). Null simulation diagnostics and 0.5x/2x dispersion
sensitivity are mandatory outputs. No SNPs or read pairs are biological replicates.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, chi2
from statsmodels.stats.multitest import multipletests

STAGES = ["seedling_leaf", "fruiting_leaf", "fruit_3d", "fruit_6d", "fruit_9d"]
SCALES = [("half", 0.5), ("main", 1.0), ("double", 2.0)]


def stable_seed(text, seed=20261013):
    return int.from_bytes(hashlib.sha256(f"{seed}:{text}".encode()).digest()[:8], "little")


def bh(values):
    x = np.asarray(values, float)
    out = np.full(len(x), np.nan)
    ok = np.isfinite(x)
    if ok.any():
        out[ok] = multipletests(x[ok], method="fdr_bh")[1]
    return out


def dispersion_moments(y, n):
    """Unbiased moment contributions for unequal binomial depths, k >= 2.

    For x_i=y_i/n_i, T=xbar*(1-xbar)+s^2/k estimates mu*(1-mu).
    E[s^2 - mean(1/n)*T] = rho*(1-mean(1/n))*mu*(1-mu).
    Pool numerators/denominators over genes before clipping the ratio.
    """
    y, n = np.asarray(y), np.asarray(n)
    keep = n >= 5
    x, depths = y[keep] / n[keep], n[keep]
    if len(x) < 2:
        return np.nan, np.nan
    s2 = np.var(x, ddof=1)
    t = x.mean() * (1 - x.mean()) + s2 / len(x)
    a = np.mean(1 / depths)
    return s2 - a * t, (1 - a) * t


def statistic(y, n, rho, split=None):
    """Vectorized quasi-score; y can have shape (draws, libraries)."""
    y, n, rho = np.asarray(y), np.asarray(n), np.asarray(rho)
    weights = np.where(n > 0, n / (1 + np.maximum(n - 1, 0) * rho), 0)
    x = np.divide(y, n, out=np.zeros_like(y, dtype=float), where=n > 0)
    if split is None:
        mean = np.sum(weights * x, axis=-1) / weights.sum()
        return (mean - 0.5) ** 2 * weights.sum() / 0.25
    w0, w1 = weights[:split], weights[split:]
    m0 = np.sum(w0 * x[..., :split], axis=-1) / w0.sum()
    m1 = np.sum(w1 * x[..., split:], axis=-1) / w1.sum()
    pooled = np.sum(weights * x, axis=-1) / weights.sum()
    variance = np.maximum(pooled * (1 - pooled), 1e-12) * (1 / w0.sum() + 1 / w1.sum())
    return (m0 - m1) ** 2 / variance


def weighted_effect(y, n, rho):
    weights = n / (1 + (n - 1) * rho)
    mu = float(np.sum(weights * y / n) / weights.sum())
    bounded = np.clip(mu, 1e-6, 1 - 1e-6)
    ratio = float(np.log2(bounded / (1 - bounded)))
    se = float(1 / np.sqrt(weights.sum() * mu * (1 - mu)) / np.log(2)) if 0 < mu < 1 else np.nan
    return mu, ratio, se


def simulate_pvalue(y, n, rho, split, draws, seed):
    """Plus-one Monte Carlo p; conservative futility stopping returns p=1.

    Stopping never reports a small partially simulated p-value. Retained events
    use the predeclared full number of simulations. Simulation uncertainty is
    recorded; p-values are never zero.
    """
    rng = np.random.default_rng(seed)
    y, n, rho = np.asarray(y), np.asarray(n, int), np.asarray(rho)
    observed = float(statistic(y, n, rho, split))
    w = n / (1 + (n - 1) * rho)
    null_mu = 0.5 if split is None else float(np.clip(np.sum(w * y / n) / w.sum(), 1e-6, 1 - 1e-6))
    concentration = 1 / rho - 1
    exceed = simulated = 0
    while simulated < draws:
        size = min(2047 if simulated == 0 else 8192, draws - simulated)
        probs = rng.beta(null_mu * concentration, (1 - null_mu) * concentration, size=(size, len(n)))
        sampled = rng.binomial(n, probs)
        exceed += int(np.count_nonzero(statistic(sampled, n, rho, split) >= observed - 1e-12))
        simulated += size
        if simulated < draws and simulated >= 2047 and exceed / simulated >= 0.2:
            return 1.0, 1.0, 1.0, simulated, exceed, True
    p = (exceed + 1) / (draws + 1)
    lo = 0.0 if exceed == 0 else float(beta.ppf(0.025, exceed, draws - exceed + 1))
    hi = 1.0 if exceed == draws else float(beta.ppf(0.975, exceed + 1, draws - exceed))
    return p, lo, max(p, hi), simulated, exceed, False


def load_events(manifest_path, counts_dir):
    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest = manifest[manifest.genotype.isin(["MD", "MS"])].copy()
    if manifest.sample_id.duplicated().any() or not manifest.groupby(["genotype", "stage"]).size().eq(3).all():
        raise ValueError("Expected unique IDs and three declared libraries per genotype/stage")
    frames = []
    for row in manifest.itertuples(index=False):
        f = pd.read_csv(Path(counts_dir) / f"{row.sample_id}.gene_allele_counts.tsv.gz", sep="\t")
        if f.gene_id.duplicated().any():
            raise ValueError(f"Duplicate gene IDs in {row.sample_id}")
        fields = ["maternal_fragments", "paternal_fragments", "ambiguous_fragments", "n_informative_snps"]
        a = f[fields].to_numpy(float)
        if not np.isfinite(a).all() or (a < 0).any() or not (a == np.floor(a)).all():
            raise ValueError(f"Invalid counts in {row.sample_id}")
        f["sample_id"], f["genotype"], f["stage"] = row.sample_id, row.genotype, row.stage
        frames.append(f)
    counts = pd.concat(frames, ignore_index=True)
    events = []
    for (genotype, stage, gene), f in counts.groupby(["genotype", "stage", "gene_id"], sort=True):
        ids = manifest.loc[(manifest.genotype == genotype) & (manifest.stage == stage), "sample_id"]
        name = str(f.gene_name.iloc[0])
        f = f.set_index("sample_id").reindex(ids).fillna(0)
        y = f.maternal_fragments.to_numpy(int)
        n = y + f.paternal_fragments.to_numpy(int)
        covered = n >= 5
        tested = bool(covered.sum() >= 2 and n[covered].sum() >= 30)
        events.append(dict(genotype=genotype, stage=stage, gene_id=gene, gene_name=name,
                           y=y.tolist(), n=n.tolist(), testable=tested,
                           covered_replicates=int(covered.sum()), total_depth=int(n.sum()),
                           maternal_total=int(y.sum()), paternal_total=int(n.sum()-y.sum()),
                           ambiguous_total=int(f.ambiguous_fragments.sum()),
                           informative_snps_max=int(f.n_informative_snps.max()),
                           replicate_depths=",".join(map(str, n))))
    return events


def crossfit_dispersion(events):
    estimates, records = {}, []
    for genotype in ("MD", "MS"):
        for stage in STAGES:
            group = [e for e in events if e["genotype"] == genotype and e["stage"] == stage and e["testable"]]
            moments = [(stable_seed(e["gene_id"]) % 5, *dispersion_moments(e["y"], e["n"])) for e in group]
            for fold in range(5):
                train = [(a, b) for f, a, b in moments if f != fold and np.isfinite(a+b) and b > 0]
                if len(train) < 20:
                    raise ValueError(f"Too few training genes for shared dispersion: {genotype}/{stage}/{fold}")
                raw = sum(a for a, b in train) / sum(b for a, b in train)
                rho = float(np.clip(raw, 1e-4, 0.95))
                estimates[(genotype, stage, fold)] = rho
                records.append(dict(genotype=genotype, stage=stage, held_out_fold=fold,
                                    n_training_genes=len(train), rho_unclipped=raw, rho=rho,
                                    rho_clipped=bool(raw != rho)))
    return estimates, pd.DataFrame(records)


def analyze_task(task):
    key, row, y, n, rho, split, draws, seed = task
    y, n, rho = np.asarray(y), np.asarray(n, int), np.asarray(rho, float)
    for label, scale in SCALES:
        scaled = np.clip(rho * scale, 1e-4, 0.95)
        if split is None:
            effect = weighted_effect(y, n, scaled)[1]
        else:
            effect = (weighted_effect(y[split:], n[split:], scaled[split:])[1]
                      - weighted_effect(y[:split], n[:split], scaled[:split])[1])
        p, lo, hi, b, k, stopped = simulate_pvalue(y, n, scaled, split, draws, stable_seed(key + label, seed))
        row.update({f"p_{label}": p, f"p_mc_lower_{label}": lo, f"p_mc_upper_{label}": hi,
                    f"draws_{label}": b, f"exceedances_{label}": k, f"conservative_stop_{label}": stopped,
                    f"asymptotic_p_{label}": float(chi2.sf(statistic(y, n, scaled, split), 1)),
                    f"effect_{label}": effect})
    if split is None:
        mu, effect, se = weighted_effect(y, n, rho)
        signs = np.sign(2*y-n)
        row.update(maternal_fraction=mu, log2_maternal_paternal_ratio=effect, log2_ratio_se=se,
                   rho=float(rho[0]), beta_binomial_concentration=float(1/rho[0]-1),
                   replicate_direction_consistent=bool(len(set(signs[signs != 0])) == 1))
    else:
        mu0, r0, se0 = weighted_effect(y[:split], n[:split], rho[:split])
        mu1, r1, se1 = weighted_effect(y[split:], n[split:], rho[split:])
        row.update(MD_maternal_fraction=mu0, MS_maternal_fraction=mu1, MD_log2_ratio=r0, MS_log2_ratio=r1,
                   delta_log2_ratio_MS_minus_MD=r1-r0, delta_maternal_fraction=mu1-mu0,
                   delta_log2_ratio_se=float(np.sqrt(se0**2+se1**2)),
                   MD_rho=float(rho[0]), MS_rho=float(rho[-1]))
    row.update(pvalue=row["p_main"], model_converged=True, calibration_scope="conditional_shared_dispersion")
    return key, row


def apply_fdr(frame, differential=False):
    family = ["stage"] if differential else ["genotype", "stage"]
    effect = "delta_log2_ratio_MS_minus_MD" if differential else "log2_maternal_paternal_ratio"
    for scale, _ in SCALES:
        for col, source in [(f"q_{scale}", f"p_{scale}"), (f"q_mc_upper_{scale}", f"p_mc_upper_{scale}")]:
            frame[col] = np.nan
            for _, idx in frame.groupby(family).groups.items():
                frame.loc[idx, col] = bh(frame.loc[idx, source])
        effect_col = f"effect_{scale}" if f"effect_{scale}" in frame else effect
        frame[f"called_{scale}"] = (frame[f"q_{scale}"] < 0.05) & (frame[effect_col].abs() >= 1)
    frame["padj"] = frame.q_main
    frame["global_padj"] = bh(frame.p_main)
    frame["differential_ase" if differential else "ase_significant"] = frame.called_main
    frame["stable_across_dispersion_scales"] = frame[[f"called_{s}" for s, _ in SCALES]].all(axis=1)
    frame["mc_interval_supports_call"] = (frame.q_mc_upper_main < 0.05) & (frame[effect].abs() >= 1)
    if not differential:
        frame["ase_direction"] = "balanced"
        frame.loc[frame.ase_significant & (frame[effect] > 0), "ase_direction"] = "maternal_C_bias"
        for genotype, parent in [("MD", "D"), ("MS", "S")]:
            frame.loc[frame.ase_significant & (frame.genotype == genotype) & (frame[effect] < 0), "ase_direction"] = f"paternal_{parent}_bias"
    return frame


def null_diagnostics(tasks, nchecks, seed):
    """Assess conditional type-I error at the observed library-depth designs.

    This checks the implemented simulation model, not adequacy of that model for
    real biological replicates. For differential ASE also probes null means
    0.05/0.25/0.5/0.75/0.95 instead of evaluating only balanced alleles.
    """
    families = {}
    for task in tasks:
        row = task[1]
        families.setdefault((row.get("genotype", "MD_vs_MS"), row["stage"]), []).append(task)
    result = []
    for family, items in families.items():
        rng = np.random.default_rng(stable_seed(str(family)+"null", seed))
        ps = []
        for j in range(nchecks):
            task = items[int(rng.integers(len(items)))]
            _, _, _, n, rho, split, _, _ = task
            n, rho = np.asarray(n, int), np.asarray(rho)
            mu = 0.5 if split is None else [0.05, 0.25, 0.5, 0.75, 0.95][j % 5]
            c = 1 / rho - 1
            y = rng.binomial(n, rng.beta(mu*c, (1-mu)*c))
            p = simulate_pvalue(y, n, rho, split, 1999, stable_seed(f"{family}/{j}/diagnostic", seed))[0]
            ps.append(p)
        k = int(np.count_nonzero(np.array(ps) < 0.05))
        result.append(dict(genotype=family[0], stage=family[1], simulations=nchecks,
                           rejection_count_005=k, rejection_fraction_005=k/nchecks,
                           rejection_ci_low=0 if k == 0 else beta.ppf(0.025, k, nchecks-k+1),
                           rejection_ci_high=1 if k == nchecks else beta.ppf(0.975, k+1, nchecks-k),
                           diagnostic_flag="review" if k/nchecks > 0.08 else "within_tolerance",
                           scope="conditional_model_check_not_real_data_validation"))
    return pd.DataFrame(result)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--counts-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--same-sites", action="store_true")
    p.add_argument("--workers", type=int, default=12)
    p.add_argument("--draws", type=int, default=199999)
    p.add_argument("--null-checks", type=int, default=500)
    p.add_argument("--seed", type=int, default=20261013)
    args = p.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(exist_ok=True)
    events = load_events(args.manifest, args.counts_dir)
    estimates, dispersions = crossfit_dispersion(events)
    dispersions.to_csv(out/"shared_dispersion.tsv", sep="\t", index=False)
    task_list, indexed, untested = [], {}, []
    for e in events:
        row = {k: v for k, v in e.items() if k not in ("y", "n")}
        if not e["testable"]:
            untested.append(row)
            continue
        y, n = np.asarray(e["y"]), np.asarray(e["n"])
        mask = n >= 5
        y, n = y[mask], n[mask]
        fold = stable_seed(e["gene_id"]) % 5
        rho = np.repeat(estimates[(e["genotype"], e["stage"], fold)], len(n))
        key = "/".join(str(e[k]) for k in ("genotype", "stage", "gene_id"))
        task_list.append((key, row, y, n, rho, None, args.draws, args.seed))
        indexed[(e["genotype"], e["stage"], e["gene_id"])] = (row, y, n, rho)
    if args.same_sites:
        for (genotype, stage, gene), (md, y0, n0, rho0) in indexed.items():
            if genotype != "MD" or ("MS", stage, gene) not in indexed:
                continue
            ms, y1, n1, rho1 = indexed[("MS", stage, gene)]
            row = dict(stage=stage, gene_id=gene, gene_name=md["gene_name"],
                       MD_total_depth=int(n0.sum()), MS_total_depth=int(n1.sum()),
                       MD_covered_replicates=len(n0), MS_covered_replicates=len(n1))
            task_list.append((f"diff/{stage}/{gene}", row, np.r_[y0,y1], np.r_[n0,n1], np.r_[rho0,rho1], len(n0), args.draws, args.seed))
    fingerprint = hashlib.sha256(Path(__file__).read_bytes())
    fingerprint.update(json.dumps(dict(draws=args.draws, seed=args.seed, same_sites=args.same_sites), sort_keys=True).encode())
    for e in events:
        fingerprint.update(json.dumps(e, sort_keys=True).encode())
    stamp = fingerprint.hexdigest()
    stamp_file = checkpoints/"input.sha256"
    if stamp_file.exists() and stamp_file.read_text().strip() != stamp:
        raise ValueError("Checkpoint inputs or code changed; use a new output directory")
    stamp_file.write_text(stamp+"\n")
    results, todo = [], []
    for task in task_list:
        file = checkpoints / f"{hashlib.sha256(task[0].encode()).hexdigest()}.json"
        if file.exists():
            results.append((task[0], json.loads(file.read_text())))
        else:
            todo.append(task)
    print(f"events={len(task_list)} cached={len(results)} queued={len(todo)}", flush=True)
    manager = ProcessPoolExecutor(max_workers=args.workers) if args.workers > 1 else nullcontext(None)
    with manager as pool:
        iterator = pool.map(analyze_task, todo, chunksize=4) if pool else map(analyze_task, todo)
        for i, (key, row) in enumerate(iterator, 1):
            file = checkpoints/f"{hashlib.sha256(key.encode()).hexdigest()}.json"
            tmp = file.with_suffix(".tmp")
            tmp.write_text(json.dumps(row, ensure_ascii=False)+"\n")
            tmp.replace(file)
            results.append((key, row))
            if i % 100 == 0 or i == len(todo):
                print(f"calibrated={len(results)}/{len(task_list)}", flush=True)
    ase = pd.DataFrame([row for key, row in results if not key.startswith("diff/")]+untested)
    ase = apply_fdr(ase).sort_values(["genotype", "stage", "gene_id"])
    ase.to_csv(out/"ase_gene_results.tsv", sep="\t", index=False)
    ase.groupby(["genotype", "stage"]).agg(testable=("testable", "sum"), called=("ase_significant", "sum"),
            scale_stable=("stable_across_dispersion_scales", "sum")).to_csv(out/"ase_stage_summary.tsv", sep="\t")
    if args.same_sites:
        diff = pd.DataFrame([row for key, row in results if key.startswith("diff/")])
        if diff.empty:
            diff = pd.DataFrame(columns=["stage", "gene_id", "p_main", "delta_log2_ratio_MS_minus_MD"]+
                                 [f"{prefix}_{scale}" for prefix in ("p", "p_mc_upper") for scale, _ in SCALES])
        diff = apply_fdr(diff, True).sort_values(["stage", "gene_id"])
        diff.to_csv(out/"differential_ase_MD_vs_MS.tsv", sep="\t", index=False)
    diagnostics = null_diagnostics(task_list, args.null_checks, args.seed)
    diagnostics.to_csv(out/"null_calibration_diagnostics.tsv", sep="\t", index=False)
    metadata = vars(args) | dict(input_code_sha256=stamp, folds=5, rho_limits=[1e-4,0.95],
        dispersion_scales=dict(SCALES), minimum_retained_replicate_depth=5, minimum_retained_group_depth=30,
        minimum_covered_replicates=2, finite_simulation_pvalue_floor=1/(args.draws+1),
        diagnostic_review_required=bool((diagnostics.diagnostic_flag == "review").any()),
        interpretation="Shared-dispersion sensitivity, conditional on estimated rho; not experimental validation. "
                       "Differential test uses a fitted null mean. Saturated fractions have no delta-method SE. "
                       "BH families are genotype x stage for ASE and stage for same-SNP differential ASE. "
                       "Futility-stopped tests retain p=1 in full BH families. Fold fitting avoids self-training.")
    (out/"CALIBRATION.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False)+"\n")
    print(json.dumps(dict(completed=True, events=len(task_list), diagnostics_review=metadata["diagnostic_review_required"])), flush=True)


if __name__ == "__main__":
    main()
