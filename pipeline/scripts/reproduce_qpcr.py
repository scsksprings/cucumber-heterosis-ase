#!/usr/bin/env python3
"""Recompute the 27-gene manuscript qPCR statistics and stratified permutation.

Statistical functions extracted from plot_qpcr_publication.py; panel-wide BH
families now refer to the supplied 27-gene manuscript set (135 panels).
"""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
from string import ascii_lowercase
import argparse, json
import numpy as np
import pandas as pd
from scipy.stats import studentized_range, t, spearmanr, ttest_ind, rankdata
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.oneway import anova_oneway
GENOTYPES = ['C','D','MD','S','MS']
CORRELATION_PERMUTATIONS=99999
CORRELATION_SEED=20261019
OUT=None

def write_table(frame,name):
    frame.to_csv(OUT/name,sep='\t',index=False)

def letters_from_pairs(groups, p_matrix, expression_order):
    """Maximal cliques cover every nonsignificant pair without merging significant pairs."""
    adjacency = p_matrix >= .05
    candidates = []
    for size in range(1, len(groups) + 1):
        for vertices in combinations(range(len(groups)), size):
            if all(adjacency[i, j] for i, j in combinations(vertices, 2)):
                candidates.append(frozenset(vertices))
    cliques = [c for c in candidates if not any(c < bigger for bigger in candidates)]
    ranks = {i: rank for rank, i in enumerate(np.argsort(-np.asarray(expression_order)))}
    cliques.sort(key=lambda c: (min(ranks[i] for i in c), tuple(sorted(ranks[i] for i in c))))
    labels = {g: '' for g in groups}
    for k, clique in enumerate(cliques):
        for i in clique:
            labels[groups[i]] += ascii_lowercase[k]
    assert all(labels.values())
    for i, j in combinations(range(len(groups)), 2):
        shared = bool(set(labels[groups[i]]) & set(labels[groups[j]]))
        assert shared == bool(adjacency[i, j]), (groups[i], groups[j], p_matrix[i, j], labels)
    return labels

def statistics(wells):
    summaries, pairs, omnibus, midparent = [], [], [], []
    panel_groups = list(wells.groupby(['gene_id', 'stage'], sort=True))
    for index, ((gene, stage), frame) in enumerate(panel_groups):
        arrays = [frame[frame.genotype == g].delta_ct.to_numpy() for g in GENOTYPES]
        means = np.array([a.mean() for a in arrays])
        variances = np.array([a.var(ddof=1) for a in arrays])
        sizes = np.array([len(a) for a in arrays])
        assert (variances > 0).all(), (gene, stage, variances)
        test = anova_oneway(arrays, use_var='unequal', welch_correction=True)
        omnibus.append({'gene_id': gene, 'stage': stage, 'F_Welch': float(test.statistic),
            'df_numerator': float(test.df[0]), 'df_denominator': float(test.df[1]), 'p_Welch': float(test.pvalue)})
        p_matrix = np.ones((5, 5))
        for i, j in combinations(range(5), 2):
            vi, vj = variances[i] / sizes[i], variances[j] / sizes[j]
            se = np.sqrt((vi + vj) / 2)
            df = (vi + vj)**2 / (vi**2 / (sizes[i] - 1) + vj**2 / (sizes[j] - 1))
            q_stat = abs(means[i] - means[j]) / se
            p_computed = float(studentized_range.sf(q_stat, len(GENOTYPES), df))
            p = max(p_computed, np.finfo(float).eps)
            p_matrix[i, j] = p_matrix[j, i] = p
            welch = ttest_ind(arrays[i], arrays[j], equal_var=False)
            np.testing.assert_allclose(abs(welch.statistic) * np.sqrt(2), q_stat, rtol=1e-10)
            pairs.append({'gene_id': gene, 'stage': stage, 'group1': GENOTYPES[i], 'group2': GENOTYPES[j],
                'n1': int(sizes[i]), 'n2': int(sizes[j]), 'log2FC_group1_over_group2': float(means[j] - means[i]),
                'q_studentized_range': q_stat, 'df_Games_Howell': df,
                'p_Games_Howell_within_panel': p, 'p_at_numerical_floor': p_computed < np.finfo(float).eps,
                'p_Welch_unadjusted': float(welch.pvalue), 'different_at_005': p < .05})
        labels = letters_from_pairs(GENOTYPES, p_matrix, -means)
        for i, g in enumerate(GENOTYPES):
            sub = frame[frame.genotype == g]
            summaries.append({'gene_id': gene, 'stage': stage, 'genotype': g, 'n_biological': len(sub),
                'delta_ct_mean': means[i], 'delta_ct_sd': np.sqrt(variances[i]),
                'relative_expression_mean': sub.relative_expression.mean(),
                'relative_expression_sd': sub.relative_expression.std(ddof=1),
                'letters_Games_Howell_005': labels[g], 'assay_number': int(sub.assay_number.iloc[0]),
                'identity_status': sub.identity_status.iloc[0]})
        # Match the paper's formal log-parental-mean contrast; do not label this
        # as a test against an arithmetic midparent expression value.
        for hybrid, paternal in [('MD', 'D'), ('MS', 'S')]:
            indices = [GENOTYPES.index(x) for x in [hybrid, 'C', paternal]]
            weights = np.array([-1., .5, .5])  # expression is -delta Ct
            estimate = float(np.dot(weights, means[indices]))
            components = weights**2 * variances[indices] / sizes[indices]
            variance = components.sum()
            df = variance**2 / np.sum(components**2 / (sizes[indices] - 1))
            p = float(2 * t.sf(abs(estimate) / np.sqrt(variance), df))
            group_means = frame.groupby('genotype').relative_expression.mean()
            arithmetic_effect = float(np.log2(group_means[hybrid] / ((group_means['C'] + group_means[paternal]) / 2)))
            midparent.append({'gene_id': gene, 'stage': stage, 'hybrid': hybrid, 'paternal': paternal,
                'D_log_parental_mean': estimate, 'standard_error': float(np.sqrt(variance)), 'df': float(df),
                'p_Welch_linear_contrast': p, 'D_arithmetic_MPV_effect_only': arithmetic_effect,
                'identity_status': frame.identity_status.iloc[0]})
        if (index + 1) % 25 == 0:
            print(f'Computed Games–Howell comparisons for {index + 1}/{len(panel_groups)} panels', flush=True)
    summary, pairwise = pd.DataFrame(summaries), pd.DataFrame(pairs)
    omni, mpv = pd.DataFrame(omnibus), pd.DataFrame(midparent)
    omni['q_BH_across_135_panels'] = multipletests(omni.p_Welch, method='fdr_bh')[1]
    mpv['q_BH_across_270_log_midparent_contrasts'] = multipletests(mpv.p_Welch_linear_contrast, method='fdr_bh')[1]
    # This second, conservative summary is available for analyses across the
    # whole assay panel; plot letters use the within-panel Games–Howell family.
    pairwise['q_BH_across_all_1350_raw_Welch_tests'] = multipletests(pairwise.p_Welch_unadjusted, method='fdr_bh')[1]
    write_table(summary, 'qPCR_group_summary_and_letters.tsv')
    write_table(pairwise, 'qPCR_all_pairwise_Games_Howell.tsv')
    write_table(omni, 'qPCR_Welch_ANOVA.tsv')
    write_table(mpv, 'qPCR_log_midparent_deviation.tsv')
    return summary, pairwise, omni, mpv

def gene_block_permutation(data, field, n_permutations=CORRELATION_PERMUTATIONS,
                           seed=CORRELATION_SEED):
    """Test gene-profile matching, preserving batch, coverage and within-gene structure.

    Under the random-matching null, whole qPCR gene profiles are exchangeable
    within qPCR batch and identical stage/contrast coverage. The pooled Spearman
    statistic need not have a zero-centered permutation null. We therefore use
    twice the smaller inclusive Monte Carlo tail, with a +1 correction.
    This is conditional profile-matching evidence, not a population-level test
    based on 524/262 independent experimental replicates.
    """
    data = data.sort_values(['gene_id', 'stage', 'comparison']).reset_index(drop=True)
    assert not data.duplicated(['gene_id', 'stage', 'comparison']).any()
    x = rankdata(data.log2FoldChange, method='average')
    y = rankdata(data[field], method='average')
    x = x - x.mean()
    y = y - y.mean()
    x = x / np.linalg.norm(x)
    y = y / np.linalg.norm(y)
    observed = float(x @ y)
    np.testing.assert_allclose(observed, spearmanr(data.log2FoldChange, data[field]).statistic, atol=1e-14)
    strata = defaultdict(list)
    for gene, group in data.groupby('gene_id', sort=True):
        assert group.batch.nunique() == 1
        key = (str(group.batch.iloc[0]), tuple(zip(group.stage, group.comparison)))
        strata[key].append((gene, group.index.to_numpy()))
    rng = np.random.default_rng(seed)
    null = np.zeros(n_permutations)
    assignments, examples = [], []
    # Retain five explicit row mappings to verify the optimized statistic
    # independently with scipy.stats.spearmanr and audit preserved boundaries.
    example_indices = np.tile(np.arange(len(data)), (min(5, n_permutations), 1))
    movable_genes = 0
    for stratum_id, ((batch, pattern), blocks) in enumerate(sorted(strata.items()), start=1):
        rows = np.stack([rows for _, rows in blocks])
        n_genes = len(blocks)
        movable_genes += n_genes if n_genes > 1 else 0
        cross_products = x[rows] @ y[rows].T
        for gene, indices in blocks:
            assignments.append({'comparison_type': field, 'stratum_id': stratum_id,
                'gene_id': gene, 'qPCR_batch': batch, 'events_per_gene': len(indices),
                'genes_in_stratum': n_genes, 'permutable': n_genes > 1,
                'stage_comparison_pattern': ';'.join(f'{stage}:{contrast}' for stage, contrast in pattern)})
        for start in range(0, n_permutations, 4096):
            stop = min(start + 4096, n_permutations)
            permutations = rng.permuted(np.tile(np.arange(n_genes), (stop - start, 1)), axis=1)
            null[start:stop] += cross_products[np.arange(n_genes), permutations].sum(axis=1)
            if start == 0:
                for i in range(len(example_indices)):
                    example_indices[i, rows.ravel()] = rows[permutations[i]].ravel()
    assert movable_genes > 1
    for i, row_mapping in enumerate(example_indices):
        np.testing.assert_array_equal(np.sort(row_mapping), np.arange(len(data)))
        for column in ['batch', 'stage', 'comparison']:
            np.testing.assert_array_equal(data[column].to_numpy(), data[column].to_numpy()[row_mapping])
        direct = float(spearmanr(data.log2FoldChange, data[field].to_numpy()[row_mapping]).statistic)
        np.testing.assert_allclose(direct, null[i], atol=1e-14)
        for gene, group in data.groupby('gene_id', sort=True):
            source = data.iloc[row_mapping[group.index]]
            assert source.gene_id.nunique() == 1
            examples.append({'comparison_type': field, 'permutation_id': i + 1,
                'RNAseq_gene': gene, 'qPCR_gene': source.gene_id.iloc[0],
                'rho_direct_scipy': direct, 'rho_rank_matrix': float(null[i])})
    tolerance = 100 * np.finfo(float).eps
    lower = int(np.count_nonzero(null <= observed + tolerance))
    upper = int(np.count_nonzero(null >= observed - tolerance))
    p_value = min(1., 2 * (min(lower, upper) + 1) / (n_permutations + 1))
    report = {'comparison_type': field, 'events': len(data), 'genes': data.gene_id.nunique(),
        'spearman_rho_descriptive': observed, 'p_gene_block_permutation_two_sided': p_value,
        'permutations': n_permutations, 'seed': seed, 'exchangeability_strata': len(strata),
        'permutable_genes': movable_genes, 'fixed_genes': data.gene_id.nunique() - movable_genes,
        'null_lower_tail_count': lower, 'null_upper_tail_count': upper,
        'minimum_reportable_two_sided_p': 2 / (n_permutations + 1),
        'p_at_monte_carlo_floor': min(lower, upper) == 0,
        'permuted_rho_mean': float(null.mean()), 'permuted_rho_sd': float(null.std(ddof=1)),
        'permuted_rho_2_5_percentile': float(np.quantile(null, .025)),
        'permuted_rho_97_5_percentile': float(np.quantile(null, .975)),
        'permutation_scope': 'whole qPCR gene profiles within batch and identical stage/contrast coverage',
        'null_hypothesis': 'random gene-profile matching conditional on batch and observed coverage',
        'assumption': 'gene profiles exchangeable within each stratum under the null'}
    return report, null, assignments, examples


def main():
    global OUT
    root=Path(__file__).resolve().parents[2]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=root/'results/S8')
    p.add_argument('--output-dir',type=Path,default=root/'work/qpcr')
    p.add_argument('--permutations',type=int,default=CORRELATION_PERMUTATIONS)
    a=p.parse_args(); OUT=a.output_dir; OUT.mkdir(parents=True,exist_ok=True)
    if OUT.resolve()==a.results.resolve():p.error('Output must differ from the published results directory')
    wells=pd.read_csv(a.results/'qPCR_biological_replicate_values.tsv',sep='\t')
    assert wells.gene_id.nunique()==27 and len(wells)==2025
    wells['delta_ct']=wells.target_ct-wells.reference_ct
    wells['log2_normalized_expression']=-wells.delta_ct
    calibrator=wells[wells.genotype=='C'].groupby(['gene_id','stage']).delta_ct.mean()
    wells['C_mean_delta_ct']=[calibrator.loc[(g,s)] for g,s in zip(wells.gene_id,wells.stage)]
    wells['delta_delta_ct']=wells.delta_ct-wells.C_mean_delta_ct
    wells['relative_expression']=2.**(-wells.delta_delta_ct)
    write_table(wells,'qPCR_biological_replicate_values.tsv')
    summary,pairwise,omni,mp=statistics(wells)
    # Check against the independent manuscript tables, not against constants alone.
    for name,fresh,keys,cols in [
      ('qPCR_group_summary_and_letters.tsv',summary,['gene_id','stage','genotype'],['delta_ct_mean','delta_ct_sd','relative_expression_mean','relative_expression_sd']),
      ('qPCR_Welch_ANOVA.tsv',omni,['gene_id','stage'],['F_Welch','p_Welch','q_BH_across_135_panels']),
      ('qPCR_all_pairwise_Games_Howell.tsv',pairwise,['gene_id','stage','group1','group2'],['p_Games_Howell_within_panel'])]:
        old=pd.read_csv(a.results/name,sep='\t').set_index(keys).sort_index()
        new=fresh.set_index(keys).sort_index()
        assert old.index.equals(new.index)
        np.testing.assert_allclose(old[cols],new[cols],rtol=1e-8,atol=1e-10)
        if 'letters_Games_Howell_005' in old:
            assert old.letters_Games_Howell_005.equals(new.letters_Games_Howell_005)
    means=wells.groupby(['gene_id','stage','genotype']).log2_normalized_expression.mean()
    reports=[]
    for name,field in [('qPCR_RNAseq_parental_effect_comparison.tsv','qPCR_log2FC'),('qPCR_RNAseq_log_midparent_comparison.tsv','D_log_parental_mean')]:
        d=pd.read_csv(a.results/name,sep='\t')
        values=[]
        for row in d.itertuples():
            hybrid,ref=row.comparison.split('_vs_')
            get=lambda g:means.loc[(row.gene_id,row.stage,g)]
            if ref=='midparent':value=get(hybrid)-(get('C')+get('D' if hybrid=='MD' else 'S'))/2
            else:value=get(hybrid)-get(ref)
            values.append(value)
        np.testing.assert_allclose(d[field],values,rtol=1e-10,atol=1e-12)
        # Keep serialized effects after the independent Ct check: sub-ULP changes
        # from spreadsheet serialization can break exact rank ties.
        report,_,strata,examples=gene_block_permutation(d[d.included_in_concordance],field,a.permutations)
        reports.append(report)
        write_table(pd.DataFrame(strata),field+'_permutation_strata.tsv')
    holm=multipletests([x['p_gene_block_permutation_two_sided'] for x in reports],method='holm')[1]
    for row,q in zip(reports,holm):row['p_Holm_across_two_concordance_tests']=float(q)
    fresh=pd.DataFrame(reports)
    old=pd.read_csv(a.results/'qPCR_RNAseq_concordance_summary.tsv',sep='\t')
    cols=['events','genes','spearman_rho_descriptive']
    if a.permutations==99999:cols+=['p_gene_block_permutation_two_sided','p_Holm_across_two_concordance_tests']
    np.testing.assert_allclose(old[cols],fresh[cols],rtol=1e-10,atol=1e-12)
    write_table(fresh,'qPCR_RNAseq_concordance_summary.tsv')
    print(fresh[['comparison_type','events','genes','spearman_rho_descriptive','p_Holm_across_two_concordance_tests']].to_string(index=False))
    (OUT/'validation.json').write_text(json.dumps({'status':'passed','genes':27,'panels':135,'replicates':2025,'permutations':a.permutations},indent=2)+'\n')

if __name__=='__main__':main()
