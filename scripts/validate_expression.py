#!/usr/bin/env python3
"""Compare recomputed DESeq2 effects, interactions and PCA with manuscript data."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]
read=lambda f:pd.read_csv(R/f,sep='\t')
new=read('work/expression/inheritance/gene_results.tsv')
old=read('results/S3/S3_expression_inheritance_all.tsv.gz')
keys=['gene_id','genotype','stage'];a=old.set_index(keys).sort_index();b=new.set_index(keys).sort_index()
assert a.index.equals(b.index)
cols=['log2FoldChange_midparent','padj_midparent','log2FoldChange_hybrid_C','padj_hybrid_C','log2FoldChange_hybrid_paternal','padj_hybrid_paternal']
np.testing.assert_allclose(a[cols],b[cols],rtol=1e-6,atol=1e-8,equal_nan=True)
assert a.midparent_significant.equals(b.midparent_significant)
new=read('work/expression/fruit_stage_interactions.tsv');old=read('results/S4/S4_fruit_stage_interactions_all.tsv.gz')
keys=['gene_id','hybrid','comparison'];a=old.set_index(keys).sort_index();b=new.set_index(keys).sort_index()
assert a.index.equals(b.index)
np.testing.assert_allclose(a[['log2FoldChange','padj']],b[['log2FoldChange','padj']],rtol=1e-6,atol=1e-8,equal_nan=True)
for name in ['all','fruit']:
    old=read(f'results/S1/PCA_{name}_scores.tsv').set_index('sample_id')
    new=read(f'work/expression/pca/PCA_{name}_scores.tsv').set_index('sample_id').loc[old.index]
    # SVD sign is arbitrary; each entire axis may flip between BLAS implementations.
    for col in ['PC1','PC2']:
        sign=1 if np.dot(old[col],new[col])>=0 else -1
        np.testing.assert_allclose(old[col],sign*new[col],rtol=1e-5,atol=1e-6)
    x=read(f'results/S1/PCA_{name}_top500_genes.tsv')
    y=read(f'work/expression/pca/PCA_{name}_top500_genes.tsv')
    assert set(x.gene_id)==set(y.gene_id)
contrasts=read('work/expression/deseq2/deseq2_all_contrasts.tsv').set_index(['gene_id','stage','comparison'])
for name in ['parental_effect','log_midparent']:
    old=read(f'results/S8/qPCR_RNAseq_{name}_comparison.tsv').set_index(['gene_id','stage','comparison'])
    new=contrasts.reindex(old.index)
    np.testing.assert_allclose(old[['log2FoldChange','padj']],new[['log2FoldChange','padj']],rtol=1e-6,atol=1e-8,equal_nan=True)
report={'status':'passed','expression_gene_context_records':len(read('work/expression/inheritance/gene_results.tsv')),
        'fruit_interaction_records':len(read('work/expression/fruit_stage_interactions.tsv')),
        'PCA_subsets':['all','fruit'],'qPCR_RNAseq_effects_checked':810,
        'numerical_tolerance':'DESeq2 rtol=1e-6, atol=1e-8; PCA allows whole-axis sign flip'}
(R/'work/expression/validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
