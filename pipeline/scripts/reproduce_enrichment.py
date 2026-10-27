#!/usr/bin/env python3
"""Recompute all 20 GO/KEGG families from supplied gene sets and term mapping."""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def main():
    root=Path(__file__).resolve().parents[2]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=root/'results/S3')
    p.add_argument('--output-dir',type=Path,default=root/'work/enrichment')
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    if a.output_dir.resolve()==a.results.resolve():p.error('Use a separate output directory')
    read=lambda n:pd.read_csv(a.results/n,sep='\t')
    mapping=read('functional_term_to_gene.tsv.gz')
    terms={key:set(frame.gene_id) for key,frame in mapping.groupby(['database','term_id'])}
    labels=mapping.drop_duplicates(['database','term_id']).set_index(['database','term_id']).term_name
    annotated={db:set(frame.gene_id) for db,frame in mapping.groupby('database')}
    sets=read('enrichment_gene_sets_and_backgrounds.tsv.gz')
    collections=[]
    for (name,direction),frame in sets.groupby(['set','direction']):
        collections.append((name,direction,set(frame.gene_id),set(frame.loc[frame.selected,'gene_id'])))
    overlap=read('S3_leaf_overlap_all_tested.tsv.gz')
    for direction,sign in [('positive',1),('negative',-1)]:
        selected=overlap[overlap.overlap_category.eq('both_same_direction') &
            (np.sign(overlap.log2FoldChange_midparent_MD)==sign)]
        collections.append(('shared_'+direction,direction,set(overlap.gene_id),set(selected.gene_id)))
    rows=[]
    for name,direction,universe,selection in collections:
        for db in ['GO','KEGG']:
            bg=universe & annotated[db]; selected=selection & bg
            for (database,term),genes in sorted(terms.items()):
                if database!=db:continue
                members=genes & bg
                if len(members)<5:continue
                hit=members & selected
                M,K,n,k=len(bg),len(members),len(selected),len(hit)
                rows.append(dict(set=name,direction=direction,database=db,term_id=term,
                    term_name=labels.loc[(db,term)],background_n=M,selected_n=n,
                    term_background_n=K,overlap_n=k,fold_enrichment=(k/n)/(K/M) if n else 0.,
                    pvalue=float(hypergeom.sf(k-1,M,K,n)) if n else 1.,gene_ids=';'.join(sorted(hit))))
    new=pd.DataFrame(rows)
    for _,idx in new.groupby(['set','direction','database']).groups.items():
        new.loc[idx,'q_BH']=multipletests(new.loc[idx,'pvalue'],method='fdr_bh')[1]
    new['q_global_sensitivity']=multipletests(new.pvalue,method='fdr_bh')[1]
    base=read('functional_enrichment_all_tests.tsv.gz');shared=read('S3_shared_direction_enrichment.tsv')
    shared['direction']=shared['set'].str.replace('shared_','',regex=False)
    old=pd.concat([base,shared],ignore_index=True)
    keys=['set','direction','database','term_id']
    old=old.set_index(keys).sort_index();new=new.set_index(keys).sort_index()
    assert old.index.equals(new.index)
    columns=['background_n','selected_n','term_background_n','overlap_n','fold_enrichment','pvalue','q_BH','q_global_sensitivity']
    np.testing.assert_allclose(old[columns],new[columns],rtol=1e-9,atol=1e-13)
    assert old.gene_ids.fillna('').equals(new.gene_ids.fillna(''))
    new.reset_index().to_csv(a.output_dir/'all_20_families.tsv.gz',sep='\t',index=False)
    report={'status':'passed','families':len(new.reset_index().groupby(keys[:3])),
            'tests':len(new),'local_significant':int((new.q_BH<.05).sum()),
            'verified':'All term memberships, hypergeometric P, family BH and global BH match the manuscript snapshot.'}
    (a.output_dir/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
