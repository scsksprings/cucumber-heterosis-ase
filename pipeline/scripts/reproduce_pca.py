#!/usr/bin/env python3
"""PCA from one joint blind VST; select top 500 genes within each plotted subset."""
from pathlib import Path
import argparse,json
import pandas as pd
from sklearn.decomposition import PCA


def main():
    root=Path(__file__).resolve().parents[2]
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vst',type=Path,default=root/'work/expression/deseq2/vst_blind_all_samples.tsv')
    p.add_argument('--metadata',type=Path,default=root/'data/rnaseq/sample_metadata.tsv')
    p.add_argument('--output-dir',type=Path,default=root/'work/expression/pca')
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(a.metadata,sep='\t').set_index('sample_id')
    vst=pd.read_csv(a.vst,sep='\t',index_col=0);variance={}
    for name,ids in [('all',list(meta.index)),('fruit',list(meta.index[meta.stage.str.startswith('fruit_')]))]:
        matrix=vst[ids].T;genes=matrix.var().nlargest(500).index
        model=PCA(n_components=2,svd_solver='full');points=model.fit_transform(matrix[genes])
        scores=pd.DataFrame(points,index=ids,columns=['PC1','PC2']).join(meta[['genotype','stage','replicate']])
        scores.to_csv(a.output_dir/f'PCA_{name}_scores.tsv',sep='\t',index_label='sample_id')
        pd.Series(genes,name='gene_id').to_csv(a.output_dir/f'PCA_{name}_top500_genes.tsv',sep='\t',index=False)
        variance[name]=(100*model.explained_variance_ratio_).tolist()
    (a.output_dir/'variance_percent.json').write_text(json.dumps(variance,indent=2)+'\n')

if __name__=='__main__':main()
