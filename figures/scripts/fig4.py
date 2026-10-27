import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
import numpy as np, pandas as pd
from matplotlib import gridspec
T=str(RES)
p=pd.read_csv(f'{T}/S5/TableS23_parental_SNP_selection.tsv',sep='\t').set_index('cross').informative_snps
s=pd.read_csv(f'{T}/S5/TableS17_shared_SNP_selection.tsv',sep='\t').set_index('minimum_parent_support')
a=pd.read_csv(f'{T}/S6/TableS7_jointly_testable_ASE.tsv',sep='\t'); a=a[a.model=='Shared']

fig=plt.figure(figsize=(W,4.9))
gs=gridspec.GridSpec(2,2,figure=fig,height_ratios=[1,1.05],width_ratios=[0.75,1.25],hspace=0.75,wspace=0.42)
def lab_bars(ax,xs,ys,fmt='{:,}'):
    for x,y in zip(xs,ys): ax.text(x,y*1.01,fmt.format(int(y)).replace(',',' '),ha='center',va='bottom',fontsize=6.3)

ax=fig.add_subplot(gs[0,0]); aA=ax
v=[p['CD'],p['CS']]
ax.bar([0,1],v,width=0.55,color=[GCOL['MD'],GCOL['MS']],edgecolor='none')
lab_bars(ax,[0,1],v)
ax.set_xticks([0,1]); ax.set_xticklabels(['C–D（MD）','C–S（MS）'])
ax.set_ylabel('亲本差异SNP数\nParent-informative SNPs'); ax.set_ylim(0,33000); ax.set_xlim(-0.6,1.6)
ax.set_yticks([0,10000,20000,30000]); ax.set_yticklabels(['0','10 000','20 000','30 000'])
ax.set_title('亲本差异位点\nParent-informative SNPs',loc='left',pad=4)
letter(ax,'A',-0.42,1.04)

ax=fig.add_subplot(gs[0,1]); aB=ax
v=[s.loc[3,'exact_allele_intersection'],s.loc[3,'retained_same_allele_single_gene_snps'],s.loc[4,'retained_same_allele_single_gene_snps']]
ax.bar([0,1,2],v,width=0.55,color=['#9E9E9E','#5E5E5E','#2F2F2F'],edgecolor='none')
lab_bars(ax,[0,1,2],v)
ax.set_xticks([0,1,2]); ax.set_xticklabels(['共同位点\nShared sites','单基因外显子\n亲本支持≥3\nExon, support ≥3','单基因外显子\n亲本支持≥4\nExon, support ≥4'],fontsize=5.8)
ax.set_ylabel('位点数 Sites'); ax.set_ylim(0,4200); ax.set_xlim(-0.6,2.6)
ax.set_yticks([0,1000,2000,3000,4000]); ax.set_yticklabels(['0','1 000','2 000','3 000','4 000'])
ax.set_title('两组合共同位点 SNPs shared by both crosses',loc='left',pad=4)
letter(ax,'B',-0.2,1.04)

ax=fig.add_subplot(gs[1,:]); aC=ax
w=0.36
for k,(g,off) in enumerate([('MD',-w/2),('MS',w/2)]):
    d=a[a.genotype==g].set_index('stage').loc[CTX]
    xs=np.arange(5)+off
    ax.bar(xs,d.ASE_percent,width=w*0.95,color=GCOL[g],label=g,edgecolor='none')
    for x,(n,pc) in zip(xs,zip(d.ASE_calls,d.ASE_percent)): ax.text(x,pc+0.06,str(n),ha='center',va='bottom',fontsize=6.3)
den=a[a.genotype=='MD'].set_index('stage').loc[CTX].jointly_testable_genes
ax.set_xticks(range(5)); ax.set_xticklabels([CTXLAB[c].replace('\n',' ')+f'\n$n$ = {n}' for c,n in zip(CTX,den)],fontsize=6.3)
ax.set_ylabel('ASE基因比例 (%)\nASE genes (%)'); ax.set_ylim(0,5.5); ax.set_xlim(-0.6,4.6)
ax.legend(frameon=False,loc='upper left',ncol=2,handlelength=1,handletextpad=0.4,columnspacing=1)
ax.set_title('组合内ASE（共享离散度模型）Within-hybrid ASE (shared-dispersion model)',loc='left',pad=4)
ax.text(1.0,-0.30,'$n$：两组合均可检验的基因数 Genes testable in both hybrids',transform=ax.transAxes,ha='right',va='top',fontsize=5.8,color='#444')
letter(ax,'C',-0.115,1.04)
save(fig,'Fig4_SNP_selection_ASE'); save_panels(fig,'Fig4',{'A':[aA],'B':[aB],'C':[aC]})
