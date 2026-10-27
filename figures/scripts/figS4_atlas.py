import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
import numpy as np, pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator, NullLocator
S8=sys.argv[1] if len(sys.argv)>1 else str(RES/'S8'); OUTPDF=sys.argv[2] if len(sys.argv)>2 else str(OUT/'supplementary'/'FigS4_qPCR_27_gene_atlas.pdf'); os.makedirs(os.path.dirname(OUTPDF),exist_ok=True)
rep=pd.read_csv(f'{S8}/qPCR_biological_replicate_values.tsv',sep='\t')
grp=pd.read_csv(f'{S8}/qPCR_group_summary_and_letters.tsv',sep='\t')
wel=pd.read_csv(f'{S8}/qPCR_Welch_ANOVA.tsv',sep='\t').set_index(['gene_id','stage'])
order=grp.groupby('gene_id').assay_number.first().sort_values().index.tolist()
G=['C','D','MD','S','MS']
CT={'seedling_leaf':'苗期叶片 SL','fruiting_leaf':'结瓜期叶片 FL','fruit_3d':'3 d果实 F3','fruit_6d':'6 d果实 F6','fruit_9d':'9 d果实 F9'}
def fmtp(p):
    if p<0.001:
        e=int(np.floor(np.log10(p))); m=p/10**e; return rf'$P$ = {m:.2f}$\times10^{{{e}}}$'
    return f'$P$ = {p:.3f}'
NPG=9   # 27 genes on 3 full pages, 9 genes x 5 sample types per page
pages=[order[i:i+NPG] for i in range(0,len(order),NPG)]
for f in Path(OUTPDF).parent.glob('FigS4_qPCR_atlas_p*'): f.unlink()   # drop page files from earlier layouts
with PdfPages(OUTPDF) as pdf:
    for pi,genes in enumerate(pages):
        fig,axs=plt.subplots(len(genes),5,figsize=(7.0,0.2+1.02*len(genes)),squeeze=False)
        plt.subplots_adjust(left=0.075,right=0.995,top=1-0.47/(0.2+1.02*len(genes)),bottom=0.2/(0.2+1.02*len(genes)),wspace=0.3,hspace=0.62)
        for r,g in enumerate(genes):
            for c,st in enumerate(CTX):
                ax=axs[r,c]
                d=rep[(rep.gene_id==g)&(rep.stage==st)]; s=grp[(grp.gene_id==g)&(grp.stage==st)].set_index('genotype')
                vals=[]
                for i,gt in enumerate(G):
                    v=d[d.genotype==gt].relative_expression.values; vals+=list(v)
                    ax.scatter(i+np.linspace(-0.12,0.12,len(v)),v,s=5,color=GCOL[gt],lw=0,zorder=3)
                    m,sd=s.loc[gt,'relative_expression_mean'],s.loc[gt,'relative_expression_sd']
                    lo=max(m-sd,m*0.05)
                    ax.errorbar(i,m,yerr=[[m-lo],[sd]],fmt='_',ms=6,mew=0.9,color='#222',elinewidth=0.55,capsize=1.4,zorder=4)
                    vals+=[m+sd,lo]
                ax.set_yscale('log',base=2)
                lo_,hi_=np.log2(min(vals)),np.log2(max(vals)); pad=max(0.25,(hi_-lo_)*0.16)
                ax.set_ylim(2**(lo_-pad),2**(hi_+pad*3.4))
                for i,gt in enumerate(G):
                    top=max(d[d.genotype==gt].relative_expression.max(),s.loc[gt,'relative_expression_mean']+s.loc[gt,'relative_expression_sd'])
                    ax.text(i,top*2**(pad*0.3),s.loc[gt,'letters_Games_Howell_005'],ha='center',va='bottom',fontsize=5.3)
                ax.axhline(1,c='#999',lw=0.5,ls='--',zorder=0)
                ax.text(0.02,0.99,'Welch '+fmtp(wel.loc[(g,st),'p_Welch']),transform=ax.transAxes,fontsize=5,va='top')
                ax.set_xticks(range(5)); ax.set_xticklabels(G,fontsize=5.3); ax.tick_params(axis='both',labelsize=5.3,length=2,pad=1.5)
                ax.yaxis.set_minor_locator(NullLocator())
                ax.set_xlim(-0.6,4.6)
                if r==0: ax.set_title(CT[st],fontsize=6.8,pad=17)
                if c==0: ax.text(0.0,1.07,itg(g),transform=ax.transAxes,fontsize=6.3,va='bottom',ha='left')
        fig.text(0.008,0.5,'相对表达量 Relative expression',rotation=90,va='center',ha='left',fontsize=6.5)
        pdf.savefig(fig)
        # PNG and editable SVG of each atlas page, and one PDF per gene row
        fig.savefig(Path(OUTPDF).parent/f'FigS4_qPCR_atlas_p{pi+1}.png',dpi=600)   # page images for the supplementary-figure Word file
        from io import BytesIO
        _b=BytesIO(); fig.savefig(_b,format='svg'); (Path(OUTPDF).parent/f'FigS4_qPCR_atlas_p{pi+1}.svg').write_bytes(svg_font_fallback(_b.getvalue()))
        save_panels(fig,'FigS4',{g:list(axs[r,:]) for r,g in enumerate(genes)},outdir=Path(OUTPDF).parent/'panels')
        plt.close(fig)
print(len(order),'genes',len(pages),'pages')
