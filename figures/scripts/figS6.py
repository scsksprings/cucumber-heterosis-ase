import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
import numpy as np, pandas as pd
from matplotlib.colors import Normalize
import matplotlib.cm as cm
S3=str(RES/'S3')
enr=pd.read_csv(f'{S3}/functional_enrichment_all_tests.tsv.gz',sep='\t')
TERMS=[('GO:0006457','蛋白质折叠','Protein folding'),('KEGG:ko00052','半乳糖代谢','Galactose metabolism'),
       ('GO:0008137','NADH脱氢酶（泛醌）活性','NADH dehydrogenase (ubiquinone) activity'),
       ('KEGG:ko04141','内质网蛋白加工','Protein processing in ER'),
       ('GO:0010997','后期促进复合物结合','Anaphase-promoting complex binding'),('GO:0008061','几丁质结合','Chitin binding')]
COLS=[('MD_F9_vs_F6','positive','MD +'),('MD_F9_vs_F6','negative','MD −'),('MS_F9_vs_F6','positive','MS +'),('MS_F9_vs_F6','negative','MS −')]
cmap=matplotlib.colormaps['viridis_r']; norm=Normalize(1.3,5)
fig,ax=plt.subplots(figsize=(W,3.3)); plt.subplots_adjust(left=0.42,right=0.97,top=0.85,bottom=0.2)
rows=[]
for i,(tid,cn,en) in enumerate(TERMS):
    tid2=tid.split(':',1)[1] if tid.startswith('KEGG') else tid
    for j,(st,dr,lab) in enumerate(COLS):
        r=enr[(enr.set==st)&(enr.direction==dr)&(enr.term_id==tid2)]; assert len(r)==1; r=r.iloc[0]
        k,q=int(r.overlap_n),float(r.q_BH); rows.append((tid2,cn,lab,k,r.fold_enrichment,r.pvalue,q,r.q_global_sensitivity))
        if k==0: continue
        s=8+k*9
        if q<0.05: ax.scatter(j,i,s=s,c=[cmap(norm(-np.log10(q)))],edgecolors='#333',linewidths=0.3,zorder=3)
        else: ax.scatter(j,i,s=s,facecolors='none',edgecolors='#A6A6A6',linewidths=0.6,zorder=3)
pd.DataFrame(rows,columns=['term_id','term_cn','column','k','fold_enrichment','p','q_family','q_global']).to_csv(OUT/'supplementary'/'FigS6_values.tsv',sep='\t',index=False)
ax.set_xlim(-0.6,3.6); ax.set_ylim(len(TERMS)-0.5,-0.5)
ax.set_xticks(range(4)); ax.set_xticklabels([c[2] for c in COLS],fontsize=7); ax.xaxis.tick_top(); ax.tick_params(length=0)
ax.set_yticks(range(len(TERMS))); ax.set_yticklabels([f'{cn}\n{en}' for _,cn,en in TERMS],fontsize=7)
for s_ in ax.spines.values(): s_.set_visible(False)
for i in range(len(TERMS)): ax.axhline(i,c='#EEEEEE',lw=0.4,zorder=0)
ax.set_title('9 d对6 d果实的偏离变化 Departure change, F9 vs F6',loc='left',pad=16,fontsize=7)
ci=ax.inset_axes([0,-0.16,0.45,0.04])   # Axes.inset_axes keeps the colour bar in place in PDF/SVG
cb=fig.colorbar(cm.ScalarMappable(norm,cmap),cax=ci,orientation='horizontal'); cb.outline.set_linewidth(0.4)
cb.set_ticks([2,3,4,5]); ci.tick_params(labelsize=7,length=1.5,width=0.4)
ci.set_title('−log$_{10}$ $q$（检验族内 within family）',fontsize=7,pad=2,loc='left')
for kk,xx in zip((1,3,10),(0.62,0.74,0.88)):
    ax.scatter(xx,-0.14,s=8+kk*9,facecolors='#BBBBBB',edgecolors='#333',linewidths=0.3,transform=ax.transAxes,clip_on=False)
    ax.text(xx,-0.215,str(kk),transform=ax.transAxes,ha='center',va='top',fontsize=7)
ax.scatter(1.0,-0.14,s=35,facecolors='none',edgecolors='#A6A6A6',linewidths=0.6,transform=ax.transAxes,clip_on=False)
ax.text(1.0,-0.215,'$q\\geq0.05$',transform=ax.transAxes,ha='center',va='top',fontsize=7)
ax.text(0.62,-0.06,'条目命中基因数 $k$ Genes',transform=ax.transAxes,fontsize=7)
save(fig,'FigS6_fruit_stage_enrichment',outdir=OUT/'supplementary')
