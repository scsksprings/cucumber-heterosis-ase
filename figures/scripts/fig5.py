import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
import numpy as np, pandas as pd
T=str(RES)
c=pd.read_csv(f'{T}/S7/TableS9_representative_fragment_counts.tsv',sep='\t')
rho=pd.read_csv(f'{T}/S6/TableS24_shared_dispersion_estimates.tsv',sep='\t')
rho=rho[rho.profile=='common3'].set_index(['genotype','stage','held_out_fold']).rho
FOLD={'CsaV3_7G007510':2,'CsaV3_6G002590':0}   # fold identified from the rho reported for each gene in S7
# cross-hybrid within-context BH q from Table S13 (support>=3): shared model (padj_shared3) and gene-wise model (padj_MLE3)
t13=pd.read_csv(f'{T}/S6/TableS13_parent_support_sensitivity.tsv',sep='\t').set_index(['gene_id','stage'])
import math
def fq(q):
    if q<0.001:
        e=math.floor(math.log10(q)); m=q/10**e
        return r'$\mathbf{%.1f\times10^{%d}}$'%(m,e)
    if q>=0.995: return '1'
    return f'{q:.2g}' if q<0.1 else f'{q:.2f}'
Q={}; QRAW=[]
for g in ('CsaV3_7G007510','CsaV3_6G002590'):
    Q[g]={}
    for st in CTX:
        if (g,st) in t13.index:
            r=t13.loc[(g,st)]; Q[g][st]=(fq(r.padj_shared3),fq(r.padj_MLE3),r.padj_shared3<0.05,r.padj_MLE3<0.05); QRAW.append((g,st,r.padj_shared3,r.padj_MLE3))
        else: Q[g][st]=('NT','NT',False,False); QRAW.append((g,st,None,None))
OFF={'MD':-0.18,'MS':0.18}
fig,axs=plt.subplots(2,1,figsize=(W,5.15))
plt.subplots_adjust(hspace=1.35,left=0.24,right=0.98,top=0.905,bottom=0.215)
rows=[]
for ax,(L,g) in zip(axs,[('A','CsaV3_7G007510'),('B','CsaV3_6G002590')]):
    ax.axhline(0.5,c='#999',lw=0.5,ls='--',zorder=0)
    info={}
    for j,st in enumerate(CTX):
        for gt in ('MD','MS'):
            d=c[(c.gene_id==g)&(c.genotype==gt)&(c.stage==st)].sort_values('replicate')
            x0=j+OFF[gt]; jit=np.array([-0.06,0,0.06])[:len(d)]
            for xx,(_,r) in zip(x0+jit,d.iterrows()):
                if r.total_depth>=5: ax.scatter(xx,r.C_fraction,s=13,c=GCOL[gt],lw=0,zorder=3)
                elif r.total_depth>=1: ax.scatter(xx,r.C_fraction,s=13,facecolors='white',edgecolors=GCOL[gt],lw=0.7,zorder=3)
            k=d[d.total_depth>=5]
            tot=int(d.total_depth.sum()); nq=len(k); nz=int((d.total_depth==0).sum())
            est=np.nan
            if nq>=2 and k.total_depth.sum()>=30:
                rr=rho[(gt,st,FOLD[g])]; w=k.total_depth/(1+(k.total_depth-1)*rr); est=float((w*k.C_fraction).sum()/w.sum())
                ax.plot([x0-0.12,x0+0.12],[est,est],c=GCOL[gt],lw=1.6,solid_capstyle='butt',zorder=4)
            info[(st,gt)]=(tot,nq,nz); rows.append((g,st,gt,tot,nq,nz,est))
    ax.set_xlim(-0.55,4.55); ax.set_ylim(-0.03,1.03)
    ax.set_yticks([0,0.25,0.5,0.75,1]); ax.set_yticklabels(['0','0.25','0.50','0.75','1'])
    ax.set_ylabel('母本C等位基因表达比例\nMaternal C-allele fraction')
    ax.set_xticks(range(5)); ax.set_xticklabels([CTXLAB[s].replace('\n',' ') for s in CTX],fontsize=7)
    ax.tick_params(axis='x',length=0,pad=3)
    ax.set_title(itg(g),loc='left',pad=4)
    letter(ax,L,-0.10,1.03)
    # information band (positions in points below the axis, independent of axes height)
    from matplotlib.transforms import offset_copy, blended_transform_factory
    labs=['片段数 Fragments','≥5片段文库 Libraries ≥5','零覆盖文库 Zero-depth','$q_{\\mathrm{s}}$ 共享模型 Shared','$q_{\\mathrm{g}}$ 逐基因 Gene-wise']
    Y0=34; DY=10.5
    def tr(y_pt,kind):
        base=ax.transAxes if kind=='ax' else blended_transform_factory(ax.transData,ax.transAxes)
        return offset_copy(base,fig=fig,x=0,y=-y_pt,units='points')
    ax.text(-0.02,0,'MD/MS',transform=tr(Y0-DY*0.95,'ax'),ha='right',va='center',fontsize=6.5,color='#666')
    for i,lb in enumerate(labs):
        yp=Y0+i*DY
        ax.text(-0.02,0,lb,transform=tr(yp,'ax'),ha='right',va='center',fontsize=6.8,color='#333')
        for j,st in enumerate(CTX):
            md,ms=info[(st,'MD')],info[(st,'MS')]
            if i<3: t=f'{md[i]}/{ms[i]}'
            else: t=Q[g][st][i-3]
            bold=(i>=3 and Q[g][st][i-1])
            ax.text(j,0,t,transform=tr(yp,'mix'),ha='center',va='center',fontsize=7,fontweight='bold' if bold else 'normal',
                    color='#999' if t=='NT' else '#1A1A1A')
    ax.plot([0,1],[0,0],transform=tr(Y0-DY*0.6,'ax'),c='#BBB',lw=0.4,clip_on=False)
pd.DataFrame(rows,columns=['gene_id','stage','genotype','fragments_common_SNPs','libraries_depth_ge5','zero_depth_libraries','shared_model_estimate']).to_csv(OUT/'Fig5_values.tsv',sep='\t',index=False)
from matplotlib.lines import Line2D
h=[Line2D([],[],marker='o',ls='',mfc=GCOL['MD'],mec='none',ms=4,label='MD'),Line2D([],[],marker='o',ls='',mfc=GCOL['MS'],mec='none',ms=4,label='MS'),
   Line2D([],[],marker='o',ls='',mfc='white',mec='#555',ms=4,label='深度1–4 Depth 1–4'),Line2D([],[],c='#555',lw=1.6,label='共享模型估计 Shared-model estimate')]
axs[0].legend(handles=h,loc='lower left',bbox_to_anchor=(-0.02,1.10),ncol=4,frameon=False,handletextpad=0.3,columnspacing=1.0,fontsize=6.8)
pd.DataFrame(QRAW,columns=['gene_id','stage','q_shared','q_genewise']).to_csv(OUT/'Fig5_q_values.tsv',sep='\t',index=False)
save(fig,'Fig5_candidate_allele_fractions'); save_panels(fig,'Fig5',{'A':[axs[0]],'B':[axs[1]]})
