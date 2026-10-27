import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
import numpy as np, pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
T=str(RES/'S8')
pe=pd.read_csv(f'{T}/qPCR_RNAseq_parental_effect_comparison.tsv',sep='\t'); pe=pe[pe.included_in_concordance==True]
mp=pd.read_csv(f'{T}/qPCR_RNAseq_log_midparent_comparison.tsv',sep='\t'); mp=mp[mp.included_in_concordance==True]
sm=pd.read_csv(f'{T}/qPCR_RNAseq_concordance_summary.tsv',sep='\t').set_index('comparison_type')
rep=pd.read_csv(f'{T}/qPCR_biological_replicate_values.tsv',sep='\t')
grp=pd.read_csv(f'{T}/qPCR_group_summary_and_letters.tsv',sep='\t')

fig=plt.figure(figsize=(W,5.75))
gs=GridSpec(2,4,figure=fig,height_ratios=[1,1.3],hspace=0.45,wspace=0.7,left=0.09,right=0.99,top=0.91,bottom=0.06)
axA=fig.add_subplot(gs[1,0:2]); axB=fig.add_subplot(gs[1,2:4])

# E-F: matched qRT-PCR and RNA-seq effects (all 27 genes)
for ax,(L,d,yc,key,title,lim) in zip([axA,axB],[('E',pe,'qPCR_log2FC','qPCR_log2FC','F$_1$相对亲本 F$_1$ vs parent',12),
                                              ('F',mp,'D_log_parental_mean','D_log_parental_mean','对数中亲偏离 Log-midparent departure',8)]):
    for g in ('MD','MS'):
        m=d.hybrid==g
        ax.scatter(d.log2FoldChange[m],d[yc][m],s=6,c=GCOL[g],alpha=0.45,lw=0,label=g,rasterized=True,zorder=3)
    ax.plot([-lim,lim],[-lim,lim],c='#777',lw=0.6,ls='--',zorder=2)
    ax.axhline(0,c='#DDD',lw=0.5,zorder=1); ax.axvline(0,c='#DDD',lw=0.5,zorder=1)
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_aspect('equal')
    tk=np.arange(-lim,lim+1,4 if lim>8 else 2); ax.set_xticks(tk); ax.set_yticks(tk)
    ax.set_xlabel('RNA-seq log$_2$ 效应 effect'); ax.set_ylabel('qRT-PCR log$_2$ 效应 effect')
    r=sm.loc[key]; rho=spearmanr(d.log2FoldChange,d[yc]).correlation
    assert abs(r.p_Holm_across_two_concordance_tests-4e-5)<1e-12; assert abs(rho-r.spearman_rho_descriptive)<1e-9 and len(d)==r.events
    ax.text(0.04,0.96,f'Spearman $\\rho$ = {r.spearman_rho_descriptive:.2f}\n$P_{{\\mathrm{{Holm}}}}$ = $4\\times10^{{-5}}$'+f'\n{int(r.events)} 事件 events / {int(r.genes)} 基因 genes',
            transform=ax.transAxes,ha='left',va='top',fontsize=6.2,linespacing=1.35)
    ax.set_title(title,loc='left',pad=4); letter(ax,L,-0.2,1.03)
    ax.legend(frameon=False,loc='lower right',markerscale=2,handletextpad=0.2)

# A-D: qRT-PCR relative expression of four genes named in the text (same data as Fig. S4)
SEL=[('CsaV3_4G002320','PAL','fruiting_leaf'),('CsaV3_3G033350','WRKY33-like','fruiting_leaf'),
     ('CsaV3_3G042150',None,'fruit_9d'),('CsaV3_6G039530',None,'fruit_9d')]
STL={'fruiting_leaf':'结瓜期叶片 FL','fruit_9d':'9 d果实 F9'}
G=['C','D','MD','S','MS']; MPV={'MD':('C','D'),'MS':('C','S')}
axs=[]
for k,(gid,sym,st) in enumerate(SEL):
    ax=fig.add_subplot(gs[0,k]); axs.append(ax)
    d=rep[(rep.gene_id==gid)&(rep.stage==st)]; s=grp[(grp.gene_id==gid)&(grp.stage==st)].set_index('genotype')
    assert len(d)==15 and set(s.index)==set(G)
    c0=s.loc['C','delta_ct_mean']; geo={g:2**-(s.loc[g,'delta_ct_mean']-c0) for g in G}   # group geometric means, C = 1
    vals=[]
    for i,gt in enumerate(G):
        v=d[d.genotype==gt].relative_expression.values; vals+=list(v)
        ax.scatter(i+np.linspace(-0.13,0.13,len(v)),v,s=8,color=GCOL[gt],lw=0,zorder=3)
        m,sd=s.loc[gt,'relative_expression_mean'],s.loc[gt,'relative_expression_sd']; lo=max(m-sd,m*0.05)
        ax.errorbar(i,m,yerr=[[m-lo],[sd]],fmt='_',ms=8,mew=1,color='#222',elinewidth=0.6,capsize=1.8,zorder=4)
        vals+=[m+sd,lo]
    for h,(p1,p2) in MPV.items():   # geometric mid-parent value for each hybrid
        x=G.index(h); mpv=np.sqrt(geo[p1]*geo[p2])
        ax.plot([x-0.32,x+0.32],[mpv,mpv],c='#444',lw=0.8,ls=(0,(2,1.2)),zorder=2)
    ax.set_yscale('log',base=2)
    lo_,hi_=np.log2(min(vals)),np.log2(max(vals)); pad=max(0.3,(hi_-lo_)*0.12)
    ax.set_ylim(2**(lo_-pad),2**(hi_+pad*3))
    for i,gt in enumerate(G):
        top=max(d[d.genotype==gt].relative_expression.max(),s.loc[gt,'relative_expression_mean']+s.loc[gt,'relative_expression_sd'])
        ax.text(i,top*2**(pad*0.35),s.loc[gt,'letters_Games_Howell_005'],ha='center',va='bottom',fontsize=6)
    ax.set_xticks(range(5)); ax.set_xticklabels(G,fontsize=6.2); ax.set_xlim(-0.6,4.6); ax.tick_params(axis='y',labelsize=6)
    ax.minorticks_off()
    if k==0: ax.set_ylabel('qRT-PCR 相对表达量\nRelative expression')
    ax.set_title(itg(gid)+('\n'+itg(sym) if sym else '\n'),loc='left',fontsize=6.4,pad=3,linespacing=1.3)
    letter(ax,'ABCD'[k],-0.34 if k==0 else -0.3,1.09)
# tissue headers spanning A-B and C-D
fig.canvas.draw()
for (a,b),st in [((axs[0],axs[1]),'fruiting_leaf'),((axs[2],axs[3]),'fruit_9d')]:
    x0=a.get_position().x0-0.02; x1=b.get_position().x1; y=a.get_position().y1+0.058
    fig.add_artist(plt.Line2D([x0,x1],[y,y],transform=fig.transFigure,c='#888',lw=0.6))
    fig.text((x0+x1)/2,y+0.005,STL[st],ha='center',va='bottom',fontsize=6.8)
from matplotlib.lines import Line2D
fig.legend([Line2D([],[],c='#444',lw=0.8,ls=(0,(2,1.2)))],['几何中亲值 Geometric mid-parent value'],frameon=False,loc='lower right',bbox_to_anchor=(0.99,axs[0].get_position().y0-0.068),fontsize=6.2,handlelength=2)
save(fig,'Fig6_qPCR_validation')
save_panels(fig,'Fig6',{**{'ABCD'[k]:[a] for k,a in enumerate(axs)},'E':[axA],'F':[axB]})
