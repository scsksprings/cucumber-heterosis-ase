import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
import numpy as np, pandas as pd
from matplotlib.lines import Line2D
T=str(RES)
ORDER=['C','D','MD','S','MS']

# ---------------- Fig S1 cumulative harvest ----------------
x=pd.read_excel(f'{T}/S2.xlsx',sheet_name='累计采收')
fig,axs=plt.subplots(1,2,figsize=(W,2.6)); plt.subplots_adjust(wspace=0.32,bottom=0.3)
for ax,(col,lab,L) in zip(axs,[('cumulative_yield_g_per_plant','累计单株产量 (g)\nCumulative yield per plant','A'),
                               ('cumulative_fruits_per_plant','累计单株瓜数\nCumulative fruits per plant','B')]):
    for g in ORDER:
        d=x[x.genotype==g]
        for b,db in d.groupby('block'): ax.plot(db.harvest_index,db[col],c=GCOL[g],lw=0.5,alpha=0.45)
        m=d.groupby('harvest_index')[col].mean(); ax.plot(m.index,m.values,c=GCOL[g],lw=1.6,label=g)
    ax.set_xlabel('采收次序 Harvest order'); ax.set_ylabel(lab); ax.set_xticks([1,4,7,10,14]); ax.set_xlim(0.6,14.4); ax.set_ylim(bottom=0)
    letter(ax,L,-0.22,1.02)
axs[0].legend(loc='upper center',bbox_to_anchor=(1.12,-0.3),ncol=5,frameon=False,handlelength=1.5,columnspacing=1.2)
fig.text(0.98,0.01,'细线：各区组 Thin lines: blocks；粗线：区组均值 Thick lines: block means',ha='right',fontsize=6.3,color='#444')
save(fig,'FigS1_cumulative_harvest',outdir=OUT/'supplementary'); save_panels(fig,'FigS1',{'A':[axs[0]],'B':[axs[1]]},outdir=OUT/'supplementary'/'panels')

# ---------------- Fig S2 eight jointly supported events ----------------
e=pd.read_csv(f'{T}/S7/S7_core_ASE_two_model_full.tsv',sep='\t')
CT={'seedling_leaf':'SL','fruiting_leaf':'FL','fruit_3d':'F3','fruit_6d':'F6','fruit_9d':'F9'}
e['lab']=e.gene_id.map(itg)+' / '+e.stage.map(CT)
e=e.sort_values('delta_log2_ratio_MS_minus_MD_Shared').reset_index(drop=True)
fig,axs=plt.subplots(1,2,figsize=(W,2.9),gridspec_kw={'width_ratios':[1.25,1]}); plt.subplots_adjust(left=0.25,right=0.98,wspace=0.08,bottom=0.2,top=0.86)
ax=axs[0]; y=np.arange(len(e))
for i,r in e.iterrows():
    ax.plot([r.MD_maternal_fraction_Shared,r.MS_maternal_fraction_Shared],[i,i],c='#BBBBBB',lw=0.8,zorder=1)
ax.scatter(e.MD_maternal_fraction_Shared,y,c=GCOL['MD'],s=18,zorder=3,label='MD')
ax.scatter(e.MS_maternal_fraction_Shared,y,c=GCOL['MS'],s=18,zorder=3,label='MS')
ax.axvline(0.5,c='#999',lw=0.5,ls='--'); ax.set_xlim(0,1); ax.set_yticks(y); ax.set_yticklabels(e.lab,fontsize=6.8)
ax.set_xlabel('母本C等位基因表达比例（共享模型）\nMaternal C-allele fraction (shared model)'); ax.legend(loc='lower left',bbox_to_anchor=(0,1.0),ncol=2,frameon=False,handletextpad=0.2)
letter(ax,'A',-0.62,1.06)
ax=axs[1]
ax.scatter(e.delta_log2_ratio_MS_minus_MD_Shared,y,c='#333333',s=16,zorder=3)
ax.axvline(0,c='#999',lw=0.5,ls='--'); ax.set_yticks(y); ax.set_yticklabels([]); ax.set_xlabel('$\\Delta R$（MS − MD，log$_2$）\nShared-model $\\Delta R$')
letter(ax,'B',-0.08,1.06)
save(fig,'FigS2_core_ASE_shared_effects',outdir=OUT/'supplementary'); save_panels(fig,'FigS2',{'A':[axs[0]],'B':[axs[1]]},outdir=OUT/'supplementary'/'panels')
e[['gene_id','stage','MD_maternal_fraction_Shared','MS_maternal_fraction_Shared','delta_log2_ratio_MS_minus_MD_Shared']].to_csv(OUT/'supplementary'/'FigS2_values.tsv',sep='\t',index=False)

# ---------------- Fig S3 model calls and null diagnostics ----------------
m=pd.read_csv(f'{T}/S6/TableS11_ASE_model_diagnostics.tsv',sep='\t'); m=m[m.profile=='all']
cal=pd.read_csv(f'{T}/S6/S6_null_family_index_40.tsv',sep='\t')
keys=[(h,st) for h in ['MD','MS'] for st in CTX]
fig,axs=plt.subplots(1,2,figsize=(W,3.3),gridspec_kw={'width_ratios':[1,1.3]}); plt.subplots_adjust(wspace=0.3,bottom=0.38,top=0.9)
ax=axs[0]
for meth,c,mk,lab in [('legacy_gene_wise_MLE','#6E6E6E','s','逐基因 Gene-wise MLE'),('shared_dispersion_simulation','#0072B2','o','共享离散度 Shared dispersion')]:
    v=m[m.model==meth].set_index(['genotype','stage']).loc[keys,'called'].to_numpy()
    ax.plot(range(10),v,c=c,marker=mk,ms=3.5,lw=1,label=lab)
ax.set_xticks(range(10)); ax.set_xticklabels([f'{h}\n{CT[s]}' for h,s in keys],fontsize=6.3)
ax.set_ylabel('组合内ASE检出数\nWithin-hybrid ASE calls'); ax.set_title('全部亲本差异SNP All diagnostic SNPs',loc='left',pad=4)
ax.set_ylim(0,300); ax.legend(frameon=False,loc='upper left',fontsize=6.5)
letter(ax,'A',-0.25,1.04)
ax=axs[1]
pal={'all':('#6E6E6E','全部亲本差异位点 All diagnostic sites'),'common3':('#0072B2','共同位点，亲本支持文库数≥3 Shared sites, parental support ≥3'),'common4':('#D55E00','共同位点，亲本支持文库数≥4 Shared sites, parental support ≥4')}
for pf,(c,lab) in pal.items():
    d=cal[cal.profile==pf]
    ax.errorbar(d.figure_S3_x,d.rejection_fraction_005,yerr=np.stack([d.rejection_fraction_005-d.rejection_ci_low,d.rejection_ci_high-d.rejection_fraction_005]),
                fmt='o',ms=2.8,color=c,elinewidth=0.6,capsize=1,label=lab)
ax.axhline(0.05,c='#777',ls='--',lw=0.7); ax.axhline(0.08,c='#BD3035',ls=':',lw=0.8)
ax.set_xlim(0,41); ax.set_ylim(0,0.12); ticks=[1,10,20,30,40]; ax.set_xticks(ticks); ax.set_xticklabels([f'F{t:02d}' for t in ticks])
ax.set_xlabel('模拟组编号 Simulation group index（表S6）'); ax.set_ylabel('零假设拒绝比例 ($P$<0.05)\nNull rejection fraction')
ax.set_title('零假设模拟诊断（95% CI）Null diagnostics',loc='left',pad=4)
ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(-0.05,-0.25),ncol=1,fontsize=6.5,handletextpad=0.2,labelspacing=0.3)
letter(ax,'B',-0.2,1.04)
save(fig,'FigS3_model_calibration',outdir=OUT/'supplementary'); save_panels(fig,'FigS3',{'A':[axs[0]],'B':[axs[1]]},outdir=OUT/'supplementary'/'panels')

# ---------------- Fig S5 full product decomposition ----------------
d=pd.read_excel(f'{T}/S2.xlsx',sheet_name='乘积完整分解').set_index('hybrid')
terms=[('exact_parent_complement_term_g','亲本互补项\nParental\ncomplement'),('fruit_number_departure_term_g','单株瓜数偏离项\nFruit number\ndeparture'),
       ('fruit_mass_departure_term_g','单瓜重偏离项\nFruit mass\ndeparture'),('departure_product_term_g','偏离乘积项\nProduct term'),
       ('block_covariance_term_g','协方差校正\nCovariance\ncorrection')]
fig,axs=plt.subplots(2,1,figsize=(W,4.4),sharex=True); plt.subplots_adjust(hspace=0.3,bottom=0.14,left=0.13)
for ax,(h,L) in zip(axs,[('MD','A'),('MS','B')]):
    r=d.loc[h]; cum=0
    for i,(c,lab) in enumerate(terms):
        v=r[c]; col='#7F7F7F' if v>=0 else '#C9C9C9'
        ax.bar(i,v,bottom=cum,color=col,width=0.62,edgecolor='#555',lw=0.3)
        ax.text(i,cum+max(v,0)+8,f'{v:+.1f}',ha='center',va='bottom',fontsize=6.3)
        if i<len(terms)-1: ax.plot([i+0.31,i+0.69],[cum+v,cum+v],c='#888',lw=0.5)
        cum+=v
    ax.bar(len(terms),r.observed_MPH_difference_g,color=GCOL[h],width=0.62)
    ax.text(len(terms),r.observed_MPH_difference_g+8,f'{r.observed_MPH_difference_g:.1f}',ha='center',va='bottom',fontsize=6.3)
    assert abs(cum-r.observed_MPH_difference_g)<1e-6
    ax.axhline(0,c='#333',lw=0.5)
    ax.set_xticks(range(len(terms)+1)); ax.set_xticklabels([t[1] for t in terms]+['观察中亲差\nObserved\ndifference'],fontsize=6.5); ax.set_ylabel('g/株 g per plant'); ax.set_ylim(min(0,ax.get_ylim()[0]),max(ax.get_ylim()[1],r.observed_MPH_difference_g)*1.15)
    ax.set_title(h,loc='left',color=GCOL[h],fontweight='bold',pad=4); letter(ax,L,-0.1,1.04)
save(fig,'FigS5_product_decomposition',outdir=OUT/'supplementary'); save_panels(fig,'FigS5',{'A':[axs[0]],'B':[axs[1]]},outdir=OUT/'supplementary'/'panels')
print('ok')
