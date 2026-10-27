import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from style import *
import numpy as np, pandas as pd
from matplotlib import gridspec
from matplotlib.colors import Normalize, TwoSlopeNorm
import matplotlib.cm as cm
S3=str(RES/'S3')
o=pd.read_csv(f'{S3}/S3_leaf_overlap_all_tested.tsv.gz',sep='\t').set_index('gene_id')
enr=pd.read_csv(f'{S3}/functional_enrichment_all_tests.tsv.gz',sep='\t')
enr=enr[enr.set.isin(['MD_FL','MS_FL'])]
sh=pd.read_csv(f'{S3}/S3_shared_direction_enrichment.tsv',sep='\t')

fig=plt.figure(figsize=(W,8.6))
gs=gridspec.GridSpec(3,2,figure=fig,height_ratios=[0.95,0.9,1.02],width_ratios=[1,1.05],hspace=0.40,wspace=0.12,left=0.12,right=0.88,top=0.975,bottom=0.115)
FS=7

# ---------- A scatter ----------
ax=fig.add_subplot(gs[0,0]); axA=ax
LIM=7
x=o.log2FoldChange_midparent_MD.clip(-LIM+0.1,LIM-0.1); y=o.log2FoldChange_midparent_MS.clip(-LIM+0.1,LIM-0.1)
cat=o.overlap_category
sty=[('neither_called','#D9D9D9',1.2,0.6,'两组合均未达阈值 Neither (14 545)'),
     ('MD_only_called',GCOL['MD'],2.2,0.55,'仅MD MD only (1 492)'),
     ('MS_only_called',GCOL['MS'],2.2,0.6,'仅MS MS only (586)'),
     ('both_same_direction','#1A1A1A',3.5,0.85,'共同、同向 Shared, same (256)'),
     ('both_opposite_direction','#CC79A7',4.5,0.95,'共同、反向 Shared, opposite (59)')]
for k,c,s,a,lab in sty:
    m=cat==k; ax.scatter(x[m],y[m],s=s,c=c,alpha=a,lw=0,label=lab,rasterized=True,zorder=2)
for v in (0,): ax.axhline(v,c='#555',lw=0.5,zorder=1); ax.axvline(v,c='#555',lw=0.5,zorder=1)
for v in (-1,1): ax.axhline(v,c='#999',lw=0.4,ls='--',zorder=1); ax.axvline(v,c='#999',lw=0.4,ls='--',zorder=1)
ax.set_xlim(-LIM,LIM); ax.set_ylim(-LIM,LIM); ax.set_aspect('equal')
ax.set_xticks(range(-6,7,2)); ax.set_yticks(range(-6,7,2))
ax.set_xlabel('MD 表达偏离 $D_{\\mathrm{MD}}$ / Departure (log$_2$)',fontsize=FS); ax.set_ylabel('MS 表达偏离 $D_{\\mathrm{MS}}$ / Departure (log$_2$)',fontsize=FS); ax.tick_params(labelsize=6.8)
ax.set_title('结瓜期叶片 FL（$n$ = 16 938）',loc='left',pad=4,fontsize=7.5)
kw=dict(fontsize=7.5,fontweight='bold')
ax.text(6.6,6.6,'74',ha='right',va='top',color='#1A1A1A',**kw)
ax.text(-6.6,-6.6,'182',ha='left',va='bottom',color='#1A1A1A',**kw)
ax.text(6.6,-6.6,'7',ha='right',va='bottom',color='#CC79A7',**kw)
ax.text(-6.6,6.6,'52',ha='left',va='top',color='#CC79A7',**kw)
leg=ax.legend(loc='center left',bbox_to_anchor=(1.12,0.5),ncol=1,frameon=False,markerscale=2.8,handletextpad=0.3,
              labelspacing=0.7,borderaxespad=0,fontsize=FS,title='各组合分别判定 Calls per hybrid\n$q<0.05$且/and $|D|\\geq1$\n角标：共有基因数 Shared gene counts',title_fontsize=6.5)
leg._legend_box.align='left'
letter(ax,'A',-0.3,1.03)

# ---------- B dot matrix ----------
ax=fig.add_subplot(gs[1,1]); axB=ax
TERMS=[('KEGG','ko00196','光合作用天线蛋白','Antenna proteins'),
       ('KEGG','ko04141','内质网蛋白加工','Protein processing in ER'),
       ('GO','GO:0055114','氧化还原过程','Oxidation–reduction'),
       ('GO','GO:0004601','过氧化物酶活性','Peroxidase activity'),
       ('KEGG','ko00940','苯丙烷生物合成','Phenylpropanoid biosynthesis'),
       ('GO','GO:0045548','苯丙氨酸解氨酶活性','PAL activity'),
       ('KEGG','ko04626','植物–病原互作','Plant–pathogen interaction'),
       ('GO','GO:0003700','序列特异转录因子活性','TF activity (seq.-specific)'),
       ('GO','GO:0010143','角质生物合成','Cutin biosynthesis')]
COLS=[('MD_FL','positive','MD +'),('MD_FL','negative','MD −'),('MS_FL','positive','MS +'),('MS_FL','negative','MS −'),
      ('shared_positive',None,'共同 +'),('shared_negative',None,'共同 −')]
XPOS=[0,1,2,3,4.6,5.95]
cmap=matplotlib.colormaps['viridis_r']; norm=Normalize(1.3,14)   # lower end = q 0.05; max observed -log10 q = 13.1
rows=[]
for i,(db,tid,cn,en) in enumerate(TERMS):
    for j,(st,dr,_) in enumerate(COLS):
        if dr: r=enr[(enr.set==st)&(enr.direction==dr)&(enr.term_id==tid)]
        else: r=sh[(sh.set==st)&(sh.term_id==tid)]
        assert len(r)==1,(st,dr,tid); r=r.iloc[0]
        k,q=int(r.overlap_n),float(r.q_BH); rows.append((tid,cn,COLS[j][2].replace('\n',' '),k,q,float(r.pvalue)))
        if k==0: continue
        s=8+k*2.6
        if q<0.05: ax.scatter(XPOS[j],i,s=s,c=[cmap(norm(-np.log10(q)))],edgecolors='#333',linewidths=0.3,zorder=3)
        else: ax.scatter(XPOS[j],i,s=s,facecolors='none',edgecolors='#A6A6A6',linewidths=0.6,zorder=3)
pd.DataFrame(rows,columns=['term_id','term_cn','column','k','q_family','p']).to_csv(OUT/'Fig3B_values.tsv',sep='\t',index=False)
ax.set_xlim(-0.6,6.5); ax.set_ylim(len(TERMS)-0.5,-0.5)
ax.set_xticks(XPOS); ax.set_xticklabels([c[2] for c in COLS],fontsize=FS)
ax.xaxis.tick_top(); ax.tick_params(axis='x',length=0,pad=2)
ax.set_yticks(range(len(TERMS))); ax.set_yticklabels([f'{cn}  {en}' for _,_,cn,en in TERMS],fontsize=FS)
ax.tick_params(axis='y',length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
for i in range(len(TERMS)): ax.axhline(i,c='#EEEEEE',lw=0.4,zorder=0)
ax.axvline(3.85,c='#888',lw=0.5,ls=':',zorder=1)
# legends
cax=fig.add_axes([0,0,1,1]); cax.remove()
# Axes.inset_axes (not mpl_toolkits inset_axes): the toolkit locator misplaces rasterized artists in PDF/SVG saved with bbox_inches='tight'
ciB=ci=ax.inset_axes([1.06,0.58,0.05,0.38])
cb=fig.colorbar(cm.ScalarMappable(norm,cmap),cax=ci,orientation='vertical'); cb.outline.set_linewidth(0.4)
cb.set_ticks([2,6,10,14]); ci.tick_params(labelsize=6.5,length=1.5,width=0.4)
ci.set_title('−log$_{10}$ $q$',fontsize=6.5,pad=7,loc='left')
ax.text(1.06,0.39,'条目命中基因数\nGenes ($k$)',transform=ax.transAxes,fontsize=6.5,ha='left',va='bottom')
for kk,yy in zip((5,20,60),(0.33,0.22,0.08)):
    ax.scatter(1.10,yy,s=8+kk*2.6,facecolors='#BBBBBB',edgecolors='#333',linewidths=0.3,transform=ax.transAxes,clip_on=False)
    ax.text(1.17,yy,str(kk),transform=ax.transAxes,ha='left',va='center',fontsize=6.5)
ax.scatter(1.10,-0.04,s=40,facecolors='none',edgecolors='#A6A6A6',linewidths=0.6,transform=ax.transAxes,clip_on=False)
ax.text(1.17,-0.04,'$q\\geq0.05$',transform=ax.transAxes,ha='left',va='center',fontsize=6.5)
fig.canvas.draw()
r_=fig.canvas.get_renderer(); lx=min(t.get_window_extent(r_).x0 for t in ax.get_yticklabels())
bb=ax.get_window_extent(r_); lab_x=(lx-bb.x0)/bb.width-0.02
letter(ax,'B',lab_x,1.06)
for (a0,a1,t) in [(0,3,'各组合 Per hybrid'),(4.6,5.95,'共有 Shared')]:
    ax.plot([a0-0.35,a1+0.35],[-1.50,-1.50],c='#333',lw=0.6,clip_on=False)
    ax.text((a0+a1)/2,-1.65,t,ha='center',va='bottom',fontsize=FS)

# ---------- C heatmap ----------
ax=fig.add_subplot(gs[2,:]); axC=ax
s_neg=sh[sh.set=='shared_negative'].set_index('term_id')
pal=s_neg.loc['GO:0045548','gene_ids'].split(';'); ppi=s_neg.loc['ko04626','gene_ids'].split(';')
ant=enr[(enr.set=='MD_FL')&(enr.direction=='positive')&(enr.term_id=='ko00196')].gene_ids.iloc[0].split(';')
assert len(set(pal)|set(ppi)|set(ant))==23
antl=pd.read_csv(f'{S3}/S3_antenna_genes_MD_MS_FL.tsv',sep='\t').set_index('gene_id').homology_label
def sym(g):
    k=str(o.loc[g,'KEGG_annotation'])
    if g in ant: return antl.get(g,'LHC')
    for key,lab in [('K10775','PAL'),('K13424','WRKY33-like'),('K13448','CML'),('K13447','RBOH'),('K13412','CDPK'),('K13420','FLS2-like')]:
        if key in k: return lab
    return '?'
order_ppi=sorted(ppi,key=lambda g:(['WRKY33-like','RBOH','CDPK','CML','FLS2-like'].index(sym(g)),g))
ant=sorted(ant,key=lambda g:(antl[g],g))
genes=pal+order_ppi+ant
ROWS=[('midparent_MD','$D_{\\mathrm{MD}}$'),('midparent_MS','$D_{\\mathrm{MS}}$'),('hybrid_C_MD','MD / C'),('hybrid_paternal_MD','MD / D'),('hybrid_C_MS','MS / C'),('hybrid_paternal_MS','MS / S')]
YP=[0,1,2.35,3.35,4.7,5.7]
xs=[]; xx=0
for gi,g in enumerate(genes):
    if gi in (7,15): xx+=0.6
    xs.append(xx); xx+=1
M=np.array([[o.loc[g,f'log2FoldChange_{r}'] for g in genes] for r,_ in ROWS])
Q=np.array([[o.loc[g,f'padj_{r}'] for g in genes] for r,_ in ROWS])
VM=4.5; hn=TwoSlopeNorm(0,-VM,VM); hc=matplotlib.colormaps['RdBu_r']
out=[]
for i,(r,lab) in enumerate(ROWS):
    for j,g in enumerate(genes):
        v=M[i,j]; ax.add_patch(plt.Rectangle((xs[j]-0.5,YP[i]-0.5),1,1,fc=hc(hn(np.clip(v,-VM,VM))),ec='white',lw=0.5))
        sig=(Q[i,j]<0.05) and abs(v)>=1
        if sig: ax.plot(xs[j],YP[i],'o',ms=2.2,mfc='white' if abs(v)>2.6 else '#1A1A1A',mec='none')
        out.append((g,sym(g),r,v,Q[i,j],sig))
pd.DataFrame(out,columns=['gene_id','label','comparison','log2FC','padj','dot_q<0.05_and_abs>=1']).to_csv(OUT/'Fig3C_values.tsv',sep='\t',index=False)
ax.set_xlim(-0.7,xs[-1]+0.7); ax.set_ylim(YP[-1]+0.7,-2.1)
ax.set_yticks(YP); ax.set_yticklabels([l for _,l in ROWS],fontsize=FS)
ax.set_xticks(xs); ax.set_xticklabels([itg(f'{sym(g)}  {g}') for g in genes],rotation=90,fontsize=6.8)
ax.tick_params(length=0,pad=1.5)
for s_ in ax.spines.values(): s_.set_visible(False)
for (a0,a1,t) in [(0,6,'苯丙氨酸解氨酶 PAL activity'),(7,14,'植物–病原互作 Plant–pathogen'),(15,22,'光合作用天线蛋白 Antenna proteins')]:
    ax.plot([xs[a0]-0.4,xs[a1]+0.4],[-1.0,-1.0],c='#333',lw=0.6,clip_on=False)
    ax.text((xs[a0]+xs[a1])/2,-1.15,t,ha='center',va='bottom',fontsize=FS)
ax.text((xs[0]+xs[14])/2,-1.85,'共同负向偏离集合 Shared negative set',ha='center',va='bottom',fontsize=FS,color='#444')
ax.text((xs[15]+xs[22])/2,-1.85,'MD正向偏离集合 MD positive set',ha='center',va='bottom',fontsize=FS,color='#444')
ciC=ci=ax.inset_axes([1.03,0.33,0.016,0.5])
assert np.nanmax(np.abs(M))<=VM   # no value beyond the colour range, so no extension arrows
cb=fig.colorbar(cm.ScalarMappable(hn,hc),cax=ci); cb.outline.set_linewidth(0.4)
cb.set_ticks([-4,-2,0,2,4]); cb.set_ticklabels(['−4','−2','0','2','4']); ci.tick_params(labelsize=6.5,length=1.5,width=0.4)
ci.set_title('log$_2$效应\nlog$_2$ effect',fontsize=6.5,pad=3,loc='left')
ax.text(1.02,0.0,'● $q<0.05$且\n|log$_2$效应|≥1\n|log$_2$ effect|≥1',transform=ax.transAxes,fontsize=6.5,va='bottom')
ax.set_title('共同负向偏离相关基因与光合天线基因 Genes of the shared negative set and antenna genes',loc='left',pad=24,fontsize=7.5)
letter(ax,'C',-0.085,1.0)
fig.canvas.draw()
ax.title.set_position((0,1))
G={'A':[axA],'B':[axB,ciB],'C':[axC,ciC]}
print('overlap',check_no_overlap(fig,G))
save(fig,'Fig3_leaf_comparison'); save_panels(fig,'Fig3',G)
