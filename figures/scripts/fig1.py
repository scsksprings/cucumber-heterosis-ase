import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import csv
import numpy as np
from style import *
from matplotlib.gridspec import GridSpec
import pandas as pd
P=list(csv.DictReader(open(ROOT/'data'/'phenotype'/'block_observations.csv'))); M=list(csv.DictReader(open(PHENO/'c22_M.csv',encoding='utf-8')))
_an=pd.read_csv(PHENO/'c22_A.csv'); _an=_an[_an.term=='genotype'].set_index('trait').P   # ANOVA from analysis/phenotype/heterosis_rcbd_lm.R
def _sci(p):
    e=int(np.floor(np.log10(p))); return f'{p/10**e:.2f}×10$^{{{e}}}$'
order=['C','D','MD','S','MS']; mk={'G1':'o','G2':'s','G3':'^'}; off={'G1':-0.14,'G2':0,'G3':0.14}
gp={k:_sci(_an[k]) for k in ('Yield_per_plant_g','Fruits_per_plant','Fruit_mass_g')}
tr=[('Yield_per_plant_g','单株产量 Yield per plant (g)',lambda r:1000*float(r['weight_kg'])/int(r['plants'])),
    ('Fruits_per_plant','单株瓜数 Fruits per plant',lambda r:int(r['fruits'])/int(r['plants'])),
    ('Fruit_mass_g','单瓜重 Fruit mass (g)',lambda r:1000*float(r['weight_kg'])/int(r['fruits']))]
import json
from PIL import Image
AX=[]; fig=plt.figure(figsize=(W,7.2)); gs=GridSpec(3,3,height_ratios=[2.35,1.45,0.75],hspace=0.42,wspace=0.45,left=0.1,right=0.98,top=0.98,bottom=0.07)
axP=fig.add_subplot(gs[0,:]); AX.append(axP)
img=np.asarray(Image.open(INP/'fruit_canvas.png')); pos=json.load(open(INP/'fruit_pos.json'))
axP.imshow(img,interpolation='lanczos'); axP.set_xticks([]); axP.set_yticks([])
for s_ in axP.spines.values(): s_.set_visible(False)
for k,(x,y) in pos.items():
    axP.text(x,y,k,color='white',ha='center',va='bottom',fontsize=8,fontweight='bold')
    axP.plot([x-48,x+48],[y+7,y+7],color=GCOL[k] if k!='C' else '#BDBDBD',lw=2,solid_capstyle='butt')
PXCM=719/30.0   # user-supplied black-background photo: MD 12 d fruit (30 cm) = 719 px green-tissue centreline
Hc,Wc=img.shape[:2]; bx=Wc-70-5*PXCM; by=Hc-48
axP.plot([bx,bx+5*PXCM],[by,by],color='white',lw=2.2,solid_capstyle='butt')
axP.text(bx+2.5*PXCM,by-11,'5 cm',color='white',ha='center',va='bottom',fontsize=6.8)
letter(axP,'A',x=-0.02,y=1.0)
for k,((key,lab,fn),L) in enumerate(zip(tr,'BCD')):
    ax=fig.add_subplot(gs[1,k]); AX.append(ax); mm={m['genotype']:m for m in M if m['trait']==key}
    top=max(max(fn(r) for r in P),max(float(m['lsmean'])+float(m['SE']) for m in mm.values()))
    for i,g in enumerate(order):
        mu=float(mm[g]['lsmean']); se=float(mm[g]['SE'])
        ax.bar(i,mu,width=0.64,color=GCOL[g],alpha=0.30,lw=0,zorder=1)
        ax.errorbar(i,mu,yerr=se,color='black',capsize=2,lw=0.7,zorder=4)
        vals=[mu+se]
        for r in P:
            if r['genotype']==g:
                v=fn(r); vals.append(v)
                ax.scatter(i+off[r['plot']],v,marker=mk[r['plot']],s=11,color=GCOL[g],edgecolor='white',linewidth=0.3,zorder=3)
        ax.text(i,max(vals)+top*0.03,mm[g]['LSD_letter'],ha='center',va='bottom',fontsize=7)
    ax.set_ylim(0,top*1.25); ax.set_xticks(range(5)); ax.set_xticklabels(order); ax.set_ylabel(lab)
    ax.text(0.03,0.98,'基因型 Genotype $P$ = '+gp[key],transform=ax.transAxes,fontsize=5.8,va='top')
    letter(ax,L,x=-0.30)
# D: decomposition
ax=fig.add_subplot(gs[2,:]); AX.append(ax)
_m=pd.read_csv(PHENO/'c22_Mult.csv').set_index('hybrid')   # product reference from heterosis_rcbd_lm.R
_B=_m.yield_expected_from_component_MPs-_m.yield_MP; _T=_m.yield_F1-_m.yield_MP
data=[(h,round(_B[h],2),round(_T[h]-_B[h],2),round(_T[h],2),f'{100*_B[h]/_T[h]:.1f}%') for h in ('MD','MS')]
cB,cR='#7F7F7F','#C9C9C9'
for y,(h,B,R,T,pct) in zip([1,0],data):
    ax.barh(y,B,color=cB,height=0.52,edgecolor='white',lw=0.4); ax.barh(y,R,left=B,color=cR,height=0.52,edgecolor='white',lw=0.4)
    if B>60: ax.text(B/2,y,f'$B$：{B:.2f} g（占{pct}）',ha='center',va='center',fontsize=6.5,color='white')
    ax.text(T+6,y,f'{T:.2f} g'+('' if B>60 else f'　$B$：{B:.2f} g（占{pct}）'),va='center',fontsize=6.5)
ax.set_yticks([1,0]); ax.set_yticklabels(['MD','MS']); [t.set_color(GCOL[t.get_text()]) for t in ax.get_yticklabels()]
ax.set_xlim(0,560); ax.set_xlabel('F$_1$产量与产量中亲值之差 Hybrid − mid-parent yield (g per plant)')
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=cB,label='乘积参照差额 $B$  Product-reference gap'),Patch(color=cR,label='F$_1$相对乘积参照的净差额 $R$  Remainder')],
          loc='lower left',frameon=False,ncol=2,bbox_to_anchor=(0.0,1.0),columnspacing=1.2)
letter(ax,'E',x=-0.075,y=1.18)
h=[plt.Line2D([],[],marker=mk[k],ls='',color='#555555',markersize=4,label=k) for k in mk]
AX[1].legend(handles=h,title=None,loc='upper left',bbox_to_anchor=(0.0,-0.13),ncol=3,frameon=False,handletextpad=0.2,columnspacing=0.6,fontsize=6.3)
AX[1].text(-0.02,-0.205,'区组 Block',transform=AX[1].transAxes,ha='right',va='center',fontsize=6.3)
save(fig,'Fig1_yield_components'); save_panels(fig,'Fig1',dict(zip('ABCDE',[[a] for a in AX])))
