import os,sys; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import csv
from style import *
P=str(RES)
A=list(csv.DictReader(open(P+'/S1/PCA_all_scores.tsv',encoding='utf-8'),delimiter='\t'))
B=list(csv.DictReader(open(P+'/S1/PCA_fruit_scores.tsv',encoding='utf-8'),delimiter='\t'))
I=list(csv.DictReader(open(P+'/S3/S3_expression_inheritance_summary.tsv',encoding='utf-8'),delimiter='\t'))
T=list(csv.DictReader(open(P+'/S4/TableS6_fruit_interactions.tsv',encoding='utf-8'),delimiter='\t'))
ccol={'seedling_leaf':'#1B9E77','fruiting_leaf':'#7570B3','fruit_3d':'#E6AB02','fruit_6d':'#66A61E','fruit_9d':'#E7298A'}
fig,axs=plt.subplots(2,2,figsize=(W,5.2)); plt.subplots_adjust(hspace=0.62,wspace=0.32)
ax=axs[0,0]
for c in CTX:
    xs=[float(r['PC1']) for r in A if r['stage']==c]; ys=[float(r['PC2']) for r in A if r['stage']==c]
    ax.scatter(xs,ys,s=10,color=ccol[c],alpha=0.85,lw=0,label=CTXLAB[c].replace('\n',' '))
ax.set_xlabel('PC1 (60.4%)'); ax.set_ylabel('PC2 (22.6%)'); ax.legend(frameon=False,fontsize=5.8,loc='lower left',bbox_to_anchor=(0.42,0.02),handletextpad=0.2,borderaxespad=0)
ax.set_title('全部样品 All samples',loc='left'); letter(ax,'A',x=-0.2)
ax=axs[0,1]; mkr={'fruit_3d':'o','fruit_6d':'s','fruit_9d':'^'}
for g in ['C','D','MD','S','MS']:
    for st in mkr:
        xs=[float(r['PC1']) for r in B if r['genotype']==g and r['stage']==st]; ys=[float(r['PC2']) for r in B if r['genotype']==g and r['stage']==st]
        ax.scatter(xs,ys,s=12,marker=mkr[st],color=GCOL[g],lw=0)
ax.set_xlabel('PC1 (77.5%)'); ax.set_ylabel('PC2 (8.6%)'); ax.set_title('果实样品 Fruit samples',loc='left'); letter(ax,'B',x=-0.2)
hg=[plt.Line2D([],[],marker='o',ls='',color=GCOL[g],markersize=4,label=g) for g in GCOL]
hs=[plt.Line2D([],[],marker=mkr[s],ls='',color='#555555',markersize=4,label=s.replace('fruit_','F').replace('d','')) for s in mkr]
ax.legend(handles=hg+hs,frameon=False,fontsize=5.8,ncol=2,loc='center left',bbox_to_anchor=(0.36,0.55),handletextpad=0.1,columnspacing=0.6)
# C stacked
ax=axs[1,0]
cls=[('above_both_parents','超高亲 Above both','#5B2C83'),('high_parent_like','高亲偏向 High-parent-like','#B39DDB'),('low_parent_like','低亲偏向 Low-parent-like','#80CBC4'),('below_both_parents','超低亲 Below both','#00695C'),('other_departure','其他 Other','#BDBDBD')]
x=0; xt=[]; xl=[]
for c in CTX:
    for h in ['MD','MS']:
        bottom=0
        for k,lab,col in cls:
            v=sum(int(r['genes']) for r in I if r['genotype']==h and r['stage']==c and r['publication_class']==k)
            ax.bar(x,v,bottom=bottom,color=col,width=0.8,lw=0); bottom+=v
        ax.text(x,bottom+25,str(bottom),ha='center',fontsize=5.5)
        xt.append(x); xl.append(h); x+=1
    x+=0.6
ax.set_xticks(xt); ax.set_xticklabels(xl,fontsize=5.8)
for i,c in enumerate(CTX):
    ax.text((xt[2*i]+xt[2*i+1])/2,-0.13,CTXLAB[c].split('\n')[1],transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.5)
ax.set_ylabel('表达偏离基因数\nDeparture genes'); ax.set_ylim(0,2000)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=c,label=l) for k,l,c in cls],frameon=False,fontsize=5.6,loc='upper right',handlelength=1)
ax.set_title('表达偏离模式 Departure patterns',loc='left'); letter(ax,'C',x=-0.2)
# D
ax=axs[1,1]; comps=['fruit_6d_vs_fruit_3d','fruit_9d_vs_fruit_3d','fruit_9d_vs_fruit_6d']; labs=['6 d vs 3 d','9 d vs 3 d','9 d vs 6 d']
for j,h in enumerate(['MD','MS']):
    vals=[int([r for r in T if r['hybrid']==h and r['comparison']==c][0]['significant_interactions']) for c in comps]
    xs=[i+(j-0.5)*0.36 for i in range(3)]
    ax.bar(xs,vals,width=0.34,color=GCOL[h],label=h)
    for xx,v in zip(xs,vals): ax.text(xx,v+6,str(v),ha='center',fontsize=5.8)
ax.set_xticks(range(3)); ax.set_xticklabels(labs); ax.set_ylabel('偏离变化基因数\nGenes with changed departure'); ax.legend(frameon=False,loc='upper left')
ax.set_title('果实取样时期间变化\nChanges between sampling times',loc='left'); letter(ax,'D',x=-0.2)
save(fig,'Fig2_expression_patterns'); save_panels(fig,'Fig2',{'A':[axs[0,0]],'B':[axs[0,1]],'C':[axs[1,0]],'D':[axs[1,1]]})
