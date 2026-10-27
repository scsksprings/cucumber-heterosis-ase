import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import os
from pathlib import Path
ROOT=Path(os.environ.get('REPO_ROOT',Path(__file__).resolve().parents[2]))   # repository root
RES=ROOT/'results'                      # full result tables (Tables S1-S8, complete version)
INP=ROOT/'figures'/'inputs'             # photo canvas for Fig. 1A
OUT=Path(os.environ.get('FIG_OUT',ROOT/'figures'/'output')); (OUT/'supplementary').mkdir(parents=True,exist_ok=True)
PHENO=Path(os.environ.get('PHENO_OUT',OUT/'phenotype'))   # output of analysis/phenotype/heterosis_rcbd_lm.R
FONTS=ROOT/'figures'/'fonts'
FP=str(FONTS/'NotoSansSCtt-Regular.ttf')  # TrueType (glyf) conversion of Noto Sans CJK SC for clean Type 42 embedding
if os.path.exists(FP):
    fm.fontManager.addfont(FP); fm.fontManager.addfont(str(FONTS/'NotoSansSCtt-Bold.ttf'))
    FN=fm.FontProperties(fname=FP).get_name()
else:  # fallback: any installed Noto Sans CJK SC (layout may differ slightly)
    FN='Noto Sans CJK SC'
plt.rcParams.update({'font.family':FN,'font.size':7,'axes.titlesize':7.5,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,
 'legend.fontsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
 'xtick.major.size':2.5,'ytick.major.size':2.5,'axes.unicode_minus':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':600})
GCOL={'C':'#6E6E6E','D':'#56A8D0','MD':'#0072B2','S':'#009E73','MS':'#D55E00'}
CTX=['seedling_leaf','fruiting_leaf','fruit_3d','fruit_6d','fruit_9d']
CTXLAB={'seedling_leaf':'苗期叶片\nSL','fruiting_leaf':'结瓜期叶片\nFL','fruit_3d':'3 d果实\nF3','fruit_6d':'6 d果实\nF6','fruit_9d':'9 d果实\nF9'}
W=5.9  # inches = 15 cm final width
import re as _re
def itg(s):
    """italicise gene symbols/IDs with mathtext (CJK font has no italic); keeps '-like' upright"""
    return _re.sub(r'[A-Za-z0-9_.]+', lambda m: m.group(0) if m.group(0)=='like' else '$\\mathit{'+m.group(0).replace('_','\\_')+'}$', s)
def letter(ax,L,x=-0.12,y=1.06):
    ax.text(x,y,L,transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom',ha='left')
SVG_FALLBACK="'Noto Sans SC TT', 'Noto Sans CJK SC', 'Source Han Sans SC', 'Microsoft YaHei', 'PingFang SC', 'SimHei', sans-serif"
def svg_font_fallback(data):
    """SVG text stays editable (svg.fonttype='none'); add common CJK fallbacks so the text renders
    sensibly where the bundled font is not installed (e.g. when the SVG is opened in Word)."""
    return data.replace(b"'Noto Sans SC TT'", SVG_FALLBACK.encode())
def save(fig,name,outdir=None):
    outdir=OUT if outdir is None else Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    from io import BytesIO
    for ext in ('png','pdf','svg'):
        buffer=BytesIO()
        fig.savefig(buffer,format=ext,bbox_inches='tight')
        data=buffer.getvalue()
        if ext=='svg': data=svg_font_fallback(data)
        (outdir/f'{name}.{ext}').write_bytes(data)
        buffer.close()

from matplotlib.transforms import Bbox
def panel_bbox(fig,artists,pad=0.04):
    r=fig.canvas.get_renderer()
    bbs=[a.get_tightbbox(r) if hasattr(a,'get_tightbbox') else a.get_window_extent(r) for a in artists]
    bb=Bbox.union([b for b in bbs if b is not None and b.width>0])
    return bb.transformed(fig.dpi_scale_trans.inverted()).padded(pad)
def check_no_overlap(fig,groups):
    """groups: dict letter->artists; raise if panel bounding boxes intersect"""
    bbs={k:panel_bbox(fig,v,pad=0) for k,v in groups.items()}
    ks=list(bbs); bad=[]
    for i in range(len(ks)):
        for j in range(i+1,len(ks)):
            a,b=bbs[ks[i]],bbs[ks[j]]
            if a.overlaps(b):
                ix=min(a.x1,b.x1)-max(a.x0,b.x0); iy=min(a.y1,b.y1)-max(a.y0,b.y0)
                if ix>0.01 and iy>0.01: bad.append((ks[i],ks[j],round(ix,3),round(iy,3)))
    return bad
def save_panels(fig,name,groups,outdir=None):
    """Export each panel as its own PDF (and PNG) by hiding all other artists."""
    outdir=OUT/'panels' if outdir is None else Path(outdir); os.makedirs(outdir,exist_ok=True)
    fig.canvas.draw()
    allart=list(fig.axes)+list(fig.texts)+list(fig.legends)
    for L,arts in groups.items():
        bb=panel_bbox(fig,arts)
        keep=set(map(id,arts))
        vis={id(a):a.get_visible() for a in allart}
        for a in allart:
            if id(a) not in keep: a.set_visible(False)
        fig.savefig(f'{outdir}/{name}_{L}.pdf',bbox_inches=bb)
        fig.savefig(f'{outdir}/{name}_{L}.png',bbox_inches=bb,dpi=600)
        for a in allart: a.set_visible(vis[id(a)])
