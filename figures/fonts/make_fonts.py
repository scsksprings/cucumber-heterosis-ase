"""Rebuild figures/fonts/NotoSansSCtt-*.ttf from the Noto Sans CJK collection.
Usage: python make_fonts.py NotoSansCJK-Regular.ttc 2 NotoSansSCtt-Regular.ttf "Noto Sans SC TT"
       python make_fonts.py NotoSansCJK-Bold.ttc 2 NotoSansSCtt-Bold.ttf "Noto Sans SC TT"
(face index 2 = Simplified Chinese in the official TTC; check with fontTools if your build differs)."""
import sys
from fontTools.ttLib import TTCollection, newTable
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools import subset
src,idx,out,fam=sys.argv[1],int(sys.argv[2]),sys.argv[3],sys.argv[4]
f=TTCollection(src).fonts[idx]
# subset to needed ranges first
uni=list(range(0x20,0x250))+list(range(0x370,0x400))+list(range(0x2000,0x2BFF))+list(range(0x3000,0x30FF))+list(range(0x4E00,0xA000))+list(range(0xFF00,0xFFF0))
opt=subset.Options(); opt.layout_features=['*']; opt.name_IDs=['*']; opt.notdef_outline=True; opt.glyph_names=False
s=subset.Subsetter(opt); s.populate(unicodes=uni); s.subset(f)
gs=f.getGlyphSet(); order=f.getGlyphOrder()
glyf={}
for n in order:
    p=TTGlyphPen(gs); gs[n].draw(Cu2QuPen(p,1.0,reverse_direction=True)); glyf[n]=p.glyph()
del f['CFF ']
if 'VORG' in f: del f['VORG']
g=newTable('glyf'); g.glyphOrder=order; g.glyphs=glyf; f['glyf']=g
f['loca']=newTable('loca')
f['maxp'].tableVersion=0x00010000
for a in ['maxZones','maxTwilightPoints','maxStorage','maxFunctionDefs','maxInstructionDefs','maxStackElements','maxSizeOfInstructions','maxComponentElements']: setattr(f['maxp'],a,0)
f['maxp'].maxZones=1
f['head'].indexToLocFormat=0; f['head'].glyphDataFormat=0
post=f['post']; post.formatType=2.0; post.extraNames=[]; post.mapping={}; post.glyphOrder=order
f.sfntVersion='\x00\x01\x00\x00'
nm=f['name']
for rec in nm.names:
    if rec.nameID in (1,16): rec.string=fam
    if rec.nameID in (4,): rec.string=fam+(' Bold' if 'Bold' in out else ' Regular')
    if rec.nameID==6: rec.string=fam.replace(' ','')+('-Bold' if 'Bold' in out else '-Regular')
f.save(out); print('saved',out)
