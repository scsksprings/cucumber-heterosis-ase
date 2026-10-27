#!/usr/bin/env python3
"""Update the release manifest and create a source/data ZIP outside the repository."""
from pathlib import Path
import argparse,hashlib,json,os,re,zipfile

ROOT=Path(__file__).resolve().parents[1]
ROOT_FILES={'README.md','LICENSE','CITATION.cff','.gitignore','.gitattributes','run_all.sh',
            'requirements.txt','requirements-dev.txt','requirements-tested.txt','requirements-figures.txt','REPRODUCTION.md'}
FOLDERS={'analysis','data','results','docs','figures','tables','pipeline','scripts','tests','validation','.github'}
SKIP={'__pycache__','.cache','.pytest_cache','.snakemake','.git','.venv','output','work','Cleandata','reference'}


def release_files():
    paths=[ROOT/n for n in ROOT_FILES if (ROOT/n).is_file()]
    for folder in sorted(FOLDERS):
        for parent,dirs,files in os.walk(ROOT/folder):
            dirs[:]=sorted(d for d in dirs if d not in SKIP and not d.startswith('.'))
            for name in sorted(files):
                p=Path(parent)/name
                if p.suffix in {'.pyc','.pyo'} or name=='config.local.yaml':continue
                if name.startswith(('.env','.DS_Store')) or p.suffix in {'.pem','.key','.bam','.fastq','.fq','.docx'}:
                    raise ValueError(f'Unexpected private/raw file in release tree: {p.relative_to(ROOT)}')
                if p.is_symlink():raise ValueError(f'Symlink not allowed in portable release: {p}')
                if p.stat().st_size>95*1024**2:raise ValueError(f'File exceeds this release size policy: {p}')
                paths.append(p)
    return sorted(set(paths),key=lambda p:p.relative_to(ROOT).as_posix())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT.parent/'cucumber-heterosis-ase-github-ready-20261026.zip')
    a=parser.parse_args()
    if a.output.resolve().is_relative_to(ROOT):parser.error('ZIP output must be outside the repository')
    paths=release_files()
    manifest=ROOT/'checksums.sha256'
    manifest.write_text(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT).as_posix()}\n' for p in paths))
    paths.append(manifest)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for p in paths:archive.write(p,arcname='cucumber-heterosis-ase/'+p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(a.output) as archive:
        bad=archive.testzip()
        if bad:raise ValueError(f'ZIP CRC failure: {bad}')
    digest=hashlib.sha256(a.output.read_bytes()).hexdigest()
    a.output.with_suffix('.zip.sha256').write_text(f'{digest}  {a.output.name}\n')
    print(json.dumps({'archive':str(a.output),'files':len(paths),'bytes':a.output.stat().st_size,'sha256':digest},indent=2))

if __name__=='__main__':main()
