"""Fetch public inputs from the original GitHub repo without vendoring them.

Downloaded files remain local and are not bundled in this repository. Users
should verify original dataset provenance and licensing before redistribution.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parents[1]
SOURCE='https://raw.githubusercontent.com/annemcq/tcr-peptide-ranking/'
FILES={
 'original':'processed/first_dataset_pairs.csv',
 'matched':'processed/matched_composition_pairs.csv',
 'contacts':'data/contact_maps_PDB.csv',
 'potential':'data/TCRen_potential.csv',
 'reference_original':'processed/data_tcren_score.csv',
 'reference_matched':'processed/data_tcren_features_matched_composition.csv',
}


def fetch(out:Path, ref:str='main', only:str|None=None):
    if '/' in ref or not ref:
        raise ValueError('ref must be a branch name or commit SHA without slashes')
    out.mkdir(parents=True,exist_ok=True)
    selected = {only: FILES[only]} if only else FILES
    for label,source_path in selected.items():
        target=out/f'{label}.csv'
        if target.exists():
            raise FileExistsError(f'{target} exists; remove explicitly to refresh')
        url=SOURCE+ref+'/'+source_path
        try:
            with urlopen(url,timeout=45) as response:
                contents=response.read()
            if not contents.startswith(b'pdb_id,') and label in ('original','matched','reference_original','reference_matched'):
                raise ValueError(f'{source_path} does not appear to be a CSV with pdb_id')
            target.write_bytes(contents)
            print(f'{label:20s} {target}: {len(contents)} bytes')
        except Exception:
            target.unlink(missing_ok=True)
            raise


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=ROOT/'data'/'raw')
    parser.add_argument('--ref',default='main',help='Source branch or commit SHA (pin for reproducibility)')
    parser.add_argument('--only',choices=sorted(FILES),help='Fetch only one named source file')
    args=parser.parse_args()
    fetch(args.out,args.ref,args.only)
if __name__=='__main__':main()
