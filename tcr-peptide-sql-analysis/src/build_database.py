"""Build a normalized SQLite database from the PUBLIC TCRen portfolio data.

This script does not redistribute third-party input CSVs. Pass paths to locally
obtained copies of the source files from the related tcr-peptide-ranking repo.
"""
from __future__ import annotations

import argparse
import csv
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))


def norm_id(x: str) -> str:
    value = (x or '').strip().lower()
    if not value:
        raise ValueError('Empty PDB identifier')
    return value


def norm_aa(x: str) -> str:
    aa = (x or '').strip().upper()
    if len(aa) != 1 or aa not in AMINO_ACIDS:
        raise ValueError(f'Unknown amino acid: {x!r}')
    return aa


def create_db(db: Path, *, original: Path, matched: Path,
              contacts: Path, potential: Path,
              reference_original: Path | None = None,
              reference_matched: Path | None = None) -> dict[str, int]:
    """Return record counts; fail fast on incomplete/ambiguous inputs."""
    if db.exists():
        raise FileExistsError(f'{db} already exists (delete it to rebuild)')
    sets = {'original': load_csv(original),
            'composition_matched': load_csv(matched)}
    source_contacts = load_csv(contacts)
    source_potentials = load_csv(potential)

    reference = {}
    for scenario, path in [('original',reference_original),
                           ('composition_matched',reference_matched)]:
        if path:
            for row in load_csv(path):
                key = scenario, norm_id(row['pdb_id']), row['peptide'].strip().upper()
                score = float(row['tcren_score'])
                if key in reference and abs(reference[key] - score) > 1e-12:
                    raise ValueError(f'Contradictory reference score for {key}')
                reference[key] = score

    candidate_rows = []
    key_counter = Counter()
    for scenario, rows in sets.items():
        for row in rows:
            pdb = norm_id(row['pdb_id'])
            pep = row['peptide'].strip().upper()
            if not pep or not set(pep) <= AMINO_ACIDS:
                raise ValueError(f'Bad peptide {pep!r} in {pdb}')
            label = int(row['label'])
            if label not in (0,1):
                raise ValueError(f'Expected binary label for {pdb}, got {label}')
            key = scenario, pdb, pep
            key_counter[key] += 1
            candidate_rows.append((scenario,pdb,pep,label,reference.get(key)))
        if not rows:
            raise ValueError(f'Scenario {scenario} is empty')
    dup = [k for k,n in key_counter.items() if n != 1]
    if dup:
        raise ValueError(f'Duplicate candidates: {dup[:3]}')

    grouped = defaultdict(list)
    for scenario,pdb,pep,label,score in candidate_rows:
        grouped[(scenario,pdb)].append(label)
    for group, labels in grouped.items():
        if sum(labels) != 1:
            raise ValueError(f'{group}: expected exactly 1 cognate peptide')
    if set(k[1] for k in grouped if k[0]=='original') != set(k[1] for k in grouped if k[0]=='composition_matched'):
        raise ValueError('Scenarios contain different TCR systems')
    cognates = defaultdict(dict)
    for scenario, pdb, pep, label, score in candidate_rows:
        if label:
            cognates[pdb][scenario] = pep
    for pdb, by_scenario in cognates.items():
        if len(set(by_scenario.values())) != 1:
            raise ValueError(f'{pdb}: cognate sequences differ across scenarios')

    # Python reference scoring uses a dictionary (last key wins). For this
    # relational variant we reject duplicates rather than silently choosing.
    potentials = {}
    for row in source_potentials:
        k = norm_aa(row['residue.aa.from']), norm_aa(row['residue.aa.to'])
        v = float(row['TCRen'])
        if k in potentials and potentials[k] != v:
            raise ValueError(f'Conflicting potential for {k}')
        potentials[k] = v

    system_ids = {pdb for _,pdb,_,_,_ in candidate_rows}
    native_contacts = []
    for row in source_contacts:
        if row['chain.type.to'].strip().upper() != 'PEPTIDE':
            continue
        pdb = norm_id(row['pdb.id'])
        if pdb not in system_ids:
            continue
        pos = int(row['residue.index.to'])  # ORIGINAL Python scorer uses 0-based positions.
        native_contacts.append((pdb,pos,norm_aa(row['residue.aa.from'])))
    if not native_contacts:
        raise ValueError('No peptide contacts for selected TCRs')
    contact_groups = defaultdict(list)
    for pdb,pos,aa in native_contacts:
        contact_groups[pdb].append((pos,aa))
    for scenario,pdb,pep,label,score in candidate_rows:
        for pos,aa in contact_groups[pdb]:
            if not 0 <= pos < len(pep):
                raise ValueError(f'Contact index {pos} outside peptide {pep} ({pdb})')
            if (aa,pep[pos]) not in potentials:
                raise ValueError(f'No potential for ({aa},{pep[pos]}) in {pdb}')

    db.parent.mkdir(parents=True,exist_ok=True)
    try:
        with sqlite3.connect(db) as conn:
            conn.execute('PRAGMA foreign_keys=ON')
            conn.executescript((ROOT/'sql'/'00_schema.sql').read_text())
            conn.executemany('INSERT INTO tcr_systems VALUES (?)',[(p,) for p in sorted(system_ids)])
            conn.executemany('INSERT INTO scenarios VALUES (?,?)', [
                ('original','Uniformly sampled synthetic negatives (source project)'),
                ('composition_matched','Scrambled-residue negatives preserving cognate composition')])
            conn.executemany('INSERT INTO tcren_potentials VALUES (?,?,?)',
                             [(a,b,v) for (a,b),v in sorted(potentials.items())])
            conn.executemany('INSERT INTO native_contacts(pdb_id,peptide_position,tcr_aa) VALUES (?,?,?)',native_contacts)
            for scenario,pdb,pep,label,score in candidate_rows:
                cur=conn.execute('INSERT INTO candidates(scenario,pdb_id,peptide,label,python_score) VALUES (?,?,?,?,?)',
                                 (scenario,pdb,pep,label,score))
                conn.executemany('INSERT INTO candidate_residues VALUES (?,?,?)',
                                 [(cur.lastrowid,i,aa) for i,aa in enumerate(pep)])
            unmatched=conn.execute('SELECT count(*) FROM candidates c LEFT JOIN candidate_scores s ON c.candidate_id=s.candidate_id WHERE s.candidate_id IS NULL').fetchone()[0]
            if unmatched:
                raise ValueError(f'{unmatched} candidates missing reconstructed SQL scores')
    except Exception:
        db.unlink(missing_ok=True)
        raise
    return {'systems':len(system_ids),'candidates':len(candidate_rows),
            'contacts':len(native_contacts),'potentials':len(potentials),
            'reference_scores':len(reference)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ('original','matched','contacts','potential'):
        parser.add_argument('--'+flag,required=True,type=Path)
    parser.add_argument('--reference-original',type=Path)
    parser.add_argument('--reference-matched',type=Path)
    parser.add_argument('--db',type=Path,default=ROOT/'results'/'tcr_analysis.sqlite')
    args=parser.parse_args()
    info=create_db(args.db,original=args.original,matched=args.matched,
                   contacts=args.contacts,potential=args.potential,
                   reference_original=args.reference_original,
                   reference_matched=args.reference_matched)
    print('Database created:',args.db)
    for key,value in info.items():
        print(f'  {key}: {value}')

if __name__=='__main__':
    main()
