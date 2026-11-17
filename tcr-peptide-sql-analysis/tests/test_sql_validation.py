import csv
import sqlite3
from pathlib import Path

import pytest

from src.build_database import create_db
from src.run_analysis import run

ROOT=Path(__file__).resolve().parents[1]


def write_csv(path, rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fixture_sources(tmp_path):
    # Two systems, each with one positive and two negatives in two scenarios.
    data={}
    for scenario,file in [('original','original.csv'),('composition_matched','matched.csv')]:
        data[scenario]=tmp_path/file
        rows=[]
        for pdb,seqs in [('1aaa',['AA','AC','CA']),('2bbb',['CC','CA','AC'])]:
            for i,pep in enumerate(seqs):
                rows.append({'pdb_id':pdb,'peptide':pep,'label':int(i==0)})
        write_csv(data[scenario],rows)
    data['contacts']=tmp_path/'contacts.csv'
    write_csv(data['contacts'],[
        {'pdb.id':pdb,'chain.type.to':'PEPTIDE','residue.index.to':str(pos),'residue.aa.from':aa}
        for pdb,pos,aa in [('1aaa',0,'A'),('1aaa',1,'A'),('2bbb',0,'C'),('2bbb',1,'C')]])
    data['potential']=tmp_path/'potential.csv'
    write_csv(data['potential'],[
        {'residue.aa.from':a,'residue.aa.to':b,'TCRen':str(v)}
        for a,b,v in [('A','A',-2),('A','C',1),('C','A',1),('C','C',-2)]])
    return data


@pytest.fixture
def sample_db(tmp_path):
    data=fixture_sources(tmp_path)
    db=tmp_path/'analysis.sqlite'
    counts=create_db(db,original=data['original'],matched=data['composition_matched'],
                     contacts=data['contacts'],potential=data['potential'])
    assert counts['candidates']==12
    return db


def test_sql_scoring_reconstructs_contact_sum(sample_db):
    with sqlite3.connect(sample_db) as conn:
        row=conn.execute("SELECT sql_score,n_contacts FROM candidate_scores WHERE scenario='original' AND pdb_id='1aaa' AND peptide='AA'").fetchone()
        assert row==(-4.0,2)
        assert conn.execute('SELECT COUNT(*) FROM candidate_scores').fetchone()[0]==12


def test_grouped_ranks_and_mrr(sample_db):
    sql=(ROOT/'sql'/'04_peptide_ranking.sql').read_text()
    with sqlite3.connect(sample_db) as conn:
        rows=conn.execute(sql).fetchall()
        assert len(rows)==4
        assert all(row[4]==1 for row in rows)  # both cognates rank first


def test_reports_run_and_have_headers(sample_db,tmp_path):
    out=tmp_path/'reports'
    run(sample_db,out)
    files=list(out.glob('*.csv'))
    assert len(files)==12
    assert all(f.read_text().strip() for f in files)


def test_duplicate_candidates_rejected(tmp_path):
    data=fixture_sources(tmp_path)
    with data['original'].open('a') as f:
        f.write('1aaa,AA,1\n')
    with pytest.raises(ValueError,match='Duplicate candidates'):
        create_db(tmp_path/'out.sqlite',original=data['original'],matched=data['composition_matched'],
                  contacts=data['contacts'],potential=data['potential'])


def test_missing_contact_lookup_fails_loudly(tmp_path):
    data=fixture_sources(tmp_path)
    write_csv(data['potential'],[{'residue.aa.from':'A','residue.aa.to':'A','TCRen':'-2'}])
    with pytest.raises(ValueError,match='No potential'):
        create_db(tmp_path/'out.sqlite',original=data['original'],matched=data['composition_matched'],
                  contacts=data['contacts'],potential=data['potential'])


def test_tie_aware_average_rank(tmp_path):
    db=tmp_path/'tie.sqlite'
    with sqlite3.connect(db) as conn:
        conn.executescript((ROOT/'sql'/'00_schema.sql').read_text())
        conn.execute("INSERT INTO tcr_systems VALUES ('X')")
        conn.execute("INSERT INTO scenarios VALUES ('original','fixture')")
        for i,(pep,label) in enumerate([('A',1),('C',0),('G',0)],1):
            conn.execute('INSERT INTO candidates VALUES (?, ?, ?, ?, ?, ?)',(i,'original','X',pep,label,None))
            conn.execute('INSERT INTO candidate_residues VALUES (?,0,?)',(i,pep))
        conn.execute("INSERT INTO native_contacts(pdb_id,peptide_position,tcr_aa) VALUES('X',0,'A')")
        for aa in 'ACG':
            conn.execute('INSERT INTO tcren_potentials VALUES (?,?,?)',('A',aa,0.0))
        rows=conn.execute((ROOT/'sql'/'04_peptide_ranking.sql').read_text()).fetchall()
        assert rows[0][4]==2.0  # three equal scores -> ranks 1,2,3 average 2


def test_mismatched_cognate_rejected(tmp_path):
    data=fixture_sources(tmp_path)
    rows=list(csv.DictReader(data['composition_matched'].open()))
    rows[0]['peptide']='GG'  # label 1 no longer matches original cognate
    write_csv(data['composition_matched'],rows)
    with pytest.raises(ValueError,match='cognate sequences differ'):
        create_db(tmp_path/'out.sqlite',original=data['original'],matched=data['composition_matched'],
                  contacts=data['contacts'],potential=data['potential'])


def test_reference_parity(sample_db):
    with sqlite3.connect(sample_db) as conn:
        conn.execute("UPDATE candidates SET python_score=-4.0 WHERE scenario='original' AND pdb_id='1aaa' AND peptide='AA'")
        rows=conn.execute((ROOT/'sql'/'03_tcren_scoring.sql').read_text().split(';')[1]).fetchall()
        assert any(row[0]=='original' and row[1]==1 and row[2]==0.0 for row in rows)
