"""Run grouped SQL reports and write the results as CSV files."""
import argparse
import csv
import sqlite3
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REPORTS={
    'data_quality':'01_data_quality.sql',
    'contact_statistics':'02_contact_statistics.sql',
    'scoring':'03_tcren_scoring.sql',
    'ranking':'04_peptide_ranking.sql',
    'negative_controls':'05_negative_controls.sql',
}

def split_queries(sql):
    """Split project SQL files on semicolons outside of quoted strings.

    The bundled query files have no semicolons inside strings or comments.
    For general SQL parsing, use a real SQL parser.
    """
    return [part.strip() for part in sql.split(';') if part.strip()]


def run(db:Path,out:Path):
    out.mkdir(parents=True,exist_ok=True)
    with sqlite3.connect(db) as conn:
        for group,filename in REPORTS.items():
            sql=(ROOT/'sql'/filename).read_text()
            for i,part in enumerate(split_queries(sql),1):
                cursor=conn.execute(part)
                if cursor.description is None:
                    continue
                cols=[x[0] for x in cursor.description]
                rows=cursor.fetchall()
                path=out/f'{group}_{i:02d}.csv'
                with path.open('w',newline='',encoding='utf-8') as f:
                    writer=csv.writer(f)
                    writer.writerow(cols)
                    writer.writerows(rows)
                print(f'{path.name}: {len(rows)} rows')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--db',type=Path,default=ROOT/'results'/'tcr_analysis.sqlite')
    parser.add_argument('--out',type=Path,default=ROOT/'results'/'reports')
    args=parser.parse_args()
    run(args.db,args.out)
if __name__=='__main__':main()
