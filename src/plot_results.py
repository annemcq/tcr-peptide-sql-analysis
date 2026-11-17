"""Generate honest, data-derived figures after building the SQLite database."""
from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

RANK_COMPARE = """
WITH ranked AS (
 SELECT scenario,pdb_id,label,
   RANK() OVER (PARTITION BY scenario,pdb_id ORDER BY sql_score ASC) AS first_rank,
   COUNT(*) OVER (PARTITION BY scenario,pdb_id,sql_score) AS tie_count
 FROM candidate_scores
), positives AS (
 SELECT scenario,pdb_id,first_rank+(tie_count-1)/2.0 AS rank
 FROM ranked WHERE label=1
)
SELECT a.pdb_id,a.rank,b.rank
FROM positives a JOIN positives b ON a.pdb_id=b.pdb_id
WHERE a.scenario='original' AND b.scenario='composition_matched'
ORDER BY a.pdb_id
"""


def plot(db: Path, out: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    out.mkdir(parents=True,exist_ok=True)
    with sqlite3.connect(db) as con:
        rows=con.execute(RANK_COMPARE).fetchall()
        vals=con.execute("SELECT sql_score,python_score FROM candidate_scores WHERE python_score IS NOT NULL").fetchall()
    if rows:
        original=np.array([r[1] for r in rows]); matched=np.array([r[2] for r in rows])
        fig,ax=plt.subplots(figsize=(6.5,5.2))
        ax.scatter(original,matched,alpha=.75)
        lim=max(2,float(max(original.max(),matched.max())))
        ax.plot([0,lim],[0,lim],linestyle='--',linewidth=1,alpha=.5)
        ax.set(xlabel='Cognate rank: original negatives',ylabel='Cognate rank: composition-matched negatives',
               title='TCRen rank sensitivity to negative-set design')
        ax.set_xlim(0,lim*1.02);ax.set_ylim(0,lim*1.02)
        fig.tight_layout();fig.savefig(out/'rank_sensitivity.png',dpi=180);plt.close(fig)
        print('Saved rank_sensitivity.png')
    if vals:
        values=np.asarray(vals,float)
        fig,ax=plt.subplots(figsize=(6.5,5.2))
        ax.scatter(values[:,1],values[:,0],s=12,alpha=.5)
        lo=float(values.min());hi=float(values.max())
        ax.plot([lo,hi],[lo,hi],linestyle='--',linewidth=1,alpha=.5)
        err=float(np.max(np.abs(values[:,0]-values[:,1])))
        ax.set(xlabel='Python TCRen score',ylabel='SQL-reconstructed TCRen score',
               title=f'SQL versus Python scoring parity (max absolute error: {err:.2g})')
        fig.tight_layout();fig.savefig(out/'sql_python_parity.png',dpi=180);plt.close(fig)
        print('Saved sql_python_parity.png')
    else:
        print('No Python reference scores provided: parity figure not generated.')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--db',type=Path,default=ROOT/'results'/'tcr_analysis.sqlite')
    p.add_argument('--out',type=Path,default=ROOT/'figures')
    a=p.parse_args();plot(a.db,a.out)
if __name__=='__main__':main()
