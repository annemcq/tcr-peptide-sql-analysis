-- Compare cognate ranks for original versus composition-matched decoys.
WITH ranked AS (
 SELECT scenario, pdb_id, label,
        RANK() OVER (PARTITION BY scenario,pdb_id ORDER BY sql_score) AS first_rank,
        COUNT(*) OVER (PARTITION BY scenario,pdb_id,sql_score) AS tie_n
 FROM candidate_scores
), positive AS (
 SELECT scenario, pdb_id,
        first_rank + (tie_n - 1)/2.0 AS cognate_rank
 FROM ranked WHERE label=1
)
SELECT orig.pdb_id, orig.cognate_rank AS original_rank,
       matched.cognate_rank AS matched_rank,
       matched.cognate_rank - orig.cognate_rank AS change_in_rank
FROM positive orig
JOIN positive matched ON orig.pdb_id = matched.pdb_id
WHERE orig.scenario='original' AND matched.scenario='composition_matched'
ORDER BY ABS(change_in_rank) DESC, orig.pdb_id;

-- Summary: performance of TCRen baseline ONLY, not ML model performance.
WITH ranked AS (
 SELECT scenario,pdb_id,label,
        RANK() OVER (PARTITION BY scenario,pdb_id ORDER BY sql_score) AS first_rank,
        COUNT(*) OVER (PARTITION BY scenario,pdb_id,sql_score) AS tie_n
 FROM candidate_scores
), pos AS (
 SELECT scenario,pdb_id,first_rank+(tie_n-1)/2.0 AS rank
 FROM ranked WHERE label=1
)
SELECT scenario, COUNT(*) AS n_tcr,
       ROUND(AVG(1.0/rank),5) AS mean_reciprocal_rank,
       ROUND(AVG(CASE WHEN rank=1 THEN 1.0 ELSE 0.0 END),5) AS top1,
       ROUND(AVG(CASE WHEN rank<=5 THEN 1.0 ELSE 0.0 END),5) AS top5
FROM pos GROUP BY scenario;
