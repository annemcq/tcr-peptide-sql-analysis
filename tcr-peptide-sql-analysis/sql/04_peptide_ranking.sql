-- SQL RANK alone does NOT reproduce pandas rank(method='average').
-- If a tie starts at rank k and spans m rows, average occupied rank = k+(m-1)/2.
WITH ranks AS (
  SELECT scenario, pdb_id, candidate_id, peptide, label, sql_score,
         RANK() OVER (
           PARTITION BY scenario, pdb_id ORDER BY sql_score ASC
         ) AS start_rank,
         COUNT(*) OVER (
           PARTITION BY scenario, pdb_id, sql_score
         ) AS tie_count
  FROM candidate_scores
),
positive_ranks AS (
  SELECT scenario, pdb_id, candidate_id, sql_score,
         start_rank + (tie_count - 1) / 2.0 AS mean_tie_rank
  FROM ranks WHERE label = 1
)
SELECT scenario, pdb_id, candidate_id, sql_score, mean_tie_rank,
       1.0 / mean_tie_rank AS reciprocal_rank,
       CASE WHEN mean_tie_rank = 1 THEN 1 ELSE 0 END AS top1,
       CASE WHEN mean_tie_rank <= 5 THEN 1 ELSE 0 END AS top5
FROM positive_ranks
ORDER BY scenario, pdb_id;
